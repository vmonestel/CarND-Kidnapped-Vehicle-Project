/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define NUM_PARTICLES 50

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Check is_initialized flag
	if (is_initialized) {
		return;
	}

	// Create the class to generate random numbers
	default_random_engine generator;

	// Set the number of particles
	num_particles = NUM_PARTICLES;

	// Get the standard deviations from input
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	// Create normal distributions for x,y,theta based on GPS input
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	// Initialize all the particles with the normal distribution with mean on GPS values.
	for (int i = 0; i < num_particles; i++) {
		// Create a new particle instance
		Particle new_particle;
		// Assign particle values (id,x,y,theta,weight)
		new_particle.id = i;
		new_particle.x = dist_x(generator);
		new_particle.y = dist_y(generator);
		new_particle.theta = dist_theta(generator);
		new_particle.weight = 1.0;
		// Add the new particle to the vector
	    particles.push_back(new_particle);
	    // Add the particle weight to the vector
	    weights.push_back(new_particle.weight);
	}

	// Set initialized flag to true
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Create the class to generate random numbers
	default_random_engine generator;

	// Theta to predict the state
	double theta;

	// Create normal distributions for x,y,theta based on GPS input
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	// Calculate the new state (x,y,theta)
	for (int i = 0; i < num_particles; i++) {
		theta = particles[i].theta;
		if (fabs(yaw_rate) < 0.0001 ) {
			// yaw is not changing (yaw rate is very near to 0)
			particles[i].x += velocity * delta_t * cos(theta);
			particles[i].y += velocity * delta_t * sin(theta);
		} else {
			// yaw rate is different than 0
			particles[i].x += velocity / yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
			particles[i].y += velocity / yaw_rate * (-cos(theta + yaw_rate * delta_t) + cos(theta));
			particles[i].theta += yaw_rate * delta_t;
	    }

	    // Add random Noise
	    particles[i].x += dist_x(generator);
	    particles[i].y += dist_y(generator);
	    particles[i].theta += dist_theta(generator);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	/*
	 * Associate observations in map co-ordinates to predicted landmarks using nearest neighbor algorithm.
	 * Here, the number of observations may be less than the total number of landmarks as some of the landmarks
	 * may be outside the range of vehicle's sensor.
	 *
	 */

	for (unsigned int observation_idx = 0; observation_idx < observations.size(); observation_idx++) {
		// Initialize the minimum distance to a big number.
		double min_distance = numeric_limits<double>::max();

		// Initialize the landmark id to something invalid
		int landmark_id = -1;
		double observation_x = observations[observation_idx].x;
		double observation_y = observations[observation_idx].y;
		for (unsigned int prediction_idx = 0; prediction_idx < predicted.size(); prediction_idx++) {
			// Compute the distance between the oservation and the prediction
			double distance = dist(observation_x, observation_y,
					predicted[prediction_idx].x, predicted[prediction_idx].y);
			// Check if distance is less than the previous min_distance
			if (distance < min_distance) {
				min_distance = distance;
				landmark_id = predicted[prediction_idx].id;
			}
		}
		// Update the observation identifier.
		observations[observation_idx].id = landmark_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {

	double weight_normalizer = 0.0;
	for (int particle_index = 0; particle_index < num_particles; particle_index++) {
		// Get particle data
		double particle_x = particles[particle_index].x;
		double particle_y = particles[particle_index].y;
		double particle_theta = particles[particle_index].theta;

		// Map the observations from vehicle to map coordinates.
		vector<LandmarkObs> transformed_observations;
		for (unsigned int observation_idx = 0; observation_idx < observations.size(); observation_idx++) {
			LandmarkObs transformed_obs;
			transformed_obs.id = observation_idx;
			transformed_obs.x = particle_x + (cos(particle_theta) * observations[observation_idx].x) -
					(sin(particle_theta) * observations[observation_idx].y);
			transformed_obs.y = particle_y + (sin(particle_theta) * observations[observation_idx].x) +
					(cos(particle_theta) * observations[observation_idx].y);
			transformed_observations.push_back(transformed_obs);
	    }

		// Filter map landmarks
		vector<LandmarkObs> predicted_landmarks;
		for (unsigned int landmark_idx = 0; landmark_idx < map_landmarks.landmark_list.size(); landmark_idx++) {
			Map::single_landmark_s current_landmark = map_landmarks.landmark_list[landmark_idx];
			if ((fabs((particle_x - current_landmark.x_f)) <= sensor_range) &&
					(fabs((particle_y - current_landmark.y_f)) <= sensor_range)) {
				predicted_landmarks.push_back(LandmarkObs {current_landmark.id_i, current_landmark.x_f, current_landmark.y_f});
			}
		}

		// Associate data observation to predicted landmarks by calling dataAssociation
		dataAssociation(predicted_landmarks, transformed_observations);
		// Compute the new weight of each particle
		particles[particle_index].weight = 1.0;
		double sigma_x = std_landmark[0];
		double sigma_y = std_landmark[1];
		double sigma_x_square = pow(sigma_x, 2);
		double sigma_y_square = pow(sigma_y, 2);
		double normalizer = (1.0 / (2.0 * M_PI * sigma_x * sigma_y));

		// Calculate the weight of each particle based on the multivariate Gaussian probability function
		for (unsigned int observation_idx = 0; observation_idx < transformed_observations.size(); observation_idx++) {
			double trans_obs_x = transformed_observations[observation_idx].x;
			double trans_obs_y = transformed_observations[observation_idx].y;
			double trans_obs_id = transformed_observations[observation_idx].id;
			double multi_prob = 1.0;

			for (unsigned int landmark_idx = 0; landmark_idx < predicted_landmarks.size(); landmark_idx++) {
				double pred_landmark_x = predicted_landmarks[landmark_idx].x;
				double pred_landmark_y = predicted_landmarks[landmark_idx].y;
				double pred_landmark_id = predicted_landmarks[landmark_idx].id;
				if (trans_obs_id == pred_landmark_id) {
					multi_prob = normalizer * exp(-1.0 * ((pow((trans_obs_x - pred_landmark_x), 2)/
							(2.0 * sigma_x_square)) + (pow((trans_obs_y - pred_landmark_y), 2)/(2.0 * sigma_y_square))));
					particles[particle_index].weight *= multi_prob;
				}
			}
		}
		weight_normalizer += particles[particle_index].weight;
	}

	// Normalize the weights of all particles
	for (unsigned int particle_index = 0; particle_index < particles.size(); particle_index++) {
		particles[particle_index].weight /= weight_normalizer;
		weights[particle_index] = particles[particle_index].weight;
	}
}

void ParticleFilter::resample() {
	/*
	 * Re-sample particles with replacement with probability proportional to their weight.
	 */
	vector<Particle> resampled_particles;

	// Create the class to generate random numbers
	default_random_engine gen;

	// Generate random particle index
	uniform_int_distribution<int> particle_index(0, num_particles - 1);
	int current_index = particle_index(gen);

	double beta = 0.0;
	double max_weight_2 = 2.0 * *max_element(weights.begin(), weights.end());

	for (unsigned int i = 0; i < particles.size(); i++) {
		uniform_real_distribution<double> random_weight(0.0, max_weight_2);
		beta += random_weight(gen);
		while (beta > weights[current_index]) {
			beta -= weights[current_index];
			current_index = (current_index + 1) % num_particles;
		}
		resampled_particles.push_back(particles[current_index]);
	}
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
