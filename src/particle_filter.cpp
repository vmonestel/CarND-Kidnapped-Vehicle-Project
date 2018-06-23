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
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
