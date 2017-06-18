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

using namespace std;

/**
 * init Initializes particle filter by initializing particles to Gaussian
 *   distribution around first position and all the weights to 1.
 * @param x Initial x position [m] (simulated estimate from GPS)
 * @param y Initial y position [m]
 * @param theta Initial orientation [rad]
 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles.
	// Initialize all particles to first position (based on estimates of
	//         x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Set the number of particles.
	num_particles = 100;
	weights.resize(num_particles);

	default_random_engine gen;

	// Create a normal (Gaussian) distributions
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_theta;

		// TODO: Sample  and from these normal distrubtions like this:
		// sample_x = dist_x(gen);
		// where "gen" is the random engine initialized earlier (line 18).
		sample_x = (double) dist_x(gen);
		sample_y = (double) dist_y(gen);
		sample_theta = (double) dist_theta(gen);

		Particle p = {
				i,                // int id;
				sample_x,         // double x;
				sample_y,         // double y;
				sample_theta,     // double theta;
				1.0/num_particles // double weight;
		};
		particles.push_back(p);
	}

	// cout << "init there are " << particles.size() << "particles" << endl;
	is_initialized = true;

}


/**
 * prediction Predicts the state for the next time step
 *   using the process model.
 * @param delta_t Time between time step t and t+1 in measurements [s]
 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 * @param velocity Velocity of car from t to t+1 [m/s]
 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
 */
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	// Create a normal (Gaussian) distributions
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);


	for (int i = 0; i < num_particles; ++i) {
		if (fabs(yaw_rate) > 0.001) {  // when the yaw rate is greater than zero
			double theta = particles[i].theta + yaw_rate * delta_t;
			particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(theta) - sin(particles[i].theta));
			particles[i].y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(theta));
			particles[i].theta = theta;
		} else {  // when the yaw rate is close to zero
			particles[i].x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
		}
		particles[i].x = particles[i].x + dist_x(gen);
		particles[i].y = particles[i].y + dist_y(gen);
		particles[i].theta = particles[i].theta + dist_theta(gen);

	}

}

/**
 * dataAssociation Finds which observations correspond to which landmarks (likely by using
 *   a nearest-neighbors data association).
 * @param predicted Vector of predicted landmark observations
 * @param observations Vector of landmark observations
 */
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.


	for (LandmarkObs& obs: observations) {
		double closest_prediction_distance = 1000000;
		double distance;

		for (LandmarkObs pre: predicted) {
			distance = dist(obs.x, obs.y, pre.x, pre.y);
			if (distance < closest_prediction_distance) {
				closest_prediction_distance = distance;
				obs.id = pre.id;
			}
		}
	}
}


/**
 * updateWeights Updates the weights for each particle based on the likelihood of the
 *   observed measurements.
 * @param sensor_range Range [m] of sensor
 * @param std_landmark[] Array of dimension 2 [standard deviation of range [m],
 *   standard deviation of bearing [rad]]
 * @param observations Vector of landmark observations
 * @param map Map class containing map landmarks
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
								   std::vector<LandmarkObs> observations, Map map_landmarks) {
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

	for (auto &p: particles) {

		int ct = 0;
		vector<LandmarkObs> predictions;

		for (auto &landmark: map_landmarks.landmark_list) {
			if (dist(landmark.x_f, landmark.y_f, p.x, p.y) <= sensor_range) {
				LandmarkObs landmark_obs = {
						ct,                     // int id;				// Id of matching landmark in the map.
						double(landmark.x_f),          // double x;			// Local (vehicle coordinates) x position of landmark observation [m]
						double(landmark.y_f)           // double y;			// Local (vehicle coordinates) y position of landmark observation [m]
				};
				ct = ct + 1;
				predictions.push_back(landmark_obs);
			}
		}

		vector<LandmarkObs> transformed_observations;
		for (auto &observation: observations) {
			double p_x = p.x;
			double p_y = p.y;
			double p_theta = p.theta;

			//  http://planning.cs.uiuc.edu/node99.html
			double transformed_x = observation.x * cos(p_theta) - observation.y * sin(p_theta) + p_x;
			double transformed_y = observation.x * sin(p_theta) + observation.y * cos(p_theta) + p_y;
			LandmarkObs transformed_observation = {observation.id, transformed_x, transformed_y};
			transformed_observations.push_back(transformed_observation);
		}

		dataAssociation(predictions, transformed_observations);

		double total_weight = 1.0;

		for (auto const& observation: transformed_observations) {

			LandmarkObs prediction = predictions[observation.id];
			double x_diff = (observation.x - prediction.x);
			double y_diff = (observation.y - prediction.y);


			double std_sq_x = std_landmark[0] * std_landmark[0];
			double std_sq_y = std_landmark[1] * std_landmark[1];

			double exp_term = -1 * ( (x_diff * x_diff) / (2 * std_sq_x) + (y_diff * y_diff) / (2* std_sq_y));

			// cout << "individual exp_term: " << exp_term << endl;

			double weight = exp (exp_term) / (2 * M_PI * std_landmark[0] * std_landmark[1]);

			total_weight = total_weight * weight;

		}

		p.weight = total_weight;
	//	weights[p.id] = weight;
	}

	for (int i = 0; i < num_particles; i++) {
		weights[i] = particles[i].weight;
	}

}

/**
 * resample Resamples from the updated set of particles to form
 *   the new set of particles.
 */
void ParticleFilter::resample() {
// TODO: Resample particles with replacement with probability proportional to their weight.
// NOTE: You may find std::discrete_distribution helpful here.
//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;
	discrete_distribution<> d(weights.begin(), weights.end());
	vector<Particle> resampled_particles;

	for (int i = 0; i < num_particles; ++i) {
		int resampled_particle_index = d(gen);
		resampled_particles.push_back(particles[resampled_particle_index]);
		// cout << "resampled_particle_index: " << resampled_particle_index << endl;
	}

	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

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