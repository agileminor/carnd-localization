/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *      Revised and completed version Author: Ian Gilmore
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
using namespace std;

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	std::random_device rd;
    	std::mt19937 gen(rd());
	//default_random_engine gen;
	is_initialized = true;
	num_particles = 10;
	//vector<Particle> particles;
	normal_distribution<double> N_x_part(x, std[0]);
	normal_distribution<double> N_y_part(y, std[1]);
	normal_distribution<double> N_theta_part(theta, std[2]);

	for(int i=0;i<num_particles;i++) {
		double n_x = N_x_part(gen);
		double n_y = N_y_part(gen);
		double n_theta = N_theta_part(gen);
		Particle new_particle;
		new_particle.x = n_x;
		new_particle.y = n_y;
		new_particle.theta = n_theta;
		new_particle.id = i;
		new_particle.weight = 1.0;
		particles.push_back(new_particle);

	}
	weights.clear();
  	weights.resize(num_particles, 1.0);
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	//default_random_engine gen;
	std::random_device rd;
    	std::mt19937 gen(rd());
	normal_distribution<double> N_x_part(0, std_pos[0]);
	normal_distribution<double> N_y_part(0, std_pos[1]);
	normal_distribution<double> N_theta_part(0, std_pos[2]);
	//cout << "MOVE" << endl;
	for(int i=0;i<num_particles;i++) {
		double n_x = N_x_part(gen);
		double n_y = N_y_part(gen);
		double n_theta = N_theta_part(gen);
		double current_x = particles[i].x;
		double current_y = particles[i].y;
		double current_theta = particles[i].theta;
		if (fabs(yaw_rate) > 0.000001) {
			particles[i].x = current_x + velocity / yaw_rate * (sin(current_theta + yaw_rate * delta_t) - sin(current_theta));
			particles[i].y = current_y + velocity / yaw_rate * (cos(current_theta) - cos(current_theta + yaw_rate * delta_t));
			particles[i].theta = current_theta + yaw_rate * delta_t;
		} else {
			particles[i].x = current_x + velocity * delta_t * cos(current_theta);
			particles[i].y = current_y + velocity * delta_t * sin(current_theta);
		}
		particles[i].x += n_x;
		particles[i].y += n_y;
		particles[i].theta += n_theta;
		//cout << "\tvelocity= " << velocity << " yaw_rate= " << yaw_rate << endl;
		//cout << "\tnew particle x= " << particles[i].x << " y= " << particles[i].y << " yaw= " << particles[i].theta << endl;

	}
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

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
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	for(int i=0;i<num_particles;i++) {
		Particle cur_part = particles[i];
		std::vector<LandmarkObs> obs_temp = observations;
		std::vector<LandmarkObs> cur_lmks;
		//cout << "UPDATE" << endl;
	        //cout << "\tparticle x = " << cur_part.x << " y= " << cur_part.y << " theta= " << cur_part.theta << " id= " << cur_part.id << endl;	
		// filter out landmarks out of sensor range, to minimize looping
		for(int k=0;k<map_landmarks.landmark_list.size();k++) {
			if (dist(cur_part.x, cur_part.y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f) <= sensor_range) {
				LandmarkObs mark;
				mark.x =map_landmarks.landmark_list[k].x_f; 
				mark.y = map_landmarks.landmark_list[k].y_f;
				mark.id = map_landmarks.landmark_list[k].id_i;
				cur_lmks.push_back(mark);
			}	
		}
		double total_weight = 1.0;
		for(int j=0;j < obs_temp.size();j++) {
			double min_dist = 1e9;
			double cur_weight = 0.0;
			double min_x;
			double min_y;
			double min_id;
			// remap observation x/y using particle x/y/yaw
			//cout << "\tobs x = " << obs_temp[j].x << " y = " << obs_temp[j].y << endl;
			double cur_x = obs_temp[j].x;
			double cur_y = obs_temp[j].y;
			obs_temp[j].x = cur_part.x + cur_x * cos(cur_part.theta) - cur_y * sin(cur_part.theta);
			obs_temp[j].y = cur_part.y + cur_x * sin(cur_part.theta) + cur_y * cos(cur_part.theta);
			//cout << "\ttransformed obs x = " << obs_temp[j].x << " y = " << obs_temp[j].y << endl;
			// compare current remapped observation with landmarks. Pick nearest one to calculate weight
			for(int m=0;m < cur_lmks.size(); m++) {
				double cur_dist = dist(obs_temp[j].x, obs_temp[j].y, cur_lmks[m].x, cur_lmks[m].y);
				if (cur_dist < min_dist) {
					min_dist = cur_dist;
					min_x = cur_lmks[m].x;
					min_y = cur_lmks[m].y;
					min_id = cur_lmks[m].id;
				}
			}
			//cout << "\tmatched landmark x=" << min_x << " y= " << min_y << " id= " << min_id << endl;
			double term1 = pow((obs_temp[j].x - min_x), 2)/ (2 * pow(std_landmark[0], 2));
			double term2 = pow((obs_temp[j].y - min_y), 2)/ (2 * pow(std_landmark[1], 2));
			double term3 = 2 * M_PI * std_landmark[0] * std_landmark[1];
			//cout << "\t term1=" << term1 << " term2= " << term2 << " term3= " << term3 << endl;
			cur_weight = 1 / term3 * exp(-(term1 + term2));
			if(cur_weight < 1e-6) {
				cur_weight = 1e-6;
			}
			//cout << "\tadded weight= " << cur_weight << endl;
			total_weight *= cur_weight;
			//cout << "\ttotal weight= " << total_weight << endl;

		}
		particles[i].weight = total_weight;
		//cout << "\ttotal weight= " << total_weight << endl;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::random_device rd;
    	std::mt19937 gen(rd());
    	std::vector<Particle> particles_new;
	/*
	double weights_sum = 0;
    	for(int i=0; i<num_particles; i++) {
        	weights_sum += particles[i].weight;
    	}
	*/
	//cout << "weights sum= " << weights_sum << endl;
	//cout << "RESAMPLE" << endl;
    	for(int i=0; i<num_particles; i++) {
        	weights[i] = particles[i].weight;
		//cout << "weight= " << particles[i].weight << endl;
    	}
    	std::discrete_distribution<> d(weights.begin(), weights.end());
    	for(int n=0; n<num_particles; ++n) {
        	particles_new.push_back(particles[d(gen)]);
    	}
    	particles = particles_new;
	for(int i=0;i<num_particles;i++) {
		//cout << "x= " << particles[i].x << " y= " << particles[i].y << " theta= " << particles[i].theta << " id= " << particles[i].id << endl;
		particles[i].id = i; // renumber to keep track of particles in next round
	}
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
