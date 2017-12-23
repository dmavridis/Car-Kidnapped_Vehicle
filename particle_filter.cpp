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
#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 50;
    particles = std::vector<Particle>(num_particles);
    weights = std::vector<double>(num_particles);

    for (int ix = 0; ix < num_particles; ix++){
        weights[ix] = 1;
        particles[ix].id = ix;
        particles[ix].x = x;
        particles[ix].y = y;
        particles[ix].theta = theta;
        particles[ix].weight = weights[ix];
    }
	is_initialized = true;
    std::cout << "Filter initialized..." << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    // Setup of normal distribution for each dimension
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> x_noise(0,std_pos[0]);
    std::normal_distribution<> y_noise(0,std_pos[1]);
    std::normal_distribution<> theta_noise(0,std_pos[2]);

    // Predicting x,y,theta for each particle according to the vehicle control data
    for (int ix = 0; ix < num_particles; ix++){
        if(yaw_rate != 0){
            particles[ix].x += velocity/yaw_rate*(sin(particles[ix].theta+yaw_rate*delta_t) - sin(particles[ix].theta));
            particles[ix].y -= velocity/yaw_rate*(cos(particles[ix].theta+yaw_rate*delta_t) - cos(particles[ix].theta));
            particles[ix].theta += yaw_rate*delta_t;
        }
        else{
            particles[ix].x += velocity*delta_t*cos(particles[ix].theta);
            particles[ix].y -= velocity*delta_t*sin(particles[ix].theta);
        }
        // Noise addition due to sensor
        particles[ix].x += x_noise(gen);
        particles[ix].y += y_noise(gen);
        particles[ix].theta += theta_noise(gen);

        // Check that theta is within [0,2*pi]
        if(particles[ix].theta < 0){
            particles[ix].theta += 2*M_PI;
        }
        else if(particles[ix].theta > 2*M_PI){
            particles[ix].theta -= 2*M_PI;
        }
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    // Not used //
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
    // TODO: Update the weights of each particle using a multi-variate Gaussian distribution. You can read
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


    // define the observations in car coordinates to map coordinates in relation to the particle
    std::vector<LandmarkObs> transformed_landmarks;
    transformed_landmarks = std::vector<LandmarkObs>(observations.size());

    for (int ix = 0; ix < num_particles; ix++){
        // transform observations car coordinates to map coordinates in relation to the particle
        double weight_ = 1;

        // Define the closest landmark
        double min_distance, min_distance_x, min_distance_y, distance, distance_x, distance_y;
        int land_id;
        min_distance = sensor_range;

        // Transform the observed data from car to map coordinates in relation to the particle
        for (int ix_obs = 0; ix_obs < observations.size(); ix_obs++){
            transformed_landmarks[ix_obs].id = ix_obs;
            transformed_landmarks[ix_obs].x = particles[ix].x + cos(particles[ix].theta)*observations[ix_obs].x \
                                                              - sin(particles[ix].theta)*observations[ix_obs].y;
            transformed_landmarks[ix_obs].y = particles[ix].y + sin(particles[ix].theta)*observations[ix_obs].x \
                                                              + cos(particles[ix].theta)*observations[ix_obs].y;

            // For the transformed observation find the closest landmark in map coordinates
            for (int ix_map = 0; ix_map < map_landmarks.landmark_list.size(); ix_map ++){
                distance_x = map_landmarks.landmark_list[ix_map].x_f - transformed_landmarks[ix_obs].x;
                distance_y = map_landmarks.landmark_list[ix_map].y_f - transformed_landmarks[ix_obs].y;
                distance = sqrt(distance_x*distance_x + distance_y*distance_y);
                if (distance <= min_distance){
                    min_distance = distance;
                    min_distance_x = distance_x;
                    min_distance_y = distance_y;
                    land_id = map_landmarks.landmark_list[ix_map].id_i;
                }
            }
            // Get the weight of the particle based on each observation and the closent landmark
            weight_ *= 1/(2*M_PI*std_landmark[0]*std_landmark[1])* \
                       pow(M_E,-(min_distance_x*min_distance_x/(2*std_landmark[0]*std_landmark[0])))* \
                       pow(M_E,-(min_distance_y*min_distance_y/(2*std_landmark[1]*std_landmark[1])));
        }

        particles[ix].weight = weight_;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    // Define the discrete random generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Resampled particle
    Particle particle_;

    // Probability of each particle to be resampled
    std::vector<double> probs;
    probs = std::vector<double>(num_particles);

    // Clone of the particles for the resampling
    std::vector<Particle> temp_particles;
    temp_particles = std::vector<Particle>(num_particles);

    // Initializing the variables
    for (int ix = 0; ix < num_particles; ix++){
        temp_particles[ix] = particles[ix];
        probs[ix] = temp_particles[ix].weight;
    }
    std::discrete_distribution<> d(probs.begin(), probs.end());

    // Resampling using the discrete distribution
    for (int ix = 0; ix < num_particles; ix++){
        particle_ = temp_particles[d(gen)];
        particles[ix] = particle_;
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
