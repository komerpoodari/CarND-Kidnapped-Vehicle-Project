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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    num_particles = 60;  // settled after experimenting with a range [10 - 100] particles.
    //cout << "Enter number of particles:\n";
    //cin >> num_particles;
   
    std::normal_distribution<double> ND_x(x, std[0]); 
    std::normal_distribution<double> ND_y(y, std[1]); 
    std::normal_distribution<double> ND_theta(theta, std[2]); 
    std::default_random_engine gen;
    const double initial_weight = 1.0;
    for (int i = 0; i < num_particles; i++) {
        Particle particle;
        particle.id = i;
        
        // Initialize all particles to first position (based on estimates of 
	    // x, y, theta and their uncertainties from GPS) and all weights to 1. 
	    // Add random Gaussian noise to each particle.
        particle.x = ND_x(gen);
        particle.y = ND_y(gen);
        particle.theta = ND_theta(gen);
        particle.weight = initial_weight;
        
        particles.push_back(particle);
        weights.push_back(initial_weight);
    }
    debug_count = 1;
    
    /**
    *
    //cout << "Init x: " << x << "; y: " << y << "; theta: " << theta << endl;
    
    for (unsigned int i = 0; i < particles.size(); i++) {
        cout << "Particle: Id: " << particles[i].id << "\t" << particles[i].x << "\t" << particles[i].y << "\t" << particles[i].theta << endl;
    }
    //cout << "debug_count: " << debug_count << endl;
    *
    */
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    std::default_random_engine gen;
    
    debug_count++;
    //cout << "debug_count: " << debug_count << endl;
    
    for (int i = 0; i < num_particles; i++) {
        double new_x, new_y, new_theta;
        Particle particle = particles[i];
        
        if (fabs(yaw_rate) < 0.001) {  // review later changed from 0.00001
            new_x = particle.x + velocity*cos(particle.theta)*delta_t;
            new_y = particle.y + velocity*sin(particle.theta)*delta_t;
            new_theta = particle.theta;
        } else {  // Lesson 14, section 7
            new_x = particle.x + (velocity/yaw_rate)*(sin(particle.theta + yaw_rate*delta_t) - 
                                                        sin(particle.theta));
            new_y = particle.y + (velocity/yaw_rate)*(cos(particle.theta) - 
                                                        cos(particle.theta + yaw_rate * delta_t));
            new_theta = particle.theta + yaw_rate*delta_t;
        }

        // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
        std::normal_distribution<double> ND_x(new_x, std_pos[0]); 
        std::normal_distribution<double> ND_y(new_y, std_pos[1]); 
        std::normal_distribution<double> ND_theta(new_theta, std_pos[2]); 
        
        // TODO: Add measurements to each particle and add random Gaussian noise.
        particles[i].x = ND_x(gen);
        particles[i].y = ND_y(gen);
        particles[i].theta = ND_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    if (predicted.size() == 0 || observations.size() == 0) {
        //cout << "Landmarks: predicted: " << predicted.size() << " observed: " << observations.size() << endl;
        return;
    }

    for (unsigned int i = 0; i < observations.size(); i++) {
        LandmarkObs an_o = observations[i];
        
        //minimum distance initialized max numerical limit
        //double m_dist = std::numeric_limits<double>::max();
        double m_dist = 1.0e99;
        
        //cout << "initial min: " << init_dist << endl;
        int m_id = -1;
        double c_dist = 0.0;
        for (unsigned int j = 0; j < predicted.size(); j++) {
            LandmarkObs a_p = predicted[j];

            //compute the distance between landmark and prediction element
            c_dist = dist(an_o.x, an_o.y, a_p.x, a_p.y);
            if (c_dist > 1.0e99) {
                cout << "an_o.x: " << an_o.x << " an_o.y: " << an_o.y << " a_p.x: " <<  a_p.x << " a_p.y: " << a_p.y << endl;
            }
            if (c_dist < m_dist) {
                m_dist = c_dist;
                m_id = a_p.id; 
                //cout << "min_id: " << m_id << " min_dist: " << m_dist << endl;
            }
        }
        
        // check whether we got an association
        observations[i].id = m_id;
        if (m_id == -1)
            cout << "dataAss: m_dist: " << m_dist << " c_dist: " << c_dist << endl;  //should not happe
    }
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
    
    //Lesson 14, sections 15 & 18
    
    // copy standard deviation into local variables
    double std_x, std_y;
    std_x = std_landmark[0];
    std_y = std_landmark[1];
 
    // create a short list of landmarks that are in sensor range as per each particle location.
    
    //initiliaze weights sum
    double weight_sum = 0.0;
    for (int i = 0; i < num_particles; i++) {
        Particle pa = particles[i];
        
        //declare a subset of landmarks that are in sensor range
         vector<LandmarkObs> nearby_lm;
        
        for (unsigned int l = 0; l < map_landmarks.landmark_list.size(); l++) {
            
            double distance = dist(pa.x, pa.y, map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f);
            if (sensor_range >= distance) {
                LandmarkObs lm;
                lm.id = map_landmarks.landmark_list[l].id_i;
                lm.x = map_landmarks.landmark_list[l].x_f;  
                lm.y = map_landmarks.landmark_list[l].y_f;
                nearby_lm.push_back(lm);
            }
            //ignore landmarks far away from particle sensor range.
        } // for loop for shortlisting nearby landmarks
         
        /*        
        if (nearby_lm.size() == 0) {
            //cout << "Particle position x: " << pa.x << "; y: " << pa.y << endl;
        }
        */  
        
        // transform the observations to map coordinates
        vector<LandmarkObs> map_xfrm_obs;
        for (unsigned int j = 0; j < observations.size(); j++) {
            LandmarkObs xfm_obs;
            xfm_obs.id = observations[j].id;
            xfm_obs.x = pa.x + cos(pa.theta)*observations[j].x - sin(pa.theta)*observations[j].y;
            xfm_obs.y = pa.y + sin(pa.theta)*observations[j].x + cos(pa.theta)*observations[j].y;
            map_xfrm_obs.push_back(xfm_obs);
            
            // debug statement
            if (xfm_obs.x > 1.0E99 || xfm_obs.y > 1.0E99){
                cout << "pa.x: " << pa.x << "; pa.y: " << pa.y << "; pa.theta: " << pa.theta << endl;
                cout << "obs.x: " << observations[j].x << "; obs.y: " << observations[j].y << endl;
            }
        }
        
        //invoke association
        dataAssociation(nearby_lm, map_xfrm_obs);
        
        //compute the weights
        double pa_weight = 1.0;
        
        for (unsigned int j = 0; j < map_xfrm_obs.size(); j++) {
            double obs_x, obs_y, pred_x, pred_y, denom, numer, assoc_id;

            obs_x = map_xfrm_obs[j].x;
            obs_y = map_xfrm_obs[j].y;
            assoc_id = map_xfrm_obs[j].id;  //bug because of typo 'i' instead of 'j'
            
            //identify the corresponding entry in predictions
            pred_x = 0; pred_y = 0;
            //debug variable
            //bool found = false;
            for (unsigned int l = 0; l < nearby_lm.size(); l++) {
                if(assoc_id == nearby_lm[l].id) {
                    pred_x = nearby_lm[l].x;
                    pred_y = nearby_lm[l].y;
                    //found = true;
                    break;  // from loop review later
                }
            }
            double dx = obs_x - pred_x;
            double dy = obs_y - pred_y;
            
            double d2x = dx * dx;
            double d2y = dy * dy;  

            /*
            if (found == false)
                cout << "id match not found" << endl;
            */
            denom = 2.0 * M_PI * std_x * std_y;
            numer =  exp (-( d2x / (2.0 * std_x * std_x) + d2y / (2.0 * std_y * std_y)));
            pa_weight *= numer / denom;
            //cout << "pa_weight: " << pa_weight << " numerator: " << numer << " denominator: " << denom << endl;

        }
        particles[i].weight = pa_weight;
        weight_sum += pa_weight;
        //cout << "particles[" << i << "].weight: " << particles[i].weight << endl;
    } // for each particle
    
    //normalize weights
    weight_sum /= num_particles;
    for (int i=0; i< num_particles; i++) {
        particles[i].weight /= weight_sum; 
        //cout << "particles[" << i << "].weight: " << particles[i].weight << endl;
        weights[i] = particles[i].weight;
    }
    
    
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    default_random_engine gen;
    discrete_distribution <int> dist(weights.begin(), weights.end());
    vector<Particle> resampled_particles;
    
    for (int i = 0; i < num_particles; i++) {
        Particle particle = particles[dist(gen)];
        resampled_particles.push_back(particle);
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
