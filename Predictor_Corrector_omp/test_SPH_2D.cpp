#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <omp.h>
#include <cstdio>
#define CHUNKSIZE 100
#define THREAD_NUM 16

using namespace std;

SPH_main domain;

int main(int argc, char* argv[]) {

	// SPH_particle::main_data->dx = 0;

	// Initialisations
	int nwall = 0;						// Initialise number of boundary particles
	domain.set_values();		// Set problem parameters (dt, h_fac)
	domain.initialise_grid();			// Initialise simulation grid
	int iterations = 0;					// Initialise iteration count
	int print_every = 30;				// Write to file every n iterations
	int smoothening_every = 10;			// Apply smoothening every n iterations
	int max_iter = 5000;

	domain.place_points(domain.min_x, domain.max_x, nwall);

	std::cout << "Number of boundary particles: " << nwall << endl;  // Print the number of boundary particles

	// Initiate write-to-file
	string fname = "tests/test_initial_file_writer00" + to_string(iterations) + ".vtp";
	write_file(fname.c_str(), &(domain.particle_list));  // Write initial state

	double start_time = omp_get_wtime(); // Start time
	int domain_size = domain.particle_list.size();  // Creates a variable of the domain size of int type
	int k;  // Iterations counter

	// Main iterations loop
	while (iterations < max_iter) {

		domain.allocate_to_grid();  // Allocates particles to grid at the beginning of each iteration


		// ----------------------------------------------------------------------------
		// Note that all the following calculations, expect pressure are only applied on 
		// the fluid particles. Since boundary particles are stored first in the list, 
		// simply looping from nwall + 1 assures that only fluid particles are computed.
		// ----------------------------------------------------------------------------


#pragma omp parallel for num_threads(THREAD_NUM) private(k) schedule(guided, CHUNKSIZE)
		for (k = nwall + 1; k < domain_size; k++) {

			domain.particle_list[k].zeroing();  // Zeroing of the all the auxiliery summation terms used in Predictor-Corrector

		}

#pragma omp parallel for num_threads(THREAD_NUM) private(k) schedule(guided, CHUNKSIZE)
		for (k = nwall + 1; k < domain_size; k++) {

			domain.neighbour_iterate(&domain.particle_list[k]);  // Half-step

		}

#pragma omp parallel for num_threads(THREAD_NUM) private(k) schedule(guided, CHUNKSIZE)
		for (k = nwall + 1; k < domain_size; k++) {

			domain.particle_list[k].zeroing_bar();  // Intermediate zeroing of auxiliery summation terms

		}

#pragma omp parallel for num_threads(THREAD_NUM) private(k) schedule(guided, CHUNKSIZE)
		for (k = nwall + 1; k < domain_size; k++) {

			domain.neighbour_iterate2(&domain.particle_list[k]);  // Full-step

		}

#pragma omp parallel for num_threads(THREAD_NUM) private(k) schedule(guided, CHUNKSIZE)
		for (k = nwall + 1; k < domain.particle_list.size(); k++) {

			domain.particle_list[k].update_all_bar();  // Calculates the final particle position and velocity

		}

#pragma omp parallel for num_threads(THREAD_NUM) private(k) schedule(guided, CHUNKSIZE)
		for (k = nwall + 1; k < domain_size; k++) {

			domain.particle_list[k].update();  // Updates the newest calculated particle properties
			// domain.particle_list[k].calc_index();  // Re-destributes particles in the grid

		}

#pragma omp parallel for num_threads(THREAD_NUM) private(k) schedule(guided, CHUNKSIZE)
		for (k = nwall + 1; k < domain_size; k++) {

			domain.particle_list[k].calc_index();

		}


		// Pressure smoothening technique
		if (iterations % smoothening_every == 0) {

#pragma omp parallel for num_threads(THREAD_NUM) private(k) schedule(guided, CHUNKSIZE)
			for (k = 0; k < domain_size; k++) {

				domain.smooth_pressure(&domain.particle_list[k]);

			}

#pragma omp parallel for num_threads(THREAD_NUM) private(k) schedule(guided, CHUNKSIZE)
			for (k = 0; k < domain_size; k++) {

				domain.particle_list[k].set_pressure();
				domain.particle_list[k].rho = domain.particle_list[k].rho_nxt;

			}

		}

		if (iterations % print_every == 0 && iterations > 0) {

			// std::cout << "Iteration : " << iterations << endl;  // Print iteration number

			// Write to file
			fname = "tests/test_initial_file_writer00" + to_string(iterations) + ".vtp";
			write_file(fname.c_str(), &(domain.particle_list));
		}

		iterations++;

	}

	double end_time = omp_get_wtime(); // End time
	std::cout << "Time elapsed: " << (end_time - start_time) << " sec " << endl; // Total time

	return 0;
}
