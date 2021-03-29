#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include "omp.h"
#define CHUNKSIZE 100
#define THREAD_NUM 16

using namespace std;

SPH_main domain;

int main(void)
{
	double dx;
	//for users to input value
	cout << "dx: ";
	cin >> dx;
	domain.set_values(dx);										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid
	int iterations = 0;
	cout << "Iterations: ";
	int max_iter = 10000;
	cin >> max_iter;
	domain.place_points(domain.min_x, domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	string fname = "output/full_timstep_" + to_string(iterations) + ".vtp";
	write_file(fname.c_str(), &(domain.particle_list));
	cout << "Maximum " << max_iter << " iterations, Full Timestep." << endl;
	double start_time = omp_get_wtime();
	while (iterations < max_iter) {
		domain.allocate_to_grid();									//needs to be called for each time step
		#pragma omp parallel for num_threads(THREAD_NUM) schedule(guided, CHUNKSIZE)	
		for (int i = 0; i < domain.particle_list.size(); i++)    //vector::size() return the size instead of number of elements
		{
			if (!domain.particle_list[i].bounadry)
			{
				domain.neighbour_iterate(&domain.particle_list[i]);
			}
		}
		#pragma omp parallel for num_threads(THREAD_NUM) schedule(guided, CHUNKSIZE)
		for (int i = 0; i < domain.particle_list.size(); i++)
		{
			if (!domain.particle_list[i].bounadry)
			{
				domain.particle_list[i].update();
				domain.particle_list[i].calc_index();
			}
		}
		if (iterations % 10 == 0)
		{
			#pragma omp parallel for num_threads(THREAD_NUM) schedule(guided, CHUNKSIZE)
			for (int i = 0; i < domain.particle_list.size(); i++)
			{
				domain.smooth_pressure(&domain.particle_list[i]);
			}
			#pragma omp parallel for num_threads(THREAD_NUM) schedule(guided, CHUNKSIZE)
			for (int i = 0; i < domain.particle_list.size(); i++)
			{
				domain.particle_list[i].rho = domain.particle_list[i].rho_nxt;
				domain.particle_list[i].set_pressure();
			}
		}
		if (iterations % 200 == 0)
		{
			fname = "output/full_timstep_" + to_string(iterations) + ".vtp";
			write_file(fname.c_str(), &(domain.particle_list));
		}
		iterations++;
	}
	double end_time = omp_get_wtime(); // end time
	std::cout << "Time elapsed: " << (end_time - start_time) << " sec " << endl; //total time
	return 0;
}
