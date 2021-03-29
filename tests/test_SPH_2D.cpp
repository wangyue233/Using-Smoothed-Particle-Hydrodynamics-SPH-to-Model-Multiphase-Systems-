#include "SPH_2D.h"
#include "file_writer.h"
#include <string>

using namespace std;

SPH_main domain;

int main(void)
{
	domain.set_values(0.2);										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid
	int iterations = 0;
	domain.place_points(domain.min_x, domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	string fname = "tests/test_output/test_full_timstep_" + to_string(iterations) + ".vtp";
	write_file(fname.c_str(), &(domain.particle_list));
	cout << "TEST for SPH_2D starts: 1000 iterations, rectangular, Full Timesteps. VTP file generated each 200 steps in tests/test_output/" << endl;
	clock_t start = clock();
	while (iterations < 1000) {
		domain.allocate_to_grid();									//needs to be called for each time step
		for (int i = 0; i < domain.particle_list.size(); i++)    //vector::size() return the size instead of number of elements
		{
			if (!domain.particle_list[i].bounadry)
			{
				domain.neighbour_iterate(&domain.particle_list[i]);
			}
		}
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
			for (int i = 0; i < domain.particle_list.size(); i++)
			{

				domain.smooth_pressure(&domain.particle_list[i]);

			}
			for (int i = 0; i < domain.particle_list.size(); i++)
			{

				domain.particle_list[i].rho = domain.particle_list[i].rho_nxt;
				domain.particle_list[i].set_pressure();

			}
		}
		if (iterations % 200 == 0)
		{
			fname = "tests/test_output/test_full_timstep_" + to_string(iterations) + ".vtp";
			write_file(fname.c_str(), &(domain.particle_list));
		}
		iterations++;
	}
	cout << "Time Elapsed: " << (double) (clock() - start) / CLOCKS_PER_SEC << endl;
	return 0;
}
