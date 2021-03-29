#pragma once
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

class SPH_main;

// All function descriptions stated in the .cpp file of the solver

class SPH_particle {

public:
	// ----- Particle properties -----//

	double x[2], v[2];						// position and velocity for current time
	double x_next[2], v_next[2];			// position and velocity for next time step
	double x_next_bar[2], v_next_bar[2];    // position and velocity for next time step
	bool bounadry;							// state for boundary points or not
	double rho, P, P_nxt;					// density and pressure for current time
	double rho_nxt;							// density and pressure for next time step
	double rho_nxt_bar;
	double mass;

	static SPH_main* main_data;				// link to SPH_main class so that it can be used in calc_index

	int list_num[2];						// index in neighbour finding array


	// ----- Class Functions -----//

	void calc_index(void);

	void set_pressure(void);

	void set_pressure_nxt(void);

	double kernel(double q, double h);

	double grad_kernel(double q, double h);

	int dist_to_index(double r, double dr);

	void precompute_grad_W(double*& pre_grad_W, double dr);

	void set_initial(void);

	void update();

	void update_rho(double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]);

	void update_rho_bar(double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]);

	void update_vel(double dist, double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]);

	void update_vel_bar(double dist, double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]);

	void update_pos(double dt);

	void update_pos_bar(double dt);

	void update_all_bar();

	void zeroing(void);

	void zeroing_bar(void);


};

class SPH_main {

public:

	// ----- Particle properties -----//

	double h;								//smoothing length
	double h_fac;
	double dx;								//particle initial spacing

	double min_x[2], max_x[2];				//dimensions of simulation region

	int max_list[2];

	vector<SPH_particle> particle_list;		//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell


	// ----- Class Functions -----//

	SPH_main();

	void set_values(void);

	void initialise_grid(void);

	void place_points(double* min, double* max, int& nwall);

	void allocate_to_grid(void);

	void neighbour_iterate(SPH_particle* part);

	void neighbour_iterate2(SPH_particle* part);

	void smooth_pressure(SPH_particle* part);

};

