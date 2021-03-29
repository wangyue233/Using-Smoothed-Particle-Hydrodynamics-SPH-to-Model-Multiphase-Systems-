#pragma once
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

double define_beach(double x);
bool determine_beach(double x, double y);
double gaussian(double x);
bool determine_gaussian(double x, double y);

class SPH_main;

class SPH_particle
{
public:
	double x[2], v[2];				//position and velocity for current time
	double x_next[2], v_next[2];    //position and velocity for next time step
	double x_prv[2], v_prv[2];      //position and velocity for previous time step
	bool bounadry;                  // state for boundary points or not
	double rho, P;					//density and pressure for current time
	double rho_nxt;			//density for next time step
	double rho_prv;         //density for previous time step
	double mass;

	static SPH_main* main_data;		//link to SPH_main class so that it can be used in calc_index

	int list_num[2];				//index in neighbour finding array

	void calc_index(void);

	void set_pressure(void);

	double kernel(double q, double h);

	double grad_kernel(double q, double h);

	void set_initial(void);

	void update(void);

	void update_beach(void);

	double boundary_acc(double dx, double xab);


	void update_rho(double& rho2, double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]);

	void update_vel(double* v_next, double dist, double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]);

	void update_pos(double dt);


};


class SPH_main
{
public:
	SPH_main();

	void set_values(double);
	void initialise_grid(void);

	void place_points(double* min, double* max);

	void place_beach(double* min, double* max);

	void allocate_to_grid(void);			//allocates all the points to the search grid (assumes that index has been appropriately updated)

	void neighbour_iterate(SPH_particle* part);

	void neighbour_iterate2(SPH_particle* part);

	void smooth_pressure(SPH_particle* part);


	double h;								//smoothing length
	double h_fac;
	double dx;								//particle initial spacing
	double dt;

	double min_x[2], max_x[2];				//dimensions of simulation region

	int max_list[2];

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell
};
