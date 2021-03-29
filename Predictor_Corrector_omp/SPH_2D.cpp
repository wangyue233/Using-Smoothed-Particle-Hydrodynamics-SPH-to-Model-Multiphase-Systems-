#include "SPH_2D.h"
#include <math.h>  
#include <omp.h>
#include <cstdio>
#define CHUNKSIZE 100
#define THREAD_NUM 8

double M_PI = 3.14159265359;

SPH_main* SPH_particle::main_data;

SPH_main::SPH_main() {
	SPH_particle::main_data = this;
}

void SPH_particle::calc_index(void) {
	/*
	 * Description:
	 * Calculates the new particle index based on the updated position, x.
	 *
	 */

	for (int i = 0; i < 2; i++)
		list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0 * main_data->h));
}

void SPH_particle::set_pressure(void) {
	/* 
	 * Description:
	 * Calculates pressure based on rho (current density).
	 *
	 */

	double c0(20.0), rho0(1000.0), gamma(7.0);
	P = rho0 * c0 * c0 / gamma * (pow(rho / rho0, gamma) - 1);
}

double SPH_particle::kernel(double q, double h) {
	/*
	 * Description:
	 * Applies Kernel function to calculate weights for particle interaction.
	 *
	 * Inputs:
	 * q: distance of interaction
	 * h: smoothening length (proglem-dependent constant)
	 *
	 */

	if (q <= 1) {
		return 10.0 / (7.0 * M_PI * h * h) * (1.0 - 3.0 / 2.0 * q * q + 3.0 / 4.0 * q * q * q);
	}
	else if (q <= 2) {
		return 10.0 / (7.0 * M_PI * h * h) * (1.0 / 4.0 * pow(2 - q, 3));
	}
	else {
		return 0;
	}

}

double SPH_particle::grad_kernel(double q, double h) {
	/*
	 * Description:
	 * Applies the derivative of the Kernel function to calculate weights for particle interaction.
	 *
	 * Inputs:
	 * q: distance of interaction
	 * h: smoothening length (proglem-dependent constant)
	 *
	 */

	if (q <= 1) {
		return 10.0 / (7.0 * M_PI * h * h * h) * (-3 * q + 9.0 / 4 * q * q);
	}
	else if (q <= 2) {
		return 10.0 / (7.0 * M_PI * h * h * h) * (-3.0 / 4.0 * pow(2 - q, 2));
	}
	else {
		return 0;
	}
}

int SPH_particle::dist_to_index(double r, double dr) {
	/*
	 * This function converts the given radius / distance into the index 
	 * in which to access in the precomputed W and dW array.
	 */

	double k = r / dr;
	int index_prec = (int)k;
	return index_prec;
}

void SPH_particle::precompute_grad_W(double*& pre_grad_W, double dr) {
	/*
	 * This function precomputes the values for the smooth kernel gradient
	 *----------
	 *	pre_W : double pointer
	 *			Initialised pointer to array to store precomputed values
	 *	dr : double
	 *			radius spacing (precision - evaluated every (dr) distance )
	 *
	 * Returns
	 * -------
	 *	grad_W : double
	 *		Smooth kernel Weight with respect to dr.
	 */

	//create array storing W values
	int dim;
	dim = 10; // ( (2*main_data->h) / dr ) + 1 ;
	pre_grad_W[dim + 1];
	//initialise smooth kernel constants
	const double PI = 3.141592653589793238463;
	//Wight common factor
	const double fact = (10.0 / (7 * PI * pow(main_data->h, 2.0)));
	double r = 0.0;
	//STROES ARRAY of size = dim+2 ; dim +1 to access the last one that is always 0 !!! To access will be dim+1 being dim = [2h/dr] + 1;
	for (int i = 0; i <= dim + 1; i++)
	{
		if (r <= main_data->h)
		{
			pre_grad_W[i] = fact * (-(3.0 / pow(main_data->h, 2.0)) * r + (9.0 / (4.0 * pow(main_data->h, 3.0))) * pow(r, 2.0));
		}
		if (r > main_data->h&& r <= 2 * main_data->h)
		{
			pre_grad_W[i] = fact * 0.25 * (-(3.0 / pow(main_data->h, 3.0)) * pow(r, 2.0) + (12.0 / pow(main_data->h, 2.0)) * r - 12 / main_data->h);
		}
		if (r > 2)
		{
			pre_grad_W[i] = 0.0;
		}
		r += dr;
	}
	//return pre_grad_W;
}


void SPH_particle::set_pressure_nxt(void) {
	/* Predictor-Corrector exclusive function
	 *
	 * Description:
	 * Calculates pressure based on rho_nxt (half-step density).
	 *
	 */

	double c0(20.0), rho0(1000.0), gamma(7.0);
	P_nxt = rho0 * c0 * c0 / gamma * (pow(rho_nxt / rho0, gamma) - 1);

}

void SPH_particle::set_initial(void) {
	/* 
	 * Description:
	 * Initialises particle's v, v_next, v_next_bar, 
	 * rho, rho_nxt_, rho_nxt_bar and mass.
	 *
	 */

	v[0] = v[1] = 0.0;
	v_next[0] = v_next[1] = 0.0;
	v_next_bar[0] = v_next_bar[1] = 0.0;
	rho = 1000.0;
	rho_nxt = 1000.0;
	rho_nxt_bar = 1000.0;
	mass = main_data->dx * main_data->dx * rho;
	set_pressure();
}

void SPH_main::set_values(void) {
	/*
	 * Description:
	 * Sets initial dimensions and problem constants (dx, h).
	 *
	 */

	min_x[0] = 0.0;  // Set min values for grid --> min_x[0] = x minimum value ; 
	min_x[1] = 0.0;  // Set min values for grid --> min_x[1] = y minimum value ;

	max_x[0] = 20.0;
	max_x[1] = 10.0;

	// Set time-step and h factor
	dx = 0.2;
	h_fac = 1.3;

	h = dx * h_fac;    // obtain the smoothing lenght  

	cout << "Time step (dt): " << 0.1 * h / 20.0 << endl;  // Print time-step
}

void SPH_main::initialise_grid(void) {
	/*
	 * Description:
	 * Initialises grid of particles to be used in computations.
	 *
	 */

	for (int i = 0; i < 2; i++) {
		min_x[i] -= 2.0 * h;
		max_x[i] += 2.0 * h;  //add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1.0);  // number of elements in i-direction of the grid
	}

	search_grid.resize(max_list[0]);
	for (int i = 0; i < max_list[0]; i++)
		search_grid[i].resize(max_list[1]);
}

void SPH_main::place_points(double* min, double* max, int& nwall) {
	/*
	 * Description:
	 * Initialises boundaries and fluid particles.
	 *
	 * Inputs:
	 * min: Minimum dimension
	 * max: Maximum dimension
	 * nwall: Number of boundary (wall) particles
	 *
	 */

	SPH_particle particle;

	// Bottom boundary
	int num_x = int((max[0] - min[0]) / dx);
	int num_y = int((0 - min[1]) / dx);
	double length = min[0] + num_x * dx;
	for (int j = 0; j <= num_y; j++) {
		particle.x[1] = 0 - j * dx;
		particle.x_next[1] = 0 - j * dx;
		for (int i = 0; i <= num_x; i++) {
			particle.x[0] = min[0] + i * dx;
			particle.x_next[0] = min[0] + i * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
			nwall++;
		}
	}

	// Top boundary
	for (int j = 0; j <= num_y; j++) {
		particle.x[1] = max[1] - j * dx;
		particle.x_next[1] = max[1] - j * dx;
		for (int i = 0; i <= num_x; i++) {
			particle.x[0] = min[0] + i * dx;
			particle.x_next[0] = min[0] + i * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
			nwall++;
		}
	}

	// Left boundary
	num_x = int((2 * h) / dx);
	num_y = int((max[1] - 2 * h) / dx);
	for (int i = 0; i <= num_x; i++) {
		particle.x[0] = min[0] + i * dx;
		particle.x_next[0] = min[0] + i * dx;
		for (int j = 0; j <= num_y; j++) {
			particle.x[1] = j * dx;
			particle.x_next[1] = j * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
			nwall++;
		}
	}
	// Right boundary
	for (int i = 0; i <= num_x; i++) {
		particle.x[0] = length - i * dx;
		particle.x_next[0] = length - i * dx;
		for (int j = 0; j <= num_y; j++) {
			particle.x[1] = j * dx;
			particle.x_next[1] = j * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
			nwall++;
		}
	}

	// Initial set for mobile fluid particles

	// Left fluid part (column)
	num_x = int(3 / dx);
	num_y = int(5 / dx);
	for (int i = 0; i <= num_x; i++) {
		particle.x[0] = i * dx;
		particle.x_next[0] = i * dx;
		for (int j = 1; j <= num_y; j++) {
			particle.x[1] = j * dx;
			particle.x_next[1] = j * dx;
			particle.bounadry = false;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
		}

	}

	// Right fluid part
	num_x = int((max[0] - 2 * h - 3) / dx);
	num_y = int(2 / dx);
	for (int i = 1; i <= num_x; i++) {
		particle.x[0] = 3 + i * dx;
		particle.x_next[0] = 3 + i * dx;
		for (int j = 1; j < num_y; j++) {
			particle.x[1] = j * dx;
			particle.x_next[1] = j * dx;
			particle.bounadry = false;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
		}
	}
}

void SPH_particle::zeroing(void) {
	/* Predictor-Corrector exclusive function
	 *
	 * Description:
	 * Zeroes all summation terms required for 
	 * predictor corrector particle's position,
	 * velocities and densities.
	 *
	 */

	x_next[0] = x_next[1] = 0;
	v_next[0] = v_next[1] = 0;
	x_next_bar[0] = x_next_bar[0] = 0;
	v_next_bar[0] = v_next_bar[1] = 0;
	rho_nxt = 0.0;
	rho_nxt_bar = 0.0;
}

void SPH_particle::zeroing_bar(void) {
	/* Predictor-Corrector exclusive function
	 *
	 * Description:
	 * Zeroes all auxiliery (bar) summation terms required for
	 * predictor corrector particle's position,
	 * velocities and densities.
	 *
	 */

	x_next_bar[0] = x_next_bar[0] = 0;
	v_next_bar[0] = v_next_bar[1] = 0;
	rho_nxt_bar = 0.0;
}

void SPH_particle::update() {
	/*
	 * Description:
	 * Updates particle x, v, rho and P based on the
	 * corresponding newest calculated *_next values.
	 * Also applies boundary conditions on the particle.
	 *
	 */

	rho = rho_nxt;
	P = P_nxt;
	double dt = 0.1 * main_data->h / 20.0;

	// Imposed boundary condition
	//if (x_next[1] < 0.1) {
	//	v_next[1] += dt * 20.0;
	//}

	for (int i = 0; i < 2; i++) {
		v[i] = v_next[i];
		x[i] = x_next[i];
	}

	// Imposed mirror boundary conditions
	if (x[0] < 0) {
		x[0] = -x[0] / 2.0;
		v[0] = -v[0];
		// v[1] = 0.0;
		v[0] = abs(v[0]) / 2.0;
	}
	if (x[1] < 0) {
		x[1] = 0.0;
		v[0] = 0.0;
		// v[1] = abs(v[1]) / 2.0;
		v[1] = abs(v[1]);
	}
	if (x[0] > 20.0) {
		x[0] = 20.0 - x[0] + 20.0;
		// v[0] = -abs(v[0]) / 2.0;
		v[0] = -abs(v[0]);
		v[1] = 0.0;
	}
	if (x[1] > 10.0) {
		x[1] = 10.0 - x[1] + 10.0;
		v[0] = 0.0;
		v[1] = -abs(v[1]) / 2.0;
	}

	//if (x[0] >  0.5*(main_data->max_x[0] - main_data->min_x[0]) && x[1] < 0.05 * x[0] ) {
	//	x[0] = 5;
	//	x[1] = 5;
	//	v[0] = 0.0;
	//	v[1] = 0.0;
	//}
}

void SPH_particle::update_rho(double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]) {
	/*
	 * Description:
	 * Calulates summation of density-bar values.
	 *
	 * Inputs
	 * dt: time-step
	 * other_part: neighbour particle
	 * smooth: smoothing length
	 * eij: relative distance between current and neighbour particle
	 * vij: relative velocity between current and neighbour particle
	 *
	 */

	rho_nxt += dt * other_part->mass * smooth * (eij[0] * vij[0] + eij[1] * vij[1]);
}

void SPH_particle::update_rho_bar(double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]) {
	/* Predictor-Corrector exclusive function
	 *
	 * Description:
	 * Calulates summation of density-bar values for the intermediate step.
	 *
	 * Inputs
	 * dt: time-step
	 * other_part: neighbour particle
	 * smooth: smoothing length
	 * eij: (difference of x values)/(distance) term
	 * vij: relative velocity between current and neighbour particle
	 *
	 */

	rho_nxt_bar += dt * other_part->mass * smooth * (eij[0] * vij[0] + eij[1] * vij[1]);

}

void SPH_particle::update_vel(double dist, double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]) {
	/*
	 * Description:
	 * Calulates velocity summation term.
	 *
	 * Inputs
	 * dist: relative distance between current and neighbour particle
	 * dt: time-step
	 * other_part: neighbour particle
	 * smooth: smoothing length
	 * eij: (difference of x values)/(distance) term
	 * vij: relative velocity between current and neighbour particle
	 *
	 */

	v_next[0] -= dt * (other_part->mass * (other_part->P / other_part->rho / other_part->rho + P / rho / rho) * smooth * eij[0] \
		+ 0.001 * other_part->mass * (1 / other_part->rho / other_part->rho + 1 / rho / rho) * smooth * vij[0] / dist);

	v_next[1] -= dt * (other_part->mass * (other_part->P / other_part->rho / other_part->rho + P / rho / rho) * smooth * eij[1] \
		+ 0.001 * other_part->mass * (1 / other_part->rho / other_part->rho + 1 / rho / rho) * smooth * vij[1] / dist);

}

void SPH_particle::update_vel_bar(double dist, double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2]) {
	/* Predictor-Corrector exclusive function
	 *	
	 * Description:
	 * Calulates velocity-bar summation term for the intermediate step.
	 *
	 * Inputs:
	 * dist: relative distance between current and neighbour particle
	 * dt: time-step
	 * other_part: neighbour particle
	 * smooth: smoothing length
	 * eij: (difference of x values)/(distance) term
	 * vij: relative velocity between current and neighbour particle
	 *
	 */

	v_next_bar[0] -= dt * (other_part->mass * (other_part->P_nxt / other_part->rho_nxt / other_part->rho_nxt + P_nxt / rho_nxt / rho_nxt) * smooth * eij[0] \
		+ 0.001 * other_part->mass * (1 / other_part->rho_nxt / other_part->rho_nxt + 1 / rho_nxt / rho_nxt) * smooth * vij[0] / dist);

	v_next_bar[1] -= dt * (other_part->mass * (other_part->P_nxt / other_part->rho_nxt / other_part->rho_nxt + P_nxt / rho_nxt / rho_nxt) * smooth * eij[1] \
		+ 0.001 * other_part->mass * (1 / other_part->rho_nxt / other_part->rho_nxt + 1 / rho_nxt / rho_nxt) * smooth * vij[1] / dist);

}

void  SPH_particle::update_pos(double dt) {
	/*
	 * Description:
	 * Calulates position x_next.
	 * For the Predictor-Corrector, it is called in the half-step.
	 *
	 * Inputs:
	 * dt: time-step
	 *
	 */

	x_next[0] = x[0] + dt * v[0];
	x_next[1] = x[1] + dt * v[1];

}

void  SPH_particle::update_pos_bar(double dt) {
	/* Predictor-Corrector exclusive function
	 *
	 * Description:
	 * Calulates position x_next_bar, here used for the full-step.
	 *
	 * Inputs:
	 * dt: time-step
	 *
	 */

	x_next_bar[0] = x[0] + dt * v_next[0];
	x_next_bar[1] = x[1] + dt * v_next[1];
}

void  SPH_particle::update_all_bar() {
	/* Predictor-Corrector exclusive function
	 *
	 * Description:
	 * Final update of the position (x_next), velocity (v_next) and
	 * density (rho_nxt) of the examined particle in the Predictor-Corrector scheme.
	 *
	 */

	x_next[0] = 2 * x_next_bar[0] - x[0];
	x_next[1] = 2 * x_next_bar[1] - x[1];

	v_next[0] = 2 * v_next_bar[0] - v[0];
	v_next[1] = 2 * v_next_bar[1] - v[1];

	rho_nxt = 2 * rho_nxt_bar - rho;

}

void SPH_main::allocate_to_grid(void) {
	/*
	 * Description:
	 * Allocates each particle in the grid each time that 
	 * all the particles have their positions updated.
	 *
	 */

	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);

}

void SPH_main::neighbour_iterate(SPH_particle* part) {
	/*
	 * Description:
	 * Predictor - corrector half-step.
	 * Iterates over all particles within 2h of part
	 * and calculates the examined particle's position, velocity,
	 * density and pressure based on the corresponding summation terms.
	 *
	 * Could made more efficient using a stencil, considering symmetric interactions.
	 */

	vector<int> stencil_i = { part->list_num[0] -1, part->list_num[0] -1, part->list_num[0] -1, part->list_num[0],  part->list_num[0] + 1};
	vector<int> stencil_j = { part->list_num[1] -1, part->list_num[1], part->list_num[1] + 1, part->list_num[1] + 1, part->list_num[1] + 1};
	double dt_fact = 0.5;   // time step factor (0.5 for Predictor-Corrector)
	SPH_particle* other_part;
	double dist;			// distance between particles
	double dn[2];			// vector from 1st to 2nd particle
	double eij[2];
	double vij[2];			// relative particle velocites
	double dt = 0.1 * h / 20.0 * dt_fact;
	double smooth;

	part->update_pos(dt);  // updates particle position
	int i, j;
	unsigned int cnt;

	// Iterate through neighbours using stencil defined by stncil_i & stencil_j vectors
	for (int nn = 0; nn < stencil_i.size(); nn++) {

		i = stencil_i[nn]; j = stencil_j[nn];  // Find i and j indices from stencil vectors

		for (cnt = 0; cnt < search_grid[i][j].size(); cnt++) {
			other_part = search_grid[i][j][cnt];

			if (i >= 0 && i < max_list[0] && j >= 0 && j < max_list[1]) { //stops particle interacting with itself

				//Calculates the distance between potential neighbours
				for (int n = 0; n < 2; n++)
					dn[n] = part->x[n] - other_part->x[n];

				dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

				if (dist < 2. * h) { // Only particle within 2h

					// All the interactions between the particles
					eij[0] = dn[0] / dist;
					eij[1] = dn[1] / dist;
					vij[0] = part->v[0] - other_part->v[0];
					vij[1] = part->v[1] - other_part->v[1];
					smooth = part->grad_kernel(dist / h, h);  // Calculate smoothing length

					// Summation terms for velocity and density
					part->update_rho(dt, other_part, smooth, eij, vij);
					part->update_vel(dist, dt, other_part, smooth, eij, vij);

				}
			}
		}
	}

	// Necessary addition for predictor-corrector half-step.
	// Adds v, rho, and gravity to the corresponding intermediate summation terms.
	// Also updates pressure value.
	// Note that the addition is positioned outside the neighbour-scanning
	// since it occurs only once for every particle.
	part->v_next[0] += part->v[0];
	part->v_next[1] += part->v[1];
	part->v_next[1] -= dt * 9.81;
	part->rho_nxt += part->rho;
	part->set_pressure_nxt();

}

void SPH_main::neighbour_iterate2(SPH_particle* part) {
	/*
	 * Description:
	 * Predictor - corrector full-step.
	 * Iterates over all particles within 2h of part
	 * and calculates the examined particle's position, velocity
	 * and density, based on the corresponding summation terms.
	 *
	 * Input:
	 * part: current particle of SPH_particle class
	 *
	 */

	double dt_fact = 0.5;  // dt factor (here 0.5)
	double dt = 0.1 * h / 20.0 * dt_fact;  
	SPH_particle* other_part;
	double dist;			// distance between particles
	double dn[2];			// vector from 1st to 2nd particle
	double vij[2];			// relative particle velocites
	double eij[2];
	double smooth;

	part->update_pos_bar(dt);
	int i, j;
	unsigned int cnt;

	// Iterate through neighbours
	for (i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1]) {
					for (cnt = 0; cnt < search_grid[i][j].size(); cnt++) {
						other_part = search_grid[i][j][cnt];

						if (part != other_part) {  //stops particle interacting with itself

							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
								dn[n] = part->x_next[n] - other_part->x_next[n];

							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

							if (dist < 2. * h) {  // only particle within 2h

								//All the interactions between the particles
								eij[0] = dn[0] / dist;
								eij[1] = dn[1] / dist;
								vij[0] = part->v_next[0] - other_part->v_next[0];
								vij[1] = part->v_next[1] - other_part->v_next[1];
								smooth = part->grad_kernel(dist / h, h);  // Calculate smoothing length

								// Makes sure that density goes to zero for summation term of the velocity.
								if (other_part->rho_nxt < 1) { other_part->rho_nxt = 1000.0; }

								// Summation terms for velocity and density
								part->update_vel_bar(dist, dt, other_part, smooth, eij, vij);
								part->update_rho_bar(dt, other_part, smooth, eij, vij);

							}
						}
					}
				}

	// Necessary addition of v, rho and gravity terms for predictor-corrector full-step.
	part->v_next_bar[0] += part->v[0];
	part->v_next_bar[1] += part->v[1];
	part->rho_nxt_bar += part->rho;
	part->v_next_bar[1] -= dt * 9.81;

}

void SPH_main::smooth_pressure(SPH_particle* part) {
	/*
	 * Description:
	 * Pressure additional smoothening technique.
	 * Called every few time-steps in the main().
	 *
	 * Input:
	 * part: current particle of SPH_particle class
	 *
	 */


	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle
	double eij[2];
	double smooth;
	double sum_1(0), sum_2(0);
	int i, j;
	unsigned int cnt;

	for (i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1]) {
					for (cnt = 0; cnt < search_grid[i][j].size(); cnt++) {
						other_part = search_grid[i][j][cnt];

						//Calculates the distance between potential neighbours
						for (int n = 0; n < 2; n++)
							dn[n] = part->x[n] - other_part->x[n];

						dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

						if (dist < 2. * h) {  //only particle within 2h

							// All the interactions between the particles
							eij[0] = dn[0] / dist;
							eij[1] = dn[1] / dist;
							smooth = part->kernel(dist / h, h);
							sum_1 += smooth;
							sum_2 += smooth / other_part->rho;
						}
					}
				}

	// Update density
	part->rho_nxt = sum_1 / sum_2;
}
