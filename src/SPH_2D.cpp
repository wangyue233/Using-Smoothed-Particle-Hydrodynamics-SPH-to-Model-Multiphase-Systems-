#include "SPH_2D.h"
#include <math.h> 

#define G 9.8
#define H 1
#define N 6
#define F_LIMIT 0.5
#define C0 20.0
#define RHO0 1000.0
#define GAMMA 7.0


SPH_main* SPH_particle::main_data;

void SPH_particle::calc_index(void)
{
/*
 * Description:
 * Calculates the new particle index based on the updated position, x.
 *
 */
	for (int i = 0; i < 2; i++)
		list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0 * main_data->h));
}
void SPH_particle::set_pressure(void)
{
/* 
 * Description:
 * Calculates pressure based on rho (current density).
 *
 */
	P = RHO0 * C0 * C0 / GAMMA * (pow(rho / RHO0, GAMMA) - 1);
}

double SPH_particle::kernel(double q, double h)
{
/*
 * Description:
 * Applies Kernel function to calculate weights for particle interaction.
 *
 * Inputs:
 * q: distance of interaction
 * h: smoothening length (proglem-dependent constant)
 *
 */
	if (q <= 1)
	{
		return 10.0 / (7.0 * M_PI * h * h) * (1.0 - 3.0 / 2.0 * q * q + 3.0 / 4.0 * q * q * q);
	}
	else if (q <= 2)
	{
		return 10.0 / (7.0 * M_PI * h * h) * (1.0 / 4.0 * pow(2 - q, 3));
	}
	else {
		return 0;
	}

}

double SPH_particle::grad_kernel(double q, double h)
{
/*
 * Description:
 * Applies the derivative of the Kernel function to calculate weights for particle interaction.
 *
 * Inputs:
 * q: distance of interaction
 * h: smoothening length (proglem-dependent constant)
 *
 */
	if (q <= 1)
	{
		return 10.0 / (7.0 * M_PI * h * h * h) * (-3 * q + 9.0 / 4 * q * q);
	}
	else if (q <= 2)
	{
		return 10.0 / (7.0 * M_PI * h * h * h) * (-3.0 / 4.0 * pow(2 - q, 2));
	}
	else {
		return 0;
	}

}

void SPH_particle::set_initial(void) {
/* 
 * Description:
 * Initialises particle's v, v_next, v_prv, 
 * rho, rho_nxt_, rho_prv and mass.
 *
 */
	// set velocity
	v[0] = v[1] = 0;
	v_next[0] = v_next[1] = 0;
	v_prv[0] = v[0];
	v_prv[1] = v[1];
	// set positions
	x_prv[0] = x[0];
	x_prv[1] = x[1];
	//set density
	rho = 1000.0;
	rho_nxt = rho;
	rho_prv = rho;
	//set mass
	mass = main_data->dx * main_data->dx * rho;
	// set pressure
	set_pressure();
}

void SPH_particle::update_beach(void) {
/*
 * Description:
 * Sets initial dimensions and problem constants (dx, h).
 *
 */

 // rho = rho_nxt;
 // set_pressure();
 double dt = 0.1 * 0.2 * 1.3 / 20;
 // v_next[1] -= dt * 9.8;

 if (x_next[0] < 0.1)
 {
  v_next[0] += dt * 20;
 }

 if (x_next[0] > 20 - 0.1)
 {
  v_next[0] -= dt * 20;
 }
 if (x_next[1] < 0.1 && x_next[0] < 10)
 {
  v_next[1] += dt * 20;
 }
 if (x_next[1] < define_beach(x_next[0]) + 0.1)
 {
  if (x_next[0] >= 10 && x_next[0] <= 15)
  {
   v_next[0] = 0;
   v_next[1] = 1;
  }

 }
 if (x_next[0] > 15 && x_next[1] < 4.9)
 {
  v_next[1] += dt * 20;
 }
 if (x_next[1] > 10 - 0.1)
 {
  v_next[1] -= dt * 20;
 }
 for (int i = 0; i < 2; i++)
 {
  v[i] = v_next[i];
  x[i] = x_next[i];
 }
 if (x[0] < 0) {
  x[0] = -x[0] / 2;
  v[1] = 0;
  v[0] = 0;
 }
 if (x[1] < 0) {
  x[1] = -x[1] / 2;
  v[0] = 0;
  v[1] = 0;
 }
 if (x[0] > 20.0)
 {
  x[0] = 20.0 - x[0] + 20.0;
  v[0] = 0;
  v[1] = 0;
 }
 if (x[1] > 10.0)
 {
  x[1] = 10.0 - x[1] + 10.0;
  v[0] = 0;
  v[1] = 0;
 }
 if (x[0] > 15 && x[1] < 4.9)
 {
  x[1] = 10.0 - x[1] + 10.0;
  v[0] = 0;
  v[1] = 0;
 }
 if (x[0] >= 10 && x[0] <= 15){
  if (determine_beach(x[0], x[1])) {
  x[1] = define_beach(x[0]) + 0.1;
  v[0] = 0;
  v[1] = 0;
 }
 }
}

double SPH_particle::boundary_acc(double dx, double xab)
{
/*
 * Description:
 * Applies the force to let boundary go waay with boundary.
 *
 * Inputs:
 * dx: distance of initial grid
 * xab: distance from particle to boundary
 *
 */
	return G * H / dx * pow(F_LIMIT * dx / xab, N) * (1 - xab / dx / F_LIMIT);
}

void SPH_particle::update(void) {
/*
 * Description:
 * Updates particle x, v, rho and P based on the
 * corresponding newest calculated *_next values.
 * Also applies boundary conditions on the particle.
 *
 */

	double dt = main_data->dt;
	double dx = main_data->dx;

	for (int i = 0; i < 2; i++) {
		x[i] = x_next[i];
		v[i] = v_next[i];
	}
	rho = rho_nxt;

	if (x[0] < F_LIMIT * dx)
		v[0] += dt * boundary_acc(dx, x[0]);
	else if (x[0] > main_data->max_x[0] - 2 * main_data->h - F_LIMIT * dx)
		v[0] -= dt * boundary_acc(dx, main_data->max_x[0] - 2 * main_data->h - x[0]);

	if (x[1] < F_LIMIT * dx)
		v[1] += dt * boundary_acc(dx, x[1]);
	else if (x[1] > main_data->max_x[1] - 2 * main_data->h - F_LIMIT * dx)
		v[1] -= dt * boundary_acc(dx, main_data->max_x[1] - 2 * main_data->h - x[1]);

}
void SPH_particle::update_rho(double& rho2, double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2])
{
/*
 * Description:
 * Calulates summation of density values.
 *
 * Inputs
 * dt: time-step
 * other_part: neighbour particle
 * smooth: smoothing length
 * eij: relative distance between current and neighbour particle
 * vij: relative velocity between current and neighbour particle
 *
 */
	rho2 += dt * other_part->mass * smooth * (eij[0] * vij[0] + eij[1] * vij[1]);
}

void SPH_particle::update_vel(double* v_next, double dist, double dt, SPH_particle* other_part, double smooth, double eij[2], double vij[2])
{
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
	v_next[0] -= dt * other_part->mass * (other_part->P / other_part->rho / other_part->rho + P / rho / rho) * smooth * eij[0];
	v_next[0] += 0.001 * dt * other_part->mass * (1 / other_part->rho / other_part->rho + 1 / rho / rho) * smooth * vij[0] / dist;

	v_next[1] -= dt * other_part->mass * (other_part->P / other_part->rho / other_part->rho + P / rho / rho) * smooth * eij[1];
	v_next[1] += 0.001 * dt * other_part->mass * (1 / other_part->rho / other_part->rho + 1 / rho / rho) * smooth * vij[1] / dist;

}

void  SPH_particle::update_pos(double dt)
{
/*
 * Description:
 * Calulates position x_next.
 * For the Predictor-Corrector, it is called in the half-step.
 *
 * Inputs:
 * dt: time-step
 *
 */
	x_next[0] += dt * v[0];
	x_next[1] += dt * v[1];

}



SPH_main::SPH_main()
{
	SPH_particle::main_data = this;
}

void SPH_main::set_values(double idx)
{
/*
 * Description:
 * Sets initial dimensions and problem constants (dx, h).
 *
 */
	min_x[0] = 0.0;  // Set min values for grid --> min_x[0] = x minimum value ; 
	min_x[1] = 0.0;  // Set min values for grid --> min_x[1] = y minimum value ;

	max_x[0] = 20.0;
	max_x[1] = 10.0;

	dx = idx;

	h_fac = 1.3;
	h = dx * h_fac;    // obtain the smoothing lenght  

	dt = 0.1 * h / 20;
}

void SPH_main::initialise_grid(void)
{
/*
 * Description:
 * Initialises grid of particles to be used in computations.
 *
 */
	for (int i = 0; i < 2; i++)
	{
		min_x[i] -= 2.0 * h;
		max_x[i] += 2.0 * h;  //add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1.0);  // number of elements in i-direction of the grid
	}

	search_grid.resize(max_list[0]);
	for (int i = 0; i < max_list[0]; i++)
		search_grid[i].resize(max_list[1]);
}

double define_beach(double x) {
/*
 * Description:
 * Initialises grid of particles to be used in sloop boundary.
 * Inputs:
 * x: x coordianate
 */
	if (x < 10)
	{
		return 0;
	}
	else if (x < 15)
	{
		return x - 10;
	}
	else {
		return 5;
	}
}

bool determine_beach(double x, double y) {
/*
 * Description:
 * Initialises grid of particles to be used in sloop boundary.
 * Inputs:
 * x: x coordianate
 * y: y coordianate
 */

	if (y <= define_beach(x)) {
		return true;
	}
	else {
		return false;
	}
}

double gaussian(double x) {
/*
 * Description:
 * Initialises grid of particles to be used as the initial wave.
 * Inputs:
 * x: x coordianate
 */
	return 6 + 1.35 * exp(-(x - 2) * (x - 2) / 2.0);
}

bool determine_gaussian(double x, double y) {
/*
 * Description:
 * Initialises grid of particles to be used as the initial wave.
 * Inputs:
 * x: x coordianate
 * y: y coordianate
 */
	if (y <= gaussian(x)) {
		return true;
	}
	else {
		return false;
	}
}

void SPH_main::place_beach(double *min, double *max)
{
/*
 * Description:
 * Initialises boundaries and fluid particles.
 *
 * Inputs:
 * min: Minimum dimension
 * max: Maximum dimension
 * 
 *
 */
	SPH_particle particle;
	int right_index = 0;
	int top_index = 0;
	double top_max = 0;
	//bottom boundary
	for(int i = 0; i * dx + min[0] < max[0]; i++){
		double x_pos = min[0] + i * dx;
		particle.x[0] = x_pos;
		particle.x_next[0] = x_pos;
		if(right_index == 0 && x_pos > max[0] - 2 * h){
			right_index = i;
		}
		for(int j = 0; j * dx + min[1] < 0; j++){
			particle.x[1] = min[1] + j * dx;
			particle.x_next[1] = min[1] + j * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
		}
	}
	//left boundary
	for(int j = 0; j * dx + min[1] < max[1]; j++){
		double y_pos = min[1] + j * dx;
		if(y_pos < 0) continue;
		particle.x[1] = y_pos;
		particle.x_next[1] = y_pos;
		if(top_max == 0 && y_pos > max[1] - 2 * h){
			top_max = y_pos;
			top_index = j;
		}
		for(int i = 0; i * dx + min[0] < 0; i++){
			particle.x[0] = min[0] + i * dx;
			particle.x_next[0] = min[0] + i * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
		}
	}
	//top boundary
	for(int i = 0; i * dx + min[0] < max[0]; i++){
		double x_pos = min[0] + i * dx;
		if(x_pos < 0) continue;
		particle.x[0] = x_pos;
		particle.x_next[0] = x_pos;
		for(int j = top_index; j * dx + min[1] < max[1]; j++){
			particle.x[1] = min[1] + j * dx;
			particle.x_next[1] = min[1] + j * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
		}
	}
	//right boundary
	for(int j = 0; j * dx + min[1] < max[1]; j++){
		double y_pos = min[1] + j * dx;
		if(y_pos < 0 || y_pos >= top_max) continue;
		particle.x[1] = y_pos;
		particle.x_next[1] = y_pos;
		if(top_max == 0 && y_pos > max[1] - 2 * h){
			top_max = y_pos;
		}
		for(int i = right_index; i * dx + min[0] < max[0]; i++){
			particle.x[0] = min[0] + i * dx;
			particle.x_next[0] = min[0] + i * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
		}
	}

	// initial set for mobile fluid particles

	for(int j = 0; min[1] + j * dx <= 7; j++){
		double pos_y = min[1] + j * dx;
		if(pos_y < 0) continue;
		particle.x[1] = pos_y;
		particle.x_next[1] = pos_y;
		for(int i = 0; i < right_index; i++){
			double pos_x = min[0] + i * dx;
			if(pos_x < 0 || !determine_gaussian(pos_x, pos_y)) continue;
			particle.x[0] = pos_x;
			particle.x_next[0] = pos_x;
			particle.bounadry = determine_beach(pos_x, pos_y);
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
		}
	}
}

void SPH_main::place_points(double* min, double* max)
{
/*
 * Description:
 * Initialises boundaries and fluid particles.
 *
 * Inputs:
 * min: Minimum dimension
 * max: Maximum dimension
 * 
 *
 */
	SPH_particle particle;
	//bottom boundary
	int num_x = int((max[0] - min[0]) / dx);
	int num_y = int((0 - min[1]) / dx);
	double length = min[0] + num_x * dx;
	for (int j = 0; j <= num_y; j++)
	{
		particle.x[1] = 0 - j * dx;
		particle.x_next[1] = 0 - j * dx;
		for (int i = 0; i <= num_x; i++)
		{
			particle.x[0] = min[0] + i * dx;
			particle.x_next[0] = min[0] + i * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
			//number++;
		}
	}
	//top bounary
	for (int j = 0; j <= num_y; j++)
	{
		particle.x[1] = max[1] - j * dx;
		particle.x_next[1] = max[1] - j * dx;
		for (int i = 0; i <= num_x; i++)
		{
			particle.x[0] = min[0] + i * dx;
			particle.x_next[0] = min[0] + i * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
			//number++;
		}
	}

	//left boundary
	num_x = int((2 * h) / dx);
	num_y = int((max[1] - 2 * h) / dx);
	for (int i = 0; i <= num_x; i++)
	{
		particle.x[0] = min[0] + i * dx;
		particle.x_next[0] = min[0] + i * dx;
		for (int j = 0; j <= num_y; j++)
		{
			particle.x[1] = j * dx;
			particle.x_next[1] = j * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
			//number++;
		}
	}
	//right boundary
	for (int i = 0; i <= num_x; i++)
	{
		// particle.x[0] = max[0] - i * dx;
		// particle.x_next[0] = max[0] - i * dx;
		particle.x[0] = length - i * dx;
		particle.x_next[0] = length - i * dx;
		for (int j = 0; j <= num_y; j++)
		{
			particle.x[1] = j * dx;
			particle.x_next[1] = j * dx;
			particle.bounadry = true;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
			//number++;
		}
	}
	// initial set for mobile fluid particles
	// left part

	num_y = int(5 / dx);
	for (int i = 1; i <= num_y; i++) {
		particle.x[1] = i * dx;
		particle.x_next[1] = i * dx;
		for (int j = 0; j * dx < max_x[0] - 3 * dx; j++) {
			if (i * dx > 2 && j * dx > 3) break;
			particle.x[0] = min_x[0] + (j + 3) * dx;
			particle.x_next[0] = min_x[0] + (j + 3) * dx;
			particle.bounadry = false;
			particle.calc_index();
			particle.set_initial();
			particle_list.push_back(particle);
		}
	}


}


void SPH_main::allocate_to_grid(void)				//needs to be called each time that all the particles have their positions updated
{
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
	{
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
}


void SPH_main::neighbour_iterate(SPH_particle* part)					//iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
/*
 * Description:
 * Predictor - corrector full-step.
 * Iterates over all particles within 2h of part
 * and calculates the examined particle's position, velocity,
 * density and pressure based on the corresponding summation terms.
 *
 * Could made more efficient using a stencil, considering symmetric interactions.
 */
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle
	double eij[2];
	double vij[2];
	double smooth;


	// update x_next
	part->x_next[0] = part->x[0] + dt * part->v[0];
	part->x_next[1] = part->x[1] + dt * part->v[1];
	part->rho_nxt = part->rho;
	part->v_next[0] = part->v[0];
	part->v_next[1] = part->v[1];

	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1])
				{
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						other_part = search_grid[i][j][cnt];

						if (part != other_part)					//stops particle interacting with itself
						{
							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
								dn[n] = part->x[n] - other_part->x[n];

							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

							if (dist < 2. * h)					//only particle within 2h
							{
								//TODO: all the interactions between the particles
								eij[0] = dn[0] / dist;
								eij[1] = dn[1] / dist;
								vij[0] = part->v[0] - other_part->v[0];
								vij[1] = part->v[1] - other_part->v[1];
								smooth = part->grad_kernel(dist / h, h);
								part->update_rho(part->rho_nxt, dt, other_part, smooth, eij, vij);
								part->update_vel(part->v_next, dist, dt, other_part, smooth, eij, vij);

								// cout << "dn: " << dn[0] << " " << dn[1] << endl;		//Should be removed from the code - simply here for you to see that it is working
							}
						}
					}
				}
	part->v_next[1] -= G * 0.5 * dt;
	part->set_pressure();
}


void SPH_main::smooth_pressure(SPH_particle* part)					//iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
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

	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1])
				{
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						other_part = search_grid[i][j][cnt];

						//Calculates the distance between potential neighbours
						for (int n = 0; n < 2; n++)
							dn[n] = part->x[n] - other_part->x[n];

						dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

						if (dist < 2. * h)					//only particle within 2h
						{
							//TODO: all the interactions between the particles
							eij[0] = dn[0] / dist;
							eij[1] = dn[1] / dist;
							smooth = part->kernel(dist / h, h);
							sum_1 += smooth;
							sum_2 += smooth / other_part->rho;
							// cout << "dn: " << dn[0] << " " << dn[1] << endl;		//Should be removed from the code - simply here for you to see that it is working
						}

					}
				}
	part->rho_nxt = sum_1 / sum_2;
}



