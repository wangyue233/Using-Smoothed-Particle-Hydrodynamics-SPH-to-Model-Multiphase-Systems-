#pragma once

#include <vector>
#include "SPH_2D_Half.h"

int write_file(const char* filename,
	       std::vector<SPH_particle> *particle_list);
