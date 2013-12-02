#ifndef __SPHSYSTEM_H__
#define __SPHSYSTEM_H__

#include "sph_type.h"

struct Particle
{
public:
	uint id;
	float3 pos;
	float3 vel;
	float3 acc;
	float3 ev;

	float dens;
	float pres;

	float surf_norm;
	float energy;
	uint surf;
	uint color;
	float dist;
	uint root;
	uint status;

	uint level;

	Particle *next;
};

struct Cell
{
	uint num;
	Particle *start;
};

struct Point 
{
public:
	float3 pos;
	uint3 vox;
	float value;
	uint tag;
};

class SPHSystem
{
public:
	float div;

	uint max_particle;
	uint num_particle;

	float kernel[2];
	float mass[2];

	float3 world_size;
	float cell_size;
	uint3 grid_size;
	uint tot_cell;

	uint vox_div;
	float vox_step;
	uint3 vox_size;
	uint tot_vox;

	uint scalar_div;
	float scalar_step;
	uint3 scalar_size;
	uint tot_scalar;

	float3 gravity;
	float wall_damping;
	float rest_density;
	float gas_constant;
	float viscosity;
	float time_step;
	float surf_norm;
	float surf_coe;
	float dist_split;
	float dist_merge;
	float energy_split;
	float energy_merge;

	uint enable_split;
	uint enable_merge;

	float poly6_value[2];
	float spiky_value[2];
	float visco_value[2];

	float grad_poly6[2];
	float lplc_poly6[2];

	float kernel_2[2];
	float self_dens[2];
	float self_lplc_color[2];

	Particle *mem;
	Particle *buf;
	Cell *cell;
	Point *point;
	float *scalar;
	float3 *scalar_pos;

	uint propagation;
	uint tick;

	uint sys_running;
	uint sys_render;

public:
	SPHSystem();
	~SPHSystem();
	void animation();
	void init_system();
	void add_particle(float3 pos, float3 vel, uint level);

public:
	void build_table();
	void comp_dens_pres();
	void comp_force_adv();
	void comp_scalar();
	void find_surface();
	void comp_dist();
	void comp_energy();
	void advection();
	void comp_split();
	void comp_merge();
	void interp_scalar();

private:
	int3 calc_cell_pos(float3 p);
	uint calc_cell_hash(int3 cell_pos);
	int3 calc_vox_pos(float3 p);
};

#endif