#include "sph_system.h"
#include "sph_header.h"
#include "sph_math.h"

SPHSystem::SPHSystem()
{
	div=3.0f;
	max_particle=1500000;
	num_particle=0;

	kernel[LEVEL1]=0.04f;
	kernel[LEVEL2]=kernel[LEVEL1]/pow(div, 1.0f/3.0f);

	mass[LEVEL1]=0.02f;
	mass[LEVEL2]=mass[LEVEL1]/div;

	world_size=make_float3(1.28f, 0.64f, 0.64f);
	cell_size=kernel[LEVEL1];
	grid_size.x=(uint)ceil(world_size.x/cell_size);
	grid_size.y=(uint)ceil(world_size.y/cell_size);
	grid_size.z=(uint)ceil(world_size.z/cell_size);
	tot_cell=grid_size.x*grid_size.y*grid_size.z;

	vox_div=4;
	vox_step=cell_size/vox_div;
	vox_size.x=(uint)ceil(world_size.x/vox_step)+3;
	vox_size.y=(uint)ceil(world_size.y/vox_step)+3;
	vox_size.z=(uint)ceil(world_size.z/vox_step)+3;
	tot_vox=vox_size.x*vox_size.y*vox_size.z;

	scalar_div=4;
	scalar_step=cell_size/scalar_div;
	scalar_size.x=(uint)ceil(world_size.x/scalar_step)+1;
	scalar_size.y=(uint)ceil(world_size.y/scalar_step)+1;
	scalar_size.z=(uint)ceil(world_size.z/scalar_step)+1;
	tot_scalar=scalar_size.x*scalar_size.y*scalar_size.z;

	gravity.x=0.0f; 
	gravity.y=-3.2f;
	gravity.z=0.0f;
	wall_damping=-0.5f;
	rest_density=1000.0f;
	gas_constant=1.0f;
	viscosity=10.5f;
	time_step=0.001f;
	surf_norm=6.0f;
	surf_coe=0.1f;
	dist_split=0.10f;
	dist_merge=0.12f;
	energy_split=8.0f;
	energy_merge=5.0f;

	enable_split=0;
	enable_merge=0;

	poly6_value[LEVEL1]=315.0f/(64.0f * PI * pow(kernel[LEVEL1], 9));
	spiky_value[LEVEL1]=-45.0f/(PI * pow(kernel[LEVEL1], 6));
	visco_value[LEVEL1]=45.0f/(PI * pow(kernel[LEVEL1], 6));

	poly6_value[LEVEL2]=315.0f/(64.0f * PI * pow(kernel[LEVEL2], 9));
	spiky_value[LEVEL2]=-45.0f/(PI * pow(kernel[LEVEL2], 6));
	visco_value[LEVEL2]=45.0f/(PI * pow(kernel[LEVEL2], 6));

	grad_poly6[LEVEL1]=-945/(32 * PI * pow(kernel[LEVEL1], 9));
	lplc_poly6[LEVEL1]=-945/(8 * PI * pow(kernel[LEVEL1], 9));

	grad_poly6[LEVEL2]=-945/(32 * PI * pow(kernel[LEVEL2], 9));
	lplc_poly6[LEVEL2]=-945/(8 * PI * pow(kernel[LEVEL2], 9));

	kernel_2[LEVEL1]=kernel[LEVEL1]*kernel[LEVEL1];
	kernel_2[LEVEL2]=kernel[LEVEL2]*kernel[LEVEL2];

	self_dens[LEVEL1]=mass[LEVEL1]*poly6_value[LEVEL1]*pow(kernel[LEVEL1], 6);
	self_dens[LEVEL2]=mass[LEVEL2]*poly6_value[LEVEL2]*pow(kernel[LEVEL2], 6);

	self_lplc_color[LEVEL1]=lplc_poly6[LEVEL1]*mass[LEVEL1]*kernel_2[LEVEL1]*(0-3/4*kernel_2[LEVEL1]);
	self_lplc_color[LEVEL2]=lplc_poly6[LEVEL2]*mass[LEVEL2]*kernel_2[LEVEL2]*(0-3/4*kernel_2[LEVEL2]);

	mem=(Particle *)malloc(sizeof(Particle)*max_particle);
	buf=(Particle *)malloc(sizeof(Particle)*max_particle);
	cell=(Cell *)malloc(sizeof(Cell)*tot_cell);
	point=(Point *)malloc(sizeof(Point)*tot_vox);
	scalar=(float *)malloc(sizeof(float)*tot_scalar);
	scalar_pos=(float3 *)malloc(sizeof(float3)*tot_scalar);

	sys_running=0;
	sys_render=0;
	tick=0;

	printf("Initialize SPH:\n");
	printf("Large Kernel: %f\n", kernel[LEVEL1]);
	printf("Small Kernel: %f\n", kernel[LEVEL2]);
	printf("World Width : %f\n", world_size.x);
	printf("World Height: %f\n", world_size.y);
	printf("World Length: %f\n", world_size.z);
	printf("Cell Size  : %f\n", cell_size);
	printf("Grid Width : %u\n", grid_size.x);
	printf("Grid Height: %u\n", grid_size.y);
	printf("Grid Length: %u\n", grid_size.z);
	printf("Total Cell : %u\n", tot_cell);
	printf("Large Poly6: %f\n", poly6_value[LEVEL1]);
	printf("Large Spiky: %f\n", spiky_value[LEVEL1]);
	printf("Large Visco: %f\n", visco_value[LEVEL1]);
	printf("Small Poly6: %f\n", poly6_value[LEVEL2]);
	printf("Large Spiky: %f\n", spiky_value[LEVEL2]);
	printf("Large Visco: %f\n", visco_value[LEVEL2]);
	printf("Large Self Density: %f\n", self_dens[LEVEL1]);
	printf("Small Self Density: %f\n", self_dens[LEVEL2]);
	printf("Vox Size  : %f\n", vox_step);
	printf("Vox Width : %u\n", vox_size.x);
	printf("Vox Height: %u\n", vox_size.y);
	printf("Vox Length: %u\n", vox_size.z);
	printf("Total Vox : %u\n", tot_vox);
	printf("Scalar Size  : %f\n", scalar_step);
	printf("Scalar Width : %u\n", scalar_size.x);
	printf("Scalar Height: %u\n", scalar_size.y);
	printf("Scalar Length: %u\n", scalar_size.z);
	printf("Total Scalar : %u\n", tot_scalar);

	uint index;
	for(uint i=0; i<vox_size.x; i++)
	{
		for(uint j=0; j<vox_size.y; j++)
		{
			for(uint k=0; k<vox_size.z; k++)
			{
				index=k*vox_size.y*vox_size.x+j*vox_size.x+i;

				point[index].pos.x=i*vox_step;
				point[index].pos.y=j*vox_step;
				point[index].pos.z=k*vox_step;

				point[index].vox.x=i;
				point[index].vox.y=j;
				point[index].vox.z=k;

				point[index].value=0.0f;
			}
		}
	}

	for(uint i=0; i<scalar_size.x; i++)
	{
		for(uint j=0; j<scalar_size.y; j++)
		{
			for(uint k=0; k<scalar_size.z; k++)
			{
				index=k*scalar_size.y*scalar_size.x+j*scalar_size.x+i;

				scalar_pos[index].x=i*scalar_step-scalar_step;
				scalar_pos[index].y=j*scalar_step-scalar_step;
				scalar_pos[index].z=k*scalar_step-scalar_step;
			}
		}
	}
}

SPHSystem::~SPHSystem()
{
	free(mem);
	free(cell);
	free(point);
}

void SPHSystem::animation()
{
	if(sys_running == 0)
	{
		return;
	}

	build_table();
	if(sys_render == 1)
	{
		interp_scalar();
	}

	if(enable_split == 1 && tick == 0 || enable_merge == 1 && tick == 1)
	{
		comp_scalar();
		find_surface();
		comp_dist();
		//comp_energy();
	}

	comp_dens_pres();
	comp_force_adv();

	if(enable_split == 1 && tick == 0)
	{
		comp_split();
	}

	if(enable_merge == 1 && tick == 1)
	{
		comp_merge();
	}

	advection();

	tick=1-tick;
}

void SPHSystem::init_system()
{
	float3 pos;
	float3 vel=make_float3(0.0f);

	for(pos.x=world_size.x*0.0f; pos.x<world_size.x*0.6f; pos.x+=(kernel[LEVEL1]*0.5f))
	{
		for(pos.y=world_size.y*0.0f; pos.y<world_size.y*1.0f; pos.y+=(kernel[LEVEL1]*0.5f))
		{
			for(pos.z=world_size.z*0.0f; pos.z<world_size.z*1.0f; pos.z+=(kernel[LEVEL1]*0.5f))
			{
				add_particle(pos, vel, LEVEL1);
			}
		}
	}

	/*for(pos.x=world_size.x*0.7f; pos.x<world_size.x*1.0f; pos.x+=(kernel[LEVEL2]*0.5f))
	{
		for(pos.y=world_size.y*0.0f; pos.y<world_size.y*1.0f; pos.y+=(kernel[LEVEL2]*0.5f))
		{
			for(pos.z=world_size.z*0.0f; pos.z<world_size.z*1.0f; pos.z+=(kernel[LEVEL2]*0.5f))
			{
				add_particle(pos, vel, 1);
			}
		}
	}*/

	printf("Init Particle: %u\n", num_particle);
}

void SPHSystem::add_particle(float3 pos, float3 vel, uint level)
{
	Particle *p=&(mem[num_particle]);

	p->id=num_particle;

	p->pos=pos;
	p->vel=vel;

	p->acc=make_float3(0.0f);
	p->ev=make_float3(0.0f);

	p->dens=rest_density;
	p->pres=0.0f;

	p->level=level;

	p->next=NULL;

	num_particle++;
}

void SPHSystem::build_table()
{
	Particle *p;
	uint hash;

	for(uint i=0; i<tot_cell; i++)
	{
		cell[i].start=NULL;
		cell[i].num=0;
	}

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]);
		hash=calc_cell_hash(calc_cell_pos(p->pos));

		p->surf=0;
		p->color=0;
		p->status=0;

		if(cell[hash].start == NULL)
		{
			p->next=NULL;
			cell[hash].start=p;
			cell[hash].num++;
		}
		else
		{
			p->next=cell[hash].start;
			cell[hash].start=p;
			cell[hash].num++;
		}
	}
}

void SPHSystem::comp_dens_pres()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float r2;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);
		float3 transform=make_float3(0.0f);
		float3 tot_trans=make_float3(0.0f);
		float3 vel_trans=make_float3(0.0f);

		p->dens=0.0f;
		p->pres=0.0f;

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;
					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash].start;
					while(np != NULL)
					{
						rel_pos=np->pos-p->pos;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

						if(r2<INF || r2>=kernel_2[p->level])
						{
							np=np->next;
							continue;
						}

						p->dens=p->dens + mass[np->level] * poly6_value[p->level] * pow(kernel_2[p->level]-r2, 3);

						if(enable_split == 1 && tick == 0 || enable_merge == 1 && tick == 1)
						{
							transform.x=wavelet(rel_pos.x/kernel[p->level]);
							transform.y=wavelet(rel_pos.y/kernel[p->level]);
							transform.z=wavelet(rel_pos.z/kernel[p->level]);

							tot_trans=tot_trans+transform;

							vel_trans.x=vel_trans.x+(np->vel.x*transform.x/sqrt(kernel[p->level]));
							vel_trans.y=vel_trans.y+(np->vel.y*transform.y/sqrt(kernel[p->level]));
							vel_trans.z=vel_trans.z+(np->vel.z*transform.z/sqrt(kernel[p->level]));
						}

						np=np->next;
					}
				}
			}
		}

		if(enable_split == 1 && tick == 0 || enable_merge == 1 && tick == 1)
		{
			vel_trans=vel_trans/tot_trans;
			p->energy=0.5f*(vel_trans.x*vel_trans.x+vel_trans.y*vel_trans.y+vel_trans.z*vel_trans.z);
		}

		p->dens=p->dens+self_dens[p->level];
		p->pres=(pow(p->dens / rest_density, 7) - 1) *gas_constant;
	}
}

void SPHSystem::comp_force_adv()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float3 rel_vel;

	float r2;
	float r;
	float i_kernel_r;
	float j_kernel_r;
	float iV;
	float jV;

	float i_press_kernel;
	float j_press_kernel;
	float i_visco_kernel;
	float j_visco_kernel;
	float temp_force;

	float3 grad_color;
	float lplc_color;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

		p->acc=make_float3(0.0f);
		grad_color=make_float3(0.0f);
		lplc_color=0.0f;
		
		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;
					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash].start;
					while(np != NULL)
					{
						rel_pos=p->pos-np->pos;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

						float max_kernel_2=kernel_2[p->level]>kernel_2[np->level]?kernel_2[p->level]:kernel_2[np->level];

						if(r2 < max_kernel_2 && r2 > INF)
						{
							r=sqrt(r2);
							iV=mass[p->level]/p->dens;
							jV=mass[np->level]/np->dens;
				
							i_kernel_r=kernel[p->level]-r;
							j_kernel_r=kernel[np->level]-r;

							if(i_kernel_r < 0.0f)
								i_kernel_r=0.0f;
							if(j_kernel_r < 0.0f)
								j_kernel_r=0.0f;

							i_press_kernel=spiky_value[p->level] * i_kernel_r * i_kernel_r;
							i_visco_kernel=visco_value[p->level] * (i_kernel_r);

							j_press_kernel=spiky_value[np->level] * j_kernel_r * j_kernel_r;
							j_visco_kernel=visco_value[np->level] * (j_kernel_r);

							temp_force=(iV+jV)/2 * (p->pres+np->pres) * (i_press_kernel+j_press_kernel)/2;
							p->acc=p->acc-rel_pos*temp_force/r;

							rel_vel=np->ev-p->ev;
							temp_force=(iV+jV)/2 * viscosity * (i_visco_kernel+j_visco_kernel)/2;
							p->acc=p->acc + rel_vel*temp_force; 

							float temp=(-1) * grad_poly6[p->level] * jV * pow(kernel_2[p->level]-r2, 2);
							grad_color += temp * rel_pos;
							lplc_color += lplc_poly6[p->level] * jV * (kernel_2[p->level]-r2) * (r2-3/4*(kernel_2[p->level]-r2));
						}

						np=np->next;
					}
				}
			}
		}

		lplc_color+=self_lplc_color[p->level]/p->dens;
		p->surf_norm=sqrt(grad_color.x*grad_color.x+grad_color.y*grad_color.y+grad_color.z*grad_color.z);

		//if(p->surf_norm > surf_norm)
		if(p->surf_norm > 0.0f && p->surf == 1)
		{
			p->acc +=surf_coe * lplc_color * grad_color / p->surf_norm;
			grad_color=grad_color/surf_norm;
		}

		p->vel=p->vel+p->acc*time_step/p->dens+gravity*time_step;
		float dist=grad_color.x*p->vel.x+grad_color.y*p->vel.y+grad_color.z*p->vel.z;
		p->vel=p->vel-grad_color*dist*0.01f*time_step/p->dens;
	}
}

void SPHSystem::comp_scalar()
{
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float r2;
	float r;
	float temp;
	float scalar;

	for(uint i=0; i<tot_vox; i++)
	{
		scalar=0.0f;
		cell_pos=calc_cell_pos(point[i].pos);

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;

					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash].start;
					while(np != NULL)
					{
						rel_pos=np->pos-point[i].pos;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;
						r=sqrt(r2);

						temp=r-kernel[np->level];

						if(temp < scalar)
						{
							scalar=temp;
						}

						np=np->next;
					}
				}
			}
		}
		point[i].value=scalar;
	}
}

void SPHSystem::find_surface()
{
	propagation=0;

	for(uint i=0; i<tot_vox; i++)
	{
		point[i].tag=0;
		if(point[i].value >= 0.0f)
		{
			continue;
		}

		uint ix=point[i].vox.x;
		uint iy=point[i].vox.y;
		uint iz=point[i].vox.z;

		uint near;

		if(iy == vox_size.y-1)
		{
			point[i].tag=1;
			goto INIT_DIST;
		}
		
		if(iy==0 || iy==vox_size.y-1)
		{
			continue;
		}

		near=iz*vox_size.y*vox_size.x + (iy+1)*vox_size.x + ix;
		if(point[near].value == 0.0f)
		{
			point[i].tag=1;
			goto INIT_DIST;
		}
		near=iz*vox_size.y*vox_size.x + (iy-1)*vox_size.x + ix;
		if(point[near].value == 0.0f)
		{
			point[i].tag=1;
			goto INIT_DIST;
		}

		if(ix==0 || ix==vox_size.x-1 || iz==0 || iz==vox_size.z-1)
		{
			continue;
		}

		near=iz*vox_size.y*vox_size.x + iy*vox_size.x + ix-1;
		if(point[near].value == 0.0f)
		{
			point[i].tag=1;
			goto INIT_DIST;
		}
		near=iz*vox_size.y*vox_size.x + iy*vox_size.x + ix+1;
		if(point[near].value == 0.0f)
		{
			point[i].tag=1;
			goto INIT_DIST;
		}

		near=(iz+1)*vox_size.y*vox_size.x + iy*vox_size.x + ix;
		if(point[near].value == 0.0f)
		{
			point[i].tag=1;
			goto INIT_DIST;
		}
		near=(iz-1)*vox_size.y*vox_size.x + iy*vox_size.x + ix;
		if(point[near].value == 0.0f)
		{
			point[i].tag=1;
			goto INIT_DIST;
		}

		continue;

INIT_DIST:

		int3 cell_pos=calc_vox_pos(point[i].pos);
		int3 near_pos;
		uint hash;
		Particle *np;
		float3 rel_pos;
		float r2;

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;

					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash].start;
					while(np != NULL)
					{
						rel_pos=np->pos-point[i].pos;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

						if(r2<kernel_2[np->level])
						{
							if(np->color != 1)
							{
								propagation++;
							}
							np->surf=1;
							np->color=1;
							np->dist=0.0f;
							np->root=np->id;
						}

						np=np->next;
					}
				}
			}
		}
	}
}

void SPHSystem::comp_dist()
{
	while(propagation != num_particle)
	{
		for(uint i=0; i<num_particle; i++)
		{
			if(mem[i].color == 1)
			{
				continue;
			}

			Particle *p=&(mem[i]);
			Particle *np;
			Particle *root;

			int3 cell_pos=calc_cell_pos(p->pos);
			int3 near_pos;
			uint hash;
			float current_dist;

			p->dist=MAX_DIST;

			for(int x=-1; x<=1; x++)
			{
				for(int y=-1; y<=1; y++)
				{
					for(int z=-1; z<=1; z++)
					{
						near_pos.x=cell_pos.x+x;
						near_pos.y=cell_pos.y+y;
						near_pos.z=cell_pos.z+z;

						hash=calc_cell_hash(near_pos);

						if(hash == 0xffffffff)
						{
							continue;
						}

						np=cell[hash].start;
						while(np != NULL)
						{
							if(np->color != 1 || np->id == i)
							{
								np=np->next;
								continue;
							}

							current_dist=length(np->pos-p->pos);

							if(current_dist > kernel[0])
							{
								np=np->next;
								continue;
							}

							root=&(mem[np->root]);
							current_dist=length(root->pos-p->pos);

							if(current_dist < p->dist)
							{
								if(p->color != 1)
								{
									propagation++;
								}
								p->color=1;
								p->dist=current_dist;
								p->root=np->root;
							}

							np=np->next;
						}
					}
				}
			}
		}
	}
}

void SPHSystem::comp_energy()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float r2;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

		float3 transform=make_float3(0.0f);
		float3 tot_trans=make_float3(0.0f);
		float3 vel_trans=make_float3(0.0f);

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;
					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash].start;
					while(np != NULL)
					{
						rel_pos=np->pos-p->pos;
						r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

						if(r2<INF || r2>=kernel_2[p->level])
						{
							np=np->next;
							continue;
						}

						transform.x=wavelet(rel_pos.x/kernel[p->level]);
						transform.y=wavelet(rel_pos.y/kernel[p->level]);
						transform.z=wavelet(rel_pos.z/kernel[p->level]);

						tot_trans=tot_trans+transform;

						vel_trans.x=vel_trans.x+(np->vel.x*transform.x/sqrt(kernel[p->level]));
						vel_trans.y=vel_trans.y+(np->vel.y*transform.y/sqrt(kernel[p->level]));
						vel_trans.z=vel_trans.z+(np->vel.z*transform.z/sqrt(kernel[p->level]));

						np=np->next;
					}
				}
			}
		}

		vel_trans=vel_trans/tot_trans;
		p->energy=0.5f*(vel_trans.x*vel_trans.x+vel_trans.y*vel_trans.y+vel_trans.z*vel_trans.z);
	}
}

void SPHSystem::comp_split()
{
	uint tot=num_particle;
	for(uint i=0; i<tot; i++)
	{
		if(mem[i].level == LEVEL2)
		{
			continue;
		}

		if(mem[i].pos.y < kernel[LEVEL1])
		{
			continue;
		}

		if(mem[i].pos.x< kernel[LEVEL1] && mem[i].pos.z< kernel[LEVEL1] && mem[i].surf!=1)
		{
			continue;
		}

		if(mem[i].pos.x< world_size.x-kernel[LEVEL1] && mem[i].pos.z< world_size.y-kernel[LEVEL1] && mem[i].surf!=1)
		{
			continue;
		}

		/*if(mem[i].pos.x < kernel[LEVEL1] || mem[i].pos.y < kernel[LEVEL1] || mem[i].pos.z < kernel[LEVEL1])
		{
			continue;
		}

		if(mem[i].pos.x > world_size.x-kernel[LEVEL1] || mem[i].pos.y > world_size.y-kernel[LEVEL1] || mem[i].pos.z > world_size.z-kernel[LEVEL1])
		{
			continue;
		}*/
		
		//dMem[index].level==0 && ( dMem[index].dist<dParam.dist_split || dMem[index].energy>dParam.energy_split 
		if(mem[i].dist < dist_split || mem[i].energy > energy_split)
		{
			uint new1st=i;
			uint new2nd=num_particle;
			uint new3rd=num_particle+1;
			//uint new4th=num_particle+2;

			//2
			mem[new2nd].level=LEVEL2;

			mem[new2nd].pos.x=mem[new1st].pos.x+kernel[LEVEL2]/3;
			mem[new2nd].pos.y=mem[new1st].pos.y;
			mem[new2nd].pos.z=mem[new1st].pos.z;

			mem[new2nd].vel=mem[new1st].vel;
			mem[new2nd].ev=mem[new1st].ev;
			mem[new2nd].acc=mem[new1st].acc;
			mem[new2nd].energy=mem[new1st].energy;

			//3
			mem[new3rd].level=LEVEL2;

			mem[new3rd].pos.x=mem[new1st].pos.x;
			mem[new3rd].pos.y=mem[new1st].pos.y+kernel[LEVEL2]/3;
			mem[new3rd].pos.z=mem[new1st].pos.z;

			mem[new3rd].vel=mem[new1st].vel;
			mem[new3rd].ev=mem[new1st].ev;
			mem[new3rd].acc=mem[new1st].acc;
			mem[new3rd].energy=mem[new1st].energy;

			/*//4
			mem[new4th].level=LEVEL2;

			mem[new4th].pos.x=mem[new1st].pos.x+kernel[LEVEL2]/2;
			mem[new4th].pos.y=mem[new1st].pos.y;
			mem[new4th].pos.z=mem[new1st].pos.z;

			mem[new4th].vel=mem[new1st].vel;
			mem[new4th].ev=mem[new1st].ev;
			mem[new4th].acc=mem[new1st].acc;
			mem[new4th].energy=mem[new1st].energy;*/

			//1
			mem[new1st].level=1;
			mem[new1st].pos.x=mem[new1st].pos.x;
			mem[new1st].pos.y=mem[new1st].pos.y;
			mem[new1st].pos.z=mem[new1st].pos.z-kernel[LEVEL2]/3;

			num_particle=num_particle+2;
		}
	}
}

void SPHSystem::comp_merge()
{
	Particle *p;
	Particle *np;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float r2;

	uint merge1st;
	uint merge2nd;
	uint merge3rd;

	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]); 
		cell_pos=calc_cell_pos(p->pos);

		if(p->status == 1 || p->level == LEVEL1 || p->dist <= dist_merge || p->energy >=energy_merge)
		{
			continue;
		}

		merge1st=p->id;
		uint two_merge=0;

		for(int x=-1; x<=1; x++)
		{
			for(int y=-1; y<=1; y++)
			{
				for(int z=-1; z<=1; z++)
				{
					near_pos.x=cell_pos.x+x;
					near_pos.y=cell_pos.y+y;
					near_pos.z=cell_pos.z+z;
					hash=calc_cell_hash(near_pos);

					if(hash == 0xffffffff)
					{
						continue;
					}

					np=cell[hash].start;
					while(np != NULL)
					{
						if(two_merge == 0)
						{
							if(np->status == 1 || np->level == LEVEL1 || np->id == p->id || np->dist <= dist_merge || p->energy >=energy_merge)
							{
								np=np->next;
								continue;
							}

							rel_pos=np->pos-p->pos;
							r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

							if(r2 <= kernel_2[LEVEL2]/2)
							{
								merge2nd=np->id;
								np=np->next;
								two_merge=1;
								continue;
								
							}
						}

						if(two_merge == 1)
						{
							if(np->status == 1 || np->level == LEVEL1 || np->id == p->id || np->dist <= dist_merge || p->energy >=energy_merge)
							{
								np=np->next;
								continue;
							}

							rel_pos=np->pos-p->pos;
							r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

							if(r2 <= kernel_2[LEVEL2]/2)
							{
								merge3rd=np->id;
								goto MERGE;
							}
						}

						np=np->next;
					}
				}
			}
		}

		continue;

		////////////////////////////////
MERGE: 
		mem[merge1st].level=LEVEL1;

		mem[merge1st].pos=(mem[merge1st].pos+mem[merge2nd].pos+mem[merge3rd].pos)/3;
		mem[merge1st].vel=(mem[merge1st].vel+mem[merge2nd].vel+mem[merge3rd].vel)/3;
		mem[merge1st].ev=(mem[merge1st].ev+mem[merge2nd].ev+mem[merge3rd].ev)/3;
		mem[merge1st].acc=(mem[merge1st].acc+mem[merge2nd].acc+mem[merge3rd].acc)/3;

		mem[merge1st].status=0;
		mem[merge2nd].status=1;
		mem[merge3rd].status=1;
	}

	uint tot=0;
	for(uint i=0; i<num_particle; i++)
	{
		if(mem[i].status == 1)
		{
			continue;
		}

		buf[tot]=mem[i];
		buf[tot].id=tot;
		tot++;
	}

	num_particle=tot;
	Particle *mid;
	mid=mem;
	mem=buf;
	buf=mid;
}

void SPHSystem::advection()
{
	Particle *p;
	for(uint i=0; i<num_particle; i++)
	{
		p=&(mem[i]);

		p->pos=p->pos+p->vel*time_step;

		if(p->pos.x >= world_size.x-BOUNDARY)
		{
			p->vel.x=p->vel.x*wall_damping;
			p->pos.x=world_size.x-BOUNDARY;
		}

		if(p->pos.x < BOUNDARY)
		{
			p->vel.x=p->vel.x*wall_damping;
			p->pos.x=BOUNDARY;
		}

		if(p->pos.y >= world_size.y-BOUNDARY)
		{
			p->vel.y=p->vel.y*wall_damping;
			p->pos.y=world_size.y-BOUNDARY;
		}

		if(p->pos.y < BOUNDARY)
		{
			p->vel.y=p->vel.y*wall_damping;
			p->pos.y=BOUNDARY;
		}

		if(p->pos.z >= world_size.z-BOUNDARY)
		{
			p->vel.z=p->vel.z*wall_damping;
			p->pos.z=world_size.z-BOUNDARY;
		}

		if(p->pos.z < BOUNDARY)
		{
			p->vel.z=p->vel.z*wall_damping;
			p->pos.z=BOUNDARY;
		}

		p->ev=(p->ev+p->vel)/2;
	}
}

void SPHSystem::interp_scalar()
{
	uint index;
	float3 pos;
	float3 weight_pos;
	float weight_kernel;
	float weight;

	int3 cell_pos;
	int3 near_pos;
	uint hash;

	float3 rel_pos;
	float r2;
	float r;
	Particle *np;
	
	for(uint i=0; i<scalar_size.x; i++)
	{
		for(uint j=0; j<scalar_size.y; j++)
		{
			for(uint k=0; k<scalar_size.z; k++)
			{
				index=k*scalar_size.y*scalar_size.x+j*scalar_size.x+i;

				if(i==0 || j==0 || k==0 || i==scalar_size.x-1 || j==scalar_size.y-1 || k==scalar_size.z-1)
				{
					scalar[index]=0.0f;
					continue;
				}

				weight_pos=make_float3(0.0f);
				weight_kernel=0.0f;
				weight=0.0f;

				pos=scalar_pos[index];
				cell_pos=calc_vox_pos(pos);
	
				for(int x=-1; x<=1; x++)
				{
					for(int y=-1; y<=1; y++)
					{
						for(int z=-1; z<=1; z++)
						{
							near_pos.x=cell_pos.x+x;
							near_pos.y=cell_pos.y+y;
							near_pos.z=cell_pos.z+z;
							hash=calc_cell_hash(near_pos);

							if(hash == 0xffffffff)
							{
								continue;
							}

							np=cell[hash].start;
							while(np != NULL)
							{
								rel_pos=np->pos-pos;
								r2=rel_pos.x*rel_pos.x+rel_pos.y*rel_pos.y+rel_pos.z*rel_pos.z;

								if(r2>kernel_2[np->level])
								{
									np=np->next;
									continue;
								}

								r=sqrt(r2);
								float w=polyk(r/kernel[np->level]);
								weight_pos=weight_pos+w*np->pos;
								weight_kernel=weight_kernel+w*kernel[np->level];
								weight=weight+w;

								np=np->next;
							}
						}
					}
				}
				if(weight == 0.0f)
				{
					scalar[index]=0.0f;
				}
				else
				{
					weight_pos=weight_pos/weight;
					weight_kernel=weight_kernel/weight;
					scalar[index]=weight_kernel-length(pos-weight_pos);
				}
			}
		}
	}
}

int3 SPHSystem::calc_cell_pos(float3 p)
{
	int3 cell_pos;
	cell_pos.x = int(floor((p.x) / cell_size));
	cell_pos.y = int(floor((p.y) / cell_size));
	cell_pos.z = int(floor((p.z) / cell_size));

    return cell_pos;
}

uint SPHSystem::calc_cell_hash(int3 cell_pos)
{
	if(cell_pos.x<0 || cell_pos.x>=(int)grid_size.x || cell_pos.y<0 || cell_pos.y>=(int)grid_size.y || cell_pos.z<0 || cell_pos.z>=(int)grid_size.z)
	{
		return (uint)0xffffffff;
	}

	cell_pos.x = cell_pos.x & (grid_size.x-1);  
    cell_pos.y = cell_pos.y & (grid_size.y-1);  
	cell_pos.z = cell_pos.z & (grid_size.z-1);  

	return ((uint)(cell_pos.z))*grid_size.y*grid_size.x + ((uint)(cell_pos.y))*grid_size.x + (uint)(cell_pos.x);
}

int3 SPHSystem::calc_vox_pos(float3 p)
{
	int3 grid_pos;

    grid_pos.x = int(floor((p.x) / cell_size));
    grid_pos.y = int(floor((p.y) / cell_size));
	grid_pos.z = int(floor((p.z) / cell_size));

	if(p.x >= world_size.x)
	{
		grid_pos.x=grid_size.x-1;
	}

	if(p.y >= world_size.y)
	{
		grid_pos.y=grid_size.y-1;
	}

	if(p.z >= world_size.z)
	{
		grid_pos.z=grid_size.z-1;
	}

    return grid_pos;
}
