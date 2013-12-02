#ifndef __MARCHINGCUBE_H__
#define __MARCHINGCUBE_H__

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sph_header.h"

typedef unsigned int uint;

class MarchingCube
{
private:

	uint row_vox;
	uint col_vox;
	uint len_vox;
	uint tot_vox;
	float step;

	float *scalar;
	float3 *normal;
	float3 *pos;
	float3 origin;
	float3 ratio;

	float isovalue;

public:

	float3 *t1;
	float3 *t2;
	float3 *t3;

	float3 *n1;
	float3 *n2;
	float3 *n3;

	uint num_triangle;
	uint max_triangle;

	MarchingCube(uint _row_vox, uint _col_vox, uint _len_vox, float *_scalar, float3 *_pos, float3 _origin, float3 _ratio, float _step, float _isovalue);
	~MarchingCube();
	void run();
};

#endif