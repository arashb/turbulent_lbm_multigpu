#include "src/cl_programs/lbm_header.h"


/**
 * return 3D position
 * \input linear_position	linear position in 3D cube (ordering: X,Y,Z)
 * \return	Vector with 3D position in cube
 */
inline T4 getCubePosition(int linear_position)
{
	T4 pos;

	// TODO: use AND operation to speed up
	// (but this function is only used during initialization)
	pos.x = (T)((int)linear_position % (int)DOMAIN_CELLS_X);
	linear_position /= DOMAIN_CELLS_X;

	pos.y = (T)((int)linear_position % (int)DOMAIN_CELLS_Y);
	linear_position /= DOMAIN_CELLS_Y;

	pos.z = linear_position;// % CUBE_SIZE_Z;
	return pos;
}


/**
 * INIT KERNEL
 *
 * \param dd	density distributions
 * \param flags	flags of cells
 */
__kernel void init_kernel(
			__global T *global_dd,	// density distributions
			__global int *flags,	// flags
			__global T *velocity_array,	// velocity array (first all x components, then all y components, then z...)
			__global T *density,	// densities
			__global int *bc, 		///< boundary conditions
			T drivenCavityVelocity			// velocity parameters for modification of density distributions
		)
{
	const size_t gid = get_global_id(0);

	__global T *current_dds = &global_dd[gid];

	// initialize flag field
	T4 pos = getCubePosition(gid);
	pos.w = 0;

	T velocity_x = 0;
	T velocity_y = 0;
	T velocity_z = 0;

	int flag = FLAG_FLUID;

	if( pos.x == 0)
		flag = bc[0];
	else if( pos.x == DOMAIN_CELLS_X-1 )
		flag = bc[1];

	else if( pos.y == 0)
		flag = bc[2];
	else if( pos.y == DOMAIN_CELLS_Y-1 )
		flag = bc[3];

	else if( pos.z == 0)
		flag = bc[4];
	else if( pos.z == DOMAIN_CELLS_Z-1 )
		flag = bc[5];

//	else if (pos.y == DOMAIN_CELLS_Y-2)
//		flag = FLAG_VELOCITY_INJECTION;

#if 0
	if (	pos.x == 0 || pos.y == 0 || pos.z == 0 ||
		pos.x == DOMAIN_CELLS_X-1 || pos.y == DOMAIN_CELLS_Y-1 || pos.z == DOMAIN_CELLS_Z-1
	)
	{
		flag = FLAG_OBSTACLE;
	}
	else
	{
#if 1
		if (pos.y == DOMAIN_CELLS_Y-2)
			flag = FLAG_VELOCITY_INJECTION;
#endif
#if 0
		if (pos.y == 10)
			flag = FLAG_OBSTACLE;

		if (pos.y == 2)
			velocity_x = 10;
		if (pos.y == 3)
			velocity_x = 10;
		if (pos.y == 4)
			velocity_x = 10;
		if (pos.y == 5)
			velocity_x = 10;
		if (pos.y == 6)
			velocity_x = 10;
#endif
#if 0
		if ((pos.x == DOMAIN_CELLS_X/2 || pos.x == DOMAIN_CELLS_X-2) && pos.y <= DOMAIN_CELLS_Y/2)
		{
			flag = FLAG_INTERFACE;
		}
		else if ((pos.y == DOMAIN_CELLS_Y/2 || pos.y == 1) && pos.x >= DOMAIN_CELLS_X/2)
		{
			flag = FLAG_INTERFACE;
		}
		else if ((pos.z == DOMAIN_CELLS_Z-1 || pos.z == 1) && pos.x >= DOMAIN_CELLS_X/2 && pos.y <= DOMAIN_CELLS_Y/2)
		{
			flag = FLAG_INTERFACE;
		}
		else if (pos.x < DOMAIN_CELLS_X/2 || pos.y > DOMAIN_CELLS_Y/2)
		{
			flag = FLAG_GAS;
		}
#endif
	}
#endif

	// density distributions
	T dd0, dd1, dd2, dd3, dd4, dd5, dd6, dd7, dd8, dd9, dd10, dd11, dd12, dd13, dd14, dd15, dd16, dd17, dd18;

	T dd_param;
	T vela2;
	T vela_velb;
	T rho = 1.0f;

	// compute and store velocity

	vela2 = velocity_x*velocity_x;
	dd_param = rho - (T)(3.0f/2.0f)*(vela2);

	dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
	*current_dds = dd0;		current_dds += DOMAIN_CELLS;
	dd1 = eq_dd_a1(velocity_x, vela2, dd_param);
	*current_dds = dd1;		current_dds += DOMAIN_CELLS;

	vela2 = velocity_y*velocity_y;

	dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
	*current_dds = dd2;		current_dds += DOMAIN_CELLS;
	dd3 = eq_dd_a1(velocity_y, vela2, dd_param);
	*current_dds = dd3;		current_dds += DOMAIN_CELLS;


#define vela_velb_2	vela2
	/***********************
	 * DD1
	 ***********************/
	vela_velb = velocity_x+velocity_y;
	vela_velb_2 = vela_velb*vela_velb;

	dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd4;		current_dds += DOMAIN_CELLS;
	dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd5;		current_dds += DOMAIN_CELLS;

	vela_velb = velocity_x-velocity_y;
	vela_velb_2 = vela_velb*vela_velb;

	dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd6;		current_dds += DOMAIN_CELLS;
	dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd7;		current_dds += DOMAIN_CELLS;

	/***********************
	 * DD2
	 ***********************/
	vela_velb = velocity_x+velocity_z;
	vela_velb_2 = vela_velb*vela_velb;

	dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd8;		current_dds += DOMAIN_CELLS;
	dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd9;		current_dds += DOMAIN_CELLS;

	vela_velb = velocity_x-velocity_z;
	vela_velb_2 = vela_velb*vela_velb;

	dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd10;		current_dds += DOMAIN_CELLS;
	dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd11;		current_dds += DOMAIN_CELLS;

	/***********************
	 * DD3
	 ***********************/
	vela_velb = velocity_y+velocity_z;
	vela_velb_2 = vela_velb*vela_velb;


	dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd12;		current_dds += DOMAIN_CELLS;
	dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd13;		current_dds += DOMAIN_CELLS;

	vela_velb = velocity_y-velocity_z;
	vela_velb_2 = vela_velb*vela_velb;

	dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd14;		current_dds += DOMAIN_CELLS;
	dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);
	*current_dds = dd15;		current_dds += DOMAIN_CELLS;


#undef vela_velb_2
	/***********************
	 * DD4
	 ***********************/
	vela2 = velocity_z*velocity_z;

	dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
	*current_dds = dd16;		current_dds += DOMAIN_CELLS;
	dd17 = eq_dd_a1(velocity_z, vela2, dd_param);
	*current_dds = dd17;		current_dds += DOMAIN_CELLS;

	dd18 = eq_dd18(dd_param);
	*current_dds = dd18;

	// flag
	flags[gid] = flag;

#if STORE_VELOCITY
	// store velocity
	current_dds = &velocity_array[gid];
	*current_dds = velocity_x;	current_dds += DOMAIN_CELLS;
	*current_dds = velocity_y;	current_dds += DOMAIN_CELLS;
	*current_dds = velocity_z;
#endif

#if STORE_DENSITY
	// store density
	density[gid] = rho;
#endif
}
