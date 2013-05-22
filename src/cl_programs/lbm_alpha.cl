#include "src/cl_programs/lbm_header.h"

#define GRAVITATION	0


/**
 * COLLISION AND PROPAGATION KERNEL - ALPHA type
 *
 * \param dd	density distributions
 * \param flags	flags of cells
 */
__kernel void lbm_kernel_alpha(
			__global T *global_dd,		// density distributions
			__global const int *flag_array,	// flags
			__global T *velocity,		// velocities
			__global T *density,		// densities
			__const T inv_tau,
			__const T gravitation_x,
			__const T gravitation_y,
			__const T gravitation_z,
			__const T drivenCavityVelocity			// velocity parameters for modification of density distributions
		)
{
	const size_t gid = get_global_id(0);

	if (gid >= GLOBAL_WORK_GROUP_SIZE)
		return;

	// load cell type flag
	int flag = flag_array[gid];
	if 	(flag == FLAG_GHOST_LAYER)
		return;
	/**
	 * we use a pointer instead of accessing the array directly
	 * first this reduces the number of use registers (according to profiling information)
	 * secondly the program runs faster and we can use more threads
	 */

	// pointer to density distributions
	__global T *current_dds = &global_dd[gid];

	// velocity
	T velocity_x, velocity_y, velocity_z;

	// density distributions
	T dd0, dd1, dd2, dd3, dd4, dd5, dd6, dd7, dd8, dd9, dd10, dd11, dd12, dd13, dd14, dd15, dd16, dd17, dd18;

	// density
	T rho;

	/*
	 * we have to sum the densities up in a specific order.
	 * otherwise it seems that we run into numerical errors.
	 */

	// +++++++++++
	// +++ DD0 +++
	// +++++++++++
	//
	// 0-3: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
	dd0 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho = dd0;
	velocity_x = dd0;
	dd1 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd1;
	velocity_x -= dd1;

	dd2 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd2;
	velocity_y = dd2;
	dd3 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd3;
	velocity_y -= dd3;

	// +++++++++++
	// +++ DD1 +++
	// +++++++++++
	//
	// 4-7: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
	dd4 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd4;
	velocity_x += dd4;
	velocity_y += dd4;
	dd5 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd5;
	velocity_x -= dd5;
	velocity_y -= dd5;

	dd6 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd6;
	velocity_x += dd6;
	velocity_y -= dd6;
	dd7 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd7;
	velocity_x -= dd7;
	velocity_y += dd7;

	// +++++++++++
	// +++ DD2 +++
	// +++++++++++
	//
	// 8-11: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
	dd8 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd8;
	velocity_x += dd8;
	velocity_z = dd8;
	dd9 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd9;
	velocity_x -= dd9;
	velocity_z -= dd9;

	dd10 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd10;
	velocity_x += dd10;
	velocity_z -= dd10;
	dd11 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd11;
	velocity_x -= dd11;
	velocity_z += dd11;

	// +++++++++++
	// +++ DD3 +++
	// +++++++++++

	// dd3: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
	dd12 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd12;
	velocity_y += dd12;
	velocity_z += dd12;
	dd13 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd13;
	velocity_y -= dd13;
	velocity_z -= dd13;

	dd14 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd14;
	velocity_y += dd14;
	velocity_z -= dd14;
	dd15 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd15;
	velocity_y -= dd15;
	velocity_z += dd15;

	// +++++++++++
	// +++ DD4 +++
	// +++++++++++
	//
	// dd4: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
	dd16 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd16;
	velocity_z += dd16;
	dd17 = *current_dds;		current_dds += DOMAIN_CELLS;
	rho += dd17;
	velocity_z -= dd17;
	dd18 = *current_dds;
	rho += dd18;

//	rho *= (float)(flag != FLAG_OBSTACLE);


	// to add something to a pointer is faster than subtracting it.
	// thus we restart here with pointer to dd0
	current_dds = &global_dd[gid];

	/**
	 * instead of storing the density distributions after modification,
	 * we store it during the modifications to hide memory waiting stalls
	 */

	T vel2;		// vel*vel
	T vela2;
	T vela_velb;
#define tmp	rho


	T dd_param;	// modified rho as temporary variable
	switch(flag)
	{
		case FLAG_FLUID:	// this is the whole collision operator
			vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
			dd_param = rho - (T)(3.0f/2.0f)*(vel2);

			tmp = gravitation_x*(T)(1.0f/18.0f)*rho;
			vela2 = velocity_x*velocity_x;
			dd1 += inv_tau*(eq_dd_a1(velocity_x, vela2, dd_param) - dd1);
			dd1 -= tmp;
			*current_dds = dd1;		current_dds += DOMAIN_CELLS;

			dd0 += inv_tau*(eq_dd_a0(velocity_x, vela2, dd_param) - dd0);
			dd0 += tmp;
			*current_dds = dd0;		current_dds += DOMAIN_CELLS;

			tmp = gravitation_y*(T)(-1.0f/18.0f)*rho;
			vela2 = velocity_y*velocity_y;
			dd3 += inv_tau*(eq_dd_a1(velocity_y, vela2, dd_param) - dd3);
			dd3 -= tmp;
			*current_dds = dd3;		current_dds += DOMAIN_CELLS;

			dd2 += inv_tau*(eq_dd_a0(velocity_y, vela2, dd_param) - dd2);
			dd2 += tmp;
			*current_dds = dd2;		current_dds += DOMAIN_CELLS;


#define vela_velb_2	vela2
			/***********************
			 * DD1
			 ***********************/
			vela_velb = velocity_x+velocity_y;
			vela_velb_2 = vela_velb*vela_velb;

			tmp = (gravitation_x - gravitation_y)*(T)(1.0f/36.0f)*rho;

			dd5 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd5);
			dd5 -= tmp;
			*current_dds = dd5;		current_dds += DOMAIN_CELLS;

			dd4 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd4);
			dd4 += tmp;
			*current_dds = dd4;		current_dds += DOMAIN_CELLS;

			vela_velb = velocity_x-velocity_y;
			vela_velb_2 = vela_velb*vela_velb;
			tmp = (gravitation_x + gravitation_y)*(T)(1.0f/36.0f)*rho;
			dd7 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd7);
			dd7 -= tmp;
			*current_dds = dd7;		current_dds += DOMAIN_CELLS;

			dd6 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd6);
			dd6 += tmp;
			*current_dds = dd6;		current_dds += DOMAIN_CELLS;

			/***********************
			 * DD2
			 ***********************/
			vela_velb = velocity_x+velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			tmp = (gravitation_x + gravitation_z)*(T)(1.0f/36.0f)*rho;

			dd9 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd9);
			dd9 -= tmp;
			*current_dds = dd9;		current_dds += DOMAIN_CELLS;

			dd8 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd8);
			dd8 += tmp;
			*current_dds = dd8;		current_dds += DOMAIN_CELLS;

			tmp = (gravitation_x - gravitation_z)*(T)(1.0f/36.0f)*rho;
			vela_velb = velocity_x-velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			dd11 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd11);
			dd11 -= tmp;
			*current_dds = dd11;		current_dds += DOMAIN_CELLS;

			dd10 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd10);
			dd10 += tmp;
			*current_dds = dd10;		current_dds += DOMAIN_CELLS;

			/***********************
			 * DD3
			 ***********************/
			vela_velb = velocity_y+velocity_z;
			vela_velb_2 = vela_velb*vela_velb;

			tmp = (gravitation_z - gravitation_y)*(T)(1.0f/36.0f)*rho;
			dd13 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd13);
			dd13 -= tmp;
			*current_dds = dd13;		current_dds += DOMAIN_CELLS;

			dd12 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd12);
			dd12 += tmp;
			*current_dds = dd12;		current_dds += DOMAIN_CELLS;

			vela_velb = velocity_y-velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			tmp = (gravitation_z + gravitation_y)*(T)(-1.0f/36.0f)*rho;

			dd15 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd15);
			dd15 -= tmp;
			*current_dds = dd15;		current_dds += DOMAIN_CELLS;

			dd14 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd14);
			dd14 += tmp;
			*current_dds = dd14;		current_dds += DOMAIN_CELLS;

#undef vela_velb_2
			/***********************
			 * DD4
			 ***********************/
			vela2 = velocity_z*velocity_z;

			tmp = gravitation_z*(T)(1.0f/18.0f)*rho;
			dd17 += inv_tau*(eq_dd_a1(velocity_z, vela2, dd_param) - dd17);
			dd17 -= tmp;
			*current_dds = dd17;		current_dds += DOMAIN_CELLS;

			dd16 += inv_tau*(eq_dd_a0(velocity_z, vela2, dd_param) - dd16);
			dd16 += tmp;
			*current_dds = dd16;		current_dds += DOMAIN_CELLS;

			dd18 += inv_tau*(eq_dd18(dd_param) - dd18);
			*current_dds = dd18;

			break;

		case FLAG_OBSTACLE:	// in case of an obstacle, we bounce back the values

			/**
			 * if we are using only bounce back method, it's not necessary to write back the results.
			 * we even wouldn't need to read the density distribution values.
			 */
#if 0
			// use simple bounce back
			*current_dds = dd0;		current_dds += DOMAIN_CELLS;
			*current_dds = dd1;		current_dds += DOMAIN_CELLS;
			*current_dds = dd2;		current_dds += DOMAIN_CELLS;
			*current_dds = dd3;		current_dds += DOMAIN_CELLS;

			*current_dds = dd4;		current_dds += DOMAIN_CELLS;
			*current_dds = dd5;		current_dds += DOMAIN_CELLS;
			*current_dds = dd6;		current_dds += DOMAIN_CELLS;
			*current_dds = dd7;		current_dds += DOMAIN_CELLS;

			*current_dds = dd8;		current_dds += DOMAIN_CELLS;
			*current_dds = dd9;		current_dds += DOMAIN_CELLS;
			*current_dds = dd10;		current_dds += DOMAIN_CELLS;
			*current_dds = dd11;		current_dds += DOMAIN_CELLS;

			*current_dds = dd12;		current_dds += DOMAIN_CELLS;
			*current_dds = dd13;		current_dds += DOMAIN_CELLS;
			*current_dds = dd14;		current_dds += DOMAIN_CELLS;
			*current_dds = dd15;		current_dds += DOMAIN_CELLS;

			*current_dds = dd16;		current_dds += DOMAIN_CELLS;
			*current_dds = dd17;		current_dds += DOMAIN_CELLS;
			*current_dds = dd18;
#endif

#if STORE_VELOCITY
			velocity_x = 0.0f;
			velocity_y = 0.0f;
			velocity_z = 0.0f;
#endif
			break;

		case FLAG_VELOCITY_INJECTION:	// this flag specifies the injection area of the fluid
			velocity_x = drivenCavityVelocity;
			velocity_y = 0;
			velocity_z = 0;

			rho = 1.0f;

			vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
			dd_param = rho - (T)(3.0f/2.0f)*(vel2);

			/***********************
			 * DD0
			 ***********************/
			vela2 = velocity_x*velocity_x;
			tmp = gravitation_x*(T)(1.0f/18.0f)*rho;

			dd1 = eq_dd_a1(velocity_x, vela2, dd_param);
			dd1 -= tmp;
			*current_dds = dd1;		current_dds += DOMAIN_CELLS;

			dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
			dd0 += tmp;
			*current_dds = dd0;		current_dds += DOMAIN_CELLS;

			vela2 = velocity_y*velocity_y;
			tmp = gravitation_y*(T)(-1.0f/18.0f)*rho;

			dd3 = eq_dd_a1(velocity_y, vela2, dd_param);
			dd3 -= tmp;
			*current_dds = dd3;		current_dds += DOMAIN_CELLS;

			dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
			dd2 += tmp;
			*current_dds = dd2;		current_dds += DOMAIN_CELLS;


#define vela_velb_2	vela2
			/***********************
			 * DD1
			 ***********************/
			vela_velb = velocity_x+velocity_y;
			vela_velb_2 = vela_velb*vela_velb;
			tmp = (gravitation_x - gravitation_y)*(T)(1.0f/36.0f)*rho;

			dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);
			dd5 -= tmp;
			*current_dds = dd5;		current_dds += DOMAIN_CELLS;

			dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd4 += tmp;
			*current_dds = dd4;		current_dds += DOMAIN_CELLS;

			vela_velb = velocity_x-velocity_y;
			vela_velb_2 = vela_velb*vela_velb;
			tmp = (gravitation_x + gravitation_y)*(T)(1.0f/36.0f)*rho;

			dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);
			dd7 -= tmp;
			*current_dds = dd7;		current_dds += DOMAIN_CELLS;

			dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd6 += tmp;
			*current_dds = dd6;		current_dds += DOMAIN_CELLS;

			/***********************
			 * DD2
			 ***********************/
			vela_velb = velocity_x+velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			tmp = (gravitation_x + gravitation_z)*(T)(1.0f/36.0f)*rho;

			dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);
			dd9 -= tmp;
			*current_dds = dd9;		current_dds += DOMAIN_CELLS;

			dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd8 += tmp;
			*current_dds = dd8;		current_dds += DOMAIN_CELLS;

			vela_velb = velocity_x-velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			tmp = (gravitation_x - gravitation_z)*(T)(1.0f/36.0f)*rho;

			dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);
			dd11 -= tmp;
			*current_dds = dd11;		current_dds += DOMAIN_CELLS;

			dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd10 += tmp;
			*current_dds = dd10;		current_dds += DOMAIN_CELLS;

			/***********************
			 * DD3
			 ***********************/
			vela_velb = velocity_y+velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			tmp = (gravitation_z - gravitation_y)*(T)(1.0f/36.0f)*rho;

			dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);
			dd13 -= tmp;
			*current_dds = dd13;		current_dds += DOMAIN_CELLS;

			dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd12 += tmp;
			*current_dds = dd12;		current_dds += DOMAIN_CELLS;

			vela_velb = velocity_y-velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			tmp = (gravitation_z + gravitation_y)*(T)(-1.0f/36.0f)*rho;

			dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);
			dd15 -= tmp;
			*current_dds = dd15;		current_dds += DOMAIN_CELLS;

			dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd14 += tmp;
			*current_dds = dd14;		current_dds += DOMAIN_CELLS;
#undef vela_velb_2

			/***********************
			 * DD4
			 ***********************/
			vela2 = velocity_z*velocity_z;

			tmp = gravitation_z*(T)(1.0f/18.0f)*rho;
			dd17 = eq_dd_a1(velocity_z, vela2, dd_param);
			dd17 -= tmp;
			*current_dds = dd17;		current_dds += DOMAIN_CELLS;

			dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
			dd16 += tmp;
			*current_dds = dd16;		current_dds += DOMAIN_CELLS;

			dd18 = eq_dd18(dd_param);
			*current_dds = dd18;
			break;

		case (FLAG_GHOST_LAYER):
				break;
	}

#if STORE_VELOCITY
	// store velocity
	current_dds = &velocity[gid];
	*current_dds = velocity_x;	current_dds += DOMAIN_CELLS;
	*current_dds = velocity_y;	current_dds += DOMAIN_CELLS;
	*current_dds = velocity_z;
#endif

#if STORE_DENSITY
	// store density (not necessary)
	density[gid] = rho;
#endif
}
