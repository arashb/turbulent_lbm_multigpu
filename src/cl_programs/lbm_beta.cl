
#define GRAVITATION			0

#define CACHED_ACCESS		0

#define USE_SHARED_MEMORY	1

#include "src/cl_programs/lbm_header.h"

/**
 * COLLISION AND PROPAGATION KERNEL - BETA type
 */
__kernel void lbm_kernel_beta(
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
	const int flag = flag_array[gid];

	/**
	 * we use a pointer instead of accessing the array directly
	 * first this reduces the number of use registers (according to profiling information)
	 * secondly the program runs faster and we can use more threads
	 */
	__global T *current_dds = global_dd;

	// velocity
	T velocity_x, velocity_y, velocity_z;

	// density distributions
	T dd0, dd1, dd2, dd3, dd4, dd5, dd6, dd7, dd8, dd9, dd10, dd11, dd12, dd13, dd14, dd15, dd16, dd17, dd18;

	// density
	T rho;

#if !USE_SHARED_MEMORY
	/*
	 * dd 0-3: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
	 */
	dd1 = current_dds[DOMAIN_WRAP(gid + DELTA_POS_X)];	current_dds += DOMAIN_CELLS;
	dd0 = current_dds[DOMAIN_WRAP(gid + DELTA_NEG_X)];	current_dds += DOMAIN_CELLS;

	rho = dd0;
	velocity_x = dd0;

	/*
	 * we have to sum the densities up in a specific order.
	 * otherwise it seems that we run into numerical errors for fluids with zero velocity.
	 */
	rho += dd1;
	velocity_x -= dd1;

	dd3 = current_dds[DOMAIN_WRAP(gid + DELTA_POS_Y)];		current_dds += DOMAIN_CELLS;
	dd2 = current_dds[DOMAIN_WRAP(gid + DELTA_NEG_Y)];		current_dds += DOMAIN_CELLS;

	rho += dd2;
	velocity_y = dd2;

	rho += dd3;
	velocity_y -= dd3;

	/*
	 * dd 4-7: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
	 */
	dd5 = current_dds[DOMAIN_WRAP(gid + (DELTA_POS_X + DELTA_POS_Y))];		current_dds += DOMAIN_CELLS;
	dd4 = current_dds[DOMAIN_WRAP(gid + (DELTA_NEG_X + DELTA_NEG_Y))];		current_dds += DOMAIN_CELLS;

	rho += dd4;
	velocity_x += dd4;
	velocity_y += dd4;

	rho += dd5;
	velocity_x -= dd5;
	velocity_y -= dd5;

	dd7 = current_dds[DOMAIN_WRAP(gid + (DELTA_POS_X + DELTA_NEG_Y))];		current_dds += DOMAIN_CELLS;
	dd6 = current_dds[DOMAIN_WRAP(gid + (DELTA_NEG_X + DELTA_POS_Y))];		current_dds += DOMAIN_CELLS;


	rho += dd6;
	velocity_x += dd6;
	velocity_y -= dd6;

	rho += dd7;
	velocity_x -= dd7;
	velocity_y += dd7;

	/*
	 * dd 8-11: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
	 */
	dd9 = current_dds[DOMAIN_WRAP(gid + (DELTA_POS_X + DELTA_POS_Z))];		current_dds += DOMAIN_CELLS;
	dd8 = current_dds[DOMAIN_WRAP(gid + (DELTA_NEG_X + DELTA_NEG_Z))];		current_dds += DOMAIN_CELLS;

	rho += dd8;
	velocity_x += dd8;
	velocity_z = dd8;

	rho += dd9;
	velocity_x -= dd9;
	velocity_z -= dd9;

	dd11 = current_dds[DOMAIN_WRAP(gid + (DELTA_POS_X + DELTA_NEG_Z))];		current_dds += DOMAIN_CELLS;
	dd10 = current_dds[DOMAIN_WRAP(gid + (DELTA_NEG_X + DELTA_POS_Z))];		current_dds += DOMAIN_CELLS;

	rho += dd10;
	velocity_x += dd10;
	velocity_z -= dd10;

	rho += dd11;
	velocity_x -= dd11;
	velocity_z += dd11;

	/*
	 * dd 12-15: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
	 */
	dd13 = current_dds[DOMAIN_WRAP(gid + (DELTA_POS_Y + DELTA_POS_Z))];	current_dds += DOMAIN_CELLS;
	dd12 = current_dds[DOMAIN_WRAP(gid + (DELTA_NEG_Y + DELTA_NEG_Z))];	current_dds += DOMAIN_CELLS;

	rho += dd12;
	velocity_y += dd12;
	velocity_z += dd12;

	rho += dd13;
	velocity_y -= dd13;
	velocity_z -= dd13;

	dd15 = current_dds[DOMAIN_WRAP(gid + (DELTA_POS_Y + DELTA_NEG_Z))];	current_dds += DOMAIN_CELLS;
	dd14 = current_dds[DOMAIN_WRAP(gid + (DELTA_NEG_Y + DELTA_POS_Z))];	current_dds += DOMAIN_CELLS;

	rho += dd14;
	velocity_y += dd14;
	velocity_z -= dd14;

	rho += dd15;
	velocity_y -= dd15;
	velocity_z += dd15;

	/*
	 * dd 16-18: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
	 */
	dd17 = current_dds[DOMAIN_WRAP(gid + DELTA_POS_Z)];		current_dds += DOMAIN_CELLS;
	dd16 = current_dds[DOMAIN_WRAP(gid + DELTA_NEG_Z)];		current_dds += DOMAIN_CELLS;

	rho += dd16;
	velocity_z += dd16;

	rho += dd17;
	velocity_z -= dd17;

	dd18 = current_dds[gid];
	rho += dd18;

#else
	__local T dd_buf[12][LOCAL_WORK_GROUP_SIZE];

	const size_t lid = get_local_id(0);

#if CACHED_ACCESS
	size_t dd_read_delta_position_2;
	size_t dd_read_delta_position_3;
	size_t dd_read_delta_position_4;
	size_t dd_read_delta_position_5;
	size_t dd_read_delta_position_6;
	size_t dd_read_delta_position_7;
	size_t dd_read_delta_position_8;
	size_t dd_read_delta_position_9;
	size_t dd_read_delta_position_10;
	size_t dd_read_delta_position_11;
	size_t dd_read_delta_position_12;
	size_t dd_read_delta_position_13;
	size_t dd_read_delta_position_14;
	size_t dd_read_delta_position_15;
	size_t dd_read_delta_position_16;
	size_t dd_read_delta_position_17;
	size_t dd_read_delta_position_18;
#endif

	/*
	 * We have to handle "misaligned" data with a shift of -1 and +1:
	 *
	 * As an example, we handle the access to the density distributions with a shift of -1:
	 *
	 * We use "(( (lid + 1) mod LOCAL_WORK_GROUP_SIZE) + (DOMAIN_CELLS+1)) mod DOMAIN_CELLS" as the reading index
	 *          shift to right                              shift back
	 *
	 * This allows us to read almost everything (except the first thread) aligned. After that the data is stored
	 * to shared memory, a sync operation is called and finally the tread can read the originally required data which
	 * was previously read by another thread.
	 *
	 * Every float is stored to a local memory array indexed by the thread id.
	 * after the dd's are read from global memory, the local memory is accessed with
	 * "(lid + (LOCAL_WORK_GROUP_SIZE-1)) mod LOCAL_WORK_GROUP_SIZE"
	 */

	/*
	 * pos_x_wrap specifies the position in the local buffer for the dd with the displacement x=+1
	 * this is used to force the thread with the largest number to read the dd at the displacement -1 (see gid_pos below)
	 * pos_x_wrap and gid_pos have the following values for given local thread ids:
	 *  
	 * thread_id:   0 1 2 3 4 5 ... 63
	 * pos_x_wrap:  1 2 3 4 5 6 ... 0
	 * gid_pos:     0 1 2 3 4 5 ... -1   <<< !!!
	 */

	int pos_x_wrap = LOCAL_WORK_GROUP_WRAP(lid + 1);
	int neg_x_wrap = LOCAL_WORK_GROUP_WRAP(lid + (LOCAL_WORK_GROUP_SIZE - 1));

#if (LOCAL_WORK_GROUP_SIZE/DOMAIN_CELLS_X)*DOMAIN_CELLS_X == LOCAL_WORK_GROUP_SIZE
	/*
	 * handle domain x-sizes specially if LOCAL_WORK_GROUP_SIZE is a multiple of DOMAIN_CELLS_X
	 * in this case, we dont have to read unaligned data!!!
	 */
	int read_delta_neg_x = gid;
	int read_delta_pos_x = gid;
#else
	/*
	 * cache variables for speedup
	 */
	int read_delta_neg_x = DOMAIN_WRAP(gid - lid + pos_x_wrap + DELTA_NEG_X);
	int read_delta_pos_x = DOMAIN_WRAP(gid - lid + neg_x_wrap + DELTA_POS_X);
#endif

	/*
	 * +++++++++++
	 * +++ DD0 +++
	 * +++++++++++
	 *
	 * dd0: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
	 * negative displacement
	 * preload to alignment buffer
	 */

	/*
	 * read negative distribution vector (-1,0,0) from relative cell (-1,0,0) and store it to positive distribution vector (1,0,0)
	 *
	 * in the alpha kernel, the density distribution values have been stored to the oppisite density distribution storage
	 * to avoid the propagation step
	 */

	/*
	 * pointer to current dd buf entry with index lid
	 */
	__local T *dd_buf_lid = &dd_buf[1][lid];

#define dd_read_delta_position_0	read_delta_neg_x
#define dd_read_delta_position_1	read_delta_pos_x

	*dd_buf_lid = current_dds[dd_read_delta_position_1];		current_dds += DOMAIN_CELLS;	dd_buf_lid -= LOCAL_WORK_GROUP_SIZE;
	*dd_buf_lid = current_dds[dd_read_delta_position_0];		current_dds += DOMAIN_CELLS;	dd_buf_lid += 5*LOCAL_WORK_GROUP_SIZE;
#if CACHED_ACCESS
	dd_read_delta_position_3 = DOMAIN_WRAP(gid + DELTA_POS_Y);
#else
	#define dd_read_delta_position_3	DOMAIN_WRAP(gid + DELTA_POS_Y)
#endif
	dd3 = current_dds[dd_read_delta_position_3];		current_dds += DOMAIN_CELLS;
	rho = dd3;
	velocity_y = -dd3;

#if CACHED_ACCESS
	dd_read_delta_position_2 = DOMAIN_WRAP(gid + DELTA_NEG_Y);
#else
	#define dd_read_delta_position_2	DOMAIN_WRAP(gid + DELTA_NEG_Y)
#endif
	dd2 = current_dds[dd_read_delta_position_2];		current_dds += DOMAIN_CELLS;
	rho += dd2;
	velocity_y += dd2;

	// DD0 STUFF
	barrier(CLK_LOCAL_MEM_FENCE);

	dd0 = dd_buf[0][neg_x_wrap];
	rho += dd0;
	velocity_x = dd0;

	dd1 = dd_buf[1][pos_x_wrap];
	rho += dd1;
	velocity_x -= dd1;


	/* +++++++++++
	 * +++ DD1 +++
	 * +++++++++++
	 *
	 * dd1: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
	 */
#if CACHED_ACCESS
	dd_read_delta_position_5 = DOMAIN_WRAP(read_delta_pos_x + DELTA_POS_Y);
#else
	#define dd_read_delta_position_5	DOMAIN_WRAP(read_delta_pos_x + DELTA_POS_Y)
#endif
	*dd_buf_lid = current_dds[dd_read_delta_position_5];	current_dds += DOMAIN_CELLS;	dd_buf_lid -= LOCAL_WORK_GROUP_SIZE;

#if CACHED_ACCESS
	dd_read_delta_position_4 = DOMAIN_WRAP(read_delta_neg_x + DELTA_NEG_Y);
#else
	#define dd_read_delta_position_4	DOMAIN_WRAP(read_delta_neg_x + DELTA_NEG_Y)
#endif
	*dd_buf_lid = current_dds[dd_read_delta_position_4];	current_dds += DOMAIN_CELLS;	dd_buf_lid += 3*LOCAL_WORK_GROUP_SIZE;

#if CACHED_ACCESS
	dd_read_delta_position_7 = DOMAIN_WRAP(read_delta_pos_x + DELTA_NEG_Y);
#else
	#define dd_read_delta_position_7	DOMAIN_WRAP(read_delta_pos_x + DELTA_NEG_Y)
#endif
	*dd_buf_lid = current_dds[dd_read_delta_position_7];	current_dds += DOMAIN_CELLS;	dd_buf_lid -= LOCAL_WORK_GROUP_SIZE;

#if CACHED_ACCESS
	dd_read_delta_position_6 = DOMAIN_WRAP(read_delta_neg_x + DELTA_POS_Y);
#else
	#define dd_read_delta_position_6	DOMAIN_WRAP(read_delta_neg_x + DELTA_POS_Y)
#endif
	*dd_buf_lid = current_dds[dd_read_delta_position_6];	current_dds += DOMAIN_CELLS;	dd_buf_lid += 3*LOCAL_WORK_GROUP_SIZE;

	barrier(CLK_LOCAL_MEM_FENCE);
	dd4 = dd_buf[4][neg_x_wrap];
	rho += dd4;
	velocity_x += dd4;
	velocity_y += dd4;

	dd5 = dd_buf[5][pos_x_wrap];
	rho += dd5;
	velocity_x -= dd5;
	velocity_y -= dd5;


	dd6 = dd_buf[6][neg_x_wrap];
	rho += dd6;
	velocity_x += dd6;
	velocity_y -= dd6;

	dd7 = dd_buf[7][pos_x_wrap];
	rho += dd7;
	velocity_x -= dd7;
	velocity_y += dd7;


	/* +++++++++++
	 * +++ DD2 +++
	 * +++++++++++
	 *
	 * dd2: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
	 */

#if CACHED_ACCESS
	dd_read_delta_position_9 = DOMAIN_WRAP(read_delta_pos_x + DELTA_POS_Z);
#else
	#define dd_read_delta_position_9	DOMAIN_WRAP(read_delta_pos_x + DELTA_POS_Z)
#endif
	*dd_buf_lid = current_dds[dd_read_delta_position_9];	current_dds += DOMAIN_CELLS;	dd_buf_lid -= LOCAL_WORK_GROUP_SIZE;

#if CACHED_ACCESS
	dd_read_delta_position_8 = DOMAIN_WRAP(read_delta_neg_x + DELTA_NEG_Z);
#else
	#define	dd_read_delta_position_8	DOMAIN_WRAP(read_delta_neg_x + DELTA_NEG_Z)
#endif
	*dd_buf_lid = current_dds[dd_read_delta_position_8];	current_dds += DOMAIN_CELLS;	dd_buf_lid += 3*LOCAL_WORK_GROUP_SIZE;

#if CACHED_ACCESS
	dd_read_delta_position_11 = DOMAIN_WRAP(read_delta_pos_x + DELTA_NEG_Z);
#else
	#define dd_read_delta_position_11	DOMAIN_WRAP(read_delta_pos_x + DELTA_NEG_Z)
#endif
	*dd_buf_lid = current_dds[dd_read_delta_position_11];	current_dds += DOMAIN_CELLS;	dd_buf_lid -= LOCAL_WORK_GROUP_SIZE;

#if CACHED_ACCESS
	dd_read_delta_position_10 = DOMAIN_WRAP(read_delta_neg_x + DELTA_POS_Z);
#else
	#define dd_read_delta_position_10	DOMAIN_WRAP(read_delta_neg_x + DELTA_POS_Z)
#endif
	*dd_buf_lid = current_dds[dd_read_delta_position_10];	current_dds += DOMAIN_CELLS;


	barrier(CLK_LOCAL_MEM_FENCE);
	dd8 = dd_buf[8][neg_x_wrap];
	rho += dd8;
	velocity_x += dd8;
	velocity_z = dd8;

	dd9 = dd_buf[9][pos_x_wrap];
	rho += dd9;
	velocity_x -= dd9;
	velocity_z -= dd9;

	dd10 = dd_buf[10][neg_x_wrap];
	rho += dd10;
	velocity_x += dd10;
	velocity_z -= dd10;

	dd11 = dd_buf[11][pos_x_wrap];
	rho += dd11;
	velocity_x -= dd11;
	velocity_z += dd11;


	/*
	 * +++++++++++
	 * +++ DD3 +++
	 * +++++++++++
	 *
	 * dd3: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
	 */

#if CACHED_ACCESS
	dd_read_delta_position_13 = DOMAIN_WRAP(gid + DELTA_POS_Y + DELTA_POS_Z);
#else
	#define dd_read_delta_position_13	DOMAIN_WRAP(gid + DELTA_POS_Y + DELTA_POS_Z)
#endif
	dd13 = current_dds[dd_read_delta_position_13];	current_dds += DOMAIN_CELLS;
	rho += dd13;
	velocity_y -= dd13;
	velocity_z -= dd13;

#if CACHED_ACCESS
	dd_read_delta_position_12 = DOMAIN_WRAP(gid + DELTA_NEG_Y + DELTA_NEG_Z);
#else
	#define dd_read_delta_position_12	DOMAIN_WRAP(gid + DELTA_NEG_Y + DELTA_NEG_Z)
#endif
	dd12 = current_dds[dd_read_delta_position_12];	current_dds += DOMAIN_CELLS;
	rho += dd12;
	velocity_y += dd12;
	velocity_z += dd12;

#if CACHED_ACCESS
	dd_read_delta_position_15 = DOMAIN_WRAP(gid + DELTA_POS_Y + DELTA_NEG_Z);
#else
	#define dd_read_delta_position_15	DOMAIN_WRAP(gid + DELTA_POS_Y + DELTA_NEG_Z)
#endif
	dd15 = current_dds[dd_read_delta_position_15];	current_dds += DOMAIN_CELLS;
	rho += dd15;
	velocity_y -= dd15;
	velocity_z += dd15;

#if CACHED_ACCESS
	dd_read_delta_position_14 = DOMAIN_WRAP(gid + DELTA_NEG_Y + DELTA_POS_Z);
#else
	#define dd_read_delta_position_14	DOMAIN_WRAP(gid + DELTA_NEG_Y + DELTA_POS_Z)
#endif
	dd14 = current_dds[dd_read_delta_position_14];	current_dds += DOMAIN_CELLS;
	rho += dd14;
	velocity_y += dd14;
	velocity_z -= dd14;


	/*
	 * +++++++++++
	 * +++ DD4 +++
	 * +++++++++++
	 *
	 * dd4: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
	 */
#if CACHED_ACCESS
	dd_read_delta_position_17 = DOMAIN_WRAP(gid + DELTA_POS_Z);
#else
	#define dd_read_delta_position_17	DOMAIN_WRAP(gid + DELTA_POS_Z)
#endif
	dd17 = current_dds[dd_read_delta_position_17];	current_dds += DOMAIN_CELLS;
	rho += dd17;
	velocity_z -= dd17;

#if CACHED_ACCESS
	dd_read_delta_position_16 = DOMAIN_WRAP(gid + DELTA_NEG_Z);
#else
	#define dd_read_delta_position_16	DOMAIN_WRAP(gid + DELTA_NEG_Z)
#endif
	dd16 = current_dds[dd_read_delta_position_16];	current_dds += DOMAIN_CELLS;
	rho += dd16;
	velocity_z += dd16;

	dd18 = current_dds[gid];
	rho += dd18;

#endif

	T vel2;		// vel*vel
	T vela2;

#define vela_velb	vel2
#define vela_velb_2	vela2

#define dd_param	rho
//	T dd_param;	// modified rho as temporary variable
	switch(flag)
	{
		case FLAG_FLUID:	// this is the whole collision operator
			vel2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
			dd_param = rho - (T)(3.0f/2.0f)*(vel2);

			vela2 = velocity_x*velocity_x;
			dd0 += inv_tau*(eq_dd_a0(velocity_x, vela2, dd_param) - dd0);
			dd1 += inv_tau*(eq_dd_a1(velocity_x, vela2, dd_param) - dd1);

			vela2 = velocity_y*velocity_y;
			dd2 += inv_tau*(eq_dd_a0(velocity_y, vela2, dd_param) - dd2);
			dd3 += inv_tau*(eq_dd_a1(velocity_y, vela2, dd_param) - dd3);


			/***********************
			 * DD1
			 ***********************/
			vela_velb = velocity_x+velocity_y;
			vela_velb_2 = vela_velb*vela_velb;

			dd4 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd4);
			dd5 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd5);

			vela_velb = velocity_x-velocity_y;
			vela_velb_2 = vela_velb*vela_velb;

			dd6 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd6);
			dd7 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd7);

			/***********************
			 * DD2
			 ***********************/
			vela_velb = velocity_x+velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			dd8 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd8);
			dd9 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd9);

			vela_velb = velocity_x-velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			dd10 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd10);
			dd11 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd11);

			/***********************
			 * DD3
			 ***********************/
			vela_velb = velocity_y+velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			dd12 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd12);
			dd13 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd13);

			vela_velb = velocity_y-velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			dd14 += inv_tau*(eq_dd4(vela_velb, vela_velb_2, dd_param) - dd14);
			dd15 += inv_tau*(eq_dd5(vela_velb, vela_velb_2, dd_param) - dd15);

			/***********************
			 * DD4
			 ***********************/
			vela2 = velocity_z*velocity_z;
			dd16 += inv_tau*(eq_dd_a0(velocity_z, vela2, dd_param) - dd16);
			dd17 += inv_tau*(eq_dd_a1(velocity_z, vela2, dd_param) - dd17);

			dd18 += inv_tau*(eq_dd18(dd_param) - dd18);
			break;

		case FLAG_OBSTACLE:	// in case of an obstacle, we bounce back the values
			// set to zero velocity and no fluid density
#if STORE_VELOCITY
			velocity_x = 0.0f;
			velocity_y = 0.0f;
			velocity_z = 0.0f;
#endif

			// use simple bounce back
			vela2 = dd1;	dd1 = dd0;		dd0 = vela2;
			vela2 = dd3;	dd3 = dd2;		dd2 = vela2;
			vela2 = dd5;	dd5 = dd4;		dd4 = vela2;
			vela2 = dd7;	dd7 = dd6;		dd6 = vela2;
			vela2 = dd9;	dd9 = dd8;		dd8 = vela2;
			vela2 = dd11;	dd11 = dd10;	dd10 = vela2;
			vela2 = dd13;	dd13 = dd12;	dd12 = vela2;
			vela2 = dd15;	dd15 = dd14;	dd14 = vela2;
			vela2 = dd17;	dd17 = dd16;	dd16 = vela2;

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
			dd0 = eq_dd_a0(velocity_x, vela2, dd_param);
			dd1 = eq_dd_a1(velocity_x, vela2, dd_param);

			vela2 = velocity_y*velocity_y;
			dd2 = eq_dd_a0(velocity_y, vela2, dd_param);
			dd3 = eq_dd_a1(velocity_y, vela2, dd_param);

			/***********************
			 * DD1
			 ***********************/
			vela_velb = velocity_x+velocity_y;
			vela_velb_2 = vela_velb*vela_velb;

			dd4 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd5 = eq_dd5(vela_velb, vela_velb_2, dd_param);

			vela_velb = velocity_x-velocity_y;
			vela_velb_2 = vela_velb*vela_velb;
			dd6 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd7 = eq_dd5(vela_velb, vela_velb_2, dd_param);

			/***********************
			 * DD2
			 ***********************/
			vela_velb = velocity_x+velocity_z;
			vela_velb_2 = vela_velb*vela_velb;

			dd8 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd9 = eq_dd5(vela_velb, vela_velb_2, dd_param);

			vela_velb = velocity_x-velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			dd10 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd11 = eq_dd5(vela_velb, vela_velb_2, dd_param);

			/***********************
			 * DD3
			 ***********************/
			vela_velb = velocity_y+velocity_z;
			vela_velb_2 = vela_velb*vela_velb;

			dd12 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd13 = eq_dd5(vela_velb, vela_velb_2, dd_param);

			vela_velb = velocity_y-velocity_z;
			vela_velb_2 = vela_velb*vela_velb;
			dd14 = eq_dd4(vela_velb, vela_velb_2, dd_param);
			dd15 = eq_dd5(vela_velb, vela_velb_2, dd_param);

			/***********************
			 * DD4
			 ***********************/
			vela2 = velocity_z*velocity_z;
			dd16 = eq_dd_a0(velocity_z, vela2, dd_param);
			dd17 = eq_dd_a1(velocity_z, vela2, dd_param);

			dd18 = eq_dd18(dd_param);
			break;
		case ( FLAG_GHOST_LAYER):
			break;
	}

#if GRAVITATION
	/*
	velocity_x += 10.0*gravitation_x;
	velocity_y += 10.0*gravitation_y;
	velocity_z += 10.0*gravitation_z;
	*/

#define tmp	vela2
	if (flag != FLAG_OBSTACLE)
	{
		tmp = gravitation_x*(T)(1.0f/18.0f)*rho;
		dd0 += tmp;
		dd1 -= tmp;
		tmp = gravitation_y*(T)(-1.0f/18.0f)*rho;
		dd2 += tmp;
		dd3 -= tmp;

		tmp = (gravitation_x - gravitation_y)*(T)(1.0f/36.0f)*rho;
		dd4 += tmp;
		dd5 -= tmp;
		tmp = (gravitation_x + gravitation_y)*(T)(1.0f/36.0f)*rho;
		dd6 += tmp;
		dd7 -= tmp;

		tmp = (gravitation_x + gravitation_z)*(T)(1.0f/36.0f)*rho;
		dd8 += tmp;
		dd9 -= tmp;
		tmp = (gravitation_x - gravitation_z)*(T)(1.0f/36.0f)*rho;
		dd10 += tmp;
		dd11 -= tmp;

		tmp = (gravitation_z - gravitation_y)*(T)(1.0f/36.0f)*rho;
		dd12 += tmp;
		dd13 -= tmp;
		tmp = (gravitation_z + gravitation_y)*(T)(-1.0f/36.0f)*rho;
		dd14 += tmp;
		dd15 -= tmp;

		tmp = gravitation_z*(T)(1.0f/18.0f)*rho;
		dd16 += tmp;
		dd17 -= tmp;
	}
#undef tmp

#endif
	// gain little speedup
	barrier(CLK_LOCAL_MEM_FENCE);

	current_dds = global_dd;

#if USE_SHARED_MEMORY
	dd_buf_lid = &dd_buf[0][lid];

	/* f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0) */
	dd_buf[0][pos_x_wrap] = dd0;
	dd_buf[1][neg_x_wrap] = dd1;
	barrier(CLK_LOCAL_MEM_FENCE);

	current_dds[dd_read_delta_position_1] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;	dd_buf_lid += LOCAL_WORK_GROUP_SIZE;
	current_dds[dd_read_delta_position_0] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;	dd_buf_lid += 3*LOCAL_WORK_GROUP_SIZE;
	current_dds[dd_read_delta_position_3] = dd2;		current_dds += DOMAIN_CELLS;
	current_dds[dd_read_delta_position_2] = dd3;		current_dds += DOMAIN_CELLS;

	/* f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0) */
	dd_buf[4][pos_x_wrap] = dd4;
	dd_buf[5][neg_x_wrap] = dd5;
	dd_buf[6][pos_x_wrap] = dd6;
	dd_buf[7][neg_x_wrap] = dd7;
	barrier(CLK_LOCAL_MEM_FENCE);

	current_dds[dd_read_delta_position_5] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;	dd_buf_lid += LOCAL_WORK_GROUP_SIZE;
	current_dds[dd_read_delta_position_4] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;	dd_buf_lid += LOCAL_WORK_GROUP_SIZE;
	current_dds[dd_read_delta_position_7] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;	dd_buf_lid += LOCAL_WORK_GROUP_SIZE;
	current_dds[dd_read_delta_position_6] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;	dd_buf_lid += LOCAL_WORK_GROUP_SIZE;

	/* f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1) */
	dd_buf[8][pos_x_wrap] = dd8;
	dd_buf[9][neg_x_wrap] = dd9;
	dd_buf[10][pos_x_wrap] = dd10;
	dd_buf[11][neg_x_wrap] = dd11;
	barrier(CLK_LOCAL_MEM_FENCE);

	current_dds[dd_read_delta_position_9] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;	dd_buf_lid += LOCAL_WORK_GROUP_SIZE;
	current_dds[dd_read_delta_position_8] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;	dd_buf_lid += LOCAL_WORK_GROUP_SIZE;
	current_dds[dd_read_delta_position_11] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;	dd_buf_lid += LOCAL_WORK_GROUP_SIZE;
	current_dds[dd_read_delta_position_10] = *dd_buf_lid;	current_dds += DOMAIN_CELLS;

	/* f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1) */
	current_dds[dd_read_delta_position_13] = dd12;	current_dds += DOMAIN_CELLS;
	current_dds[dd_read_delta_position_12] = dd13;	current_dds += DOMAIN_CELLS;
	current_dds[dd_read_delta_position_15] = dd14;	current_dds += DOMAIN_CELLS;
	current_dds[dd_read_delta_position_14] = dd15;	current_dds += DOMAIN_CELLS;

	/* f(0,0,1), f(0,0,-1),  f(0,0,0) */
	current_dds[dd_read_delta_position_17] = dd16;	current_dds += DOMAIN_CELLS;
	current_dds[dd_read_delta_position_16] = dd17;	current_dds += DOMAIN_CELLS;
	current_dds[gid] = dd18;
#else

	/* f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0) */
	current_dds[DOMAIN_WRAP(gid + DELTA_POS_X)] = dd0;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_NEG_X)] = dd1;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_POS_Y)] = dd2;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_NEG_Y)] = dd3;	current_dds += DOMAIN_CELLS;

	/* f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0) */
	current_dds[DOMAIN_WRAP(gid + DELTA_POS_X + DELTA_POS_Y)] = dd4;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_NEG_X + DELTA_NEG_Y)] = dd5;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_POS_X + DELTA_NEG_Y)] = dd6;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_NEG_X + DELTA_POS_Y)] = dd7;	current_dds += DOMAIN_CELLS;

	/* f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1) */
	current_dds[DOMAIN_WRAP(gid + DELTA_POS_X + DELTA_POS_Z)] = dd8;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_NEG_X + DELTA_NEG_Z)] = dd9;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_POS_X + DELTA_NEG_Z)] = dd10;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_NEG_X + DELTA_POS_Z)] = dd11;	current_dds += DOMAIN_CELLS;

	/* f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1) */
	current_dds[DOMAIN_WRAP(gid + DELTA_POS_Y + DELTA_POS_Z)] = dd12;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_NEG_Y + DELTA_NEG_Z)] = dd13;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_POS_Y + DELTA_NEG_Z)] = dd14;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_NEG_Y + DELTA_POS_Z)] = dd15;	current_dds += DOMAIN_CELLS;

	/* f(0,0,1), f(0,0,-1),  f(0,0,0) */
	current_dds[DOMAIN_WRAP(gid + DELTA_POS_Z)] = dd16;	current_dds += DOMAIN_CELLS;
	current_dds[DOMAIN_WRAP(gid + DELTA_NEG_Z)] = dd17;	current_dds += DOMAIN_CELLS;
	current_dds[gid] = dd18;
#endif

	if ( flag == FLAG_GHOST_LAYER)
		return;

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
