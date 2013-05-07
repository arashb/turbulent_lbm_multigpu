/**
 * <<< prefix created from cpu program >>>
 *
 * the cpu program defines:
 *  - DOMAIN_CELLS
 *  - LOCAL_WORK_GROUP_SIZE
 *  - DOMAIN_CELLS_X
 *  - DOMAIN_CELLS_Y
 *  - DOMAIN_CELLS_Z
 *
 *
 * velocity components for density distributions:
 *
 *                  <-- X -->
 * dd0: f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
 * dd1: f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
 * dd2: f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
 * dd3: f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
 * dd4: f(0,0,1), f(0,0,-1),  f(0,0,0),  (not used)
 *
 *
 * MEMORY_LAYOUT:
 * dd0[0],dd0[1],dd0[2],dd0[3],...
 * dd1[0],dd1[1],dd1[2],dd1[3],...
 * ...
 * thus, each velocity components can be accessed linearly
 *
 */

#include "src/cl_programs/lbm_defaults.h"

#define USE_SWIZZLING	0

#define DOMAIN_CELLS		(DOMAIN_CELLS_X*DOMAIN_CELLS_Y*DOMAIN_CELLS_Z)
#define DOMAIN_SLICE_CELLS	(DOMAIN_CELLS_X*DOMAIN_CELLS_Y)

/**
 * Next we define the delta values for the neighbor cells.
 * Because we use modulo computations (in reality we use a bitmask in the
 * case that the domain resolutions are pow2), e.g. we dont subtract 1 to access the left
 * neighbor, but use a additive components of (DOMAIN_CELLS-1). Otherwise we would run
 * into more computations due to the modulo computations of negative numbers.
 */
#define DELTA_POS_X		(1)
#define DELTA_NEG_X		(DOMAIN_CELLS-1)

#define DELTA_POS_Y		(DOMAIN_CELLS_X)
#define DELTA_NEG_Y		(DOMAIN_CELLS-DOMAIN_CELLS_X)

#define DELTA_POS_Z		(DOMAIN_SLICE_CELLS)
#define DELTA_NEG_Z		(DOMAIN_CELLS-DOMAIN_SLICE_CELLS)

#include "src/cl_programs/wrap.h"

/**
 * equilibrium distributions f_eq for incompressible fluids
 * (for more information, have a look at Thuerey_070313.pdf, page 14)
 * f_eq = w_i * (	\rho + 3*(e_i*u) + (9/2)*(e_i*u)^2 - (3/2)*(u)^2	)
 * w_i = 1/3	for (0,0,0)
 * w_i = 1/18   for (+-1,0,0), (0,+-1, 0), (0,0,+-1)
 * w_i = 1/36   for (+-1,+-1,+-1)
 */

/*
 * we can reuse the following function because of its symmetry
 * f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0) f(0,0,1) f(0,0,-1)
 */
#define eq_dd_a0(vela, vela2, rho_alpha)	((T)(1.0f/18.0f)*((rho_alpha) + (T)(3.0f)*(vela) + (T)(9.0f/2.0f)*(vela2)))
#define eq_dd_a1(vela, vela2, rho_alpha)	((T)(1.0f/18.0f)*((rho_alpha) + (T)(-3.0f)*(vela) + (T)(9.0f/2.0f)*(vela2)))

/*
 * we can reuse the following functions because of the symmetry of the density distributions!
 *
 * f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
 * f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
 * f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
 */
#define eq_dd4(velx_add_vely, velx_add_vely_2, rho_alpha)		\
			((T)(1.0f/36.0f)*((rho_alpha) + (T)(3.0f)*(velx_add_vely) + (T)(9.0f/2.0f)*(velx_add_vely_2)))

#define eq_dd5(velx_add_vely, velx_add_vely_2, rho_alpha)		\
			((T)(1.0f/36.0f)*((rho_alpha) + (T)(-3.0f)*(velx_add_vely) + (T)(9.0f/2.0f)*(velx_add_vely_2)))

#define eq_dd6(velx_sub_vely, velx_sub_vely_2, rho_alpha)		\
			((T)(1.0f/36.0f)*((rho_alpha) + (T)(3.0f)*(velx_sub_vely) + (T)(9.0f/2.0f)*(velx_sub_vely_2)))

#define eq_dd7(velx_sub_vely, velx_sub_vely_2, rho_alpha)		\
			((T)(1.0f/36.0f)*((rho_alpha) + (T)(-3.0f)*(velx_sub_vely) + (T)(9.0f/2.0f)*(velx_sub_vely_2)))

/*
 * f(0,0,0)
 */
#define eq_dd18(rho_alpha)				\
			((T)(1.0f/3.0f)*(rho_alpha))


