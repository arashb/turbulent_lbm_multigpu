/*
 * Copyright 2010 Martin Schreiber
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef CLBMOPENCL_HH
#define CLBMOPENCL_HH

#include "CLbmSkeleton.hpp"
#include "libcl/CCL.hpp"
#include "lib/CError.hpp"
#include "libmath/CVector.hpp"
#include <typeinfo>
#include <iomanip>
#include <list>

#define LBM_FLAG_OBSTACLE			(1<<0)
#define LBM_FLAG_FLUID				(1<<1)
#define LBM_FLAG_VELOCITY_INJECTION	(1<<2)

/**
 * the density distributions are packed cell wise to optimize the collision operation and stored
 * linearly (first x, then y, then z) in the buffer cMemDensityDistributions
 *
 * this implementation is using the D3Q19 implementation
 *
 *
 * D3Q19 Vector enumerations + directions
 *
 * f(1,0,0), f(-1,0,0),  f(0,1,0),  f(0,-1,0)
 * f(1,1,0), f(-1,-1,0), f(1,-1,0), f(-1,1,0)
 * f(1,0,1), f(-1,0,-1), f(1,0,-1), f(-1,0,1)
 * f(0,1,1), f(0,-1,-1), f(0,1,-1), f(0,-1,1)
 * f(0,0,1), f(0,0,-1),  f(0,0,0),  X
 *
 * this enumeration makes the computation of the equilibrium distribution more efficient
 *
 * velocity(vx, vy, vz, phi)
 *
 *
 *       14,12
 *     15  |  18
 *       \ | /
 *        \|/
 * 7,9----- ----- 8,10
 *        /|\
 *       / | \
 *     17  |  16
 *       11,13
 *
 */

template <typename T>
class CLbmOpenCl
	:	public CLbmSkeleton<T>
{
public:
	CVector<4,T> drivenCavityVelocity;
private:
	static const size_t SIZE_DD_HOST = 19;
	static const size_t SIZE_DD_HOST_BYTES = SIZE_DD_HOST*sizeof(T);

	// opencl handlers
	CCL::CCommandQueue &cCommandQueue;
	CCL::CContext &cContext;
	CCL::CDevice &cDevice;
	CCL::CDeviceInfo cDeviceInfo; 

	// INITIALIZATION KERNEL
	CCL::CKernel cKernelInit;
	size_t cKernelInit_WorkGroupSize;
	size_t cKernelInit_GlobalWorkGroupSize;
	size_t cKernelInit_MaxRegisters;

	// COLLISION KERNELS
	CCL::CKernel cLbmKernelAlpha;
	size_t cLbmKernelAlpha_WorkGroupSize;
	size_t cLbmKernelAlpha_GlobalWorkGroupSize;
	size_t cLbmKernelAlpha_MaxRegisters;

	CCL::CKernel cLbmKernelBeta;
	size_t cLbmKernelBeta_WorkGroupSize;
	size_t cLbmKernelBeta_GlobalWorkGroupSize;
	size_t cLbmKernelBeta_MaxRegisters;

	// density distributions (D3Q19 model)
	CCL::CMem cMemDensityDistributions;

	// flags giving e.g. obstacles, (maybe gas) or other properties
	CCL::CMem cMemCellFlags;

	// Velocity (3 components) and Density (1 component) of cell
	CCL::CMem cMemVelocity;
	CCL::CMem cMemDensity;

	size_t computation_kernel_count;
	size_t domain_cells_count;

	bool store_velocity;	// store velocity for visualization
	bool store_density;

	bool debug;		// output debug informations

public:
	size_t simulation_step_counter;
	CError error;

	std::list<int> lbm_opencl_number_of_work_items_list;		///< list with number of threads for each successively created kernel
	std::list<int> lbm_opencl_number_of_registers_list;		///< list with number of registers for each thread threads for each successively created kernel


	// viscosity parameters:
	// http://en.wikipedia.org/wiki/Viscosity
	//
	CLbmOpenCl(	CCL::CCommandQueue &p_cCommandQueue,
				CCL::CContext &p_cContext,
				CCL::CDevice &p_cDevice,

				CVector<3,int> &p_domain_cells,
				T p_d_domain_x_length,
				CVector<3,T> &p_d_gravitation,
				T p_d_viscosity,
				size_t p_computation_kernel_count,
				bool p_debug,
				bool p_store_velocity,
				bool p_store_density,
				T p_d_timestep,

				std::list<int> &p_lbm_opencl_number_of_work_items_list,		///< list with number of threads for each successively created kernel
				std::list<int> &p_lbm_opencl_number_of_registers_list		///< list with number of registers for each thread threads for each successively created kernel
		) :
		drivenCavityVelocity(100.0, 0, 0, 1),

		cCommandQueue(p_cCommandQueue),
		cContext(p_cContext),
		cDevice(p_cDevice),
		cDeviceInfo(p_cDevice),
		computation_kernel_count(p_computation_kernel_count),
		lbm_opencl_number_of_work_items_list(p_lbm_opencl_number_of_work_items_list),
		lbm_opencl_number_of_registers_list(p_lbm_opencl_number_of_registers_list)
	{
		CLbmSkeleton<T>::init(p_domain_cells, p_d_domain_x_length, p_d_gravitation, p_d_viscosity, 1.0);

		if (CLbmSkeleton<T>::error())
			error << CLbmSkeleton<T>::error.getString();

		store_velocity = p_store_velocity;
		store_density = p_store_density;
		debug = p_debug;

		reload();
	}

	void addDrivenCavityValue(T value)
	{
		drivenCavityVelocity[0] += value;
		CVector<4,T> paramDrivenCavityVelocity = drivenCavityVelocity;
		paramDrivenCavityVelocity *= CLbmSkeleton<T>::d_timestep;

		cKernelInit.setArg(4, paramDrivenCavityVelocity[0]);
		cLbmKernelAlpha.setArg(8, paramDrivenCavityVelocity[0]);
		cLbmKernelBeta.setArg(8, paramDrivenCavityVelocity[0]);
	}

	void reload()
	{
		/**
		 * WORK GROUP SIZE
		 *
		 * initialize the variable (postfixed appropriately with _WorkGroupSize and _MaxRegisters)
		 * with either the standard value (max_local_work_group_size) or with the value from the list
		 */
#define INIT_WORK_GROUP_SIZE(variable)										\
		if (it != this->lbm_opencl_number_of_work_items_list.end())		\
		{																\
			variable##_WorkGroupSize = (*it != 0 ? *it : computation_kernel_count);	\
			it++;														\
		}																\
		else															\
		{																\
			variable##_WorkGroupSize = computation_kernel_count;		\
		}																\
																		\
		/* enlarge global work group size to be a multiple of collprop_work_group_size */	\
		if (domain_cells_count % variable##_WorkGroupSize != 0)			\
			variable##_GlobalWorkGroupSize = (domain_cells_count / variable##_WorkGroupSize + 1) * variable##_WorkGroupSize;	\
																		\
		if (ir != this->lbm_opencl_number_of_registers_list.end())		\
		{																\
			variable##_MaxRegisters = (*ir != 0 ? *ir : 0);				\
			ir++;														\
		}																\
		else															\
		{																\
			variable##_MaxRegisters  = 0;								\
		}

		std::list<int>::iterator it = this->lbm_opencl_number_of_work_items_list.begin();
		std::list<int>::iterator ir = this->lbm_opencl_number_of_registers_list.begin();

		INIT_WORK_GROUP_SIZE(cKernelInit);
		INIT_WORK_GROUP_SIZE(cLbmKernelAlpha);
		INIT_WORK_GROUP_SIZE(cLbmKernelBeta);

		/**
		 * program defines
		 */
		domain_cells_count = this->domain_cells.elements();

		std::ostringstream cl_program_defines;

		cl_program_defines << "#define SIZE_DD_HOST_BYTES (" << SIZE_DD_HOST_BYTES << ")" << std::endl;

		if (typeid(T) == typeid(float))
		{
			cl_program_defines << "typedef float T;" << std::endl;
			cl_program_defines << "typedef float2 T2;" << std::endl;
			cl_program_defines << "typedef float4 T4;" << std::endl;
			cl_program_defines << "typedef float8 T8;" << std::endl;
			cl_program_defines << "typedef float16 T16;" << std::endl;
		}
		else if (typeid(T) == typeid(double))
		{
			cl_program_defines << "#pragma OPENCL EXTENSION cl_khr_fp64 : enable" << std::endl;
			cl_program_defines << "#ifndef cl_khr_fp64" << std::endl;
			cl_program_defines << "	#error cl_khr_fp64 not supported - please switch to single precision to run the simulation" << std::endl;
			cl_program_defines << "#endif cl_khr_fp64" << std::endl;

			cl_program_defines << "typedef double T;" << std::endl;
			cl_program_defines << "typedef double2 T2;" << std::endl;
			cl_program_defines << "typedef double4 T4;" << std::endl;
			cl_program_defines << "typedef double8 T8;" << std::endl;
			cl_program_defines << "typedef double16 T16;" << std::endl;
		}
		else
		{
			error << "unsupported class type T" << std::endl;
			return;
		}

		cl_program_defines << "#define DOMAIN_CELLS_X	(" << this->domain_cells[0] << ")" << std::endl;
		cl_program_defines << "#define DOMAIN_CELLS_Y	(" << this->domain_cells[1] << ")" << std::endl;
		cl_program_defines << "#define DOMAIN_CELLS_Z	(" << this->domain_cells[2] << ")" << std::endl;
		cl_program_defines << "#define GLOBAL_WORK_GROUP_SIZE	(DOMAIN_CELLS_X*DOMAIN_CELLS_Y*DOMAIN_CELLS_Z)" << std::endl;

		cl_program_defines << "#define FLAG_OBSTACLE	(" << LBM_FLAG_OBSTACLE << ")" << std::endl;
		cl_program_defines << "#define FLAG_FLUID	(" << LBM_FLAG_FLUID << ")" << std::endl;
		cl_program_defines << "#define FLAG_VELOCITY_INJECTION	(" << LBM_FLAG_VELOCITY_INJECTION << ")" << std::endl;

		if (store_velocity)
			cl_program_defines << "#define STORE_VELOCITY 1" << std::endl;

		if (store_density)
			cl_program_defines << "#define STORE_DENSITY 1" << std::endl;

		if (debug)
			std::cout << cl_program_defines.str() << std::endl;

		/*
		 * ALLOCATE BUFFERS
		 */
		cMemDensityDistributions.create(cContext, CL_MEM_READ_WRITE, sizeof(T)*domain_cells_count*SIZE_DD_HOST, NULL);
		//cMemCellFlags.create(cContext, CL_MEM_READ_WRITE, sizeof(cl_char)*domain_cells_count, NULL);
		cMemCellFlags.create(cContext, CL_MEM_READ_WRITE, sizeof(cl_int)*domain_cells_count, NULL);
		cMemVelocity.create(cContext, CL_MEM_READ_WRITE, sizeof(T)*domain_cells_count*3, NULL);
		cMemDensity.create(cContext, CL_MEM_READ_WRITE, sizeof(T)*domain_cells_count, NULL);



		/**
		 * prepare INIT KERNEL DATA
		 * nothing to do because we use COLLISION KERNEL DATA
		 */

		cl_program_defines << std::endl;

		std::string cProgramDefinesPostfixString;
		std::string cProgramCompileOptionsString;
		char charbuf[255];

		/*
		 * INIT kernel
		 */
		sprintf(charbuf, "%i", (int)cKernelInit_WorkGroupSize);
		cProgramDefinesPostfixString = "#define LOCAL_WORK_GROUP_SIZE	(";
		cProgramDefinesPostfixString += charbuf;
		cProgramDefinesPostfixString +=  ")";

		/*cProgramCompileOptionsString = "-Werror -I./";*/
		cProgramCompileOptionsString = "-I./";
		if (cKernelInit_MaxRegisters != 0)
		{
			/* TODO: check for cl_nv_compiler_options extension */
			cProgramCompileOptionsString += " -cl-nv-maxrregcount=";
			cProgramCompileOptionsString += cKernelInit_MaxRegisters;
		}

		cKernelInit_GlobalWorkGroupSize = domain_cells_count;
		if (cKernelInit_GlobalWorkGroupSize % cKernelInit_WorkGroupSize != 0)
			cKernelInit_GlobalWorkGroupSize = (cKernelInit_GlobalWorkGroupSize / cKernelInit_WorkGroupSize + 1) * cKernelInit_WorkGroupSize;

		CCL::CProgram cProgramInit;
		cProgramInit.load(cContext, cl_program_defines.str()+cProgramDefinesPostfixString, "data/cl_programs/lbm_init.cl");
		cProgramInit.build(cDevice, cProgramCompileOptionsString.c_str());
		if (cProgramInit.error())
		{
			error << "failed to compile lbm_init.cl" << CError::endl;
			error << cProgramInit.error.getString() << CError::endl;
			return;
		}
		if (debug)
			std::cout << "KernelInit:	local_work_group_size: " << cKernelInit_WorkGroupSize << "		max_registers: " << cKernelInit_MaxRegisters << std::endl;


		/*
		 * ALPHA
		 */
		sprintf(charbuf, "%i", (int)cLbmKernelAlpha_WorkGroupSize);
		cProgramDefinesPostfixString = "#define LOCAL_WORK_GROUP_SIZE	(";
		cProgramDefinesPostfixString += charbuf;
		cProgramDefinesPostfixString +=  ")";

		// cProgramCompileOptionsString = "-Werror -I./";
		cProgramCompileOptionsString = "-I./";
		if (cLbmKernelAlpha_MaxRegisters != 0)
		{
			/* TODO: check for cl_nv_compiler_options extension */
			cProgramCompileOptionsString += " -cl-nv-maxrregcount=";
			cProgramCompileOptionsString += cLbmKernelAlpha_MaxRegisters;
		}

		cLbmKernelAlpha_GlobalWorkGroupSize = domain_cells_count;
		if (cLbmKernelAlpha_GlobalWorkGroupSize % cLbmKernelAlpha_WorkGroupSize != 0)
			cLbmKernelAlpha_GlobalWorkGroupSize = (cKernelInit_GlobalWorkGroupSize / cLbmKernelAlpha_WorkGroupSize + 1) * cLbmKernelAlpha_WorkGroupSize;

		CCL::CProgram cProgramAlpha;
		cProgramAlpha.load(cContext, cl_program_defines.str()+cProgramDefinesPostfixString, "data/cl_programs/lbm_alpha.cl");
		cProgramAlpha.build(cDevice, cProgramCompileOptionsString.c_str());
		if (cProgramAlpha.error())
		{
			error << "failed to compile lbm_alpha.cl" << CError::endl;
			error << cProgramAlpha.error.getString() << CError::endl;
			return;
		}
		if (debug)
			std::cout << "KernelAlpha:	local_work_group_size: " << cLbmKernelAlpha_WorkGroupSize << "		max_registers: " << cLbmKernelAlpha_MaxRegisters << std::endl;

		/*
		 * BETA
		 */
		sprintf(charbuf, "%i", (int)cLbmKernelBeta_WorkGroupSize);
		cProgramDefinesPostfixString = "#define LOCAL_WORK_GROUP_SIZE	(";
		cProgramDefinesPostfixString += charbuf;
		cProgramDefinesPostfixString +=  ")";

		cProgramCompileOptionsString = "-Werror -I./";
		if (cLbmKernelBeta_MaxRegisters != 0)
		{
			/* TODO: check for cl_nv_compiler_options extension */
			cProgramCompileOptionsString += " -cl-nv-maxrregcount=";
			cProgramCompileOptionsString += cLbmKernelBeta_MaxRegisters;
		}

		cLbmKernelBeta_GlobalWorkGroupSize = domain_cells_count;
		if (cLbmKernelBeta_GlobalWorkGroupSize % cLbmKernelBeta_WorkGroupSize != 0)
			cLbmKernelBeta_GlobalWorkGroupSize = (cKernelInit_GlobalWorkGroupSize / cLbmKernelBeta_WorkGroupSize + 1) * cLbmKernelBeta_WorkGroupSize;

		CCL::CProgram cProgramBeta;
		cProgramBeta.load(cContext, cl_program_defines.str()+cProgramDefinesPostfixString, "data/cl_programs/lbm_beta.cl");
		cProgramBeta.build(cDevice, cProgramCompileOptionsString.c_str());
		if (cProgramBeta.error())
		{
			error << "failed to compile lbm_beta.cl" << CError::endl;
			error << cProgramBeta.error.getString() << CError::endl;

			return;
		}
		if (debug)
			std::cout << "KernelBeta:	local_work_group_size: " << cLbmKernelBeta_WorkGroupSize << "		max_registers: " << cLbmKernelBeta_MaxRegisters << std::endl;

		/**
		 * create kernels and setup arguments
		 */
		CVector<4,T> paramDrivenCavityVelocity = drivenCavityVelocity;
		paramDrivenCavityVelocity *= CLbmSkeleton<T>::d_timestep;

		if (debug)
		{
			std::cout << "driven cavity velocity: " << paramDrivenCavityVelocity << std::endl;
			std::cout << "inverse tau: " << this->inv_tau << std::endl;
			std::cout << "d_timestep: " << this->d_timestep << std::endl;
			std::cout << "gravitaton: " << this->gravitation << std::endl;
		}
		/**
		 * SETUP ARGUMENTS
		 */
		// initialization kernel
		cKernelInit.create(cProgramInit, "init_kernel");
		cKernelInit.setArg(0, cMemDensityDistributions);
		cKernelInit.setArg(1, cMemCellFlags);
		cKernelInit.setArg(2, cMemVelocity);
		cKernelInit.setArg(3, cMemDensity);
		cKernelInit.setArg(4, paramDrivenCavityVelocity[0]);

		// collision and propagation kernels (alpha and beta)
		cLbmKernelAlpha.create(cProgramAlpha, "lbm_kernel_alpha");
		cLbmKernelAlpha.setArg(0, cMemDensityDistributions);
		cLbmKernelAlpha.setArg(1, cMemCellFlags);
		cLbmKernelAlpha.setArg(2, cMemVelocity);
		cLbmKernelAlpha.setArg(3, cMemDensity);
		cLbmKernelAlpha.setArg(4, this->inv_tau);
		cLbmKernelAlpha.setArg(5, this->gravitation[0]);
		cLbmKernelAlpha.setArg(6, this->gravitation[1]);
		cLbmKernelAlpha.setArg(7, this->gravitation[2]);
		cLbmKernelAlpha.setArg(8, paramDrivenCavityVelocity[0]);

		cLbmKernelBeta.create(cProgramBeta, "lbm_kernel_beta");
		cLbmKernelBeta.setArg(0, cMemDensityDistributions);
		cLbmKernelBeta.setArg(1, cMemCellFlags);
		cLbmKernelBeta.setArg(2, cMemVelocity);
		cLbmKernelBeta.setArg(3, cMemDensity);
		cLbmKernelBeta.setArg(4, this->inv_tau);
		cLbmKernelBeta.setArg(5, this->gravitation[0]);
		cLbmKernelBeta.setArg(6, this->gravitation[1]);
		cLbmKernelBeta.setArg(7, this->gravitation[2]);
		cLbmKernelBeta.setArg(8, paramDrivenCavityVelocity[0]);


		reset();
	}

	void reset()
	{
		simulation_step_counter = 0;

		if (debug)
			std::cout << "Init Simulation: " << std::flush;

		cCommandQueue.enqueueNDRangeKernel(	cKernelInit,	// kernel
											1,				// dimensions
											NULL,			// global work offset
											&cKernelInit_GlobalWorkGroupSize,
											&cKernelInit_WorkGroupSize
						);

		cCommandQueue.enqueueBarrier();

		if (debug)
			std::cout << "OK" << std::endl;
	}

	/**
	 * start one simulation step (enqueue kernels)
	 */
	void simulationStep()
	{
		/*
		 * collision kernels are inserted as 1d kernels because they work cell wise without neighboring information
		 */

		if (simulation_step_counter & 1)
		{
			cCommandQueue.enqueueNDRangeKernel(	cLbmKernelAlpha,	// kernel
												1,						// dimensions
												NULL,					// global work offset
												&cLbmKernelAlpha_GlobalWorkGroupSize,
												&cLbmKernelAlpha_WorkGroupSize
							);
		}
		else
		{
			cCommandQueue.enqueueNDRangeKernel(	cLbmKernelBeta,	// kernel
												1,						// dimensions
												NULL,					// global work offset
												&cLbmKernelBeta_GlobalWorkGroupSize,
												&cLbmKernelBeta_WorkGroupSize
							);
		}
		cCommandQueue.enqueueBarrier();

		simulation_step_counter++;
	}

	/**
	 * wait until all kernels have finished
	 */
	void wait()
	{
		cCommandQueue.finish();
	}

	/**
	 * store velocity and density values to host memory
	 * this is useful for a host memory based visualization
	 */
	void storeVelocity(T *dst)
	{
		size_t byte_size = cMemVelocity.getSize();

		cCommandQueue.enqueueReadBuffer(	cMemVelocity,
							CL_TRUE,	// sync reading
							0,
							byte_size,
							dst);
	}

	void storeDensity(T *dst)
	{
		size_t byte_size = cMemDensity.getSize();

		cCommandQueue.enqueueReadBuffer(	cMemDensity,
							CL_TRUE,	// sync reading
							0,
							byte_size,
							dst);
	}

	void storeFlags(int *dst)
	{
		size_t byte_size = cMemDensity.getSize();

		cCommandQueue.enqueueReadBuffer(	cMemCellFlags,
							CL_TRUE,	// sync reading
							0,
							byte_size,
							dst);
	}

private:
	void debugChar(CCL::CMem &cMem, size_t wrap_size = 20)
	{
		size_t char_size = cMem.getSize();
		char *buffer = new char[char_size];

		cCommandQueue.enqueueReadBuffer(	cMem,
							CL_TRUE,	// sync reading
							0,
							char_size,
							buffer);


		for (size_t i = 0; i < char_size; i++)
		{
			if (i % wrap_size == 0)
			{
				std::cout << std::endl;
				std::cout << (i/wrap_size) << ": ";
			}

			std::cout << (int)buffer[i] << " ";
		}

		delete [] buffer;
	}

	void debugFloat(CCL::CMem &cMem, size_t wrap_size = 20)
	{
		size_t byte_size = cMem.getSize();
		size_t T_size = byte_size / sizeof(T);

		T *buffer = new T[T_size];

		cCommandQueue.enqueueReadBuffer(	cMem,
							CL_TRUE,	// sync reading
							0,
							byte_size,
							buffer);

		std::streamsize ss = std::cout.precision();
		std::cout.precision(4);
		std::cout.setf(std::ios::fixed,std::ios::floatfield);

		for (size_t i = 0; i < T_size; i++)
		{
			if (i % wrap_size == 0)
			{
				std::cout << std::endl;
				std::cout << (i/wrap_size) << ": ";
			}

			std::cout << buffer[i] << " ";
		}

		std::cout.precision(ss);
		std::cout << std::resetiosflags(std::ios::fixed);

		delete [] buffer;
	}

	/**
	 * debug handler to output some useful information
	 */
public:
	void debug_print()
	{
		// read out DensityDistributions
		std::cout << "DENSITY DISTRIBUTIONS:";
		//debugFloat(cMemDensityDistributions, SIZE_DD_HOST);
		debugFloat(cMemDensityDistributions, 16);
		std::cout << std::endl;

		// read out Velocity
		std::cout << std::endl;
		std::cout << "VELOCITY:";
		debugFloat(cMemVelocity, 4*3);
		std::cout << std::endl;

		// read out VelocityDensity
		std::cout << std::endl;
		std::cout << "DENSITY:";
		debugFloat(cMemDensity, 4);
		std::cout << std::endl;

		// read out Flags
		std::cout << std::endl;
		std::cout << "FLAGS:";
		debugChar(cMemCellFlags, 4*4);
		std::cout << std::endl;
	}

	/**
	 * debug density distributions (works only in linear mode!)
	 * \param	dd_id	specifies the dd number
	 */
	void debugDD(size_t dd_id = 0, size_t wrap_size = 16, size_t empty_line = 16)
	{
		CCL::CMem &cMem = cMemDensityDistributions;

		size_t byte_size = cMem.getSize();
		size_t T_size = byte_size / sizeof(T);

		T *buffer = new T[T_size];

		cCommandQueue.enqueueReadBuffer(	cMem,
							CL_TRUE,	// sync reading
							0,
							byte_size,
							buffer);

		std::streamsize ss = std::cout.precision();
		std::cout.precision(4);
		std::cout.setf(std::ios::fixed,std::ios::floatfield);

		size_t start = this->domain_cells.elements()*dd_id;
		size_t end = this->domain_cells.elements()*(dd_id+1);

		for (size_t i = start; i < end; i++)
		{
			if (empty_line != wrap_size)
				if (i % empty_line == 0 && i != start)
				{
					std::cout << std::endl;
				}

			if (i % wrap_size == 0)
			{
				if (i != start)
					std::cout << std::endl;
				std::cout << (i/wrap_size) << ": ";
			}


			std::cout << buffer[i] << " ";
		}

		std::cout.precision(ss);
		std::cout << std::resetiosflags(std::ios::fixed);
		std::cout << std::endl;

		delete [] buffer;
	}


	float getVelocityChecksum()
	{
		T *velocity = new T[this->domain_cells.elements()*3];
		int *flags = new int[this->domain_cells.elements()*3];

		storeVelocity(velocity);
		storeFlags(flags);

		T *velx = velocity;
		T *vely = velocity+this->domain_cells.elements();
		T *velz = velocity+this->domain_cells.elements()*2;

		float checksum = 0;
		for (int a = 0; a < this->domain_cells.elements(); a++)
			if (flags[a] == LBM_FLAG_FLUID)
				checksum += velx[a] + vely[a] + velz[a];

		delete [] flags;
		delete [] velocity;

		return checksum;
	}
};

#endif
