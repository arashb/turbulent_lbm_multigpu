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

#include "common.h"

#define LBM_FLAG_OBSTACLE			(1<<0)
#define LBM_FLAG_FLUID				(1<<1)
#define LBM_FLAG_VELOCITY_INJECTION	(1<<2)
#define LBM_FLAG_GHOST_LAYER		(1<<3)

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

template<typename T>
class CLbmSolver: public CLbmSkeleton<T>
{
public:
	CVector<4, T> drivenCavityVelocity;
	static const size_t SIZE_DD_HOST = 19;
private:
	int _UID;
	static const size_t SIZE_DD_HOST_BYTES = SIZE_DD_HOST * sizeof(T);
	int _BC[3][2]; ///< Boundary conditions. First index specifys the dimension and second the upper or the lower boundary.
	// opencl handlers
	CCL::CCommandQueue &cCommandQueue;
	CCL::CContext &cContext;
	CCL::CDevice &cDevice;
	CCL::CDeviceInfo cDeviceInfo;
	OPENCL_VERSION _cl_version;

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

	// INITIALIZATION KERNEL
	CCL::CKernel cKernelCopyRect;
	size_t cKernelCopyRect_WorkGroupSize;
	size_t cKernelCopyRect_GlobalWorkGroupSize;
	size_t cKernelCopyRect_MaxRegisters;

	// density distributions (D3Q19 model)
	CCL::CMem cMemDensityDistributions;

	// flags giving e.g. obstacles, (maybe gas) or other properties
	CCL::CMem cMemCellFlags;
	CCL::CMem cMemBC;

	// Velocity (3 components) and Density (1 component) of cell
	CCL::CMem cMemVelocity;
	CCL::CMem cMemDensity;

	size_t computation_kernel_count;
	size_t domain_cells_count;

	bool store_velocity; // store velocity for visualization
	bool store_density;
private:

	std::vector<std::string> split(const std::string& str,
			const std::string& delimiter = " ") {
		std::vector<std::string> tokens;

		std::string::size_type lastPos = 0;
		std::string::size_type pos = str.find(delimiter, lastPos);

		while (std::string::npos != pos) {
			// Found a token, add it to the vector.
			//std::cout << str.substr(lastPos, pos - lastPos) << std::endl;
			tokens.push_back(str.substr(lastPos, pos - lastPos));
			lastPos = pos + delimiter.size();
			pos = str.find(delimiter, lastPos);
		}

		tokens.push_back(str.substr(lastPos, str.size() - lastPos));
		return tokens;
	}

	OPENCL_VERSION getOpenCLVersion() {
		OPENCL_VERSION cl_version;
		// format:
		// OpenCL<space><major_version.minor_version><space><vendor-specific information>
		std::string major_minor_version = split(
				std::string(cDeviceInfo.version))[1];
		std::vector<std::string> major_minor_vector = split(major_minor_version,
				".");
		int major_version = atoi(major_minor_vector[0].c_str());
		int minor_version = atoi(major_minor_vector[1].c_str());
		int version = major_version * 100 + minor_version * 10;
		switch (version) {
		case OPENCL_VERSION_1_0_0:
			cl_version = OPENCL_VERSION_1_0_0;
			break;
		case OPENCL_VERSION_1_1_0:
			cl_version = OPENCL_VERSION_1_1_0;
			break;
		case OPENCL_VERSION_1_2_0:
			cl_version = OPENCL_VERSION_1_2_0;
			break;
		default:
			cl_version = OPENCL_VERSION_UNKNOWN;
			break;
		}
		return cl_version;
	}

	inline void enqueueCopyRectKernel(CCL::CMem& src, CCL::CMem& dst,
			const int src_offset, CVector<3, int> src_origin,
			CVector<3, int> src_size, const int dst_offset,
			CVector<3, int> dst_origin, CVector<3, int> dst_size,
			CVector<3, int> block_size, bool withBarrier = true) {
		// set kernel args
		// source args
		cKernelCopyRect.setArg(0, src);
		cKernelCopyRect.setArg(1, src_offset);
		cKernelCopyRect.setArg(2, src_origin[0]);
		cKernelCopyRect.setArg(3, src_origin[1]);
		cKernelCopyRect.setArg(4, src_origin[2]);
		cKernelCopyRect.setArg(5, src_size[0]);
		cKernelCopyRect.setArg(6, src_size[1]);
		cKernelCopyRect.setArg(7, src_size[2]);

		// destination args
		cKernelCopyRect.setArg(8, dst);
		cKernelCopyRect.setArg(9, dst_offset);
		cKernelCopyRect.setArg(10, dst_origin[0]);
		cKernelCopyRect.setArg(11, dst_origin[1]);
		cKernelCopyRect.setArg(12, dst_origin[2]);
		cKernelCopyRect.setArg(13, dst_size[0]);
		cKernelCopyRect.setArg(14, dst_size[1]);
		cKernelCopyRect.setArg(15, dst_size[2]);
		cKernelCopyRect.setArg(16, block_size[0]);

		size_t lGlobalSize[2];
		lGlobalSize[0] = block_size[1];
		lGlobalSize[1] = block_size[2];
		// enqueue the CopyRect kernel
		cCommandQueue.enqueueNDRangeKernel(cKernelCopyRect, // kernel
				2, // dimensions
				NULL, // global work offset
				lGlobalSize, NULL);
		if (withBarrier)
			cCommandQueue.enqueueBarrier();
	}

public:
	size_t simulation_step_counter;
	CError error;

	std::list<int> lbm_opencl_number_of_work_items_list; ///< list with number of threads for each successively created kernel
	std::list<int> lbm_opencl_number_of_registers_list; ///< list with number of registers for each thread threads for each successively created kernel

	// viscosity parameters:
	// http://en.wikipedia.org/wiki/Viscosity
	//
	CLbmSolver(int UID, CCL::CCommandQueue &p_cCommandQueue,
			CCL::CContext &p_cContext, CCL::CDevice &p_cDevice, int BC[3][2],
			CDomain<T> &domain, CVector<3, T> &p_d_gravitation, T p_d_viscosity,
			size_t p_computation_kernel_count,
			//bool p_debug,
			bool p_store_velocity, bool p_store_density, T p_d_timestep,
			CVector<4, T>& _drivenCavityVelocity,
			std::list<int> &p_lbm_opencl_number_of_work_items_list, ///< list with number of threads for each successively created kernel
			std::list<int> &p_lbm_opencl_number_of_registers_list ///< list with number of registers for each thread threads for each successively created kernel
			) :
			CLbmSkeleton<T>(CDomain<T>(domain), _drivenCavityVelocity), //drivenCavityVelocity(100.0, 0,0, 1),
			drivenCavityVelocity(_drivenCavityVelocity), _UID(UID), cCommandQueue(
					p_cCommandQueue), cContext(p_cContext), cDevice(p_cDevice), cDeviceInfo(
					p_cDevice), computation_kernel_count(
					p_computation_kernel_count), lbm_opencl_number_of_work_items_list(
					p_lbm_opencl_number_of_work_items_list), lbm_opencl_number_of_registers_list(
					p_lbm_opencl_number_of_registers_list) {
		CLbmSkeleton<T>::init(p_d_gravitation, p_d_viscosity, 1.0);

		if (CLbmSkeleton<T>::error())
			error << CLbmSkeleton<T>::error.getString();

		store_velocity = p_store_velocity;
		store_density = p_store_density;
		//debug = p_debug;

		_cl_version = OPENCL_VERSION_1_0_0; //getOpenCLVersion();
		if (_cl_version == OPENCL_VERSION_UNKNOWN)
			throw "OpenCL Version is unknown!";

#if DEBUG
		std::cout << "CL_VERSION " << _cl_version << std::endl;
#endif
		// setting the boundary conditions
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 2; j++)
				_BC[i][j] = BC[i][j];

		reload();
	}

	void addDrivenCavityValue(T value) {
		drivenCavityVelocity[0] += value;
		CVector<4, T> paramDrivenCavityVelocity = drivenCavityVelocity;
		paramDrivenCavityVelocity *= CLbmSkeleton<T>::d_timestep;

		cKernelInit.setArg(4, paramDrivenCavityVelocity[0]);
		cLbmKernelAlpha.setArg(8, paramDrivenCavityVelocity[0]);
		cLbmKernelBeta.setArg(8, paramDrivenCavityVelocity[0]);
	}

	void reload() {
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

		std::list<int>::iterator it =
				this->lbm_opencl_number_of_work_items_list.begin();
		std::list<int>::iterator ir =
				this->lbm_opencl_number_of_registers_list.begin();

		INIT_WORK_GROUP_SIZE(cKernelInit);
		INIT_WORK_GROUP_SIZE(cLbmKernelAlpha);
		INIT_WORK_GROUP_SIZE(cLbmKernelBeta);
		INIT_WORK_GROUP_SIZE(cKernelCopyRect);

		/**
		 * program defines
		 */
		domain_cells_count = this->domain_cells.elements();

		std::ostringstream cl_program_defines;

		cl_program_defines << "#define SIZE_DD_HOST_BYTES ("
				<< SIZE_DD_HOST_BYTES << ")" << std::endl;

		if (typeid(T) == typeid(float)) {
			cl_program_defines << "typedef float T;" << std::endl;
			cl_program_defines << "typedef float2 T2;" << std::endl;
			cl_program_defines << "typedef float4 T4;" << std::endl;
			cl_program_defines << "typedef float8 T8;" << std::endl;
			cl_program_defines << "typedef float16 T16;" << std::endl;
		} else if (typeid(T) == typeid(double)) {
			cl_program_defines
					<< "#pragma OPENCL EXTENSION cl_khr_fp64 : enable"
					<< std::endl;
			cl_program_defines << "#ifndef cl_khr_fp64" << std::endl;
			cl_program_defines
					<< "	#error cl_khr_fp64 not supported - please switch to single precision to run the simulation"
					<< std::endl;
			cl_program_defines << "#endif cl_khr_fp64" << std::endl;

			cl_program_defines << "typedef double T;" << std::endl;
			cl_program_defines << "typedef double2 T2;" << std::endl;
			cl_program_defines << "typedef double4 T4;" << std::endl;
			cl_program_defines << "typedef double8 T8;" << std::endl;
			cl_program_defines << "typedef double16 T16;" << std::endl;
		} else {
			error << "unsupported class type T" << std::endl;
			return;
		}

		cl_program_defines << "#define DOMAIN_CELLS_X	("
				<< this->domain_cells[0] << ")" << std::endl;
		cl_program_defines << "#define DOMAIN_CELLS_Y	("
				<< this->domain_cells[1] << ")" << std::endl;
		cl_program_defines << "#define DOMAIN_CELLS_Z	("
				<< this->domain_cells[2] << ")" << std::endl;
		cl_program_defines
				<< "#define GLOBAL_WORK_GROUP_SIZE	(DOMAIN_CELLS_X*DOMAIN_CELLS_Y*DOMAIN_CELLS_Z)"
				<< std::endl;

		cl_program_defines << "#define FLAG_OBSTACLE	(" << LBM_FLAG_OBSTACLE
				<< ")" << std::endl;
		cl_program_defines << "#define FLAG_FLUID	(" << LBM_FLAG_FLUID << ")"
				<< std::endl;
		cl_program_defines << "#define FLAG_VELOCITY_INJECTION	("
				<< LBM_FLAG_VELOCITY_INJECTION << ")" << std::endl;
		cl_program_defines << "#define FLAG_GHOST_LAYER	("
				<< LBM_FLAG_GHOST_LAYER << ")" << std::endl;

		if (store_velocity)
			cl_program_defines << "#define STORE_VELOCITY 1" << std::endl;

		if (store_density)
			cl_program_defines << "#define STORE_DENSITY 1" << std::endl;

#if DEBUG
		std::cout << cl_program_defines.str() << std::endl;
#endif
		/*
		 * ALLOCATE BUFFERS
		 */
		cMemDensityDistributions.create(cContext, CL_MEM_READ_WRITE,
				sizeof(T) * domain_cells_count * SIZE_DD_HOST, NULL);
		//cMemCellFlags.create(cContext, CL_MEM_READ_WRITE, sizeof(cl_char)*domain_cells_count, NULL);
		cMemCellFlags.create(cContext, CL_MEM_READ_WRITE,
				sizeof(cl_int) * domain_cells_count, NULL);
		cMemVelocity.create(cContext, CL_MEM_READ_WRITE,
				sizeof(T) * domain_cells_count * 3, NULL);
		cMemDensity.create(cContext, CL_MEM_READ_WRITE,
				sizeof(T) * domain_cells_count, NULL);

		int bc_linear[6];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 2; j++)
				bc_linear[i * 2 + j] = _BC[i][j];

#if DEBUG
		std::cout << "BOUNDARY CONDITION: "<< std::endl;
		for(int i = 0; i < 6; i++)
		std::cout << " " << bc_linear[i];
		std::cout << std::endl;
#endif
		cMemBC.create(cContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
				sizeof(cl_int) * 6, bc_linear);

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
		sprintf(charbuf, "%i", (int) cKernelInit_WorkGroupSize);
		cProgramDefinesPostfixString = "#define LOCAL_WORK_GROUP_SIZE	(";
		cProgramDefinesPostfixString += charbuf;
		cProgramDefinesPostfixString += ")";

		/*cProgramCompileOptionsString = "-Werror -I./";*/
		cProgramCompileOptionsString = "-I./";
		if (cKernelInit_MaxRegisters != 0) {
			/* TODO: check for cl_nv_compiler_options extension */
			cProgramCompileOptionsString += " -cl-nv-maxrregcount=";
			cProgramCompileOptionsString += cKernelInit_MaxRegisters;
		}

		cKernelInit_GlobalWorkGroupSize = domain_cells_count;
		if (cKernelInit_GlobalWorkGroupSize % cKernelInit_WorkGroupSize != 0)
			cKernelInit_GlobalWorkGroupSize = (cKernelInit_GlobalWorkGroupSize
					/ cKernelInit_WorkGroupSize + 1)
					* cKernelInit_WorkGroupSize;

		CCL::CProgram cProgramInit;
		cProgramInit.load(cContext,
				cl_program_defines.str() + cProgramDefinesPostfixString,
				"src/cl_programs/lbm_init.cl");
		cProgramInit.build(cDevice, cProgramCompileOptionsString.c_str());
		if (cProgramInit.error()) {
			error << "failed to compile lbm_init.cl" << CError::endl;
			error << cProgramInit.error.getString() << CError::endl;
			return;
		}
#if DEBUG
		std::cout << "KernelInit:	local_work_group_size: " << cKernelInit_WorkGroupSize << "		max_registers: " << cKernelInit_MaxRegisters << std::endl;
#endif

		/*
		 * ALPHA
		 */
		sprintf(charbuf, "%i", (int) cLbmKernelAlpha_WorkGroupSize);
		cProgramDefinesPostfixString = "#define LOCAL_WORK_GROUP_SIZE	(";
		cProgramDefinesPostfixString += charbuf;
		cProgramDefinesPostfixString += ")";

		// cProgramCompileOptionsString = "-Werror -I./";
		cProgramCompileOptionsString = "-I./";
		if (cLbmKernelAlpha_MaxRegisters != 0) {
			/* TODO: check for cl_nv_compiler_options extension */
			cProgramCompileOptionsString += " -cl-nv-maxrregcount=";
			cProgramCompileOptionsString += cLbmKernelAlpha_MaxRegisters;
		}

		cLbmKernelAlpha_GlobalWorkGroupSize = domain_cells_count;
		if (cLbmKernelAlpha_GlobalWorkGroupSize % cLbmKernelAlpha_WorkGroupSize
				!= 0)
			cLbmKernelAlpha_GlobalWorkGroupSize =
					(cKernelInit_GlobalWorkGroupSize
							/ cLbmKernelAlpha_WorkGroupSize + 1)
							* cLbmKernelAlpha_WorkGroupSize;

		CCL::CProgram cProgramAlpha;
		cProgramAlpha.load(cContext,
				cl_program_defines.str() + cProgramDefinesPostfixString,
				"src/cl_programs/lbm_alpha.cl");
		cProgramAlpha.build(cDevice, cProgramCompileOptionsString.c_str());
		if (cProgramAlpha.error()) {
			error << "failed to compile lbm_alpha.cl" << CError::endl;
			error << cProgramAlpha.error.getString() << CError::endl;
			return;
		}
#if DEBUG
		std::cout << "KernelAlpha:	local_work_group_size: " << cLbmKernelAlpha_WorkGroupSize << "		max_registers: " << cLbmKernelAlpha_MaxRegisters << std::endl;
#endif
		/*
		 * BETA
		 */
		sprintf(charbuf, "%i", (int) cLbmKernelBeta_WorkGroupSize);
		cProgramDefinesPostfixString = "#define LOCAL_WORK_GROUP_SIZE	(";
		cProgramDefinesPostfixString += charbuf;
		cProgramDefinesPostfixString += ")";

		//cProgramCompileOptionsString = "-Werror -I./";
		cProgramCompileOptionsString = "-I./";
		if (cLbmKernelBeta_MaxRegisters != 0) {
			/* TODO: check for cl_nv_compiler_options extension */
			cProgramCompileOptionsString += " -cl-nv-maxrregcount=";
			cProgramCompileOptionsString += cLbmKernelBeta_MaxRegisters;
		}

		cLbmKernelBeta_GlobalWorkGroupSize = domain_cells_count;
		if (cLbmKernelBeta_GlobalWorkGroupSize % cLbmKernelBeta_WorkGroupSize
				!= 0)
			cLbmKernelBeta_GlobalWorkGroupSize =
					(cKernelInit_GlobalWorkGroupSize
							/ cLbmKernelBeta_WorkGroupSize + 1)
							* cLbmKernelBeta_WorkGroupSize;

		CCL::CProgram cProgramBeta;
		cProgramBeta.load(cContext,
				cl_program_defines.str() + cProgramDefinesPostfixString,
				"src/cl_programs/lbm_beta.cl");
		cProgramBeta.build(cDevice, cProgramCompileOptionsString.c_str());
		if (cProgramBeta.error()) {
			error << "failed to compile lbm_beta.cl" << CError::endl;
			error << cProgramBeta.error.getString() << CError::endl;

			return;
		}
#if DEBUG
		std::cout << "KernelBeta:	local_work_group_size: " << cLbmKernelBeta_WorkGroupSize << "		max_registers: " << cLbmKernelBeta_MaxRegisters << std::endl;
#endif

		/*
		 * INIT CopyBufferRect
		 */
		sprintf(charbuf, "%i", (int) cKernelCopyRect_WorkGroupSize);
		cProgramDefinesPostfixString = "#define LOCAL_WORK_GROUP_SIZE	(";
		cProgramDefinesPostfixString += charbuf;
		cProgramDefinesPostfixString += ")";

		cProgramCompileOptionsString = "-Werror -I./";
		//cProgramCompileOptionsString = "-I./";
		if (cKernelCopyRect_MaxRegisters != 0) {
			/* TODO: check for cl_nv_compiler_options extension */
			cProgramCompileOptionsString += " -cl-nv-maxrregcount=";
			cProgramCompileOptionsString += cKernelCopyRect_MaxRegisters;
		}

		cKernelCopyRect_GlobalWorkGroupSize = domain_cells_count;
		if (cKernelInit_GlobalWorkGroupSize % cKernelInit_WorkGroupSize != 0)
			cKernelInit_GlobalWorkGroupSize = (cKernelInit_GlobalWorkGroupSize
					/ cKernelInit_WorkGroupSize + 1)
					* cKernelInit_WorkGroupSize;

		CCL::CProgram cProgramCopyRect;
		cProgramCopyRect.load(cContext,
				cl_program_defines.str() + cProgramDefinesPostfixString,
				"src/cl_programs/copy_buffer_rect.cl");
		cProgramCopyRect.build(cDevice, cProgramCompileOptionsString.c_str());
		if (cProgramCopyRect.error()) {
			error << "failed to compile copy_buffer_rect.cl" << CError::endl;
			error << cProgramCopyRect.error.getString() << CError::endl;
			return;
		}
#if DEBUG
		std::cout << "KernelCopyRect:	local_work_group_size: " << cKernelCopyRect_WorkGroupSize << "		max_registers: " << cKernelCopyRect_MaxRegisters << std::endl;
#endif

		/**
		 * create kernels and setup arguments
		 */
		CVector<4, T> paramDrivenCavityVelocity = CLbmSkeleton<T>::drivenCavityVelocity;
//		CVector<4, T> paramDrivenCavityVelocity = drivenCavityVelocity;
//		paramDrivenCavityVelocity *= CLbmSkeleton<T>::d_timestep;

#if DEBUG
		{
			std::cout << "driven cavity velocity: " << paramDrivenCavityVelocity << std::endl;
			std::cout << "inverse tau: " << this->inv_tau << std::endl;
			std::cout << "d_timestep: " << this->d_timestep << std::endl;
			std::cout << "gravitaton: " << this->gravitation << std::endl;
		}
#endif
		/**
		 * SETUP ARGUMENTS
		 */
		// initialization kernel
		cKernelInit.create(cProgramInit, "init_kernel");
		cKernelInit.setArg(0, cMemDensityDistributions);
		cKernelInit.setArg(1, cMemCellFlags);
		cKernelInit.setArg(2, cMemVelocity);
		cKernelInit.setArg(3, cMemDensity);
		cKernelInit.setArg(4, cMemBC);
		cKernelInit.setArg(5, paramDrivenCavityVelocity[0]);

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

		cKernelCopyRect.create(cProgramCopyRect, "copy_buffer_rect");

		reset();
	}

	void reset() {
		simulation_step_counter = 0;

#if DEBUG
		std::cout << "Init Simulation: " << std::flush;
#endif
		cCommandQueue.enqueueNDRangeKernel(cKernelInit, // kernel
				1, // dimensions
				NULL, // global work offset
				&cKernelInit_GlobalWorkGroupSize, &cKernelInit_WorkGroupSize);

		cCommandQueue.enqueueBarrier();

#if DEBUG
		std::cout << "OK" << std::endl;
#endif
	}

	void simulationStepAlpha() {
#if DEBUG
		std::cout << "--> Running Alpha kernel" << std::endl;
#endif
		cCommandQueue.enqueueNDRangeKernel(
				cLbmKernelAlpha, // kernel
				1, // dimensions
				NULL, // global work offset
				&cLbmKernelAlpha_GlobalWorkGroupSize,
				&cLbmKernelAlpha_WorkGroupSize);
	}

	void simulationStepBeta() {
#if DEBUG
		std::cout << "--> Running BETA kernel" << std::endl;
#endif
		cCommandQueue.enqueueNDRangeKernel(
				cLbmKernelBeta, // kernel
				1, // dimensions
				NULL, // global work offset
				&cLbmKernelBeta_GlobalWorkGroupSize,
				&cLbmKernelBeta_WorkGroupSize);
	}

	/**
	 * start one simulation step (enqueue kernels)
	 */
	void simulationStep() {
		/*
		 * collision kernels are inserted as 1d kernels because they work cell wise without neighboring information
		 */
		if (simulation_step_counter & 1) {
			simulationStepAlpha();
		} else {
			simulationStepBeta();
		}
		cCommandQueue.enqueueBarrier();

		simulation_step_counter++;
	}

	/**
	 * wait until all kernels have finished
	 */
	void wait() {
		cCommandQueue.finish();
	}

	/**
	 * store density distribution values to host memory
	 */
	void storeDensityDistribution(T *dst) {
		size_t byte_size = cMemDensityDistributions.getSize();

		cCommandQueue.enqueueReadBuffer(cMemDensityDistributions, CL_TRUE, // sync reading
				0, byte_size, dst);
	}

	void storeDensityDistribution(T *dst, CVector<3, int> &origin,
			CVector<3, int> &size) {
		CCL::CMem cBuffer;
		cBuffer.create(cContext, CL_MEM_READ_WRITE,
				sizeof(T) * size.elements() * SIZE_DD_HOST, NULL);

		if (_cl_version >= OPENCL_VERSION_1_1_0) // OpenCL 1.1 and later
				{
			// TODO: implement the clEnqueueReadBufferRect
		} else if (_cl_version >= OPENCL_VERSION_1_0_0) // OpenCL 1.0 and later
				{
			for (int f = 0; f < SIZE_DD_HOST; f++) {
				enqueueCopyRectKernel(cMemDensityDistributions, cBuffer,
						f * this->domain_cells_count, origin,
						this->domain_cells, f * size.elements(),
						CVector<3, int>(0, 0, 0), size, size, false);
			}
			cCommandQueue.enqueueBarrier();
		}
		clEnqueueBarrier(cCommandQueue.command_queue);
		cCommandQueue.enqueueReadBuffer(cBuffer, CL_TRUE, // sync reading
				0, cBuffer.getSize(), dst);
	}

	void setDensityDistribution(T *src, CVector<3, int> &origin,
			CVector<3, int> &size) {
		CCL::CMem cBuffer;
		cBuffer.create(cContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
				sizeof(T) * size.elements() * SIZE_DD_HOST, src);

		if (_cl_version >= OPENCL_VERSION_1_1_0) {
			// TODO: implement the clEnqueueWriteBufferRect
		} else if (_cl_version >= OPENCL_VERSION_1_0_0) {
			for (int f = 0; f < SIZE_DD_HOST; f++) {
				enqueueCopyRectKernel(cBuffer, cMemDensityDistributions,
						f * size.elements(), CVector<3, int>(0, 0, 0), size,
						f * this->domain_cells_count, origin,
						this->domain_cells, size, false);
			}
			cCommandQueue.enqueueBarrier();
		}
	}

	void setDensityDistribution(T *src, CVector<3, int> &origin,
			CVector<3, int> &size, CVector<3, int> norm) {
		CCL::CMem cBuffer;
		cBuffer.create(cContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
				sizeof(T) * size.elements() * SIZE_DD_HOST, src);

		if (_cl_version >= OPENCL_VERSION_1_1_0) {
			// TODO: implement the clEnqueueWriteBufferRect
		} else if (_cl_version >= OPENCL_VERSION_1_0_0) {
			for (int f = 0; f < SIZE_DD_HOST; f++) {
				if (norm.dotProd(lbm_units[f]) > 0) {
					enqueueCopyRectKernel(cBuffer, cMemDensityDistributions,
							f * size.elements(), CVector<3, int>(0, 0, 0), size,
							f * this->domain_cells_count, origin,
							this->domain_cells, size, false);
				}
			}
			cCommandQueue.enqueueBarrier();
		}
	}
	/**
	 * store velocity and density values to host memory
	 * this is useful for a host memory based visualization
	 *
	 * @param dst The buffer that will contain the return values
	 */
	void storeVelocity(T *dst) {
		size_t byte_size = cMemVelocity.getSize();
		cCommandQueue.enqueueReadBuffer(cMemVelocity, CL_TRUE, // sync reading
				0, byte_size, dst);
	}

	/**
	 * Store a block of velocity data from device to host
	 *
	 * @param dst The buffer that will contain the return values
	 * @param origin The origin point of data block
	 * @param size The size of data block
	 */
	void storeVelocity(T* dst, CVector<3, int> &origin, CVector<3, int> &size) {
		CCL::CMem cBuffer;
		cBuffer.create(cContext, CL_MEM_READ_WRITE,
				sizeof(T) * size.elements() * 3, NULL);

		if (_cl_version >= OPENCL_VERSION_1_1_0) {
			// TODO: implement the clEnqueueWriteBufferRect
		} else if (_cl_version >= OPENCL_VERSION_1_0_0) {
			for (int dim = 0; dim < 3; dim++) {
				enqueueCopyRectKernel(cMemVelocity, cBuffer,
						dim * this->domain_cells_count, origin,
						this->domain_cells, dim * size.elements(),
						CVector<3, int>(0, 0, 0), size, size, false);
			}
			cCommandQueue.enqueueBarrier();
		}
		cCommandQueue.enqueueReadBuffer(cBuffer, CL_TRUE, // sync reading
				0, cBuffer.getSize(), dst);
	}

	/**
	 * Set a block of velocity data from host to device
	 *
	 * @param src The buffer that contain the input values
	 * @param origin The origin point of data block
	 * @param size The size of data block
	 */
	void setVelocity(T* src, CVector<3, int> &origin, CVector<3, int> &size) {
		CCL::CMem cBuffer;
		cBuffer.create(cContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
				sizeof(T) * size.elements() * 3, src);

		if (_cl_version >= OPENCL_VERSION_1_1_0) {
			// TODO: implement the clEnqueueWriteBufferRect

		} else if (_cl_version >= OPENCL_VERSION_1_0_0) {
			for (int dim = 0; dim < 3; dim++) {
				enqueueCopyRectKernel(cBuffer, cMemVelocity,
						dim * size.elements(), CVector<3, int>(0, 0, 0), size,
						dim * this->domain_cells_count, origin,
						this->domain_cells, size, false);
			}
			cCommandQueue.enqueueBarrier();
		}
	}
	/**
	 * Store velocity data from device to host
	 *
	 * @param dst The buffer that will contain the return values
	 */
	void storeDensity(T *dst) {
		size_t byte_size = cMemDensity.getSize();

		cCommandQueue.enqueueReadBuffer(cMemDensity, CL_TRUE, // sync reading
				0, byte_size, dst);
	}

	/**
	 * Store a block of density data from device to host
	 *
	 * @param dst The buffer that will contain the return values
	 * @param origin The origin point of data block
	 * @param size The size of data block
	 */
	void storeDensity(T *dst, CVector<3, int> &origin, CVector<3, int> &size) {
		CCL::CMem cBuffer;
		cBuffer.create(cContext, CL_MEM_READ_WRITE, sizeof(T) * size.elements(),
				NULL);

		if (_cl_version >= OPENCL_VERSION_1_1_0) {
//			// TODO: test this part
//			size_t buffer_origin[3] = {origin[0], origin[1], origin[2]};
//			size_t host_origin[3] = {0,0,0};
//			size_t region[3] = {size[0], size[1], size[2]};
//			cCommandQueue.enqueueReadBufferRect(
//					cMemDensity,
//					CL_TRUE,
//					buffer_origin,
//					host_origin,
//					region,
//					this->domain_cells[0],
//					this->domain_cells[0]*this->domain_cells[1],
//					0,
//					0,
//					dst
//			);
		} else if (_cl_version >= OPENCL_VERSION_1_0_0) {
			enqueueCopyRectKernel(cMemDensity, cBuffer, 0, origin,
					this->domain_cells, 0, CVector<3, int>(0, 0, 0), size,
					size);
		}
		cCommandQueue.enqueueReadBuffer(cBuffer, CL_TRUE, // sync reading
				0, cBuffer.getSize(), dst);
	}

	/**
	 * Set a block of density data from host to device
	 *
	 * @param src The buffer that will contain the source values
	 * @param origin The origin point of data block
	 * @param size The size of data block
	 */
	void setDensity(T *src, CVector<3, int> &origin, CVector<3, int> &size) {
		if (_cl_version >= OPENCL_VERSION_1_1_0) {
//			// TODO: test this part
//			size_t buffer_origin[3] = {origin[0], origin[1], origin[2]};
//			size_t host_origin[3] = {0,0,0};
//			size_t region[3] = {size[0], size[1], size[2]};
//			cCommandQueue.enqueueWriteBufferRect(
//					cMemDensity,
//					CL_TRUE,
//					buffer_origin,
//					host_origin,
//					region,
//					this->domain_cells[0],
//					this->domain_cells[0]*this->domain_cells[1],
//					0,
//					0,
//					src
//			);
		} else if (_cl_version >= OPENCL_VERSION_1_0_0) {
			CCL::CMem cBuffer;
			cBuffer.create(cContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
					sizeof(T) * size.elements(), src);
			enqueueCopyRectKernel(cBuffer, cMemDensity, 0,
					CVector<3, int>(0, 0, 0), size, 0, origin,
					this->domain_cells, size);
		}
	}

	void storeFlags(int *dst) {
		size_t byte_size = cMemDensity.getSize();

		cCommandQueue.enqueueReadBuffer(cMemCellFlags, CL_TRUE, // sync reading
				0, byte_size, dst);
	}

	void storeFlags(int *dst, CVector<3, int> &origin, CVector<3, int> &size) {
		if (_cl_version >= OPENCL_VERSION_1_1_0) {
//			// TODO: test this part
//			size_t buffer_origin[3] = {origin[0], origin[1], origin[2]};
//			size_t host_origin[3] = {0,0,0};
//			size_t region[3] = {size[0], size[1], size[2]};
//			cCommandQueue.enqueueReadBufferRect(
//					cMemCellFlags,
//					CL_TRUE,
//					buffer_origin,
//					host_origin,
//					region,
//					this->domain_cells[0],
//					this->domain_cells[0]*this->domain_cells[1],
//					0,
//					0,
//					dst
//			);
		} else {
			CCL::CMem cBuffer;
			cBuffer.create(cContext, CL_MEM_READ_WRITE,
					sizeof(int) * size.elements(), NULL);

			if (_cl_version >= OPENCL_VERSION_1_0_0) {
				enqueueCopyRectKernel(cMemCellFlags, cBuffer, 0, origin,
						this->domain_cells, 0, CVector<3, int>(0, 0, 0), size,
						size);
			}
			cCommandQueue.enqueueReadBuffer(cBuffer, CL_TRUE, // sync reading
					0, cBuffer.getSize(), dst);
		}
	}

	void setFlags(int *src, CVector<3, int> &origin, CVector<3, int> &size) {
		if (_cl_version >= OPENCL_VERSION_1_1_0) // OpenCL 1.2 and later
				{
//			// TODO: test this part
//			size_t buffer_origin[3] = {origin[0], origin[1], origin[2]};
//			size_t host_origin[3] = {0,0,0};
//			size_t region[3] = {size[0], size[1], size[2]};
//			cCommandQueue.enqueueWriteBufferRect(
//					cMemCellFlags,
//					CL_TRUE,
//					buffer_origin,
//					host_origin,
//					region,
//					this->domain_cells[0],
//					this->domain_cells[0]*this->domain_cells[1],
//					0,
//					0,
//					src
//			);
		} else {
			CCL::CMem cBuffer;
			cBuffer.create(cContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
					sizeof(int) * size.elements(), src);

			if (_cl_version >= OPENCL_VERSION_1_0_0) // OpenCL 1.0 and later
					{
				enqueueCopyRectKernel(cBuffer, cMemCellFlags, 0,
						CVector<3, int>(0, 0, 0), size, 0, origin,
						this->domain_cells, size);
			}
		}
	}

private:
	void debugChar(CCL::CMem &cMem, size_t wrap_size = 20) {
		size_t char_size = cMem.getSize();
		char *buffer = new char[char_size];

		cCommandQueue.enqueueReadBuffer(cMem, CL_TRUE, // sync reading
				0, char_size, buffer);

		for (size_t i = 0; i < char_size; i++) {
			if (i % wrap_size == 0) {
				std::cout << std::endl;
				std::cout << (i / wrap_size) << ": ";
			}

			std::cout << (int) buffer[i] << " ";
		}

		delete[] buffer;
	}

	void debugFloat(CCL::CMem &cMem, size_t wrap_size = 20) {
		size_t byte_size = cMem.getSize();
		size_t T_size = byte_size / sizeof(T);

		T *buffer = new T[T_size];

		cCommandQueue.enqueueReadBuffer(cMem, CL_TRUE, // sync reading
				0, byte_size, buffer);

		std::streamsize ss = std::cout.precision();
		std::cout.precision(4);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);

		for (size_t i = 0; i < T_size; i++) {
			if (i % wrap_size == 0) {
				std::cout << std::endl;
				std::cout << (i / wrap_size) << ": ";
			}

			std::cout << buffer[i] << " ";
		}

		std::cout.precision(ss);
		std::cout << std::resetiosflags(std::ios::fixed);

		delete[] buffer;
	}

	/**
	 * debug handler to output some useful information
	 */
public:
	void debug_print() {
		// read out DensityDistributions
		std::cout << "DENSITY DISTRIBUTIONS:";
		//debugFloat(cMemDensityDistributions, SIZE_DD_HOST);
		debugFloat(cMemDensityDistributions, 16);
		std::cout << std::endl;

		// read out Velocity
		std::cout << std::endl;
		std::cout << "VELOCITY:";
		debugFloat(cMemVelocity, 4 * 3);
		std::cout << std::endl;

		// read out VelocityDensity
		std::cout << std::endl;
		std::cout << "DENSITY:";
		debugFloat(cMemDensity, 4);
		std::cout << std::endl;

		// read out Flags
		std::cout << std::endl;
		std::cout << "FLAGS:";
		debugChar(cMemCellFlags, 4 * 4);
		std::cout << std::endl;
	}

	/**
	 * debug density distributions (works only in linear mode!)
	 * \param	dd_id	specifies the dd number
	 */
	void debugDD(size_t dd_id = 0, size_t wrap_size = 16,
			size_t empty_line = 16) {
		CCL::CMem &cMem = cMemDensityDistributions;

		size_t byte_size = cMem.getSize();
		size_t T_size = byte_size / sizeof(T);

		T *buffer = new T[T_size];

		cCommandQueue.enqueueReadBuffer(cMem, CL_TRUE, // sync reading
				0, byte_size, buffer);

		std::streamsize ss = std::cout.precision();
		std::cout.precision(4);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);

		size_t start = this->domain_cells.elements() * dd_id;
		size_t end = this->domain_cells.elements() * (dd_id + 1);

		for (size_t i = start; i < end; i++) {
			if (empty_line != wrap_size)
				if (i % empty_line == 0 && i != start) {
					std::cout << std::endl;
				}

			if (i % wrap_size == 0) {
				if (i != start)
					std::cout << std::endl;
				std::cout << (i / wrap_size) << ": ";
			}

			std::cout << buffer[i] << " ";
		}

		std::cout.precision(ss);
		std::cout << std::resetiosflags(std::ios::fixed);
		std::cout << std::endl;

		delete[] buffer;
	}

	float getVelocityChecksum() {
		T *velocity = new T[this->domain_cells.elements() * 3];
		int *flags = new int[this->domain_cells.elements() * 3];

		storeVelocity(velocity);
		storeFlags(flags);

		T *velx = velocity;
		T *vely = velocity + this->domain_cells.elements();
		T *velz = velocity + this->domain_cells.elements() * 2;

		float checksum = 0;
		for (int a = 0; a < this->domain_cells.elements(); a++)
			if (flags[a] == LBM_FLAG_FLUID)
				checksum += velx[a] + vely[a] + velz[a];

		delete[] flags;
		delete[] velocity;

		return checksum;
	}
};

#endif
