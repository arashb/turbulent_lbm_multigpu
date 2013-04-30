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


#ifndef CLATTICE_BOLTZMANN_HPP
#define CLATTICE_BOLTZMANN_HPP

#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include "libcl/CCL.hpp"
#include "libtools/CStopwatch.hpp"
#include "CLbmOpenCl.hpp"
#include "libvis/ILbmVisualization.hpp"
#include "libvis/CLbmVisualizationVTK.hpp"

#include <list>

// simulation type
typedef float T;
//typedef double T;

template <typename T>
class CMain
{
	bool reset;
	ILbmVisualization<T>* cLbmVisualization;
	int next_simulation_steps_count;
	CLbmOpenCl<T> *cLbmPtr;

	T vector_checksum;

	// amount of vectors to omit while visualization of velicities
	//int visualization_increment;

	// true if debug mode is active
	bool debug_mode;

	void outputDD(int dd_i)
	{
		std::cout << "DD " << dd_i << std::endl;
		//					int gcd = CMath<int>::gcd(cLbmPtr->domain_cells[0],wrap_max_line);
		//					if (wrap_max_line % gcd == 0)
		//						gcd = wrap_max_line;
		int wrap_max_line = 16;
		int gcd = wrap_max_line;
		cLbmPtr->debugDD(dd_i, gcd, cLbmPtr->domain_cells[0]*cLbmPtr->domain_cells[1]);
	}

public:

	CMain()	:
		next_simulation_steps_count(-1),
		cLbmVisualization(NULL)
	{
	}

	int run(
			bool p_debug_mode,
			CVector<3,int> domain_size,
			CVector<3,T> gravitation,
			T viscosity,
			size_t computation_kernel_count,
			int device_nr,
			bool do_visualization,
			bool pause,
			T timestep,
			bool take_frame_screenshots,
			int loops,

			std::list<int> &p_lbm_opencl_number_of_work_items_list,		///< list with number of threads for each successively created kernel
			std::list<int> &p_lbm_opencl_number_of_registers_list		///< list with number of registers for each thread threads for each successively created kernel
	)		
	{
		if (loops < 0)
			loops = 100;

		debug_mode = p_debug_mode;

		if (debug_mode)
			std::cout << "domain size: " << domain_size << std::endl;

		if (pause)
			next_simulation_steps_count = 0;

		// load platform information
		if (debug_mode)	std::cout << "loading platforms" << std::endl;
		CCL::CPlatforms cPlatforms;
		cPlatforms.load();

		if (cPlatforms.platform_ids_count == 0)
		{
			std::cerr << "no platform found!" << std::endl;
			return -1;
		}

		int platform_id_nr = -1;
		for (size_t i = 0; i < cPlatforms.platform_ids_count; i++)
		{
			CCL::CPlatform cPlatform(cPlatforms.platform_ids[i]);
			cPlatform.loadPlatformInfo();

			if (platform_id_nr == -1)
				if (strcmp(cPlatform.profile, "FULL_PROFILE") == 0)
				{
					platform_id_nr = i;
					if (debug_mode)
						std::cout << "Using Platform " << (i+1) << " for computation" << std::endl;
				}

			if (debug_mode)
			{
				std::cout << "Platform " << (i) << ":" << std::endl;
				std::cout << "        Name: " << cPlatform.name << std::endl;
				std::cout << "     Profile: " << cPlatform.profile << std::endl;
				std::cout << "     Version: " << cPlatform.version << std::endl;
				std::cout << "      Vendor: " << cPlatform.vendor << std::endl;
				std::cout << "  Extensions: " << cPlatform.extensions << std::endl;
				std::cout << std::endl;
			}
		}

		if (platform_id_nr == -1)
		{
			std::cout << "no usable platform found" << std::endl;
			return -1;
		}

		CCL::CPlatform cPlatform(cPlatforms.platform_ids[platform_id_nr]);

		// load standard context for GPU devices
		if (debug_mode)	std::cout << "loading gpu context" << std::endl;
		CCL::CContext cContext(cPlatform, CL_DEVICE_TYPE_GPU);

		// load devices belonging to cContext
		if (debug_mode)	std::cout << "loading devices" << std::endl;
		CCL::CDevices cDevices(cContext);

		if (cDevices.size() == 0)
		{
			std::cerr << "no device found - aborting" << std::endl;
			return -1;
		}

		if (device_nr == -1)
		{
			// list available devices
			for (int i = 0; i < (int)cDevices.size(); i++)
			{
				CCL::CDeviceInfo cDeviceInfo(cDevices[i]);
				std::cout << "Device " << (i) << ":" << std::endl;
				std::cout << "        Name: " << cDeviceInfo.name << std::endl;
				std::cout << "     Profile: " << cDeviceInfo.profile << std::endl;
				std::cout << "     Version: " << cDeviceInfo.version << std::endl;
				std::cout << "      Vendor: " << cDeviceInfo.vendor << std::endl;
				std::cout << "  Extensions: " << cDeviceInfo.extensions << std::endl;
				std::cout << std::endl;
			}
			return -1;
		}

		if (device_nr < 0 || device_nr >= (int)cDevices.size())
		{
			std::cerr << "invalid device number - use option \"-d -1\" to list all devices" << std::endl;
			return -1;
		}

		CCL::CDevice &cDevice = cDevices[device_nr];

		// load information about first device - e.g. max_work_group_size
		if (debug_mode)	std::cout << "loading device information" << std::endl;
		CCL::CDeviceInfo cDeviceInfo(cDevice);

		if (debug_mode)
		{
			std::cout << "Device " << (device_nr) << ":" << std::endl;
			std::cout << "        Name: " << cDeviceInfo.name << std::endl;
			std::cout << "     Profile: " << cDeviceInfo.profile << std::endl;
			std::cout << "     Version: " << cDeviceInfo.version << std::endl;
			std::cout << "      Vendor: " << cDeviceInfo.vendor << std::endl;
			std::cout << "  Extensions: " << cDeviceInfo.extensions << std::endl;
			std::cout << std::endl;
		}

		// initialize queue
		if (debug_mode)	std::cout << "creating command queue" << std::endl;
		CCL::CCommandQueue cCommandQueue(cContext, cDevice);

		T domain_length = 0.05;


		vector_checksum = 0;

		// approximate bandwidth
		double floats_per_cell = 0.0;

		// 19 density distribution which are read and written
		floats_per_cell += 19.0*2.0;

		// flag (obstacle, injection and fluid) is read
		floats_per_cell += 1.0;

		// velocity vector is also stored
		if (do_visualization || debug_mode)
			floats_per_cell += 3;

		reset = true;

		// INIT LATTICE BOLTZMANN!
		CLbmOpenCl<T> cLbm(	cCommandQueue, cContext, cDevice,
				domain_size,		// domain size
				domain_length,		// length of domain size in x direction
				gravitation,	// gravitation vector
				viscosity,
				computation_kernel_count,
				debug_mode,
				do_visualization || debug_mode,
				false,
				timestep,

				p_lbm_opencl_number_of_work_items_list,
				p_lbm_opencl_number_of_registers_list
		);

		if (cLbm.error())
		{
			std::cout << cLbm.error.getString();
			return -1;
		}

		cLbm.wait();


		cLbmPtr = &cLbm;

		CStopwatch cStopwatch;

		// setting up the visualization
		std::string outputfilename = "./vtkOutput/OUTPUT";
		if (do_visualization)
		{
			cLbmVisualization = new CLbmVisualizationVTK<T>(outputfilename);
			cLbmVisualization->setup(cLbm);
		}

		for (int i = 0; i < loops/10; i++)
		{
			// simulation
			cLbm.simulationStep();
			std::cout << "." << std::flush;
//			if (do_visualization)
//				cLbmVisualization->render();

		}
		std::cout << "|" << std::flush;
		cLbm.wait();

		cStopwatch.start();
		for (int i = 0; i < loops; i++)
		{
			// simulation
			cLbm.simulationStep();
			std::cout << "." << std::flush;
			if (do_visualization)
				cLbmVisualization->render(i);
		}
		cLbm.wait();
		cStopwatch.stop();

		std::cout << std::endl;

		if (domain_size.elements() <= 512)
			if (debug_mode)
				cLbm.debug_print();

		std::cout << std::endl;

		std::cout << "Cube: " << domain_size << std::endl;
		std::cout << "Seconds: " << cStopwatch.time << std::endl;
		double fps = (((double)loops) / cStopwatch.time);
		std::cout << "FPS: " << fps << std::endl;

		double mlups = ((double)fps*(double)cLbm.domain_cells.elements())*(double)0.000001;
		std::cout << "MLUPS: " << mlups << std::endl;

		std::cout << "Bandwidth: " << (mlups*floats_per_cell*(double)sizeof(T)) << " MB/s (RW, bidirectional)" << std::endl;

		std::streamsize ss = std::cout.precision();
		std::cout.precision(8);
		std::cout.setf(std::ios::fixed,std::ios::floatfield);

		if (debug_mode)
		{
			// The velocity checksum is only stored in debug mode!
			vector_checksum = cLbm.getVelocityChecksum();
			std::cout << "Checksum: " << (vector_checksum*1000.0f) << std::endl;
		}

		std::cout.precision(ss);
		std::cout << std::resetiosflags(std::ios::fixed);

		std::cout << std::endl;
		std::cout << "exit" << std::endl;

		return EXIT_SUCCESS;
	}
};

#endif
