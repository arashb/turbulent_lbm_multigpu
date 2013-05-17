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

#include "CDomain.hpp"
#include "CLbmSolver.hpp"
#include "libvis/ILbmVisualization.hpp"
#include "libvis/CLbmVisualizationVTK.hpp"
#include "common.h"

#include <list>

// simulation type
typedef float T;
//typedef double T;

#define MPI_TAG_ALPHA_SYNC 0
#define MPI_TAG_BETA_SYNC 1

/*
 * Class CConroller is responsible for controlling and managing of simulation and visualization
 * of a subdomain from the whole grid.
 *
 */
template <typename T>
class CController
{
	int _UID;						///< Unique ID of each controller
	CDomain<T> _domain;				///< Domain data
	ILbmVisualization<T>* cLbmVisualization; ///< Visualization class
	CLbmSolver<T> *cLbmPtr;
	int _BC[3][2]; 		///< Boundary conditions. First index specifies the dimension and second the upper or the lower boundary.

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

	CController(int UID, CDomain<T> domain)	:
		cLbmVisualization(NULL),
		_UID(UID),
		_domain(domain)
	{
		// default boundary conditions is obstacle
		for(int i = 0; i < 3; i++)
			for (int j = 0; j < 2; j++)
				_BC[i][j] = FLAG_OBSTACLE;
	}

	void syncAlpha() {

		// TODO: the data related to communication is hardcoded here.
		// This should change in a way to get this data from a Communication Type Object
		CVector<3,int> send_size(1,_domain.getSize()[1],_domain.getSize()[2]);
		CVector<3,int> recv_size(1,_domain.getSize()[1],_domain.getSize()[2]);
		CVector<3,int> send_origin;
		CVector<3,int> recv_origin;
		if (_UID == 0) {
			send_origin[0] = _domain.getSize()[0] - 2;
			send_origin[1] = 0;
			send_origin[2] = 0;
			recv_origin[0] = _domain.getSize()[0] - 1;
			recv_origin[1] = 0;
			recv_origin[2] = 0;
		}else if ( _UID == 1) {
			send_origin[0] = 1;
			send_origin[1] = 0;
			send_origin[2] = 0;
			recv_origin[0] = 0;
			recv_origin[1] = 0;
			recv_origin[2] = 0;
		}

		// send buffer
		int send_buffer_size = send_size.elements()*cLbmPtr->SIZE_DD_HOST;
		int recv_buffer_size = recv_size.elements()*cLbmPtr->SIZE_DD_HOST;
		T* send_buffer = new T[send_buffer_size];
		T* recv_buffer = new T[recv_buffer_size];

		MPI_Request req[2];
		MPI_Status status[2];

		// Download data from device to host
		cLbmPtr->storeDensityDistribution(send_buffer, send_origin, send_size);
		//cLbmPtr->wait();
		int my_rank, num_procs;
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
		MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processes

		int dst_rank;
		if (my_rank == 0) {
			dst_rank = 1;
		} else if ( my_rank == 1)
			dst_rank = 0;

		// TODO: check the MPI_TYPE
		MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[0]);
		MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[1]);
		MPI_Waitall(2, req, status );

		cLbmPtr->setDensityDistribution(recv_buffer, recv_origin, recv_size);
		cLbmPtr->wait();

		delete send_buffer;
		delete recv_buffer;
	}

	void syncBeta() {

		// TODO: the data related to communication is hardcoded here.
		// This should change in a way to get this data from a Communication Type Object
		CVector<3,int> send_size(1,_domain.getSize()[1],_domain.getSize()[2]);
		CVector<3,int> recv_size(1,_domain.getSize()[1],_domain.getSize()[2]);
		CVector<3,int> send_origin;
		CVector<3,int> recv_origin;
		if (_UID == 0) {
			send_origin[0] = _domain.getSize()[0] - 2;
			send_origin[1] = 0;
			send_origin[2] = 0;
			recv_origin[0] = _domain.getSize()[0] - 1;
			recv_origin[1] = 0;
			recv_origin[2] = 0;
		}else if ( _UID == 1) {
			send_origin[0] = 1;
			send_origin[1] = 0;
			send_origin[2] = 0;
			recv_origin[0] = 0;
			recv_origin[1] = 0;
			recv_origin[2] = 0;
		}

		// send buffer
		int send_buffer_size = send_size.elements()*cLbmPtr->SIZE_DD_HOST;
		int recv_buffer_size = recv_size.elements()*cLbmPtr->SIZE_DD_HOST;
		T* send_buffer = new T[send_buffer_size];
		T* recv_buffer = new T[recv_buffer_size];

		MPI_Request req[2];
		MPI_Status status[2];

		// Download data from device to host
		cLbmPtr->storeDensityDistribution(send_buffer, send_origin, send_size);
		//cLbmPtr->wait();
		int my_rank, num_procs;
		MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
		MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processes

		int dst_rank;
		if (my_rank == 0) {
			dst_rank = 1;
		} else if ( my_rank == 1)
			dst_rank = 0;

		// TODO: check the MPI_TYPE
		MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[0]);
		MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[1]);
		MPI_Waitall(2, req, status );

		// TODO: OPTI: you need to wait only for receiving to execute following command
		cLbmPtr->setDensityDistribution(recv_buffer, recv_origin, recv_size);
		cLbmPtr->wait();

		delete send_buffer;
		delete recv_buffer;
	}

	void computeNextStep(){
		cLbmPtr->simulationStep();
		if (cLbmPtr->simulation_step_counter & 1)
			syncBeta();
		else
			syncAlpha();
	}
/*
 * This function starts the simulation for the particular subdomain corresponded to
 * this class.
 */
	int run(
			bool p_debug_mode, 				///< Set this variable to true to have a verbose output of simulation process.
			//CDomain<T> domain, 				///< Specify domain properties
			CVector<3,T> gravitation,		///< Specify the gravitation vector
			T viscosity,
			size_t computation_kernel_count,
			int device_nr,
			bool do_visualization,
			bool pause,
			T timestep,
			bool take_frame_screenshots,
			int loops,

			std::list<int> &p_lbm_opencl_number_of_work_items_list,		///< List with number of threads for each successively created kernel
			std::list<int> &p_lbm_opencl_number_of_registers_list		///< List with number of registers for each thread threads for each successively created kernel
	)		
	{
		CVector<3,int> domain_size = _domain.getSize();
		T domain_length = _domain.getLength()[0];

		if (loops < 0)
			loops = 100;

		debug_mode = p_debug_mode;

		if (debug_mode)
			std::cout << "domain size: " << domain_size << std::endl;

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

		// INIT LATTICE BOLTZMANN!
		CLbmSolver<T> cLbm(	cCommandQueue, cContext, cDevice,
				_BC,
				_domain,
				gravitation,	// gravitation vector
				viscosity,
				computation_kernel_count,
				debug_mode,
				do_visualization || debug_mode,
				do_visualization || debug_mode,
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
			cLbmVisualization = new CLbmVisualizationVTK<T>(_UID,outputfilename);
			cLbmVisualization->setup(cLbm);
		}

//		for (int i = 0; i < loops/10; i++)
//		{
//			// simulation
//			computeNextStep();
//			std::cout << "." << std::flush;
//			if (do_visualization)
//				cLbmVisualization->render(i);
//
//		}
//		std::cout << "|" << std::flush;
//		cLbm.wait();

		cStopwatch.start();
		for (int i = 0; i < loops; i++)
		{
			// simulation
			computeNextStep();
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

	void setBC(int BC[3][2]) {
		for(int i = 0; i < 3; i++)
			for (int j = 0; j < 2; j++)
				_BC[i][j] = BC[i][j];
	}
};

#endif
