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

#include "mpi.h"
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <string>
#include <sstream>
#include <fstream>

#include "libcl/CCL.hpp"
#include "libtools/CStopwatch.hpp"

#include "CDomain.hpp"
#include "CLbmSolver.hpp"
#include "libvis/ILbmVisualization.hpp"
#include "libvis/CLbmVisualizationVTK.hpp"
#include "CComm.hpp"
#include "common.h"
#include "CConfiguration.hpp"
#include "Singleton.hpp"

#include <list>
#include <map>
#include <vector>
// simulation type
typedef float T;
//typedef double T;

#define MPI_TAG_ALPHA_SYNC 0
#define MPI_TAG_BETA_SYNC 1

typedef std::multimap<MPI_COMM_DIRECTION, CComm<T>* > comm_map;
typedef std::multimap<MPI_COMM_DIRECTION, CComm<T>* >::value_type comm_map_pair;
typedef std::multimap<MPI_COMM_DIRECTION, CComm<T>* >::iterator comm_map_ptr;

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
	comm_map _comm_container; ///< A std::multimap containing all the communication objects for the subdomain

	T vector_checksum;
	CCL::CPlatforms* cPlatforms;
	CCL::CPlatform* cPlatform;
	CCL::CContext* cContext;
	CCL::CDevices* cDevices;
	CCL::CDevice* cDevice;
	CCL::CCommandQueue* cCommandQueue;
	size_t _simulation_step_counter;

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

	int initLBMSolver() {

#if DEBUG
			//std::cout << "domain size: " << _domain.getSize() << std::endl;
		// load platform information
			std::cout << "loading platforms" << std::endl;
#endif
		cPlatforms = new CCL::CPlatforms();
		cPlatforms->load();



		if (cPlatforms->platform_ids_count == 0)
		{
			std::cerr << "no platform found!" << std::endl;
			return -1;
		}

		int platform_id_nr = -1;
		for (size_t i = 0; i < cPlatforms->platform_ids_count; i++)
		{
			CCL::CPlatform cPlatform(cPlatforms->platform_ids[i]);
			cPlatform.loadPlatformInfo();

			if (platform_id_nr == -1)
				if (strcmp(cPlatform.profile, "FULL_PROFILE") == 0)
				{
					platform_id_nr = i;
#if DEBUG
						std::cout << "Using Platform " << (i+1) << " for computation" << std::endl;
#endif
				}

#if DEBUG
				std::cout << "Platform " << (i) << ":" << std::endl;
				std::cout << "        Name: " << cPlatform.name << std::endl;
				std::cout << "     Profile: " << cPlatform.profile << std::endl;
				std::cout << "     Version: " << cPlatform.version << std::endl;
				std::cout << "      Vendor: " << cPlatform.vendor << std::endl;
				std::cout << "  Extensions: " << cPlatform.extensions << std::endl;
				std::cout << std::endl;
#endif
		}

		if (platform_id_nr == -1)
		{
			std::cout << "no usable platform found" << std::endl;
			return -1;
		}

		cPlatform = new CCL::CPlatform(cPlatforms->platform_ids[platform_id_nr]);

		// load standard context for GPU devices
#if DEBUG
			std::cout << "loading gpu context" << std::endl;
#endif

		cContext = new CCL::CContext(*cPlatform, CL_DEVICE_TYPE_GPU);

		// load devices belonging to cContext
#if DEBUG
			std::cout << "loading devices" << std::endl;
#endif
		cDevices = new CCL::CDevices(*cContext);

		if (cDevices->size() == 0)
		{
			std::cerr << "no device found - aborting" << std::endl;
			return -1;
		}

		if (ConfigSingleton::Instance()->device_nr == -1)
		{
			// list available devices
			for (int i = 0; i < (int)cDevices->size(); i++)
			{
				CCL::CDeviceInfo cDeviceInfo((*cDevices)[i]);
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

		if (ConfigSingleton::Instance()->device_nr < 0 || ConfigSingleton::Instance()->device_nr >= (int)cDevices->size())
		{
			std::cerr << "invalid device number - use option \"-d -1\" to list all devices" << std::endl;
			return -1;
		}

		cDevice = &((*cDevices)[ConfigSingleton::Instance()->device_nr]);
		//cDeviceInfo = new CCL::CDeviceInfo(*cDevice);

		// load information about first device - e.g. max_work_group_size
#if DEBUG
		std::cout << "loading device information" << std::endl;
		CCL::CDeviceInfo cDeviceInfo(*cDevice);

		std::cout << "Device " << (ConfigSingleton::Instance()->device_nr) << ":" << std::endl;
		std::cout << "        Name: " << cDeviceInfo.name << std::endl;
		std::cout << "     Profile: " << cDeviceInfo.profile << std::endl;
		std::cout << "     Version: " << cDeviceInfo.version << std::endl;
		std::cout << "      Vendor: " << cDeviceInfo.vendor << std::endl;
		std::cout << "  Extensions: " << cDeviceInfo.extensions << std::endl;
		std::cout << std::endl;

		std::cout << "creating command queue" << std::endl;
#endif
		// initialize queue
		cCommandQueue = new CCL::CCommandQueue(*cContext, *cDevice);

		//CVector<3,int> halo_domain(_domain[0]+2,_domain[1]+2,_domain[2]+2);
		// INIT LATTICE BOLTZMANN!
		cLbmPtr = new CLbmSolver<T>(	_UID, *cCommandQueue, *cContext, *cDevice,
				_BC,
				_domain,
				ConfigSingleton::Instance()->gravitation,	// gravitation vector
				ConfigSingleton::Instance()->viscosity,
				ConfigSingleton::Instance()->computation_kernel_count,
				//ConfigSingleton::Instance()->debug_mode,
				ConfigSingleton::Instance()->do_visualization || ConfigSingleton::Instance()->debug_mode,
				ConfigSingleton::Instance()->do_visualization || ConfigSingleton::Instance()->debug_mode,
				ConfigSingleton::Instance()->timestep,
				ConfigSingleton::Instance()->lbm_opencl_number_of_threads_list, // p_lbm_opencl_number_of_work_items_list,
				ConfigSingleton::Instance()->lbm_opencl_number_of_threads_list
		);



		if (cLbmPtr->error())
		{
			std::cout << cLbmPtr->error.getString();
			return -1;
		}

		cLbmPtr->wait();
		CStopwatch cStopwatch;
		return 0;
	}

public:

	CController(int UID, CDomain<T> domain, int BC[3][2])	:
		_UID(UID),
		_domain(domain),
		cLbmVisualization(NULL)
	{
		_simulation_step_counter = 0;
		for(int i = 0; i < 3; i++)
			for (int j = 0; j < 2; j++)
				_BC[i][j] = BC[i][j];

		// initialize the LBMSolver
		if ( -1 == initLBMSolver() )
			throw "Initialization of LBM Solver failed!";
	}
	~CController(){
		if (cPlatforms)
			delete cPlatform;

		if (cContext)
			delete cContext;

		if (cDevices)
			delete cDevices;

		if ( cCommandQueue )
			delete cCommandQueue;


		if (cLbmVisualization)
			delete cLbmVisualization;

		if (cLbmPtr)
			delete cLbmPtr;

		comm_map_ptr it = _comm_container.begin();
		for( ;it != _comm_container.end(); it++){
			delete (*it).second;
		}
	}

	void storeDataAlpha(MPI_COMM_DIRECTION direction) {
		std::pair<comm_map_ptr, comm_map_ptr> ppp;
		ppp = _comm_container.equal_range(direction);
		comm_map_ptr it = ppp.first;
		for( ;it != ppp.second; it++){
			// Store data from device to host
			cLbmPtr->storeDensityDistribution(	(*it).second->getSendBuffer(),
					(*it).second->getSendOrigin(),
					(*it).second->getSendSize());
		}
	}

	void setDataAlpha(MPI_COMM_DIRECTION direction) {
		std::pair<comm_map_ptr, comm_map_ptr> ppp;
		ppp = _comm_container.equal_range(direction);
		comm_map_ptr it = ppp.first;
		for( ;it != ppp.second; it++){
			// Store data from host to device
			cLbmPtr->setDensityDistribution((*it).second->getRecvBuffer(),
					(*it).second->getRecvOrigin(),
					(*it).second->getRecvSize());
		}
		cLbmPtr->wait();
	}

	void storeDataBeta(MPI_COMM_DIRECTION direction) {
		std::pair<comm_map_ptr, comm_map_ptr> ppp;
		ppp = _comm_container.equal_range(direction);
		comm_map_ptr it = ppp.first;
		for( ;it != ppp.second; it++){
			// Store data from device to host
			// the send and receive origin and size values in beta sync is the opposite values of
			// CComm object of current communication, since the ghost layer data
			// will be sent back to their origin
			cLbmPtr->storeDensityDistribution(	(*it).second->getSendBuffer(),
					(*it).second->getRecvOrigin(),
					(*it).second->getRecvSize());
		}
	}

	void setDataBeta(MPI_COMM_DIRECTION direction) {
		std::pair<comm_map_ptr, comm_map_ptr> ppp;
		ppp = _comm_container.equal_range(direction);
		comm_map_ptr it = ppp.first;
		for( ;it != ppp.second; it++){
			// Store data from host to device
			// the send and receive origin and size values in beta sync is the opposite values of
			// CComm object of current communication, since the ghost layer data
			// will be sent back to their origin
			cLbmPtr->setDensityDistribution((*it).second->getRecvBuffer(),
					(*it).second->getSendOrigin(),
					(*it).second->getSendSize(),
					(*it).second->getCommDirection());
		}
		cLbmPtr->wait();
	}

	void syncAlpha() {
#if DEBUG
			std::cout << "--> Sync alpha" << std::endl;
#endif
		// TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking way.
		comm_map_ptr it = _comm_container.begin();
		for( ;it != _comm_container.end(); it++){
			int dst_rank = (*it).second->getDstId();
			int send_buffer_size = (*it).second->getSendBufferSize();
			int recv_buffer_size = (*it).second->getRecvBufferSize();
			T* send_buffer = (*it).second->getSendBuffer();
			T* recv_buffer = (*it).second->getRecvBuffer();

			MPI_Request req[2];
			MPI_Status status[2];
			int my_rank, num_procs;
			MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
			MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processes

			if (typeid(T) == typeid(float)){
				MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[1]);
			}else if (typeid(T) == typeid(double)) {
				MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else {
				throw "Type id of MPI send/receive buffer is unknown!";
			}
			MPI_Waitall(2, req, status );
		}
	}

	void syncAlpha(MPI_COMM_DIRECTION direction) {
#if DEBUG
			std::cout << "--> Sync alpha: " << direction <<  std::endl;
#endif
		// TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking way.
		// iterating over all communication objects for a specific direction
		std::pair<comm_map_ptr, comm_map_ptr> ppp;
		ppp = _comm_container.equal_range(direction);
		comm_map_ptr it = ppp.first;
		for( ;it != ppp.second; it++){
			int dst_rank = (*it).second->getDstId();
			int send_buffer_size = (*it).second->getSendBufferSize();
			int recv_buffer_size = (*it).second->getRecvBufferSize();
			T* send_buffer = (*it).second->getSendBuffer();
			T* recv_buffer = (*it).second->getRecvBuffer();

			MPI_Request req[2];
			MPI_Status status[2];
			if (typeid(T) == typeid(float)){
				MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[1]);
			}else if (typeid(T) == typeid(double)) {
				MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else {
				throw "Type id of MPI send/receive buffer is unknown!";
			}
			MPI_Waitall(2, req, status );
		}
	}

	void syncBeta() {
#if DEBUG
			std::cout << "--> Sync beta" << std::endl;
#endif
		// TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking form.
		comm_map_ptr it = _comm_container.begin();
		for( ;it != _comm_container.end(); it++) {
			int dst_rank = (*it).second->getDstId();
			int send_buffer_size = (*it).second->getSendBufferSize();
			int recv_buffer_size = (*it).second->getRecvBufferSize();
			T* send_buffer = (*it).second->getSendBuffer();
			T* recv_buffer = (*it).second->getRecvBuffer();

			MPI_Request req[2];
			MPI_Status status[2];
			if (typeid(T) == typeid(float)){
				MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else if (typeid(T) == typeid(double)) {
				MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else {
				throw "Type id of MPI send/receive buffer is unknown!";
			}
			MPI_Waitall(2, req, status );
		}
	}

	void syncBeta(MPI_COMM_DIRECTION direction) {
#if DEBUG
		std::cout << "--> Sync beta: " << direction <<  std::endl;
#endif
		// TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking form.
		// iterating over all communication objects for a specific direction
		std::pair<comm_map_ptr, comm_map_ptr> ppp;
		ppp = _comm_container.equal_range(direction);
		comm_map_ptr it = ppp.first;
		for( ;it != ppp.second; it++) {
			int dst_rank = (*it).second->getDstId();
			int send_buffer_size = (*it).second->getSendBufferSize();//send_size.elements()*cLbmPtr->SIZE_DD_HOST;
			int recv_buffer_size = (*it).second->getRecvBufferSize();
			T* send_buffer = (*it).second->getSendBuffer();
			T* recv_buffer = (*it).second->getRecvBuffer();

			MPI_Request req[2];
			MPI_Status status[2];
			if (typeid(T) == typeid(float)){
				MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else if (typeid(T) == typeid(double)) {
				MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else {
				throw "Type id of MPI send/receive buffer is unknown!";
			}
			MPI_Waitall(2, req, status );
		}
	}

	void simulationStepAlpha() {
		// TODO: Optimization: some cells are computed twice. Change the size of computations in each direction to avoid it.
		// SIMULATION_STEP_ALPHA
		CVector<3,int> x_size(1						, _domain.getSize()[1]	, _domain.getSize()[2]	);
		CVector<3,int> y_size(_domain.getSize()[0]	, 1						, _domain.getSize()[2]	);
		CVector<3,int> z_size(_domain.getSize()[0]	, _domain.getSize()[1]	, 1						);

		// --> Simulation step alpha x boundary
		CVector<3,int> x0_origin(1, 0, 0);
		CVector<3,int> x1_origin(_domain.getSize()[0]-2, 0, 0);
		cLbmPtr->simulationStepAlphaRect(x0_origin, x_size);
		cLbmPtr->simulationStepAlphaRect(x1_origin, x_size);

		// --> Store x boundary
		storeDataAlpha(MPI_COMM_DIRECTION_X);
		// --> Communication x boundary
		syncAlpha(MPI_COMM_DIRECTION_X);
		// --> Set x boundary
		setDataAlpha(MPI_COMM_DIRECTION_X);

		// --> Simulation step alpha y boundary
		CVector<3,int> y0_origin(0, 1, 0);
		CVector<3,int> y1_origin(0,_domain.getSize()[1]-2, 0);
		cLbmPtr->simulationStepAlphaRect(y0_origin, y_size);
		cLbmPtr->simulationStepAlphaRect(y1_origin, y_size);

		// --> Store y boundary
		storeDataAlpha(MPI_COMM_DIRECTION_Y);
		// --> Communication y boundary
		syncAlpha(MPI_COMM_DIRECTION_Y);
		// --> Set y boundary
		setDataAlpha(MPI_COMM_DIRECTION_Y);

		// --> Simulation step alpha z boundary
		CVector<3,int> z0_origin(0, 0, 1);
		CVector<3,int> z1_origin(0, 0, _domain.getSize()[2]-2);
		cLbmPtr->simulationStepAlphaRect(z0_origin, z_size);
		cLbmPtr->simulationStepAlphaRect(z1_origin, z_size);

		// --> Store z boundary
		storeDataAlpha(MPI_COMM_DIRECTION_Z);
		// --> Communication z boundary
		syncAlpha(MPI_COMM_DIRECTION_Z);
		// --> Set z boundary
		setDataAlpha(MPI_COMM_DIRECTION_Z);

		// --> Computation of inner part
		CVector<3,int> inner_origin(2, 2, 2);
		CVector<3,int> inner_size(_domain.getSize()[0]	- 4, _domain.getSize()[1] - 4, _domain.getSize()[2] - 4);
		cLbmPtr->simulationStepAlphaRect(inner_origin, inner_size);
	}

	void simulationStepBeta() {
		// TODO: Optimization: some cells are computed twice. Change the size of computations in each direction to avoid it.
		// SIMULATION_STEP_ALPHA
		CVector<3,int> x_size(1						, _domain.getSize()[1]	, _domain.getSize()[2]	);
		CVector<3,int> y_size(_domain.getSize()[0]	- 4, 1						, _domain.getSize()[2]	);
		CVector<3,int> z_size(_domain.getSize()[0]	- 4, _domain.getSize()[1] - 4	, 1						);

		// --> Simulation step alpha x boundary
		CVector<3,int> x0_origin(1, 0, 0);
		CVector<3,int> x1_origin(_domain.getSize()[0]-2, 0, 0);
		cLbmPtr->simulationStepBetaRect(x0_origin, x_size);
		cLbmPtr->simulationStepBetaRect(x1_origin, x_size);

		// --> Store x boundary
		storeDataBeta(MPI_COMM_DIRECTION_X);
		// --> Communication x boundary
		syncBeta(MPI_COMM_DIRECTION_X);
		// --> Set x boundary
		setDataBeta(MPI_COMM_DIRECTION_X);

		// --> Simulation step alpha y boundary
		CVector<3,int> y0_origin(2, 1, 0);
		CVector<3,int> y1_origin(2,_domain.getSize()[1]-2, 0);
		cLbmPtr->simulationStepBetaRect(y0_origin, y_size);
		cLbmPtr->simulationStepBetaRect(y1_origin, y_size);

		// --> Store y boundary
		storeDataBeta(MPI_COMM_DIRECTION_Y);
		// --> Communication y boundary
		syncBeta(MPI_COMM_DIRECTION_Y);
		// --> Set y boundary
		setDataBeta(MPI_COMM_DIRECTION_Y);

		// --> Simulation step alpha z boundary
		CVector<3,int> z0_origin(2, 2, 1);
		CVector<3,int> z1_origin(2, 2, _domain.getSize()[2]-2);
		cLbmPtr->simulationStepBetaRect(z0_origin, z_size);
		cLbmPtr->simulationStepBetaRect(z1_origin, z_size);

		// --> Store z boundary
		storeDataBeta(MPI_COMM_DIRECTION_Z);
		// --> Communication z boundary
		syncBeta(MPI_COMM_DIRECTION_Z);
		// --> Set z boundary
		setDataBeta(MPI_COMM_DIRECTION_Z);

		// --> Computation of inner part
		CVector<3,int> inner_origin(2, 2, 2);
		CVector<3,int> inner_size(_domain.getSize()[0]	- 4, _domain.getSize()[1] - 4, _domain.getSize()[2] - 4);
		cLbmPtr->simulationStepBetaRect(inner_origin, inner_size);
	}

	void computeNextStep(){
		if (_simulation_step_counter & 1)
			simulationStepAlpha();
		else
			simulationStepBeta();
		cCommandQueue->enqueueBarrier();
		_simulation_step_counter++;

	}
/*
 * This function starts the simulation for the particular subdomain corresponded to
 * this class.
 */
	int run()
	{
		CVector<3,int> domain_size = _domain.getSize();
		int loops = ConfigSingleton::Instance()->loops;
		if (loops < 0)
			loops = 100;

		vector_checksum = 0;

		// approximate bandwidth
		double floats_per_cell = 0.0;

		// 19 density distribution which are read and written
		floats_per_cell += 19.0*2.0;

		// flag (obstacle, injection and fluid) is read
		floats_per_cell += 1.0;

		// velocity vector is also stored
		if (ConfigSingleton::Instance()->do_visualization || ConfigSingleton::Instance()->debug_mode)
			floats_per_cell += 3;
		CStopwatch cStopwatch;

		// setting up the visualization
		std::string outputfilename = "OUTPUT";
		std::stringstream ss_file;
		ss_file << "./" << VTK_OUTPUT_DIR << "/" << outputfilename ;
		std::string outputfile = ss_file.str();
		if (ConfigSingleton::Instance()->do_visualization)
		{
			cLbmVisualization = new CLbmVisualizationVTK<T>(_UID,outputfile);
			cLbmVisualization->setup(cLbmPtr);
		}

		cStopwatch.start();
		for (int i = 0; i < loops; i++)
		{
			computeNextStep();
			//simulation
			if (ConfigSingleton::Instance()->do_visualization)
				cLbmVisualization->render(i);
		}
		cLbmPtr->wait();
		cStopwatch.stop();
#if DEBUG
		if (domain_size.elements() <= 512) {
		  cLbmPtr->debug_print();
		}
#endif

#if BENCHMARK
		double ltime = cStopwatch.time;
		double gtime;
		//int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
		MPI_Reduce(&ltime, &gtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		if (_UID == 0) {
		  double gfps = (((double)loops) / gtime);
		  double gmlups = ((double)gfps*(double)ConfigSingleton::Instance()->domain_size.elements())*(double)0.000001;
		  double gbandwidth = (gmlups*floats_per_cell*(double)sizeof(T));
		  std::stringstream benchmark_file_name;
		  benchmark_file_name << "./" << BENCHMARK_OUTPUT_DIR << "/" << 
		    "benchmark_" << ConfigSingleton::Instance()->subdomain_num.elements()  //<< "_" << _UID
				 << ".ini";
		  const std::string& tmp = benchmark_file_name.str();
		  const char* cstr = tmp.c_str();
		  std::ofstream benchmark_file (cstr, std::ios::out | std::ios::app );
		  if (benchmark_file.is_open())
		    {
		      //benchmark_file << "[RESULTS]" << std::endl; 
		      benchmark_file << "CUBE : " << ConfigSingleton::Instance()->domain_size << std::endl;
		      benchmark_file << "SECONDS : " << gtime << std::endl;
		      benchmark_file << "FPS : " << gfps << std::endl;
		      benchmark_file << "MLUPS : " << gmlups << std::endl;
		      benchmark_file << "BANDWIDTH : " << gbandwidth // << " MB/s (RW, bidirectional)"
				<< std::endl;
		      benchmark_file << std::endl;
		
		    }
		  else std::cout << "Unable to open file";
		}
#endif // end of BENCHMARK
		std::cout << std::endl;
		std::cout << "Cube: " << domain_size << std::endl;
		std::cout << "Seconds: " << cStopwatch.time << std::endl;
		double fps = (((double)loops) / cStopwatch.time);
		std::cout << "FPS: " << fps << std::endl;
		double mlups = ((double)fps*(double)cLbmPtr->domain_cells.elements())*(double)0.000001;
		std::cout << "MLUPS: " << mlups << std::endl;
		std::cout << "Bandwidth: " << (mlups*floats_per_cell*(double)sizeof(T)) << " MB/s (RW, bidirectional)" << std::endl;
		std::streamsize ss = std::cout.precision();
		std::cout.precision(8);
		std::cout.setf(std::ios::fixed,std::ios::floatfield);
#if DEBUG
		// The velocity checksum is only stored in debug mode!
		vector_checksum = cLbmPtr->getVelocityChecksum();
		std::cout << "Checksum: " << (vector_checksum*1000.0f) << std::endl;
#endif // end of DEBUG
		std::cout.precision(ss);
		std::cout << std::resetiosflags(std::ios::fixed);

		std::cout << "done." << std::endl;

#if PROFILE
		std::stringstream profile_file_name;
		profile_file_name << "./" << PROFILE_OUTPUT_DIR << "/" << 
		  "profile_" << ConfigSingleton::Instance()->subdomain_num.elements()  << "_" << _UID
				    << ".ini";
		const std::string& tmp = profile_file_name.str();
		const char* pcstr = tmp.c_str();
		std::ofstream prof_file (pcstr, std::ios::out | std::ios::app );
		if (prof_file.is_open()) {
		  prof_file << "[METADATA]" << std::endl;
		  prof_file << "TOTAL_NUM_PROC : " << ConfigSingleton::Instance()->subdomain_num.elements() << std::endl;
		  prof_file << "CURRENT_PROC_ID : " << _UID << std::endl;
		  prof_file << std::endl;
		} else std::cout << "Unable to open file: " << pcstr << std::endl;
		// const std::string& tmp = profile_file_name.str();
		// const char* cstr = tmp.c_str();
		ProfilerSingleton::Instance()->saveEvents(profile_file_name.str());
#endif
		return EXIT_SUCCESS;
	}

	void addCommunication( MPI_COMM_DIRECTION key, CComm<T>* comm) {
		_comm_container.insert(comm_map_pair(key, comm));
	}

	/*
	 * This Function is used the set the geometry (e.g obstacles, velocity injections, ...) of corresponding domain
	 */
	// TODO: implement this in a general form
	void setGeometry() {
#if DEBUG
		std::cout << "Setting Geometry for Domain " << _UID << std::endl;
#endif
		CVector<3,int> origin(1,_domain.getSize()[1] - 2,1);
		CVector<3,int> size(_domain.getSize()[0] - 2 ,1, _domain.getSize()[2] - 2);
#if DEBUG
		std::cout << "GEOMETRY: " << size << std::endl;
#endif
		int * src = new int[size.elements()];
		for( int i = 0; i < size.elements(); i++)
			src[i] = FLAG_VELOCITY_INJECTION;
		cLbmPtr->setFlags(src,origin,size);
		delete[] src;
	}

	CLbmSolver<T>* getSolver() const {
		return cLbmPtr;
	}

	void setSolver(CLbmSolver<T>* lbmPtr) {
		cLbmPtr = lbmPtr;
	}

	CDomain<T> getDomain() const {
		return _domain;
	}

	int getUid() const {
		return _UID;
	}
};

#endif
