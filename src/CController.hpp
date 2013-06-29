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
#include <stdint.h>

#include <list>
#include <map>
#include <vector>
// simulation type
typedef float T;
//typedef double T;

#define MPI_TAG_ALPHA_SYNC 0
#define MPI_TAG_BETA_SYNC 1

typedef std::map<MPI_COMM_DIRECTION, CComm<T>* > comm_map;
typedef std::map<MPI_COMM_DIRECTION, CComm<T>* >::value_type comm_map_pair;
typedef std::map<MPI_COMM_DIRECTION, CComm<T>* >::iterator comm_map_ptr;

template <typename T>
class CController;

struct arg_block_syncbeta
{
  CController<T>* controller;
  MPI_COMM_DIRECTION direction;
  MPI_Request** req_send;
  MPI_Request** req_recv;
  cl_event ue_sync;
};

struct arg_block{
  int arg0; 
};

void CL_CALLBACK callback_profile(cl_event ev, cl_int event_status,
                                  void * user_data);
void CL_CALLBACK callback_syncbeta(cl_event ev, cl_int event_status,
                                   void * user_data);
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

	inline void storeDataAlpha(MPI_COMM_DIRECTION direction) {
    if ( _comm_container.find(direction) == _comm_container.end() ) 
      return;
#if DEBUG
    DEBUGPRINT("--> store data alpha: %s\n", get_string_direction(direction) )
#endif
    CComm<T>* comm = _comm_container[direction];
    cLbmPtr->storeDensityDistribution(	comm->getSendBuffer(),
                                          comm->getSendOrigin(),
                                          comm->getSendSize()
                                          );
	}

	inline void setDataAlpha(MPI_COMM_DIRECTION direction) {
    if ( _comm_container.find(direction) == _comm_container.end() ) 
      return;
#if DEBUG
    DEBUGPRINT("--> set data alpha: %s\n", get_string_direction(direction) )
#endif
    CComm<T>* comm = _comm_container[direction];
			cLbmPtr->setDensityDistribution(comm->getRecvBuffer(),
					comm->getRecvOrigin(),
					comm->getRecvSize());
		cLbmPtr->wait();
	}

	inline void storeDataBeta(MPI_COMM_DIRECTION direction) {
    // Store data from device to host
    if ( _comm_container.find(direction) == _comm_container.end() ) 
      return;
#if DEBUG
    DEBUGPRINT("--> store data beta: %s\n", get_string_direction(direction) )
#endif
    CComm<T>* comm = _comm_container[direction];
    // the send and receive origin and size values in beta sync is the opposite values of
    // CComm object of current communication, since the ghost layer data
    // will be sent back to their origin
    cLbmPtr->storeDensityDistribution( comm->getSendBuffer(),
        comm->getRecvOrigin(),
        comm->getRecvSize());
	}

	inline bool storeDataBeta(MPI_COMM_DIRECTION direction, 
                            cl_uint num_events_in_wait_list,
                            const cl_event *event_wait_list,
                            cl_event *event              
                            ) 
  {
    if ( _comm_container.find(direction) == _comm_container.end() ) 
      return false;
#if DEBUG
    DEBUGPRINT("--> store data beta: %s\n", get_string_direction(direction) )
#endif
    CComm<T>* comm = _comm_container[direction];
    cLbmPtr->storeDensityDistribution(	comm->getSendBuffer(),
                                        comm->getRecvOrigin(),
                                        comm->getRecvSize(),
                                        num_events_in_wait_list,
                                        event_wait_list,
                                        event
                                        );
    return true;
	}

  inline	void setDataBeta(MPI_COMM_DIRECTION direction) {
    if ( _comm_container.find(direction) == _comm_container.end() ) 
      return;
#if DEBUG
    DEBUGPRINT("--> set data beta: %s\n", get_string_direction(direction) )
#endif
    CComm<T>* comm = _comm_container[direction];
    // the send and receive origin and size values in beta sync is the opposite values of
    // CComm object of current communication, since the ghost layer data
    // will be sent back to their origin
    cLbmPtr->setDensityDistribution(comm->getRecvBuffer(),
        comm->getSendOrigin(),
        comm->getSendSize(),
        comm->getCommDirection());
		cLbmPtr->wait();
	}
    
    inline	void setDataBeta(MPI_COMM_DIRECTION direction,
        cl_uint num_events_in_wait_list,
        const cl_event *event_wait_list,
        cl_event *event              
        ) {
    if ( _comm_container.find(direction) == _comm_container.end() ) 
      return;
#if DEBUG
    DEBUGPRINT("--> set data beta: %s\n", get_string_direction(direction) )
#endif
        CComm<T>* comm = _comm_container[direction];
    // the send and receive origin and size values in beta sync is the opposite values of
    // CComm object of current communication, since the ghost layer data
    // will be sent back to their origin
    cLbmPtr->setDensityDistribution(comm->getRecvBuffer(),
        comm->getSendOrigin(),
        comm->getSendSize(),
        comm->getCommDirection(),
        num_events_in_wait_list,
        event_wait_list,
        event);
    //		cLbmPtr->wait();
	}

	void syncAlpha() {
#if DEBUG
    DEBUGPRINT( "--> Sync alpha\n")
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
    DEBUGPRINT("--> Sync alpha: %s\n", get_string_direction(direction))
#endif
		// TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking way.
		// iterating over all communication objects for a specific direction
    MPI_COMM_DIRECTION lower_bound;
    MPI_COMM_DIRECTION upper_bound;
    switch(direction) {
      case MPI_COMM_DIRECTION_X: 
        lower_bound = MPI_COMM_DIRECTION_X_0;
        upper_bound = MPI_COMM_DIRECTION_X_1;
        break;
      case MPI_COMM_DIRECTION_Y: 
        lower_bound = MPI_COMM_DIRECTION_Y_0;
        upper_bound = MPI_COMM_DIRECTION_Y_1;
        break;
      case MPI_COMM_DIRECTION_Z: 
        lower_bound = MPI_COMM_DIRECTION_Z_0;
        upper_bound = MPI_COMM_DIRECTION_Z_1;
        break;
      default:
        lower_bound = direction;
        upper_bound = direction;
        break;
    }
    comm_map_ptr it,itlow,itup;
    itlow = _comm_container.lower_bound(lower_bound);
    itup =  _comm_container.upper_bound(upper_bound);
		for( it= itlow;it != itup; it++){
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
    DEBUGPRINT( "--> Sync beta\n")
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
		DEBUGPRINT( "--> Sync beta: %s\n", get_string_direction(direction))
#endif
		// TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking form.
		// iterating over all communication objects for a specific direction
    MPI_COMM_DIRECTION lower_bound;
    MPI_COMM_DIRECTION upper_bound;
    switch(direction) {
      case MPI_COMM_DIRECTION_X: 
        lower_bound = MPI_COMM_DIRECTION_X_0;
        upper_bound = MPI_COMM_DIRECTION_X_1;
        break;
      case MPI_COMM_DIRECTION_Y: 
        lower_bound = MPI_COMM_DIRECTION_Y_0;
        upper_bound = MPI_COMM_DIRECTION_Y_1;
        break;
      case MPI_COMM_DIRECTION_Z: 
        lower_bound = MPI_COMM_DIRECTION_Z_0;
        upper_bound = MPI_COMM_DIRECTION_Z_1;
        break;
      default:
        lower_bound = direction;
        upper_bound = direction;
    }
    comm_map_ptr it,itlow,itup;
    itlow = _comm_container.lower_bound(lower_bound);
    itup =  _comm_container.upper_bound(upper_bound);
		// std::pair<comm_map_ptr, comm_map_ptr> ppp;
		// ppp = _comm_container.equal_range(direction);
		// comm_map_ptr it = ppp.first;
		for( it = itlow ;it != itup; it++) {
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

  inline void syncAlpha(MPI_COMM_DIRECTION direction, MPI_Request** req_send, MPI_Request** req_recv) {
    if ( _comm_container.find(direction) == _comm_container.end() ) 
      return;
#if DEBUG
    DEBUGPRINT( "--> Sync alpha: %s\n", get_string_direction(direction))
#endif
    *req_send = (MPI_Request*) malloc(1*sizeof(MPI_Request));
    *req_recv = (MPI_Request*) malloc(1*sizeof(MPI_Request));
    CComm<T>* comm = _comm_container[direction];
			int dst_rank = comm->getDstId();
			int send_buffer_size = comm->getSendBufferSize();
			int recv_buffer_size = comm->getRecvBufferSize();
			T* send_buffer = comm->getSendBuffer();
			T* recv_buffer = comm->getRecvBuffer();
			if (typeid(T) == typeid(float)){
				MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, *req_send);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, *req_recv);
			}else if (typeid(T) == typeid(double)) {
    MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, *req_send);
    MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, *req_recv);
			} else {
				throw "Type id of MPI send/receive buffer is unknown!";
			}
	}

	inline void syncBeta(MPI_COMM_DIRECTION direction, MPI_Request** req_send, MPI_Request** req_recv) {
    if ( _comm_container.find(direction) == _comm_container.end() ) 
      return;
#if DEBUG
		DEBUGPRINT("--> Sync beta: %s\n", get_string_direction(direction))
#endif

    *req_send = (MPI_Request*) malloc(1*sizeof(MPI_Request));
    *req_recv = (MPI_Request*) malloc(1*sizeof(MPI_Request));
    CComm<T>* comm = _comm_container[direction];
    int dst_rank = comm->getDstId();
    int send_buffer_size = comm->getSendBufferSize();//send_size.elements()*cLbmPtr->SIZE_DD_HOST;
    int recv_buffer_size = comm->getRecvBufferSize();
    T* send_buffer = comm->getSendBuffer();
    T* recv_buffer = comm->getRecvBuffer();
    if (typeid(T) == typeid(float)) {
       MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, *req_send);
       MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, *req_recv);
    } else if (typeid(T) == typeid(double)) {
       MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, *req_send);
       MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank, MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, *req_recv);
    } else {
       throw "Type id of MPI send/receive buffer is unknown!";
    }
	}

	void simulationStepAlpha() {
		// SIMULATION_STEP_ALPHA
		CVector<3,int> x_size(1						, _domain.getSize()[1]	, _domain.getSize()[2]	);
		CVector<3,int> y_size(_domain.getSize()[0]	- 4, 1						, _domain.getSize()[2]	);
		CVector<3,int> z_size(_domain.getSize()[0]	- 4, _domain.getSize()[1] - 4	, 1						);

		// --> Simulation step alpha x boundary
		CVector<3,int> x0_origin(1, 0, 0);
		CVector<3,int> x1_origin(_domain.getSize()[0]-2, 0, 0);

    cl_event cl_ss_x_events[2];   // events of simulation step alpha in x direction
		cLbmPtr->simulationStepAlphaRect(x0_origin, x_size, 0, NULL, &cl_ss_x_events[0]);
		cLbmPtr->simulationStepAlphaRect(x1_origin, x_size, 0, NULL, &cl_ss_x_events[1]);

		// --> Store x boundary
		storeDataAlpha(MPI_COMM_DIRECTION_X_0);
		storeDataAlpha(MPI_COMM_DIRECTION_X_1);

		// --> Communication x boundary
    MPI_Request* req_send_x0 = NULL;
    MPI_Request* req_recv_x0 = NULL;

		syncAlpha(MPI_COMM_DIRECTION_X_0, &req_send_x0, &req_recv_x0);
    if( req_send_x0 ) {
      MPI_Status stat_send;
      MPI_Wait( req_send_x0, &stat_send);
    }
    if( req_recv_x0 ) {
      MPI_Status stat_recv;
      MPI_Wait( req_recv_x0, &stat_recv);
    }

    MPI_Request* req_send_x1 = NULL;
    MPI_Request* req_recv_x1 = NULL;

		syncAlpha(MPI_COMM_DIRECTION_X_1, &req_send_x1, &req_recv_x1);
    if( req_send_x1 ) {
      MPI_Status stat_send;
      MPI_Wait( req_send_x1, &stat_send);
    }
    if( req_recv_x1 ) {
      MPI_Status stat_recv;
      MPI_Wait( req_recv_x1, &stat_recv);
    }

    // --> Set x boundary
		setDataAlpha(MPI_COMM_DIRECTION_X_0);
		setDataAlpha(MPI_COMM_DIRECTION_X_1);

		// --> Simulation step alpha y boundary
		CVector<3,int> y0_origin(2, 1, 0);
		CVector<3,int> y1_origin(2,_domain.getSize()[1]-2, 0);
		cLbmPtr->simulationStepAlphaRect(y0_origin, y_size);
		cLbmPtr->simulationStepAlphaRect(y1_origin, y_size);

		// --> Store y boundary
		storeDataAlpha(MPI_COMM_DIRECTION_Y_0);
		storeDataAlpha(MPI_COMM_DIRECTION_Y_1);
		// --> Communication y boundary
		syncAlpha(MPI_COMM_DIRECTION_Y);
    // --> Set y boundary
		setDataAlpha(MPI_COMM_DIRECTION_Y_0);
		setDataAlpha(MPI_COMM_DIRECTION_Y_1);

		// --> Simulation step alpha z boundary
		CVector<3,int> z0_origin(2, 2, 1);
		CVector<3,int> z1_origin(2, 2, _domain.getSize()[2]-2);
		cLbmPtr->simulationStepAlphaRect(z0_origin, z_size);
		cLbmPtr->simulationStepAlphaRect(z1_origin, z_size);

		// --> Store z boundary
		storeDataAlpha(MPI_COMM_DIRECTION_Z_0);
		storeDataAlpha(MPI_COMM_DIRECTION_Z_1);
		// --> Communication z boundary
		syncAlpha(MPI_COMM_DIRECTION_Z);
    // --> Set z boundary
		setDataAlpha(MPI_COMM_DIRECTION_Z_0);
		setDataAlpha(MPI_COMM_DIRECTION_Z_1);

		// --> Computation of inner part
		CVector<3,int> inner_origin(2, 2, 2);
		CVector<3,int> inner_size(_domain.getSize()[0]	- 4, _domain.getSize()[1] - 4, _domain.getSize()[2] - 4);
		cLbmPtr->simulationStepAlphaRect(inner_origin, inner_size);
	}


	void simulationStepBeta() {
		// SIMULATION_STEP_ALPHA
		CVector<3,int> x_size(1						, _domain.getSize()[1]	, _domain.getSize()[2]	);
		CVector<3,int> y_size(_domain.getSize()[0]	- 4, 1						, _domain.getSize()[2]	);
		CVector<3,int> z_size(_domain.getSize()[0]	- 4, _domain.getSize()[1] - 4	, 1						);

		// --> Simulation step alpha x boundary
		CVector<3,int> x0_origin(1, 0, 0);
		CVector<3,int> x1_origin(_domain.getSize()[0]-2, 0, 0);
    cl_int err;
    cl_event ue_ss_trigger = clCreateUserEvent(cContext->context, &err);

    cl_event ev_ss_x0;
    cl_event ev_ss_x1;
    cLbmPtr->simulationStepBetaRect(x0_origin, x_size, 1, &ue_ss_trigger, &ev_ss_x0);
		cLbmPtr->simulationStepBetaRect(x1_origin, x_size, 1, &ue_ss_trigger, &ev_ss_x1);

		// --> Store x boundary
    cl_event ev_store_x0;
    cl_event ev_store_x1;

    bool b_ss_x0 = storeDataBeta(MPI_COMM_DIRECTION_X_0, 1, &ev_ss_x0, &ev_store_x0);
    bool b_ss_x1 = storeDataBeta(MPI_COMM_DIRECTION_X_1, 1, &ev_ss_x1, &ev_store_x1);

    MPI_Request* req_send_x0 = NULL;
    MPI_Request* req_recv_x0 = NULL;
    MPI_Request* req_send_x1 = NULL;
    MPI_Request* req_recv_x1 = NULL;

    struct arg_block_syncbeta ipargs_syncbeta_x0;
    ipargs_syncbeta_x0.controller = this;
    ipargs_syncbeta_x0.direction = MPI_COMM_DIRECTION_X_0;
    ipargs_syncbeta_x0.req_send = &req_send_x0;
    ipargs_syncbeta_x0.req_recv = &req_recv_x0;
    ipargs_syncbeta_x0.ue_sync =  clCreateUserEvent(cContext->context, &err);

    if (b_ss_x0)
      clSetEventCallback(ev_store_x0, CL_COMPLETE,
                         &callback_syncbeta,(void *)&ipargs_syncbeta_x0);

    struct arg_block_syncbeta ipargs_syncbeta_x1;
    ipargs_syncbeta_x1.controller = this;
    ipargs_syncbeta_x1.direction = MPI_COMM_DIRECTION_X_1;
    ipargs_syncbeta_x1.req_send = &req_send_x1;
    ipargs_syncbeta_x1.req_recv = &req_recv_x1;
    ipargs_syncbeta_x1.ue_sync =  clCreateUserEvent(cContext->context, &err);

    if (b_ss_x1) 
      clSetEventCallback(ev_store_x1, CL_COMPLETE,
                         &callback_syncbeta,(void *)&ipargs_syncbeta_x1);

    err = clSetUserEventStatus(ue_ss_trigger, CL_SUCCESS);
    cLbmPtr->wait();

    if( _BC[0][0] == FLAG_GHOST_LAYER  ) {
      DEBUGPRINT("enter_waitx0_send\n");
      MPI_Status stat_send;
      MPI_CHECK_ERROR(MPI_Wait( req_send_x0, &stat_send));
      DEBUGPRINT("waitx0_send\n");

      DEBUGPRINT("enter_waitx0_recv\n");
      MPI_Status stat_recv;
      MPI_CHECK_ERROR(MPI_Wait( req_recv_x0, &stat_recv));
      DEBUGPRINT("waitx0_recv\n");
    }
    if( _BC[0][1] == FLAG_GHOST_LAYER) {
      DEBUGPRINT("enter_waitx1_send\n")
      MPI_Status stat_send;
      MPI_CHECK_ERROR(MPI_Wait( req_send_x1, &stat_send));
      DEBUGPRINT("waitx1_send\n")

      DEBUGPRINT("enter_waitx1_recv\n")
      MPI_Status stat_recv;
      MPI_CHECK_ERROR(MPI_Wait( req_recv_x1, &stat_recv));
      DEBUGPRINT("waitx1_recv\n")
    }

		// // --> Set x boundary
		setDataBeta(MPI_COMM_DIRECTION_X_0, 1, &ipargs_syncbeta_x0.ue_sync, NULL);
		setDataBeta(MPI_COMM_DIRECTION_X_1, 1, &ipargs_syncbeta_x1.ue_sync, NULL);

		// --> Simulation step alpha y boundary
		CVector<3,int> y0_origin(2, 1, 0);
		CVector<3,int> y1_origin(2,_domain.getSize()[1]-2, 0);
		cLbmPtr->simulationStepBetaRect(y0_origin, y_size);
		cLbmPtr->simulationStepBetaRect(y1_origin, y_size);

		// --> Store y boundary
		storeDataBeta(MPI_COMM_DIRECTION_Y_0);
		storeDataBeta(MPI_COMM_DIRECTION_Y_1);
		// --> Communication y boundary
		syncBeta(MPI_COMM_DIRECTION_Y);
		// --> Set y boundary
		setDataBeta(MPI_COMM_DIRECTION_Y_0);
		setDataBeta(MPI_COMM_DIRECTION_Y_1);

		// --> Simulation step alpha z boundary
		CVector<3,int> z0_origin(2, 2, 1);
		CVector<3,int> z1_origin(2, 2, _domain.getSize()[2]-2);
		cLbmPtr->simulationStepBetaRect(z0_origin, z_size);
		cLbmPtr->simulationStepBetaRect(z1_origin, z_size);

		// --> Store z boundary
		storeDataBeta(MPI_COMM_DIRECTION_Z_0);
		storeDataBeta(MPI_COMM_DIRECTION_Z_1);
		// --> Communication z boundary
		syncBeta(MPI_COMM_DIRECTION_Z);
		// --> Set z boundary
		setDataBeta(MPI_COMM_DIRECTION_Z_0);
		setDataBeta(MPI_COMM_DIRECTION_Z_1);

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

void CL_CALLBACK callback_syncbeta(cl_event ev, cl_int event_status,
                                    void * user_data)
{
  arg_block_syncbeta* ipargs = (arg_block_syncbeta*) user_data;
#if DEBUG
  DEBUGPRINT( "--> callback DEBUG beta: %d\n", ipargs->direction)
#endif
#if DEBUG
  DEBUGPRINT( "--> callback sync beta: %s\n", get_string_direction(ipargs->direction))
#endif
  ipargs->controller->syncBeta(ipargs->direction, ipargs->req_send, ipargs->req_recv);
  CL_CHECK_ERROR(clSetUserEventStatus(ipargs->ue_sync, CL_SUCCESS));
  std::cout << "DEBUGER" << std::endl;
}

void CL_CALLBACK callback_profile(cl_event ev, cl_int event_status,
                                  void * user_data) {

  arg_block* ipargs = (arg_block*) user_data;
  int evID = ipargs->arg0;
  std::cout << "CALLBACK FUNCITON " << std::endl;
  cl_ulong ev_start_time=(cl_ulong)0;
  cl_ulong ev_end_time=(cl_ulong)0;
  size_t return_bytes;  cl_ulong run_time;
  printf(" Event callback %d %d \n",(int)event_status, evID);
  CL_CHECK_ERROR( clGetEventProfilingInfo( ev, CL_PROFILING_COMMAND_QUEUED,
                                           sizeof(cl_ulong),  &ev_start_time, &return_bytes));
  CL_CHECK_ERROR( clGetEventProfilingInfo( ev, CL_PROFILING_COMMAND_END,
                                           sizeof(cl_ulong),  &ev_end_time,   &return_bytes));
  run_time = (ev_end_time - ev_start_time);
  std::cout << "kernel start "<<  ev_start_time <<" nanaosecs\n";
  std::cout << "kernel end "<<  ev_end_time <<" nanaosecs\n";
  std::cout << "kernel runtime "<<  run_time <<" nanaosecs\n";
}

#endif
