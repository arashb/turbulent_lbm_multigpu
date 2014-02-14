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
#include <vector>
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
template<typename T>
class CController
{
	int _UID; ///< Unique ID of each controller
	CDomain<T> _domain; ///< Domain data
	ILbmVisualization<T>* cLbmVisualization; ///< Visualization class
	CLbmSolver<T> *cLbmPtr;
	int _BC[3][2]; ///< Boundary conditions. First index specifies the dimension and second the upper or the lower boundary.
	std::vector<CComm<T>*> _comm_container; ///< A std::Vector containing all the communcation objects for the subdomain

	T vector_checksum;
	CCL::CPlatforms* cPlatforms;
	CCL::CPlatform* cPlatform;
	CCL::CContext* cContext;
	CCL::CDevices* cDevices;
	CCL::CDevice* cDevice;
	CCL::CCommandQueue* cCommandQueue;

	void outputDD(int dd_i) {
		std::cout << "DD " << dd_i << std::endl;
		//					int gcd = CMath<int>::gcd(cLbmPtr->domain_cells[0],wrap_max_line);
		//					if (wrap_max_line % gcd == 0)
		//						gcd = wrap_max_line;
		int wrap_max_line = 16;
		int gcd = wrap_max_line;
		cLbmPtr->debugDD(dd_i, gcd,
				cLbmPtr->domain_cells[0] * cLbmPtr->domain_cells[1]);
	}

	int initLBMSolver() {

#if DEBUG
		std::cout << "loading platforms" << std::endl;
#endif
		cPlatforms = new CCL::CPlatforms();
		cPlatforms->load();

		if (cPlatforms->platform_ids_count == 0) {
			std::cerr << "no platform found!" << std::endl;
			return -1;
		}

		int platform_id_nr = -1;
		for (size_t i = 0; i < cPlatforms->platform_ids_count; i++) {
			CCL::CPlatform cPlatform(cPlatforms->platform_ids[i]);
			cPlatform.loadPlatformInfo();

			if (platform_id_nr == -1)
				if (strcmp(cPlatform.profile, "FULL_PROFILE") == 0) {
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

		if (platform_id_nr == -1) {
			std::cout << "no usable platform found" << std::endl;
			return -1;
		}

		cPlatform = new CCL::CPlatform(
				cPlatforms->platform_ids[platform_id_nr]);

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

		if (cDevices->size() == 0) {
			std::cerr << "no device found - aborting" << std::endl;
			return -1;
		}

		if (ConfigSingleton::Instance()->device_nr == -1) {
			// list available devices
			for (int i = 0; i < (int) cDevices->size(); i++) {
				CCL::CDeviceInfo cDeviceInfo((*cDevices)[i]);
				std::cout << "Device " << (i) << ":" << std::endl;
				std::cout << "        Name: " << cDeviceInfo.name << std::endl;
				std::cout << "     Profile: " << cDeviceInfo.profile
						<< std::endl;
				std::cout << "     Version: " << cDeviceInfo.version
						<< std::endl;
				std::cout << "      Vendor: " << cDeviceInfo.vendor
						<< std::endl;
				std::cout << "  Extensions: " << cDeviceInfo.extensions
						<< std::endl;
				std::cout << std::endl;
			}
			return -1;
		}

		if (ConfigSingleton::Instance()->device_nr < 0
				|| ConfigSingleton::Instance()->device_nr
						>= (int) cDevices->size()) {
			std::cerr
					<< "invalid device number - use option \"-d -1\" to list all devices"
					<< std::endl;
			return -1;
		}
		int dev_nr = 0;
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		int name_len;
		MPI_Get_processor_name(processor_name, &name_len);
		if (strstr(processor_name, "mac-nvd") != NULL) {
			if (_UID & 1)
				dev_nr = 1;
		}
		cDevice =
				&((*cDevices)[/* ConfigSingleton::Instance()->device_nrd*/dev_nr]);

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

		// INIT LATTICE BOLTZMANN!
		cLbmPtr = new CLbmSolver<T>(_UID, *cCommandQueue, *cContext, *cDevice,
				_BC, _domain,
				ConfigSingleton::Instance()->gravitation, // gravitation vector
				ConfigSingleton::Instance()->viscosity,
				ConfigSingleton::Instance()->computation_kernel_count,
				//ConfigSingleton::Instance()->debug_mode,
				ConfigSingleton::Instance()->do_visualization
						|| ConfigSingleton::Instance()->debug_mode,
				ConfigSingleton::Instance()->do_visualization
						|| ConfigSingleton::Instance()->debug_mode,
				ConfigSingleton::Instance()->timestep,
				 ConfigSingleton::Instance()->drivenCavityVelocity,
				ConfigSingleton::Instance()->lbm_opencl_number_of_threads_list, // p_lbm_opencl_number_of_work_items_list,
				ConfigSingleton::Instance()->lbm_opencl_number_of_threads_list);

		if (cLbmPtr->error()) {
			std::cout << cLbmPtr->error.getString();
			return -1;
		}

		cLbmPtr->wait();
		CStopwatch cStopwatch;
		return 0;
	}

public:

	CController(int UID, CDomain<T> domain, int BC[3][2]) :
			_UID(UID), _domain(domain), cLbmVisualization(NULL) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 2; j++)
				_BC[i][j] = BC[i][j];

		// initialize the LBMSolver
		if (-1 == initLBMSolver())
			throw "Initialization of LBM Solver failed!";
	}
	~CController() {
		if (cPlatforms)
			delete cPlatform;

		if (cContext)
			delete cContext;

		if (cDevices)
			delete cDevices;

		if (cCommandQueue)
			delete cCommandQueue;

		if (cLbmVisualization) ///< Visualization class
			delete cLbmVisualization;

		if (cLbmPtr)
			delete cLbmPtr;

		typename std::vector<CComm<T>*>::iterator it = _comm_container.begin();
		for (; it != _comm_container.end(); it++) {
			delete *it;
		}
	}

	void syncAlpha() {
#if DEBUG
		std::cout << "--> Sync alpha" << std::endl;
#endif
		// TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking way.
		typename std::vector<CComm<T>*>::iterator it = _comm_container.begin();
		for (; it != _comm_container.end(); it++) {
			CVector<3, int> send_size = (*it)->getSendSize();
			CVector<3, int> recv_size = (*it)->getRecvSize();
			CVector<3, int> send_origin = (*it)->getSendOrigin();
			CVector<3, int> recv_origin = (*it)->getRecvOrigin();
			int dst_rank = (*it)->getDstId();
#if DEBUG
			std::cout << "ALPHA RANK: " << dst_rank << std::endl;
#endif

			// send buffer
			int send_buffer_size = send_size.elements() * cLbmPtr->SIZE_DD_HOST;
			int recv_buffer_size = recv_size.elements() * cLbmPtr->SIZE_DD_HOST;
			T* send_buffer = new T[send_buffer_size];
			T* recv_buffer = new T[recv_buffer_size];

			MPI_Request req[2];
			MPI_Status status[2];

			// Download data from device to host
			cLbmPtr->storeDensityDistribution(send_buffer, send_origin,
					send_size);
			//cLbmPtr->wait();
			int my_rank, num_procs;
			MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); /// Get current process id
			MPI_Comm_size(MPI_COMM_WORLD, &num_procs); /// get number of processes

			if (typeid(T) == typeid(float)) {
				MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank,
						MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank,
						MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else if (typeid(T) == typeid(double)) {
				MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank,
						MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank,
						MPI_TAG_ALPHA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else {
				throw "Type id of MPI send/receive buffer is unknown!";
			}
			MPI_Waitall(2, req, status);

			cLbmPtr->setDensityDistribution(recv_buffer, recv_origin,
					recv_size);
			cLbmPtr->wait();

			delete[] send_buffer;
			delete[] recv_buffer;
		}
	}

	void syncBeta() {
#if DEBUG
		std::cout << "--> Sync beta" << std::endl;
#endif

		// TODO: OPTIMIZATION: communication of different neighbors can be done in Non-blocking form.
		typename std::vector<CComm<T>*>::iterator it = _comm_container.begin();
		for (; it != _comm_container.end(); it++) {
			// the send and receive values in beta sync is the opposite values of
			// Comm instance related to current communication, since the ghost layer data
			// will be sent back to their origin
			CVector<3, int> send_size = (*it)->getRecvSize();
			CVector<3, int> recv_size = (*it)->getSendSize();
			;
			CVector<3, int> send_origin = (*it)->getRecvOrigin();
			CVector<3, int> recv_origin = (*it)->getSendOrigin();
			CVector<3, int> normal = (*it)->getCommDirection();
			int dst_rank = (*it)->getDstId();
#if DEBUG
			std::cout << "BETA RANK: " << dst_rank << std::endl;
#endif
			// send buffer
			int send_buffer_size = send_size.elements() * cLbmPtr->SIZE_DD_HOST;
			int recv_buffer_size = recv_size.elements() * cLbmPtr->SIZE_DD_HOST;
			T* send_buffer = new T[send_buffer_size];
			T* recv_buffer = new T[recv_buffer_size];

			MPI_Request req[2];
			MPI_Status status[2];

			// Download data from device to host
			cLbmPtr->storeDensityDistribution(send_buffer, send_origin,
					send_size);
			//cLbmPtr->wait();
			int my_rank, num_procs;
			MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); /// Get current process id
			MPI_Comm_size(MPI_COMM_WORLD, &num_procs); /// get number of processes

			if (typeid(T) == typeid(float)) {
				MPI_Isend(send_buffer, send_buffer_size, MPI_FLOAT, dst_rank,
						MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_FLOAT, dst_rank,
						MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else if (typeid(T) == typeid(double)) {
				MPI_Isend(send_buffer, send_buffer_size, MPI_DOUBLE, dst_rank,
						MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[0]);
				MPI_Irecv(recv_buffer, recv_buffer_size, MPI_DOUBLE, dst_rank,
						MPI_TAG_BETA_SYNC, MPI_COMM_WORLD, &req[1]);
			} else {
				throw "Type id of MPI send/receive buffer is unknown!";
			}
			MPI_Waitall(2, req, status);

			// TODO: OPTIMIZATION: you need to wait only for receiving to execute following command
			cLbmPtr->setDensityDistribution(recv_buffer, recv_origin, recv_size,
					normal);
			cLbmPtr->wait();

			delete[] send_buffer;
			delete[] recv_buffer;
		}
	}

	void computeNextStep() {
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
	int run() {
		CVector<3, int> domain_size = _domain.getSize();
		int loops = ConfigSingleton::Instance()->loops;
		if (loops < 0)
			loops = 100;

		vector_checksum = 0;

		// approximate bandwidth
		double floats_per_cell = 0.0;

		// 19 density distribution which are read and written
		floats_per_cell += 19.0 * 2.0;

		// flag (obstacle, injection and fluid) is read
		floats_per_cell += 1.0;

		// velocity vector is also stored
		if (ConfigSingleton::Instance()->do_visualization
				|| ConfigSingleton::Instance()->debug_mode)
			floats_per_cell += 3;
		CStopwatch cStopwatch;

		// setting up the visualization
		std::string outputfilename = "OUTPUT";
		std::stringstream ss_file;
		ss_file << "./" << VTK_OUTPUT_DIR << "/" << outputfilename;
		std::string outputfile = ss_file.str();
		if (ConfigSingleton::Instance()->do_visualization) {
			cLbmVisualization = new CLbmVisualizationVTK<T>(_UID, outputfile);
			cLbmVisualization->setup(cLbmPtr);
		}

		cStopwatch.start();
		for (int i = 0; i < loops; i++) {
			computeNextStep();
			//simulation
			if (ConfigSingleton::Instance()->do_visualization // && (i %  100 == 0)
			    )
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
			"benchmark_" << ConfigSingleton::Instance()->subdomain_num.elements() //<< "_" << _UID
			<< ".ini";
			const std::string& tmp = benchmark_file_name.str();
			const char* cstr = tmp.c_str();
			std::ofstream benchmark_file (cstr, std::ios::out | std::ios::app );
			if (benchmark_file.is_open())
			{
				//benchmark_file << "[RESULTS]" << std::endl;
				benchmark_file << "CUBE_X : " << ConfigSingleton::Instance()->domain_size[0] << std::endl;
				benchmark_file << "CUBE_Y : " << ConfigSingleton::Instance()->domain_size[1] << std::endl;
				benchmark_file << "CUBE_Z : " << ConfigSingleton::Instance()->domain_size[2] << std::endl;
				benchmark_file << "SECONDS : " << gtime << std::endl;
				benchmark_file << "FPS : " << gfps << std::endl;
				benchmark_file << "MLUPS : " << gmlups << std::endl;
				benchmark_file << "BANDWIDTH : " << gbandwidth// << " MB/s (RW, bidirectional)"
				<< std::endl;
				benchmark_file << std::endl;

			}
			else std::cout << "Unable to open file";
		}
#endif // end of BENCHMARK
		std::cout << std::endl;
		std::cout << "Cube: " << domain_size << std::endl;
		std::cout << "Seconds: " << cStopwatch.time << std::endl;
		double fps = (((double) loops) / cStopwatch.time);
		std::cout << "FPS: " << fps << std::endl;
		double mlups =
				((double) fps * (double) cLbmPtr->domain_cells.elements())
						* (double) 0.000001;
		std::cout << "MLUPS: " << mlups << std::endl;
		std::cout << "Bandwidth: "
				<< (mlups * floats_per_cell * (double) sizeof(T))
				<< " MB/s (RW, bidirectional)" << std::endl;
		std::streamsize ss = std::cout.precision();
		std::cout.precision(8);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
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
		"profile_" << ConfigSingleton::Instance()->subdomain_num.elements() << "_" << _UID
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

	void addCommunication(CComm<T>* comm) {
		_comm_container.push_back(comm);
	}

	/*
	 * This Function is used the set the geometry (e.g obstacles, velocity injections, ...) of corresponding domain
	 */
	void setGeometry() {
#if DEBUG
		std::cout << "Setting Geometry for Domain " << _UID << std::endl;
#endif
		CVector<3, int> origin(1, _domain.getSize()[1] - 2, 1);
		CVector<3, int> size(_domain.getSize()[0] - 2, 1,
				_domain.getSize()[2] - 2);
#if DEBUG
		std::cout << "GEOMETRY: " << size << std::endl;
#endif
		int * src = new int[size.elements()];
		for (int i = 0; i < size.elements(); i++)
			src[i] = FLAG_VELOCITY_INJECTION;
		cLbmPtr->setFlags(src, origin, size);
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
