#include <UnitTest++.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include "libcl/CCL.hpp"
#include "libtools/CStopwatch.hpp"

#include "CDomain.hpp"
#include "CLbmSolver.hpp"

typedef float T;

struct CLbmSolverFixture
{

	CLbmSolverFixture() {

		domain_size[0] = 16;
		domain_size[1] = 16;
		domain_size[2] = 16;

		CVector<3,int> origin_1(0,0,0);
		CVector<3,T> length(0.05,0.05,0.05);
		CDomain<T> domain(0, domain_size, origin_1, length);

		CVector<3, T> gravitation(0, -9.81, 0);
		T viscosity = 0.001308;
		T timestep = -1.0;
		int steps = -1;
		bool gui = false;
		size_t computation_kernel_count = 128;
		int device_nr = 0;

		std::list<int> p_lbm_opencl_number_of_work_items_list; ///< list with number of threads for each successively created kernel
		std::list<int> p_lbm_opencl_number_of_registers_list;

		int loops = steps;
		if (loops < 0)
			loops = 100;

		bool debug_mode = false;

		if (debug_mode)
			std::cout << "domain size: " << domain_size << std::endl;

//		if (pause)
//			next_simulation_steps_count = 0;

		// load platform information
		if (debug_mode)
			std::cout << "loading platforms" << std::endl;
		CCL::CPlatforms cPlatforms;
		cPlatforms.load();

		if (cPlatforms.platform_ids_count == 0) {
			std::cerr << "no platform found!" << std::endl;
			CHECK(false);
		}

		int platform_id_nr = -1;
		for (size_t i = 0; i < cPlatforms.platform_ids_count; i++) {
			CCL::CPlatform cPlatform(cPlatforms.platform_ids[i]);
			cPlatform.loadPlatformInfo();

			if (platform_id_nr == -1)
				if (strcmp(cPlatform.profile, "FULL_PROFILE") == 0) {
					platform_id_nr = i;
					if (debug_mode)
						std::cout << "Using Platform " << (i + 1)
								<< " for computation" << std::endl;
				}

			if (debug_mode) {
				std::cout << "Platform " << (i) << ":" << std::endl;
				std::cout << "        Name: " << cPlatform.name << std::endl;
				std::cout << "     Profile: " << cPlatform.profile << std::endl;
				std::cout << "     Version: " << cPlatform.version << std::endl;
				std::cout << "      Vendor: " << cPlatform.vendor << std::endl;
				std::cout << "  Extensions: " << cPlatform.extensions
						<< std::endl;
				std::cout << std::endl;
			}
		}

		if (platform_id_nr == -1) {
			std::cout << "no usable platform found" << std::endl;
			CHECK(false);
		}

		CCL::CPlatform cPlatform(cPlatforms.platform_ids[platform_id_nr]);

		// load standard context for GPU devices
		if (debug_mode)
			std::cout << "loading gpu context" << std::endl;
		CCL::CContext cContext(cPlatform, CL_DEVICE_TYPE_GPU);

		// load devices belonging to cContext
		if (debug_mode)
			std::cout << "loading devices" << std::endl;
		CCL::CDevices cDevices(cContext);

		if (cDevices.size() == 0) {
			std::cerr << "no device found - aborting" << std::endl;
			CHECK(false);
		}

		if (device_nr == -1) {
			// list available devices
			for (int i = 0; i < (int) cDevices.size(); i++) {
				CCL::CDeviceInfo cDeviceInfo(cDevices[i]);
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
			CHECK(false);
		}

		if (device_nr < 0 || device_nr >= (int) cDevices.size()) {
			std::cerr
					<< "invalid device number - use option \"-d -1\" to list all devices"
					<< std::endl;
			CHECK(false);
		}

		CCL::CDevice &cDevice = cDevices[device_nr];

		// load information about first device - e.g. max_work_group_size
		if (debug_mode)
			std::cout << "loading device information" << std::endl;
		CCL::CDeviceInfo cDeviceInfo(cDevice);

		if (debug_mode) {
			std::cout << "Device " << (device_nr) << ":" << std::endl;
			std::cout << "        Name: " << cDeviceInfo.name << std::endl;
			std::cout << "     Profile: " << cDeviceInfo.profile << std::endl;
			std::cout << "     Version: " << cDeviceInfo.version << std::endl;
			std::cout << "      Vendor: " << cDeviceInfo.vendor << std::endl;
			std::cout << "  Extensions: " << cDeviceInfo.extensions
					<< std::endl;
			std::cout << std::endl;
		}

		// initialize queue
		if (debug_mode)
			std::cout << "creating command queue" << std::endl;
		CCL::CCommandQueue cCommandQueue(cContext, cDevice);

		int BC[3][2];
		// INIT LATTICE BOLTZMANN!
		 cLbm = new CLbmSolver<T>(	0, cCommandQueue,
							cContext,
							cDevice,
							BC,
							domain,
							gravitation, // gravitation vector
							viscosity,
							computation_kernel_count,
							gui || debug_mode,
							gui || debug_mode,
							timestep,
							p_lbm_opencl_number_of_work_items_list,
							p_lbm_opencl_number_of_registers_list);
			if (cLbm->error())
			{
				std::cout << cLbm->error.getString();
				CHECK(false);
			}
		 cLbm->wait();
		for (int i = 0; i < 10; i++) {
			// simulation
		  //cLbm->simulationStep();
			std::cout << "." << std::flush;
		}
		std::cout << std::endl;
		cLbm->wait();

		CVector<3,int> origin(0,0,0);
		velocity = new T[domain_size.elements()*3];
		velocity_blockwise = new T[domain_size.elements()*3];

		cLbm->storeVelocity(velocity);
		cLbm->storeVelocity(velocity_blockwise, origin, domain_size);

		density = new T[domain_size.elements()];
		density_blockwise = new T[domain_size.elements()];

		cLbm->storeDensity(density);
		cLbm->storeDensity(density_blockwise, origin, domain_size);

		dd = new T[domain_size.elements()*SIZE_DD_HOST];
		dd_blockwise = new T[domain_size.elements()*SIZE_DD_HOST];

		cLbm->storeDensityDistribution(dd);
		cLbm->storeDensityDistribution(dd_blockwise, origin, domain_size);

		flags = new int[domain_size.elements()];
		flags_blockwise = new int[domain_size.elements()];

		cLbm->storeFlags(flags);
		cLbm->storeFlags(flags_blockwise, origin, domain_size);


		part_size[0] = 10;
		part_size[1] = 1;
		part_size[2] = 1;
		dd_part_blockwise = new T[part_size.elements()*SIZE_DD_HOST];

		cLbm->storeDensityDistribution(dd_part_blockwise, origin, part_size);

	}

	~CLbmSolverFixture() {
		delete cLbm;
		delete velocity;
		delete velocity_blockwise;
		delete density;
		delete density_blockwise;
		delete dd;
		delete dd_blockwise;
		delete flags;
		delete flags_blockwise;
	}

	static const size_t SIZE_DD_HOST = 19;
	CVector<3, int> domain_size;
	CVector<3, int> part_size;
	CLbmSolver<T> *cLbm;

	T* velocity;
	T* velocity_blockwise;

	T* density;
	T* density_blockwise;

	T* dd;
	T* dd_blockwise;

	int* flags;
	int* flags_blockwise;

	T* dd_part;
	T* dd_part_blockwise;

};

struct CLbmSolverFixtureWrite {

	CLbmSolverFixtureWrite(){
		domain_size[0] = 16;
		domain_size[1] = 16;
		domain_size[2] = 16;

		CVector<3,int> origin_1(0,0,0);
		CVector<3,T> length(0.05,0.05,0.05);
		CDomain<T> domain(0, domain_size, origin_1, length);

		CVector<3, T> gravitation(0, -9.81, 0);
		T viscosity = 0.001308;
		T timestep = -1.0;
		int steps = -1;

		bool gui = false;
		size_t computation_kernel_count = 128;
		int device_nr = 0;

		std::list<int> p_lbm_opencl_number_of_work_items_list; ///< list with number of threads for each successively created kernel
		std::list<int> p_lbm_opencl_number_of_registers_list;

		int loops = steps;
		if (loops < 0)
			loops = 100;

		bool debug_mode = false;

		if (debug_mode)
			std::cout << "domain size: " << domain_size << std::endl;

//		if (pause)
//			next_simulation_steps_count = 0;

		// load platform information
		if (debug_mode)
			std::cout << "loading platforms" << std::endl;
		CCL::CPlatforms cPlatforms;
		cPlatforms.load();

		if (cPlatforms.platform_ids_count == 0) {
			std::cerr << "no platform found!" << std::endl;
			CHECK(false);
		}

		int platform_id_nr = -1;
		for (size_t i = 0; i < cPlatforms.platform_ids_count; i++) {
			CCL::CPlatform cPlatform(cPlatforms.platform_ids[i]);
			cPlatform.loadPlatformInfo();

			if (platform_id_nr == -1)
				if (strcmp(cPlatform.profile, "FULL_PROFILE") == 0) {
					platform_id_nr = i;
					if (debug_mode)
						std::cout << "Using Platform " << (i + 1)
								<< " for computation" << std::endl;
				}

			if (debug_mode) {
				std::cout << "Platform " << (i) << ":" << std::endl;
				std::cout << "        Name: " << cPlatform.name << std::endl;
				std::cout << "     Profile: " << cPlatform.profile << std::endl;
				std::cout << "     Version: " << cPlatform.version << std::endl;
				std::cout << "      Vendor: " << cPlatform.vendor << std::endl;
				std::cout << "  Extensions: " << cPlatform.extensions
						<< std::endl;
				std::cout << std::endl;
			}
		}

		if (platform_id_nr == -1) {
			std::cout << "no usable platform found" << std::endl;
			CHECK(false);
		}

		CCL::CPlatform cPlatform(cPlatforms.platform_ids[platform_id_nr]);

		// load standard context for GPU devices
		if (debug_mode)
			std::cout << "loading gpu context" << std::endl;
		CCL::CContext cContext(cPlatform, CL_DEVICE_TYPE_GPU);

		// load devices belonging to cContext
		if (debug_mode)
			std::cout << "loading devices" << std::endl;
		CCL::CDevices cDevices(cContext);

		if (cDevices.size() == 0) {
			std::cerr << "no device found - aborting" << std::endl;
			CHECK(false);
		}

		if (device_nr == -1) {
			// list available devices
			for (int i = 0; i < (int) cDevices.size(); i++) {
				CCL::CDeviceInfo cDeviceInfo(cDevices[i]);
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
			CHECK(false);
		}

		if (device_nr < 0 || device_nr >= (int) cDevices.size()) {
			std::cerr
					<< "invalid device number - use option \"-d -1\" to list all devices"
					<< std::endl;
			CHECK(false);
		}

		CCL::CDevice &cDevice = cDevices[device_nr];

		// load information about first device - e.g. max_work_group_size
		if (debug_mode)
			std::cout << "loading device information" << std::endl;
		CCL::CDeviceInfo cDeviceInfo(cDevice);

		if (debug_mode) {
			std::cout << "Device " << (device_nr) << ":" << std::endl;
			std::cout << "        Name: " << cDeviceInfo.name << std::endl;
			std::cout << "     Profile: " << cDeviceInfo.profile << std::endl;
			std::cout << "     Version: " << cDeviceInfo.version << std::endl;
			std::cout << "      Vendor: " << cDeviceInfo.vendor << std::endl;
			std::cout << "  Extensions: " << cDeviceInfo.extensions
					<< std::endl;
			std::cout << std::endl;
		}

		// initialize queue
		if (debug_mode)
			std::cout << "creating command queue" << std::endl;
		CCL::CCommandQueue cCommandQueue(cContext, cDevice);

		int BC[3][2];
		// INIT LATTICE BOLTZMANN!
		 cLbm = new CLbmSolver<T>(	0, cCommandQueue,
							cContext,
							cDevice,
							BC,
							domain,
							gravitation, // gravitation vector
							viscosity,
							computation_kernel_count,
							gui || debug_mode,
							gui || debug_mode,
							timestep,
							p_lbm_opencl_number_of_work_items_list,
							p_lbm_opencl_number_of_registers_list);
			if (cLbm->error())
			{
				std::cout << cLbm->error.getString();
				CHECK(false);
			}
		 cLbm->wait();
		for (int i = 0; i < 10; i++) {
			// simulation
		  //cLbm->simulationStep();
			std::cout << "." << std::flush;
		}
		std::cout << std::endl;
		cLbm->wait();

		CVector<3,int> origin(0,0,0);

		// density setup
		wdensity =  new T[domain_size.elements()];
		wdensity_blockwise = new T[domain_size.elements()];

		for (int i = 0 ; i < domain_size.elements(); i++){
			wdensity[i] = i;
		}

		cLbm->setDensity(wdensity, origin, domain_size);
		cLbm->storeDensity(wdensity_blockwise);

		// velocity setup
		wvelocity =  new T[domain_size.elements()*3];
		wvelocity_blockwise = new T[domain_size.elements()*3];

		for (int i = 0 ; i < domain_size.elements()*3; i++){
			wvelocity[i] = i;
		}

		this->cLbm->setVelocity(wvelocity, origin, domain_size);
		this->cLbm->storeVelocity(wvelocity_blockwise);

		// density distribution setup
		wdd = new T[domain_size.elements()*SIZE_DD_HOST];
		wdd_blockwise = new T[domain_size.elements()*SIZE_DD_HOST];

		for (int i = 0 ; i < domain_size.elements()*SIZE_DD_HOST; i++){
			wdd[i] = i;
		}

		cLbm->setDensityDistribution(wdd, origin, domain_size);
		cLbm->storeDensityDistribution(wdd_blockwise);

		// flags setup
		wflags = new int[domain_size.elements()];
		wflags_blockwise = new int[domain_size.elements()];

		for (int i = 0 ; i < domain_size.elements(); i++){
			wflags[i] = i;
		}

		cLbm->setFlags(wflags, origin, domain_size);
		cLbm->storeFlags(wflags_blockwise);

	}

	~CLbmSolverFixtureWrite(){
		delete cLbm;
		delete wdensity;
		delete wdensity_blockwise;
		delete wvelocity;
		delete wvelocity_blockwise;
	}
	static const size_t SIZE_DD_HOST = 19;
	CVector<3, int> domain_size;
	CLbmSolver<T> *cLbm;

	T* wdensity;
	T* wdensity_blockwise;

	T* wvelocity;
	T* wvelocity_blockwise;

	T* wdd;
	T* wdd_blockwise;

	int* wflags;
	int* wflags_blockwise;
};



SUITE(CLBMSolverSuite)
{

TEST_FIXTURE(CLbmSolverFixture, StoreVelociyBlockwise) {
	T* tmpvelocity = velocity;
	T* tmpvelocity_blockwise = velocity_blockwise;

	CHECK_ARRAY_EQUAL(tmpvelocity, tmpvelocity_blockwise, domain_size.elements()*3); // succeeds
}

TEST_FIXTURE(CLbmSolverFixture, StoreDensityBlockwise) {
	T* tmpdensity = density;
	T* tmpdensity_blockwise = density_blockwise;

	CHECK_ARRAY_EQUAL(tmpdensity, tmpdensity_blockwise, domain_size.elements()); // succeeds
}

TEST_FIXTURE(CLbmSolverFixture, StoreDensityDistributionBlockwise) {
	T* tmpdd = dd;
	T* tmpdd_blockwise = dd_blockwise;

	CHECK_ARRAY_EQUAL(tmpdd, tmpdd_blockwise, domain_size.elements()*SIZE_DD_HOST); // succeeds
}

TEST_FIXTURE(CLbmSolverFixture, StoreDensityDistributionPartBlockwise) {
	T* tmpdd = dd;
	T* tmpdd_blockwise = dd_part_blockwise;

	CHECK_ARRAY_EQUAL(tmpdd, tmpdd_blockwise, part_size[0]); // succeeds
}

TEST_FIXTURE(CLbmSolverFixture, StoreFlagsBlockwise) {
	int* tmpflags = flags;
	int* tmpflags_blockwise = flags_blockwise;

	CHECK_ARRAY_EQUAL(tmpflags, tmpflags_blockwise, domain_size.elements()); // succeeds
}




TEST_FIXTURE(CLbmSolverFixtureWrite, SetDensityBlockwise) {

	T* tmpdensity = wdensity;
	T* tmpdensity_blockwise = wdensity_blockwise;

	CHECK_ARRAY_EQUAL(tmpdensity, tmpdensity_blockwise, domain_size.elements()); // succeeds
}

TEST_FIXTURE(CLbmSolverFixtureWrite, SetVelocityBlockwise) {

	T* tmpwvelociy = wvelocity;
	T* tmpwvelocity_blockwise = wvelocity_blockwise;

	CHECK_ARRAY_EQUAL(tmpwvelociy, tmpwvelocity_blockwise, domain_size.elements()*3); // succeeds
}

TEST_FIXTURE(CLbmSolverFixtureWrite, SetDensityDistributionBlockwise) {

	T* tmpdd = wdd;
	T* tmpdd_blockwise = wdd_blockwise;

	CHECK_ARRAY_EQUAL(tmpdd, tmpdd_blockwise, domain_size.elements()*SIZE_DD_HOST); // succeeds
}

TEST_FIXTURE(CLbmSolverFixtureWrite, SetFlagsBlockwise) {

	int* tmpwflags = wflags;
	int* tmpwflags_blockwise = wflags_blockwise;

	CHECK_ARRAY_EQUAL(tmpwflags, tmpwflags_blockwise, domain_size.elements()); // succeeds
}
}
