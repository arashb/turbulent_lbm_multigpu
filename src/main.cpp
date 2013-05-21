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


// standard
#include <stdlib.h>
#include <iostream>
#include <list>
#include <string>

// externals
#include <UnitTest++.h>
#include "TestReporterStdout.h"

#include "mpi.h"

// internals
#include "CDomain.hpp"
#include "CController.hpp"
#include "common.h"

CVector<3,int> E0(1,0,0) 	;
CVector<3,int> E1(-1,0,0)	;
CVector<3,int> E2(0,1,0)	;
CVector<3,int> E3(0,-1,0)	;

CVector<3,int> E4(1,1,0)	;
CVector<3,int> E5(-1,-1,0)	;
CVector<3,int> E6(1,-1,0)	;
CVector<3,int> E7(-1,1,0)	;

CVector<3,int> E8(1,0,1)	;
CVector<3,int> E9(-1,0,-1)	;
CVector<3,int> E10(1,0,-1)	;
CVector<3,int> E11(-1,0,1)	;

CVector<3,int> E12(0,1,1)	;
CVector<3,int> E13(0,-1,-1)	;
CVector<3,int> E14(0,1,-1)	;
CVector<3,int> E15(0,-1,1)	;

CVector<3,int> E16(0,0,1)	;
CVector<3,int> E17(0,0,-1)	;
CVector<3,int> E18(0,0,0)	;
CVector<3,int> lbm_units[] = {	E0,E1,E2,E3,
		E4,E5,E6,E7,
		E8,E9,E10,E11,
		E12,E13,E14,E15,
		E16,E17,E18
};

// simulation type
typedef float T;
//typedef double T;

void extract_comma_separated_integers(std::list<int> &int_list, std::string &int_string)
{
	size_t start_pos = 0;

	size_t comma_pos;

	comma_pos = int_string.find_first_of(',', start_pos);

	while (comma_pos != std::string::npos)
	{
		std::string substring = int_string.substr(start_pos, comma_pos-start_pos);

		int number;

		if (substring.empty())	number = 0;
		else					number = atoi(substring.c_str());

		int_list.push_back(number);
		start_pos = comma_pos+1;
		comma_pos = int_string.find_first_of(',', start_pos);
	}
	int_list.push_back(atoi(int_string.substr(start_pos).c_str()));
}

int main(int argc, char** argv)
{
	bool debug = false;

	CVector<3,int> domain_size(64,32,32);
	CVector<3,T> gravitation(0,-9.81,0);
	T viscosity = 0.001308;
	T timestep = -1.0;
	int steps = -1;

	bool gui = false;
	bool pause = false;
	bool take_frame_screenshots = false;
	bool unit_test = false;
	size_t computation_kernel_count = 128;

	std::string number_of_registers_string;	///< string storing the number of registers for opencl threads separated with comma
	std::string number_of_threads_string;		///< string storing the number of threads for opencl separated with comma
	std::string test_suite;

	int device_nr = 0;

	char optchar;
	while ((optchar = getopt(argc, argv, "x:y:z:d:vr:k:gG:pt:sl:R:T:X:u:")) > 0)
	{
		switch(optchar)
		{
			case 'd':
				device_nr = atoi(optarg);
				break;

			case 'v':
				debug = true;
				break;

			case 'R':
				number_of_registers_string = optarg;
				break;

			case 'T':
				number_of_threads_string = optarg;
				break;

			case 'p':
				pause = true;
				break;

			case 'x':
				domain_size[0] = atoi(optarg);
				break;

			case 'X':
				domain_size[0] = atoi(optarg);
				domain_size[1] = domain_size[0];
				domain_size[2] = domain_size[0];
				break;

			case 'l':
				steps = atoi(optarg);
				break;

			case 'y':
				domain_size[1] = atoi(optarg);
				break;

			case 'k':
				computation_kernel_count = atoi(optarg);
				break;

			case 'z':
				domain_size[2] = atoi(optarg);
				break;

			case 'r':
				viscosity = atof(optarg);
				break;

			case 'g':
				gui = true;
				break;

			case 'G':
				gravitation[1] = atof(optarg);
				break;

			case 't':
				timestep = atof(optarg);
				break;

			case 's':
				take_frame_screenshots = true;
				break;

			case 'u':
				unit_test = true;
				test_suite = optarg;
				break;
			default:
				goto parameter_error;
		}
	}

	goto parameter_error_ok;
parameter_error:
	std::cout << "usage: " << argv[0] << std::endl;
	std::cout << "		[-x resolution_x, default: 32]" << std::endl;
	std::cout << "		[-y resolution_y, default: 32]" << std::endl;
	std::cout << "		[-z resolution_z, default: 32]" << std::endl;
	std::cout << std::endl;
	std::cout << "		[-G gravitation in down direction, default: -9.81]" << std::endl;
	std::cout << std::endl;
	std::cout << "		[-r viscosity, default: 0.001308]" << std::endl;
	std::cout << "		[-v]	(debug mode on, be verbose)" << std::endl;
	std::cout << "		[-d device_num]	(-1: list available devices, or select device)" << std::endl;
	std::cout << "		[-k number of kernels to run ]	(default: 0 - autodetect maximum)" << std::endl;
	std::cout << "		[-g]	(activate gui, default:disabled)" << std::endl;
	std::cout << "		[-p]	(pause simulation at start, default:disabled)" << std::endl;
	std::cout << std::endl;
	std::cout << "		[-R registers for OpenCL]	(comma separated list of number of registers per threads in OpenCL)" << std::endl;
	std::cout << "		[-T work group threads]		(comma separated list of number of threads in OpenCL)" << std::endl;
	std::cout << std::endl;
	std::cout << "		[-t timestep]	(default: -1 for automatic detection)" << std::endl;
	std::cout << "		[-s]	(take a screenshot every frame - default: disabled)" << std::endl;
	std::cout << "    	[-u] run unit tests" << std::endl;
	return -1;

parameter_error_ok:


	  if (unit_test) {

		  if ( strcmp(  "all", test_suite.c_str() ) == 0 ) {
			  if( debug)
				  std::cout << "running all test." << std::endl;
			  return UnitTest::RunAllTests();
		  }

		  else {
			  const UnitTest::TestList& allTests( UnitTest::Test::GetTestList() );
			  UnitTest::TestList selectedTests;
			  UnitTest::Test* p = allTests.GetHead();
			  while( p )
			  {
				  //			  for( int i = 1 ; i < argc ; ++i )
				  if( strcmp(  p->m_details.suiteName , test_suite.c_str() ) == 0 ) {
					  selectedTests.Add( p );
					  if( debug)
						  std::cout << "Added test " << p->m_details.testName << "from suite " <<  p->m_details.suiteName << " to tes list." << std::endl;
				  }
				  p = p->next;
			  }

			  //run selected test(s) only
			  UnitTest::TestReporterStdout reporter;
			  UnitTest::TestRunner runner( reporter );
			  return runner.RunTestsIf( selectedTests, 0, UnitTest::True(), 0 );
		  }

	  } else {

		  CVector<3,int> domain_size_1(32,32,32);
		  CVector<3,int> domain_size_2(32,32,32);
		  CVector<3,int> origin_1(0,0,0);
		  CVector<3,int> origin_2(32,0,0);

		  CVector<3,T> length(0.05,0.05,0.05);
		  CDomain<T> domain_1(0, domain_size_1, origin_1, length);
		  CDomain<T> domain_2(1, domain_size_2, origin_2, length);

		  CController<T> lbmController1(0,domain_1);
		  CController<T> lbmController2(1,domain_2);

		  int BC1[3][2] = { /* x BC */FLAG_OBSTACLE,FLAG_GHOST_LAYER,
				  	  	  	/* y BC */FLAG_OBSTACLE,FLAG_OBSTACLE,
				  	  	  	/* z BC */FLAG_OBSTACLE,FLAG_OBSTACLE};

		  int BC2[3][2] = { /* x BC */FLAG_GHOST_LAYER, FLAG_OBSTACLE,
				  	  	  	/* y BC */FLAG_OBSTACLE,FLAG_OBSTACLE,
				  	  	  	/* z BC */FLAG_OBSTACLE,FLAG_OBSTACLE};

		  lbmController1.setBC(BC1);
		  lbmController2.setBC(BC2);

		  std::list<int> lbm_opencl_number_of_registers_list;
		  std::list<int> lbm_opencl_number_of_threads_list;

		  if (!number_of_threads_string.empty())
			  extract_comma_separated_integers(lbm_opencl_number_of_threads_list, number_of_threads_string);

		  if (!number_of_registers_string.empty())
			  extract_comma_separated_integers(lbm_opencl_number_of_registers_list, number_of_registers_string);
		  int my_rank, num_procs;
		  MPI_Init(&argc, &argv);    /// Start MPI
		  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// Get current process id
		  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processes

		  if ( my_rank == 0)
		  lbmController1.run(	debug,
				  //domain_1,
				  gravitation,
				  viscosity,
				  computation_kernel_count,
				  device_nr,
				  gui,
				  pause,
				  timestep,
				  take_frame_screenshots,
				  steps,

				  lbm_opencl_number_of_threads_list,
				  lbm_opencl_number_of_registers_list
		  );

		  if ( my_rank == 1 )
		  lbmController2.run(	debug,
				  gravitation,
				  viscosity,
				  computation_kernel_count,
				  device_nr,
				  gui,
				  pause,
				  timestep,
				  take_frame_screenshots,
				  steps,

				  lbm_opencl_number_of_threads_list,
				  lbm_opencl_number_of_registers_list
		  );
		  MPI_Finalize();    /// Cleanup MPI
	  }


}

