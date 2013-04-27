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



#include <stdlib.h>
#include <iostream>
#include "CMain.hpp"
#include <list>
#include <string>

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
	CVector<3,int> domain_size(32,32,32);

	CVector<3,T> gravitation(0,-9.81,0);
	T viscosity = 0.001308;
	T timestep = -1.0;
	int steps = -1;

	bool gui = false;
	bool pause = false;
	bool take_frame_screenshots = false;

	size_t computation_kernel_count = 128;

	std::string number_of_registers_string;	///< string storing the number of registers for opencl threads separated with comma
	std::string number_of_threads_string;		///< string storing the number of threads for opencl separated with comma

	int device_nr = 0;

	char optchar;
	while ((optchar = getopt(argc, argv, "x:y:z:d:vr:k:gG:pt:sl:R:T:X:")) > 0)
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
	return -1;

parameter_error_ok:

	CMain<T> lbmain;

	std::list<int> lbm_opencl_number_of_registers_list;
	std::list<int> lbm_opencl_number_of_threads_list;

	if (!number_of_threads_string.empty())
		extract_comma_separated_integers(lbm_opencl_number_of_threads_list, number_of_threads_string);

	if (!number_of_registers_string.empty())
		extract_comma_separated_integers(lbm_opencl_number_of_registers_list, number_of_registers_string);

	return lbmain.run(	debug,
						domain_size,
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
}

