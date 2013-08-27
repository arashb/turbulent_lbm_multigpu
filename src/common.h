/*
 * common.h
 *
 *  Created on: May 16, 2013
 *      Author: Arash Bakhtiari
 */

#ifndef COMMON_H_
#define COMMON_H_
#include "mpi.h"
#include <stdio.h>
#include <string>
#include <sstream>
#include <fstream>

#include "CConfiguration.hpp"
#include "Singleton.hpp"

#define FLAG_OBSTACLE	(1 << 0)
#define FLAG_FLUID	(1 << 1)
#define FLAG_VELOCITY_INJECTION	(1 << 2)
#define FLAG_GHOST_LAYER (1 << 3)
#define FLAG_GHOST_LAYER_BETA (FLAG_GHOST_LAYER | (1 << 4))

//#define HALO_SIZE 0
#define BENCHMARK_OUTPUT_DIR "output/benchmark"
#define PROFILE_OUTPUT_DIR "output/profile"
#define VTK_OUTPUT_DIR "output/vtk"
#define LOG_OUTPUT_DIR "output/log"

#define LOG_OUTPUT_FILE_PREFIX "log"

extern CVector<3, int> E0;
extern CVector<3, int> E1;
extern CVector<3, int> E2;
extern CVector<3, int> E3;

extern CVector<3, int> E4;
extern CVector<3, int> E5;
extern CVector<3, int> E6;
extern CVector<3, int> E7;

extern CVector<3, int> E8;
extern CVector<3, int> E9;
extern CVector<3, int> E10;
extern CVector<3, int> E11;

extern CVector<3, int> E12;
extern CVector<3, int> E13;
extern CVector<3, int> E14;
extern CVector<3, int> E15;

extern CVector<3, int> E16;
extern CVector<3, int> E17;
extern CVector<3, int> E18;

extern CVector<3, int> lbm_units[];

typedef Singleton<CConfiguration<T> > ConfigSingleton; // Global declaration
typedef Singleton<CProfiler> ProfilerSingleton; // Global declaration

typedef enum {
	OPENCL_VERSION_UNKNOWN = 0,
	OPENCL_VERSION_1_0_0 = 100,
	OPENCL_VERSION_1_1_0 = 110,
	OPENCL_VERSION_1_2_0 = 120,
} OPENCL_VERSION;

typedef enum {
	MPI_COMM_DIRECTION_UNKNOWN = 0,
	MPI_COMM_DIRECTION_X,
	MPI_COMM_DIRECTION_X_0,
	MPI_COMM_DIRECTION_X_1,
	MPI_COMM_DIRECTION_Y,
	MPI_COMM_DIRECTION_Y_0,
	MPI_COMM_DIRECTION_Y_1,
	MPI_COMM_DIRECTION_Z,
	MPI_COMM_DIRECTION_Z_0,
	MPI_COMM_DIRECTION_Z_1,
	MPI_COMM_DIRECTION_ALL,
} MPI_COMM_DIRECTION;

const char* get_string_direction(MPI_COMM_DIRECTION direction) {
	std::string dir;
	switch (direction) {
	case MPI_COMM_DIRECTION_X:
		dir = "X";
		break;
	case MPI_COMM_DIRECTION_X_0:
		dir = "X0";
		break;
	case MPI_COMM_DIRECTION_X_1:
		dir = "X1";
		break;
	case MPI_COMM_DIRECTION_Y:
		dir = "Y";
		break;
	case MPI_COMM_DIRECTION_Y_0:
		dir = "Y0";
		break;
	case MPI_COMM_DIRECTION_Y_1:
		dir = "Y1";
		break;
	case MPI_COMM_DIRECTION_Z:
		dir = "Z";
		break;
	case MPI_COMM_DIRECTION_Z_0:
		dir = "Z0";
		break;
	case MPI_COMM_DIRECTION_Z_1:
		dir = "Z1";
		break;
	default:
		dir = "UNKNOWN";
		break;
	}
	return dir.c_str();
}

#ifndef LOG_TO_FILE
#define  DEBUGPRINT(...) \
  { \
  int my_rank; \
  int num_procs;                           \
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs); \
  fprintf(stderr, "P %d: ", my_rank); \
  fprintf(stderr, __VA_ARGS__); \
  }
#else
#define  DEBUGPRINT(...)			\
  { \
  int my_rank; \
  int num_procs;                           \
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); \
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs); \
  std::string outputfilename = LOG_OUTPUT_DIR; \
  std::stringstream ss_file; \
  ss_file << "./"  << LOG_OUTPUT_DIR  << "/" << LOG_OUTPUT_FILE_PREFIX << "_" << my_rank << ".log";	\
  std::string outputfile = ss_file.str(); \
  FILE* file = fopen(outputfile.c_str(),"a");	\
  fprintf(file, "P %d: ", my_rank); \
  fprintf(file, __VA_ARGS__); \
file.close() \
  }
#endif

static const char * mpiGetErrorString(cl_int errorNum) {
	switch (errorNum) {
	case MPI_SUCCESS:
		return "CL_SUCCESS";
	case MPI_ERR_REQUEST:
		return "MPI_ERR_REQUEST";
	case MPI_ERR_ARG:
		return "MPI_ERR_ARG";
	}
	return "UNKNOWN";
}
#define MPI_CHECK_ERROR(val_a)	if ((val_a) != MPI_SUCCESS)		\
	{								\
		std::cerr	<< "MPI_ERROR: file: '" << __FILE__	\
	       			<< "', line: " << __LINE__		\
	       			<< " - " << mpiGetErrorString(val_a) << std::endl;	\
	}

#endif /* COMMON_H_ */
