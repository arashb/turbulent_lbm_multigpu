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


#ifndef CCOMPLANGERRORS_HH__
#define CCOMPLANGERRORS_HH__

#include <CL/cl.h>

static const char *cclGetErrorString(cl_int errorNum)
{
	switch(errorNum)
	{
		case CL_SUCCESS:		return "CL_SUCCESS";
		case CL_DEVICE_NOT_FOUND:		return "CL_DEVICE_NOT_FOUND";
		case CL_DEVICE_NOT_AVAILABLE:		return "CL_DEVICE_NOT_AVAILABLE";
		case CL_COMPILER_NOT_AVAILABLE:		return "CL_COMPILER_NOT_AVAILABLE";
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:		return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
		case CL_OUT_OF_RESOURCES:		return "CL_OUT_OF_RESOURCES";
		case CL_OUT_OF_HOST_MEMORY:		return "CL_OUT_OF_HOST_MEMORY";
		case CL_PROFILING_INFO_NOT_AVAILABLE:		return "CL_PROFILING_INFO_NOT_AVAILABLE";
		case CL_MEM_COPY_OVERLAP:		return "CL_MEM_COPY_OVERLAP";
		case CL_IMAGE_FORMAT_MISMATCH:		return "CL_IMAGE_FORMAT_MISMATCH";
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:		return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
		case CL_BUILD_PROGRAM_FAILURE:		return "CL_BUILD_PROGRAM_FAILURE";
		case CL_MAP_FAILURE:		return "CL_MAP_FAILURE";

		case CL_INVALID_VALUE:		return "CL_INVALID_VALUE";
		case CL_INVALID_DEVICE_TYPE:		return "CL_INVALID_DEVICE_TYPE";
		case CL_INVALID_PLATFORM:		return "CL_INVALID_PLATFORM";
		case CL_INVALID_DEVICE:		return "CL_INVALID_DEVICE";
		case CL_INVALID_CONTEXT:		return "CL_INVALID_CONTEXT";
		case CL_INVALID_QUEUE_PROPERTIES:		return "CL_INVALID_QUEUE_PROPERTIES";
		case CL_INVALID_COMMAND_QUEUE:		return "CL_INVALID_COMMAND_QUEUE";
		case CL_INVALID_HOST_PTR:		return "CL_INVALID_HOST_PTR";
		case CL_INVALID_MEM_OBJECT:		return "CL_INVALID_MEM_OBJECT";
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:		return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
		case CL_INVALID_IMAGE_SIZE:		return "CL_INVALID_IMAGE_SIZE";
		case CL_INVALID_SAMPLER:		return "CL_INVALID_SAMPLER";
		case CL_INVALID_BINARY:		return "CL_INVALID_BINARY";
		case CL_INVALID_BUILD_OPTIONS:		return "CL_INVALID_BUILD_OPTIONS";
		case CL_INVALID_PROGRAM:		return "CL_INVALID_PROGRAM";
		case CL_INVALID_PROGRAM_EXECUTABLE:		return "CL_INVALID_PROGRAM_EXECUTABLE";
		case CL_INVALID_KERNEL_NAME:		return "CL_INVALID_KERNEL_NAME";
		case CL_INVALID_KERNEL_DEFINITION:		return "CL_INVALID_KERNEL_DEFINITION";
		case CL_INVALID_KERNEL:		return "CL_INVALID_KERNEL";
		case CL_INVALID_ARG_INDEX:		return "CL_INVALID_ARG_INDEX";
		case CL_INVALID_ARG_VALUE:		return "CL_INVALID_ARG_VALUE";
		case CL_INVALID_ARG_SIZE:		return "CL_INVALID_ARG_SIZE";
		case CL_INVALID_KERNEL_ARGS:		return "CL_INVALID_KERNEL_ARGS";
		case CL_INVALID_WORK_DIMENSION:		return "CL_INVALID_WORK_DIMENSION";
		case CL_INVALID_WORK_GROUP_SIZE:		return "CL_INVALID_WORK_GROUP_SIZE";
		case CL_INVALID_WORK_ITEM_SIZE:		return "CL_INVALID_WORK_ITEM_SIZE";
		case CL_INVALID_GLOBAL_OFFSET:		return "CL_INVALID_GLOBAL_OFFSET";
		case CL_INVALID_EVENT_WAIT_LIST:		return "CL_INVALID_EVENT_WAIT_LIST";
		case CL_INVALID_EVENT:		return "CL_INVALID_EVENT";
		case CL_INVALID_OPERATION:		return "CL_INVALID_OPERATION";
		case CL_INVALID_GL_OBJECT:		return "CL_INVALID_GL_OBJECT";
		case CL_INVALID_BUFFER_SIZE:		return "CL_INVALID_BUFFER_SIZE";
		case CL_INVALID_MIP_LEVEL:		return "CL_INVALID_MIP_LEVEL";
	}
	return "UNKNOWN";
}

#define CL_CHECK_ERROR(val_a)	if ((val_a) != CL_SUCCESS)		\
	{								\
		std::cerr	<< "CL_ERROR: file: '" << __FILE__	\
	       			<< "', line: " << __LINE__		\
	       			<< " - " << cclGetErrorString(val_a) << std::endl;	\
	}
//		exit(val_a);

/*
#define CL_CHECK_ERROR(val_a, val_b)	CL_CHECK_ERROR_((val_a), (val_b), __FILE__, __LINE__)
inline void CL_CHECK_ERROR_(cl_int fun_ret_val)
{
	CL_CHECK_ERROR_(fun_ret_val, CL_SUCCESS);
}
*/

#endif
