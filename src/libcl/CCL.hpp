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


#ifndef CCOMPLANG_HPP
#define CCOMPLANG_HPP

#include <CL/cl.h>
#include <iostream>
#include <strings.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include "CCLErrors.hpp"
#include "lib/CError.hpp"
#include "lib/CFile.hpp"
#include "libtools/CProfiler.hpp"
#include "libmath/CVector.hpp"
#include "../common.h"
#include <sys/types.h>

#ifdef C_GL_TEXTURE_HPP
	#if GL3_PROTOTYPES
		#define __gl_h_
	#endif
	#include <CL/cl_gl.h>

	#ifndef CL_GL_CONTEXT_KHR
		#define		CL_GL_CONTEXT_KHR		0x2008
		#define		CL_EGL_DISPLAY_KHR		0x2009
		#define		CL_GLX_DISPLAY_KHR		0x200A
		#define		CL_WGL_HDC_KHR			0x200B
		#define		CL_CGL_SHAREGROUP_KHR	0x200C
	#endif

	#include <GL/glx.h>
#endif


/**
 * \brief OpenCL C++ abstraction class
 */
class CCL
{
public:

	/**
	 * \brief Load and get information about an OpenCL platform
	 */
	class CPlatform
	{
	public:
		cl_platform_id platform_id;	///< OpenCL platform handler

		char *profile;	///< profile string
		char *version;	///< version string
		char *name;		///< name string
		char *vendor;	///< vendor string
		char *extensions;	///< extensions string

private:
		/**
		 * initialize data structure with default values (NULL)
		 */
		inline void initCPlatform()
		{
			profile = NULL;
			version = NULL;
			name = NULL;
			vendor = NULL;
			extensions = NULL;
		}

		/**
		 * cleanup platform (free strings)
		 */
		inline void cleanupCPlatform()
		{
			delete[] profile;
			delete[] version;
			delete[] name;
			delete[] vendor;
			delete[] extensions;
		}

public:
		/**
		 * load information about platform and store to class
		 */
		inline void loadPlatformInfo()
		{
			size_t size_retval;
#define loadP(flag, variable)	\
			CL_CHECK_ERROR(clGetPlatformInfo(	platform_id,	flag,	0,	NULL,	&size_retval));		\
			variable = new char[size_retval];											\
			CL_CHECK_ERROR(clGetPlatformInfo(	platform_id,	flag,	size_retval,	variable,	NULL));	\

			loadP(CL_PLATFORM_PROFILE, profile)
			loadP(CL_PLATFORM_VERSION, version)
			loadP(CL_PLATFORM_NAME, name)
			loadP(CL_PLATFORM_VENDOR, vendor)
			loadP(CL_PLATFORM_EXTENSIONS, extensions)
#undef load
		}

		/**
		 * setup platform id to load information for platform_id
		 */
		inline void load(cl_platform_id p_platform_id)
		{
			cleanupCPlatform();
			initCPlatform();
			platform_id = p_platform_id;
		}

		/**
		 * load information about Platform given by platform_id
		 */
		inline CPlatform(cl_platform_id p_platform_id)
		{
			initCPlatform();
			load(p_platform_id);
		}

		inline CPlatform()
		{
			initCPlatform();
		}

		inline ~CPlatform()
		{
			cleanupCPlatform();
		}
	};


	/**
	 * \brief Handle an OpenCL device
	 */
	class CDevice
	{
	public:
		cl_device_id	device_id;		///< OpenCL Device ID

		/**
		 * initialize CDevice class with existing device id
		 */
		inline CDevice(cl_device_id p_device_id)
		{
			device_id = p_device_id;
		}
		/**
		 * initialize CDevice class with existing class CDevice
		 */
		inline CDevice(const CDevice &cDevice)
		{
			device_id = cDevice.device_id;
		}

		/**
		 * default constructor
		 */
		inline CDevice()
		{
			device_id = 0;
		}

		/**
		 * initialize CDevice class with existing class CDevice
		 */
		inline void initWithExistingDevice(const CDevice &cDevice)
		{
			device_id = cDevice.device_id;
		}

		/**
		 * set device_id to p_device_id
		 */
		inline void set(cl_device_id p_device_id)
		{
			device_id = p_device_id;
		}

		/**
		 * set device_id to device_id in cDevice
		 */
		inline void set(const CDevice &cDevice)
		{
			device_id = cDevice.device_id;
		}
	};

	class CContext;


	/**
	 * \brief Load a list of devices for different OpenCL contexts, platforms and/or types
	 */
	class CDevices	: public std::vector<CDevice>
	{
	public:
		cl_context context;			///< OpenCL context
		cl_device_id *device_ids;	///< array with devices available for the context
		uint device_ids_count;		///< number of contexts in device_ids[]

		/**
		 * load device list belonging to context cContext
		 * \param cContext	existing context
		 */
		inline void load(	CContext &cContext
		)
		{
			// load device information
			size_t value_size_ret;
			CL_CHECK_ERROR(clGetContextInfo(cContext.context, CL_CONTEXT_DEVICES, 0, NULL, &value_size_ret));
			device_ids_count = value_size_ret / sizeof(cl_device_id);
			if (device_ids_count == 0)	std::cerr << "Warning: no device found!" << std::endl;
			delete device_ids;
			device_ids = new cl_device_id[device_ids_count];

			CL_CHECK_ERROR(clGetContextInfo(cContext.context, CL_CONTEXT_DEVICES, value_size_ret, device_ids, NULL));

			initDeviceVector();
		}


		/**
		 * load device list of type 'device type' belonging to platform cPlatform
		 */
		inline void load(	CPlatform &cPlatform,
							cl_device_type device_type = CL_DEVICE_TYPE_ALL
		)
		{
			// read out device ids for given platform
			CL_CHECK_ERROR(clGetDeviceIDs(	cPlatform.platform_id, device_type, 0, NULL, &device_ids_count));

			delete device_ids;
			device_ids = new cl_device_id[device_ids_count];

			CL_CHECK_ERROR(clGetDeviceIDs(	cPlatform.platform_id, device_type, device_ids_count, device_ids, NULL));

			initDeviceVector();
		}

		/**
		 * initialize device list with NULL
		 */
		inline void initCDevices()
		{
			device_ids = NULL;
			device_ids_count = 0;
			clear();
		}

		/**
		 * initialize devices
		 */
		inline CDevices()
		{
			initCDevices();
		}

		/**
		 * load device list belonging to cContext
		 */
		inline CDevices(CContext &cContext)
		{
			initCDevices();
			load(cContext);
		}

		/**
		 * load device list of type 'device type' belonging to platform cPlatform
		 */
		inline CDevices(CPlatform &cPlatform, cl_device_type device_type = CL_DEVICE_TYPE_ALL)
		{
			initCDevices();
			load(cPlatform, device_type);
		}

		/**
		 * deconstructor: free unnecessary device data
		 */
		inline ~CDevices()
		{
			if (device_ids != NULL)
				delete[] device_ids;
			clear();
		}

	private:
		inline void initDeviceVector()
		{
			//clear();
			(*this).resize(device_ids_count);
			for (cl_uint i = 0; i < device_ids_count; i++)
			{
				(*this)[i].set(device_ids[i]);
			}
		}
/*
		inline CDevices operator=(const CDevices& c)
		{
			return c;
		}
*/
	};



	/**
	 * \brief Create and manage an OpenCL context
	 */
	class CContext
	{
	public:
		CError error;		///< error handler
		cl_context context;	///< OpenCL context id

		/**
		 * load context and device list by platform and type
		 */
		inline void load(	const CPlatform &cPlatform,		///< platform for parameters
							cl_device_type device_type		///< type of context and device list
		)
		{
			release();

			cl_int err_ret;
			cl_context_properties context_properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)cPlatform.platform_id, 0, 0};
			context = clCreateContextFromType(context_properties, device_type, NULL, NULL, &err_ret);
			CL_CHECK_ERROR(err_ret);

			retain();
		}

		/**
		 * load context by platform and device
		 */
		inline void load(	const CPlatform &cPlatform,	///< platform for parameters
							const CDevice &cDevice		///< OpenCL device
		)
		{
			release();

			cl_int err_ret;
			cl_context_properties context_properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)cPlatform.platform_id, 0, 0};
			context = clCreateContext(context_properties, 1, &cDevice.device_id, NULL, NULL, &err_ret);
			CL_CHECK_ERROR(err_ret);

			retain();
		}

		/**
		 * load context by platform and device list
		 */
		inline void load(	const CPlatform &cPlatform,	///< platform for parameters
							const CDevices &cDevices	///< device list available for context
		)
		{
			release();

			cl_int err_ret;
			cl_context_properties context_properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)cPlatform.platform_id, 0, 0};
			context = clCreateContext(context_properties, cDevices.device_ids_count, cDevices.device_ids, NULL, NULL, &err_ret);
			CL_CHECK_ERROR(err_ret);

			retain();
		}

#ifdef C_GL_TEXTURE_HPP

		/**
		 * load context from currently running GL context
		 * \param cPlatform	desired CL Platform information
		 */
		inline bool loadCurrentGlContext(CPlatform &cPlatform)
		{
			release();

			cl_int err_ret;

#ifdef WIN32
	#error TODO
#elifdef MACOSX
	#error TODO
#else
			GLXContext glxContext = glXGetCurrentContext();
			if (glxContext == NULL)
			{
				error << "no current glx context found" << std::endl;
				return false;
			}

			Display *display = glXGetCurrentDisplay();
			if (display == NULL)
			{
				error << "no current glx display found" << std::endl;
				return false;
			}

			cl_context_properties context_properties[] = {
							CL_CONTEXT_PLATFORM, (cl_context_properties)cPlatform.platform_id,
							CL_GL_CONTEXT_KHR,	(cl_context_properties)glxContext,
							CL_GLX_DISPLAY_KHR,	(cl_context_properties)display,
							NULL, NULL
					};
#endif
			context = clCreateContextFromType(context_properties, CL_DEVICE_TYPE_GPU, NULL, NULL, &err_ret);
//			context = clCreateContext(context_properties, 0, NULL, NULL, NULL, &err_ret);
			CL_CHECK_ERROR(err_ret);

			retain();
			return true;
		}
#endif

		/**
		 * initialize with NULL context
		 */
		inline void initCContext()
		{
			context = NULL;
		}

		/**
		 * initialize with NULL context
		 */
		inline CContext()
		{
			context = NULL;
		}

		/**
		 * initialize and load context with device_type and platform
		 */
		inline CContext(	const CPlatform &cPlatform,	///< platform for parameters
							cl_device_type device_type	///< type of context and device list
		)
		{
			initCContext();
			load(cPlatform, device_type);
		}

		/**
		 * initialize and load context with device_type and platform
		 */
		inline CContext(	CPlatform &cPlatform,	///< platform for parameters
							CDevice &cDevice		///< OpenCL device
		)
		{
			load(cPlatform, cDevice);
		}


		/**
		 * load context by platform and device list
		 */
		inline CContext(	CPlatform &cPlatform,	///< platform for parameters
							CDevices &cDevices		///< device list available for context
		)
		{
			initCContext();
			load(cPlatform, cDevices);
		}

		inline ~CContext()
		{
			release();
		}


		/**
		 * initialize context handler with existing context
		 */
		inline void initWithExistingContext(const CContext &cContext)
		{
			release();
			context = cContext.context;
			retain();
		}


private:
		/*
		inline CContext(const CContext &)
		{
		}
		*/
/*
		inline const CContext& operator=(const CContext &c)
		{
			return c;
		}
*/
		/**
		 * increment reference counter to context to avoid deletion of context if
		 * multiple classes use the same context
		 */
		inline void retain()
		{
			CL_CHECK_ERROR(clRetainContext(context));
		}

		/**
		 * decrement reference counter to context if context is valid
		 */
		inline void release()
		{
			if (context != NULL)
			{
				clReleaseContext(context);
				context = NULL;
			}
		}
	};

	/**
	 * \brief	handle multiple platforms available for computations
	 */
	class CPlatforms
	{
	public:
		cl_platform_id *platform_ids;	///< array with platform ids
		cl_uint platform_ids_count;		///< number of platforms in platform array

		/**
		 * initialize with NULL values
		 */
		inline void initCPlatforms()
		{
			platform_ids = NULL;
			platform_ids_count = 0;
		}

		/**
		 * load all available platforms
		 */
		inline void load()
		{
			initCPlatforms();

			CL_CHECK_ERROR(clGetPlatformIDs(
						0,
						0,
						&platform_ids_count
					));

			delete platform_ids;
			platform_ids = new cl_platform_id[platform_ids_count];

			CL_CHECK_ERROR(clGetPlatformIDs(
						platform_ids_count,
						platform_ids,
						NULL
					));
		}

		/**
		 * default constructor: load all available platforms
		 */
		inline CPlatforms()
		{
			initCPlatforms();
		}

		inline ~CPlatforms()
		{
			if (platform_ids != NULL)
				delete platform_ids;
		}
	};

	/**
	 * \brief handler for OpenCL buffer objects
	 *
	 * a buffer object is initialized via the CContextDevices class
	 */
	class CMem
	{
	public:
		cl_mem memobj;	///< OpenCL memory object handler

		inline CMem()
		{
			memobj = NULL;
		}

		/**
		 * initialize reference of memory object with existing reference
		 */
		inline CMem(const CMem &bo)
		{
			memobj = bo.memobj;
		}

		/**
		 * create OpenCL memory buffer
		 */
		inline CMem(	CContext&		cContext,	///< context for buffer
						cl_mem_flags	flags,		///< OpenCL flags
						size_t			size,		///< Size of memory buffer
						void*			host_ptr	///< Host pointer to initialize buffer
			)
		{
			memobj = NULL;
			create(cContext, flags, size, host_ptr);
		}

		/**
		 * release memory object (decrements OpenCL reference counter)
		 */
		void release()
		{
			if (memobj != NULL)
				clReleaseMemObject(memobj);
			memobj = NULL;
		}

		/**
		 * deconstructor (release)
		 */
		inline ~CMem()
		{
			release();
		}

		/**
		 * return the size of the memory object
		 */
		inline size_t getSize()
		{
			size_t mem_size;
			clGetMemObjectInfo(memobj, CL_MEM_SIZE, sizeof(size_t), &mem_size, NULL);
			return mem_size;
		}

		/**
		 * create OpenCL 2D image memory object
		 */
		inline void createImage2D(
									CContext&		cContext,	///< context for buffer
									cl_mem_flags	flags,		///< OpenCL flags
									const cl_image_format	*image_format,	///< image format
									size_t	image_width,		///< width of image
									size_t	image_height,		///< height of image
									size_t	image_row_pitch = 0,	///< row pitch to improve aligned access
									void 	*host_ptr = NULL	///< host pointer to initalize with existing data
			)
		{
			release();

			cl_int errcode_ret;
			memobj = clCreateImage2D(	cContext.context, flags, image_format,
										image_width, image_height, image_row_pitch,
										host_ptr, &errcode_ret);

			CL_CHECK_ERROR(errcode_ret);
		}

		/**
		 * create new memory object
		 */
		inline void create(
						CContext&		cContext,	///< context for buffer
						cl_mem_flags	flags,		///< OpenCL flags
						size_t			size,		///< Size of memory buffer
						void*			host_ptr	///< Host pointer to initialize buffer
			)
		{
			release();

			cl_int errcode_ret;
			memobj = clCreateBuffer(cContext.context, flags, size, host_ptr, &errcode_ret);
			CL_CHECK_ERROR(errcode_ret);
		}
/*
 * CGlTexture has to be included before CCL.hpp to enable the creation of
 * OpenCL memory handlers from OpenGL textures
 */
#ifdef C_GL_TEXTURE_HPP
		/**
		 * create OpenCL memory object from existing 3D OpenGL texture
		 */
		inline void createFromGLTexture3D(
				CContext&		cContext,	///< context for buffer
				cl_mem_flags	flags,		///< OpenCL flags
				CGlTexture 		&texture	///< reference to 3d texture
		)
		{
			cl_int errcode_ret;
			memobj = clCreateFromGLTexture3D(cContext.context, flags, GL_TEXTURE_3D, 0, texture.textureid, &errcode_ret);
			CL_CHECK_ERROR(errcode_ret);
		}

		/**
		 * create OpenCL memory object from existing 2D OpenGL texture
		 */
		inline void createFromGLTexture2D(
				CContext&		cContext,	///< context for buffer
				cl_mem_flags	flags,		///< OpenCL flags
				CGlTexture 		&texture	///< reference to 3d texture
		)
		{
			cl_int errcode_ret;
			memobj = clCreateFromGLTexture2D(cContext.context, flags, texture.target, 0, texture.textureid, &errcode_ret);
			CL_CHECK_ERROR(errcode_ret);
		}

#endif


	};


	/**
	 * \brief Load, use and manage OpenCL programs to create kernels
	 */
	class CProgram
	{
public:
		CError error;			///< error handler
		cl_program program;		///< OpenCL program id
		std::string filepath;	///< string with filename (for debugging purposes)

		inline CProgram()
		{
			program = 0;
		}

		/**
		 * load program given by p_filepath
		 */
		inline CProgram(	CContext &cContext,	// context
					const std::string &p_filepath	// source file
		)
		{
			CProgram();
			load(cContext, p_filepath);
		}

		/**
		 * load program given by p_filepath and place prefix_string before the filecontent
		 */
		inline CProgram(	CContext &cContext,	// context
					const std::string &p_filepath,	// source file
					const std::string &prefix_string	// prefix placed before source file
		)
		{
			CProgram();
			load(cContext, p_filepath, prefix_string);
		}

		/**
		 * initialize program with source code
		 */
		inline CProgram(	CContext &cContext,		// context
							cl_uint count,			// number of strings
							const char **strings,	// source code
							const size_t *lengths	// length of source code strings. if NULL, the strings are \0 terminated
		)
		{
			CProgram();
			load(cContext, count, strings, lengths);
		}

		/**
		 * load kernel using source given in strings
		 */
		inline void load(	CContext &cContext,	// context
							cl_uint count,		// number of strings
							const char **strings,	// source code
							const size_t *lengths	// length of source code strings. if NULL, the strings are \0 terminated
		)
		{
			cl_int errcode_ret;
			program = clCreateProgramWithSource(cContext.context, count, strings, lengths, &errcode_ret);
			if (errcode_ret != CL_SUCCESS)
				error << cclGetErrorString(errcode_ret) << std::endl;
		}


		/**
		 * load kernel source from file and prepent prefix_string
		 */
		inline void load(	CContext &cContext,					// context
							const std::string &sPrefixString,	// prefix placed before source file
							const std::string &p_filepath		// source file
		)
		{
			filepath = p_filepath;

			std::string fileContent;
			std::string errorLog;
#if 1
			std::string source = sPrefixString;
			source += "\n";
			source += "#include \"";
			source += filepath;
			source += "\"";

			const char* strings[] = {source.c_str()};
			size_t lengths[] = {source.size()};
			load(cContext, 1, strings, lengths);
#else
			if (!CFile::fileContents(filepath, fileContent, errorLog))
			{
				error << "CProgram: Failed to open file '" << filepath << "' - " << strerror(errno) << std::endl;
				return;
			}

			const char* strings[] = {sPrefixString.c_str(), fileContent.c_str()};
			size_t lengths[] = {sPrefixString.size(), fileContent.size()};
			load(cContext, 2, strings, lengths);
#endif
		}

		/**
		 * load kernel source from file
		 */
		inline void load(	CContext &cContext,				// context
							const std::string &p_filepath	// source file
		)
		{
			filepath = p_filepath;

			std::string fileContent;
			std::string errorLog;

			CFile cFile;
			cFile.fileContents(filepath, fileContent, errorLog);

			if (cFile.error())
			{
				std::cerr << "CProgram: Failed to open file '" << filepath << "' - " << cFile.error.getString() << std::endl;
				return;
			}

			const char* strings[] = {fileContent.c_str()};
			size_t lengths[] = {fileContent.size()};
			load(cContext, 1, strings, lengths);
		}

		/**
		 * return build information string
		 */
		std::string getBuildInfo(	CDevice &device		///< device on which the code was build for
								)
		{
			char *build_log;
			size_t ret_val_size;
//			CL_CHECK_ERROR(clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size));
			clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);
			build_log = new char[ret_val_size+1];
//			CL_CHECK_ERROR(clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL));
			clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);

			// to be carefully, terminate with \0
			// there's no information in the reference wheter the string is 0 terminated or not
			build_log[ret_val_size] = '\0';

			std::string ret_string = build_log;
			delete build_log;
			return ret_string;
		}

		/**
		 * build code for device
		 */
		inline void build(	CDevice &device,			///< device to build code for
							const char* options = ""	///< options for build operation
		)
		{
			if (error())
				return;

			// build
			cl_int ret_val = clBuildProgram(program, 1, &device.device_id, options, NULL, NULL);
			//cl_int ret_val = clBuildProgram(program, 0, NULL, options, NULL, NULL);

			// avoid abortion due to CL_BILD_PROGRAM_FAILURE
			if (ret_val != CL_SUCCESS)
			{
				error << cclGetErrorString(ret_val) << std::endl;
				error << getBuildInfo(device) << std::endl;
				return;
			}

			cl_build_status build_status;
//			CL_CHECK_ERROR(clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL));
			clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &build_status, NULL);
			if (build_status == CL_SUCCESS)
				return;

			char *build_log;
			size_t ret_val_size;
//			CL_CHECK_ERROR(clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size));
			clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &ret_val_size);
			build_log = new char[ret_val_size+1];
//			CL_CHECK_ERROR(clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL));
			clGetProgramBuildInfo(program, device.device_id, CL_PROGRAM_BUILD_LOG, ret_val_size, build_log, NULL);

			// to be carefully, terminate with \0
			// there's no information in the reference wheter the string is 0 terminated or not
			build_log[ret_val_size] = '\0';

			std::cout << "BUILD LOG: '" << filepath << "'" << std::endl;
			std::cout << build_log << std::endl;

			delete[] build_log;
		}

		/**
		 * output OpenCL binary code
		 */
		void printBinaries()
		{
			cl_uint program_num_devices;
			CL_CHECK_ERROR(clGetProgramInfo(	program,
								CL_PROGRAM_NUM_DEVICES,
								sizeof(cl_uint),
								&program_num_devices,
								NULL
					));

			if (program_num_devices == 0)
			{
				std::cerr << "no valid binary was found" << std::endl;
				return;
			}

			size_t *binaries_sizes = new size_t[program_num_devices];

			CL_CHECK_ERROR(clGetProgramInfo(	program,
								CL_PROGRAM_BINARY_SIZES,
								program_num_devices*sizeof(size_t),
								binaries_sizes,
								NULL
					));

			char **binaries = new char*[program_num_devices];

			for (size_t i = 0; i < program_num_devices; i++)
				binaries[i] = new char[binaries_sizes[i]+1];

			CL_CHECK_ERROR(clGetProgramInfo(program, CL_PROGRAM_BINARIES, program_num_devices*sizeof(size_t), binaries, NULL));

			for (size_t i = 0; i < program_num_devices; i++)
			{
				binaries[i][binaries_sizes[i]] = '\0';
				std::cout << "Program " << i << ":" << std::endl;
				std::cout << binaries[i];
			}

			for (size_t i = 0; i < program_num_devices; i++)
				delete [] binaries[i];

			delete [] binaries;
			delete [] binaries_sizes;
		}
	};


	/**
	 * \brief Create kernels from CProgram, use and set arguments for them
	 */
	class CKernel
	{
	public:
		cl_kernel kernel;	///< OpenCL kernel handler
		CError error;		///< error handler

		/**
		 * create kernel from OpenCL program
		 */
		inline void create(	CProgram &program,			///< OpenCL program
							const char *kernel_name		///< name of kernel function to create kernel from
			)
		{
			if (error())
				std::cout << "FUCK" << std::endl;
			cl_int errcode_ret;
			kernel = clCreateKernel(program.program, kernel_name, &errcode_ret);
			if (errcode_ret != CL_SUCCESS)
				error << cclGetErrorString(errcode_ret) << std::endl;
		}

		/**
		 * create kernel from OpenCL program
		 */
		inline CKernel(	CProgram &program,			///< OpenCL program
						const char *kernel_name	///< name of kernel function to create kernel from
		)
		{
			create(program, kernel_name);
		}

		/**
		 * create NULL kernel
		 */
		inline CKernel()
		{
			kernel = NULL;
		}

		/**
		 * create kernel from existing kernel
		 */
		inline CKernel(const CKernel &k)
		{
			kernel = k.kernel;
		}

		/**
		 * default deconstructor:
		 * release kernel
		 */
		inline ~CKernel()
		{
			if (kernel != NULL)
				CL_CHECK_ERROR(clReleaseKernel(kernel));
		}

		/**
		 * set memory object kernel argument
		 */
		inline void setArg(	cl_uint arg_index,	///< argument number
							CMem &cmem			///< OpenCL memory object
		)
		{
			cl_int ret_val = clSetKernelArg(kernel, arg_index, sizeof(cl_mem), &(cmem.memobj));

			if (ret_val != CL_SUCCESS)
				error << cclGetErrorString(ret_val) << std::endl;
		}

		/**
		 * set cl_float kernel argument
		 */
		inline void setArg(	cl_uint arg_index,	///< argument number
							cl_float &arg		///< reference to float argument
		)
		{
			cl_int ret_val = clSetKernelArg(kernel, arg_index, sizeof(float), &arg);

			if (ret_val != CL_SUCCESS)
				error << cclGetErrorString(ret_val) << std::endl;
		}

		/**
		 * set cl_double kernel argument
		 */
		inline void setArg(	cl_uint arg_index,	///< argument number
							cl_double &arg		///< reference to float argument
		)
		{
			cl_int ret_val = clSetKernelArg(kernel, arg_index, sizeof(double), &arg);

			if (ret_val != CL_SUCCESS)
				error << cclGetErrorString(ret_val) << std::endl;
		}

		/**
		 * set cl_int kernel argument
		 */
		inline void setArg(cl_uint arg_index, cl_int size)
		{
			cl_int ret_val = clSetKernelArg(kernel, arg_index, sizeof(cl_int), &size);

			if (ret_val != CL_SUCCESS)
				error << cclGetErrorString(ret_val) << std::endl;
		}

		/**
		 * return the maximum work group size
		 * \param cDevice	device where the kernel should be executed
		 */
		inline size_t getMaxWorkGroupSize(CDevice &cDevice)
		{
			size_t param_value;
			clGetKernelWorkGroupInfo(kernel, cDevice.device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &param_value, NULL);
			return param_value;
		}
	};

	/**
	 * \brief OpenCL event handler
	 */
	class CEvent
	{
	public:
		cl_event event;	///< openCL event

		inline CEvent()
		{
			event = 0;
		}

		/**
		 * wait for event to complete
		 */
		inline void waitAndRelease()
		{
			if (event != 0)
			{
				CL_CHECK_ERROR(clWaitForEvents(1, &event));
				release();
				event = 0;
			}
		}

		/**
		 * return the execution status belonging to the event
		 */
		inline cl_int getExecutionStatus()
		{
			cl_int param_value;
			CL_CHECK_ERROR(clGetEventInfo(event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &param_value, NULL));
			return param_value;
		}

		inline ~CEvent()
		{
			if (event != 0)
				clReleaseEvent(event);
		}

		/**
		 * release the event
		 */
		inline void release()
		{
			clReleaseEvent(event);
			event = 0;
		}
	};

	/**
	 * \brief OpenCL event list management (up to 10 events)
	 */
	class CEvents
	{
	public:
		cl_event events[10];	///< events (max: 10)
		cl_uint events_count;	///< number of events in events list

		/**
		 * initialize empty event list
		 */
		inline CEvents()
		{
			events_count = 0;
		}

		/**
		 * initialize with 1 event
		 */
		inline CEvents(CEvent &event0)
		{
			set(event0);
		}

		/**
		 * initialize with 2 events
		 */
		inline CEvents(CEvent &event0, CEvent &event1)
		{
			set(event0, event1);
		}

		/**
		 * initialize with 1 event
		 */
		inline void set(CEvent &event0)
		{
			events_count = 1;
			events[0] = event0.event;
		}

		/**
		 * initialize with 2 events
		 */
		inline void set(CEvent &event0, CEvent &event1)
		{
			events_count = 2;
			events[0] = event0.event;
			events[1] = event1.event;
		}
	};


	/**
	 * \brief OpenCL command queue handler
	 */
	class CCommandQueue
	{
	public:
		CError error;		///< error handler

		cl_command_queue command_queue;	///< OpenCL command queue handler

		/**
		 * increment OpenCL reference counter to command queue
		 */
		inline void retain()
		{
			CL_CHECK_ERROR(clRetainCommandQueue(command_queue));
		}

		/**
		 * decrement OpenCL reference counter to command queue
		 */
		inline void release()
		{
			if (command_queue != 0)
			{
				CL_CHECK_ERROR(clReleaseCommandQueue(command_queue));
			}
		}

		/**
		 * initialize existing command queue
		 */
		inline void init(	const CContext &cContext,	///< existing context handler
							const CDevice &cDevice		///< existing device handler
		)
		{
		  cl_command_queue_properties properties = 0;
#if PROFILE
		  properties |= CL_QUEUE_PROFILING_ENABLE;
		  CCL::CDeviceInfo cDeviceInfo(cDevice);
		  printf("device timer resolution: %zu nanoseconds\n", cDeviceInfo.profiling_timer_resolution);
#endif
			cl_int errcode_ret;
			command_queue = clCreateCommandQueue(
									cContext.context,
									cDevice.device_id,
									properties,
									&errcode_ret
							);

			CL_CHECK_ERROR(errcode_ret);
			retain();
		}


		/**
		 * initialize from existing command queue and increment reference counter
		 */
		void initWithExistingCommandQueue(const CCommandQueue &cCommandQueue)
		{
			command_queue = cCommandQueue.command_queue;
			retain();
		}

		/**
		 * initialize command queue from existing context and device
		 */
		inline CCommandQueue(	const CContext &cContext,	///< existing context handler
								const CDevice &cDevice		///< existing device handler
		)
		{
			init(cContext, cDevice);
		}

		inline CCommandQueue()
		{
			command_queue = 0;
		}

		inline ~CCommandQueue()
		{
			release();
		}

		/**
		 * enqueue barrier to command queue
		 */
		inline void enqueueBarrier()
		{
			cl_int errcode_ret;
			errcode_ret = clEnqueueBarrier(command_queue);
			if (errcode_ret != CL_SUCCESS)
				error << cclGetErrorString(errcode_ret) << std::endl;
		}
	

		/**
		 * copy data from buffer in host memory to CL device
		 */
		inline void enqueueWriteBuffer(	CMem &cMem,
						cl_bool block_write,
						size_t offset,
						size_t buffer_size,
						const void *buffer_ptr
		)
		{
			CL_CHECK_ERROR(	clEnqueueWriteBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
		}

		/**
		 * enqueue writing a buffer
		 */
		inline void enqueueWriteBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										const void *buffer_ptr,	///< host memory pointer
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list		///< list of events
		)
		{
			CL_CHECK_ERROR(	clEnqueueWriteBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								num_events_in_wait_list,
								event_wait_list,
								NULL
					));
		}

		/**
		 * enqueue writing a buffer
		 */
		inline void enqueueWriteBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										const void *buffer_ptr,	///< host memory pointer
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list,	///< list of events
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueWriteBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								num_events_in_wait_list,
								event_wait_list,
								&event.event
					));
		}

		/**
		 * enqueue writing a buffer
		 */
		inline void enqueueWriteBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										const void *buffer_ptr,	///< host memory pointer
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueWriteBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								0,
								NULL,
								&event.event
					));
		}

		/**
		 * copy data from buffer in host memory to CL device
		 */
		inline void enqueueWriteBufferRect(	CMem &cMem,
						cl_bool block_write,
						const size_t buffer_origin[3],
						const size_t host_origin[3],
						const size_t region[3],
						const void *buffer_ptr
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueWriteBufferRect(	command_queue,
								cMem.memobj,
								block_write,
								buffer_origin,
								host_origin,
								region,
								0,
								0,
								0,
								0,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
#endif
		}
		inline void enqueueWriteBufferRect(	CMem &cMem,
						cl_bool block_write,
						const size_t buffer_origin[3],
						const size_t host_origin[3],
						const size_t region[3],
						const size_t buffer_row_pitch,
						const size_t buffer_slice_pitch,
						const size_t host_row_pitch,
						const size_t host_slice_pitch,
						const void *buffer_ptr
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueWriteBufferRect(	command_queue,
								cMem.memobj,
								block_write,
								buffer_origin,
								host_origin,
								region,
								buffer_row_pitch,
								buffer_slice_pitch,
								host_row_pitch,
								host_slice_pitch,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
#endif
		}
		/**
		 * read from CL device to buffer located in host memory
		 */
		inline void enqueueReadBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										void *buffer_ptr		///< host memory pointer
		)
		{
			CL_CHECK_ERROR(	clEnqueueReadBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
		}


		/**
		 * read from CL device to buffer located in host memory
		 */
		inline void enqueueReadBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										void *buffer_ptr,		///< host memory pointer
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list		///< list of events
		)
		{
			CL_CHECK_ERROR(	clEnqueueReadBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								num_events_in_wait_list,
								event_wait_list,
								NULL
					));
		}


		/**
		 * read from CL device to buffer located in host memory
		 */
		inline void enqueueReadBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										void *buffer_ptr,		///< host memory pointer
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list,	///< list of events
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueReadBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								num_events_in_wait_list,
								event_wait_list,
								&event.event
					));
		}


		/**
		 * read from CL device to buffer located in host memory
		 */
		inline void enqueueReadBuffer(	CMem &cMem,				///< memory object
										cl_bool block_write,	///< don't return until data is written
										size_t offset,			///< offset to start writing from
										size_t buffer_size,		///< size of buffer to write
										void *buffer_ptr,		///< host memory pointer
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueReadBuffer(	command_queue,
								cMem.memobj,
								block_write,
								offset,
								buffer_size,
								buffer_ptr,
								0,
								NULL,
								&event.event
					));
		}

		inline void enqueueReadBufferRect(	CMem &cMem,
				cl_bool block_write,
				const size_t buffer_origin[3],
				const size_t host_origin[3],
				const size_t region[3],
				void *buffer_ptr
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueReadBufferRect(	command_queue,
								cMem.memobj,
								block_write,
								buffer_origin,
								host_origin,
								region,
								0,
								0,
								0,
								0,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
#endif

		}

		inline void enqueueReadBufferRect(	CMem &cMem,
				cl_bool block_write,
				const size_t buffer_origin[3],
				const size_t host_origin[3],
				const size_t region[3],
				const size_t buffer_row_pitch,
				const size_t buffer_slice_pitch,
				const size_t host_row_pitch,
				const size_t host_slice_pitch,
				void *buffer_ptr
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueReadBufferRect(	command_queue,
								cMem.memobj,
								block_write,
								buffer_origin,
								host_origin,
								region,
								buffer_row_pitch,
								buffer_slice_pitch,
								host_row_pitch0,
								host_slice_pitch,
								buffer_ptr,
								0,
								NULL,
								NULL
					));
#endif

		}

		/**
		 * copy buffer to buffer on device
		 */
		inline void enqueueCopyBuffer(	CMem &cSrcMem,		///< source memory object
										CMem &cDstMem,		///< destination memory object
										size_t src_offset,	///< offset in source memory
										size_t dst_offset,	///< offset in destination memory
										size_t cb			///< number of bytes to be copied
		)
		{
			CL_CHECK_ERROR(	clEnqueueCopyBuffer(
								command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_offset,
								dst_offset,
								cb,
								0,
								NULL,
								NULL
					));
		}

		/**
		 * copy buffer to buffer on device
		 */
		inline void enqueueCopyBuffer(	CMem &cSrcMem,		///< source memory object
										CMem &cDstMem,		///< destination memory object
										size_t src_offset,	///< offset in source memory
										size_t dst_offset,	///< offset in destination memory
										size_t cb,			///< number of bytes to be copied
										cl_uint num_events_in_wait_list,	///< number of events in waiting list
										const cl_event *event_wait_list,	///< list of events
										CEvent &event			///< event for enqueueWriteBuffer
		)
		{
			CL_CHECK_ERROR(	clEnqueueCopyBuffer(
								command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_offset,
								dst_offset,
								cb,
								num_events_in_wait_list,
								event_wait_list,
								&event.event
					));
		}

		/**
		 * copy buffer to buffer on device
		 */
		inline void enqueueCopyBufferRect(	CMem &cSrcMem,		///< source memory object
										CMem &cDstMem,		///< destination memory object
										const size_t src_origin[3],
										const size_t dst_origin[3],
										const size_t region[3]
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueCopyBufferRect(	command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_origin,
								dst_origin,
								region,
								0,
								0,
								0,
								0,
								0,
								NULL,
								NULL
					));
#endif
		}

		inline void enqueueCopyBufferRect(	CMem &cSrcMem,		///< source memory object
										CMem &cDstMem,		///< destination memory object
										const size_t src_origin[3],
										const size_t dst_origin[3],
										const size_t region[3],
										const size_t src_row_pitch,
										const size_t src_slice_pitch,
										const size_t dst_row_pitch,
										const size_t dst_slice_pitch
		)
		{
#ifdef CL_VERSION_1_1
			CL_CHECK_ERROR(	clEnqueueCopyBufferRect(	command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_origin,
								dst_origin,
								region,
								src_row_pitch,
								src_slice_pitch,
								dst_row_pitch,
								dst_slice_pitch,
								0,
								NULL,
								NULL
					));
#endif
		}
		/**
		 * copy buffer to buffer on device
		 */
		inline void enqueueCopyBufferToImage(	CMem &cSrcMem,		///< source memory object
												CMem &cDstMem,		///< destination image object
												size_t src_offset,	///< offset in source memory
												const size_t dst_origin[3],	///< coordinates in source image
												const size_t region[3]		///< area to be copied
		)
		{
			CL_CHECK_ERROR(	clEnqueueCopyBufferToImage(
								command_queue,
								cSrcMem.memobj,
								cDstMem.memobj,
								src_offset,
								dst_origin,
								region,
								0,
								NULL,
								NULL
					));
		}

		/**
		 * enqueue nd range kernel
		 */
		inline void enqueueNDRangeKernel(	CKernel &cKernel,		///< enqueue a OpenCL kernel
											cl_uint work_dim,		///< number of work dimensions (0, 1 or 2)
											const size_t *global_work_offset,	///< global work offset
											const size_t *global_work_size,		///< global work size
											const size_t *local_work_size		///< local work size
		)
		{
		  cl_event* event = NULL;
#if PROFILE
		  event = new cl_event();
#endif
			CL_CHECK_ERROR(	clEnqueueNDRangeKernel(	command_queue,
								cKernel.kernel,
								work_dim,
								global_work_offset,
								global_work_size,
								local_work_size,
								0,
								NULL,
								event
					));
#if PROFILE
			clWaitForEvents(1, event);
			char kernel_name[128];
			CL_CHECK_ERROR( clGetKernelInfo (cKernel.kernel,
							 CL_KERNEL_FUNCTION_NAME,
							 128,
							 kernel_name,
							 NULL
							 ));
			CProfilerEvent* profEvent = new CProfilerEvent(kernel_name, event);
			ProfilerSingleton::Instance()->addProfilerEvent(profEvent);
			delete event;
#endif
		}

		/**
		 * enqueue nd range kernel
		 */
		inline void enqueueNDRangeKernel(	CKernel &cKernel,		///< enqueue a OpenCL kernel
											cl_uint work_dim,		///< number of work dimensions (0, 1 or 2)
											const size_t *global_work_offset,	///< global work offset
											const size_t *global_work_size,		///< global work size
											const size_t *local_work_size,		///< local work size
											CEvents &events			///< events to wait for
		)
		{
			CL_CHECK_ERROR(	clEnqueueNDRangeKernel(	command_queue,
								cKernel.kernel,
								work_dim,
								global_work_offset,
								global_work_size,
								local_work_size,
								events.events_count,
								events.events,
								NULL
					));
		}

		/**
		 * enqueue nd range kernel
		 */
		inline void enqueueNDRangeKernel(	CKernel &cKernel,		///< enqueue a OpenCL kernel
											cl_uint work_dim,		///< number of work dimensions (0, 1 or 2)
											const size_t *global_work_offset,	///< global work offset
											const size_t *global_work_size,		///< global work size
											const size_t *local_work_size,		///< local work size
											CEvent &event			///< event to wait for
		)
		{
			CL_CHECK_ERROR(	clEnqueueNDRangeKernel(	command_queue,
								cKernel.kernel,
								work_dim,
								global_work_offset,
								global_work_size,
								local_work_size,
								0,
								NULL,
								&event.event
					));
		}

		/**
		 * wait until all enqueued object from command queue are finished
		 */
		inline void finish()
		{
			CL_CHECK_ERROR(	clFinish(command_queue));
		}

#ifdef C_GL_TEXTURE_HPP
		/**
		 * acquire memory object from OpenGL context
		 */
		inline void enqueueAcquireGLObject(CMem &cMem)
		{
			CL_CHECK_ERROR(	clEnqueueAcquireGLObjects(command_queue, 1, &(cMem.memobj), 0, NULL, NULL));
		}

		/**
		 * release memory object from OpenGL context
		 */
		inline void enqueueReleaseGLObject(CMem &cMem)
		{
			CL_CHECK_ERROR(	clEnqueueReleaseGLObjects(command_queue, 1, &(cMem.memobj), 0, NULL, NULL));
		}
#endif
	};


	/**
	 * output informations about all profiles and devices
	 */
	inline void printPlatformInfo(cl_device_type device_type = CL_DEVICE_TYPE_ALL)
	{
		CPlatforms platforms;
		std::cout << "platforms_ids: " << platforms.platform_ids_count << std::endl;

		for (cl_uint i = 0; i < platforms.platform_ids_count; i++)
		{
			CPlatform platform(platforms.platform_ids[i]);
			platform.loadPlatformInfo();

			std::cout << "  [" << i << "] PLATFORM_PROFILE: " << platform.profile << std::endl;
			std::cout << "  [" << i << "] PLATFORM_VERSION: " << platform.version << std::endl;

			CDevices cDevices(platform, device_type);

			for (cl_uint d = 0; d < cDevices.device_ids_count; d++)
			{
				CDeviceInfo deviceInfo(cDevices[d]);

				std::ostringstream sPrefix;
			       	sPrefix << "    [" << d << "] ";

				deviceInfo.printDeviceInfo(sPrefix.str());
			}
		}
	}


	/**
	 * \brief load information about a specific device
	 */
	class CDeviceInfo : public CDevice
	{
	public:
		cl_device_type device_type;		///< OpenCL device type
		cl_uint vendor_id;				///< OpenCL vendor id
		cl_uint max_compute_units;		///< maximum compute units available
		cl_uint max_work_item_dimensions;	///< maximum number of dimensions

		size_t *max_work_item_sizes;	///< maximum amount of work items
		size_t max_work_group_size;		///< maximum group size

		cl_uint preferred_vector_width_char;	///< preferred vector width for type char
		cl_uint preferred_vector_width_short;	///< preferred vector width for type short
		cl_uint preferred_vector_width_int;		///< preferred vector width for type int
		cl_uint preferred_vector_width_long;	///< preferred vector width for type long
		cl_uint preferred_vector_width_float;	///< preferred vector width for type float
		cl_uint preferred_vector_width_double;	///< preferred vector width for type float

		cl_uint max_clock_frequency;			///< maximum clock frequency
		cl_uint address_bits;					///< address bits for device

		cl_ulong max_mem_alloc_size;			///< maximum number of allocatable bytes

		cl_bool image_support;					///< image support available
		cl_uint max_read_image_args;			///< maximum number of images as read kernel arguments
		cl_uint max_write_image_args;			///< maximum number of images as write kernel arguments

		size_t image2d_max_width;				///< maximum 2d image width
		size_t image2d_max_height;				///< maximum 2d image height

		size_t image3d_max_width;				///< maximum 3d image width
		size_t image3d_max_height;				///< maximum 3d image height
		size_t image3d_max_depth;				///< maximum 3d image depth

		cl_uint max_samplers;					///< maximum number of samplers
		size_t max_parameter_size;				///< maximum number of kernel parameters
		cl_uint mem_base_addr_align;			///< alignment of device base memory
		cl_uint min_data_type_align_size;		///< minimum alignment needed for device memory

		cl_device_fp_config single_fp_config;	///< single precision floating point capabilities
		cl_device_mem_cache_type global_mem_cache_type;	///< cache type of global memory
		cl_uint global_mem_cacheline_size;		///< size of a line of global memory cache
		cl_ulong global_mem_cache_size;			///< size of global memory cache
		cl_ulong global_mem_size;				///< size of global memory

		cl_ulong max_constant_buffer_size;		///< maximum bytes for constant buffer
		cl_uint max_constant_args;				///< maximum number of constant arguments
		cl_device_local_mem_type local_mem_type;	///< type of local memory
		cl_ulong local_mem_size;				///< size of local memory
		cl_bool error_correction;				///< error correction available
		size_t profiling_timer_resolution;		///< resolution of profiling timer
		cl_bool endian_little;					///< little endian device
		cl_bool available;						///< true, if device available
		cl_bool compiler_available;				///< true, if compiler for device is available
		cl_device_exec_capabilities execution_capabilities;	///< kernel execution capabilities
		cl_command_queue_properties queue_properties;	///< queue properties

		char *name;				///< name of device
		char *vendor;			///< vendor of device
		char *driver_version;	///< driver version of device
		char *profile;			///< profile of device
		char *version;			///< version of device
		char *extensions;		///< extensions available for device

		/**
		 * initialize device information with NULL data
		 */
		inline void initCDeviceInfo()
		{
			max_work_item_sizes = NULL;

			name = NULL;
			vendor = NULL;
			driver_version = NULL;
			profile = NULL;
			version = NULL;
			extensions = NULL;
		}

		inline CDeviceInfo()
			: CDevice()
		{
			initCDeviceInfo();
		}

		/**
		 * initialize device information from existing device
		 */
		inline CDeviceInfo(const CDevice &cDevice)
			: CDevice(cDevice)
		{
			initCDeviceInfo();
			loadDeviceInfo(cDevice);
		}

		inline ~CDeviceInfo()
		{
			delete[] max_work_item_sizes;
			delete[] name;
			delete[] vendor;
			delete[] driver_version;
			delete[] profile;
			delete[] version;
			delete[] extensions;

		}

		/**
		 * return type of the device
		 */
		inline const char* getTypeString()
		{
			switch(device_type)
			{
				case CL_DEVICE_TYPE_CPU:	return "CPU";
				case CL_DEVICE_TYPE_GPU:	return "GPU";
				case CL_DEVICE_TYPE_ACCELERATOR:	return "ACCELERATOR";
				case CL_DEVICE_TYPE_DEFAULT:	return "DEFAULT";
				case CL_DEVICE_TYPE_ALL:	return "ALL";
				default:			return "unknown";
			}
		}

		/**
		 * load device information given by device_id
		 */
		inline void loadDeviceInfo(const CDevice &cDevice)
		{
			set(cDevice.device_id);	// set device id

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_TYPE,			sizeof(cl_device_type),	&device_type, NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_VENDOR_ID,			sizeof(cl_uint),	&vendor_id,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_COMPUTE_UNITS,		sizeof(cl_uint),	&max_compute_units,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,	sizeof(cl_uint),	&max_work_item_dimensions,	NULL));

			delete max_work_item_sizes;
			max_work_item_sizes = new size_t[max_work_item_dimensions];

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_work_group_size),	&max_work_group_size,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*max_work_item_dimensions,	max_work_item_sizes,	NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,	sizeof(size_t),	&preferred_vector_width_char,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,	sizeof(size_t),	&preferred_vector_width_short,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,	sizeof(size_t),	&preferred_vector_width_int,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,	sizeof(size_t),	&preferred_vector_width_long,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,	sizeof(size_t),	&preferred_vector_width_float,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,	sizeof(size_t),	&preferred_vector_width_double,		NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_CLOCK_FREQUENCY,		sizeof(size_t),	&max_clock_frequency,		NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_ADDRESS_BITS,		sizeof(size_t),		&address_bits,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_MEM_ALLOC_SIZE,	sizeof(cl_ulong),	&max_mem_alloc_size,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE_SUPPORT,		sizeof(cl_bool),	&image_support,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_READ_IMAGE_ARGS,	sizeof(cl_uint),	&max_read_image_args,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_WRITE_IMAGE_ARGS,	sizeof(cl_uint),	&max_write_image_args,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE2D_MAX_WIDTH,		sizeof(size_t),		&image2d_max_width,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE2D_MAX_HEIGHT,	sizeof(size_t),		&image2d_max_height,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE3D_MAX_WIDTH,		sizeof(size_t),		&image3d_max_width,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE3D_MAX_HEIGHT,	sizeof(size_t),		&image3d_max_height,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_IMAGE3D_MAX_DEPTH,		sizeof(size_t),		&image3d_max_depth,		NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_SAMPLERS,		sizeof(cl_uint),	&max_samplers,			NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_PARAMETER_SIZE,	sizeof(size_t),		&max_parameter_size,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MEM_BASE_ADDR_ALIGN,	sizeof(cl_uint),	&mem_base_addr_align,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE,	sizeof(cl_uint),	&min_data_type_align_size,	NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_SINGLE_FP_CONFIG,		sizeof(cl_device_fp_config),	&single_fp_config,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_GLOBAL_MEM_CACHE_TYPE,	sizeof(cl_device_mem_cache_type),	&global_mem_cache_type,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE,	sizeof(cl_uint),		&global_mem_cacheline_size,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,	sizeof(cl_ulong),		&global_mem_cache_size,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_GLOBAL_MEM_SIZE,		sizeof(cl_ulong),		&global_mem_size,		NULL));

			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,	sizeof(cl_ulong),	&max_constant_buffer_size,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_MAX_CONSTANT_ARGS,		sizeof(cl_uint),	&max_constant_args,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_LOCAL_MEM_TYPE,		sizeof(cl_device_local_mem_type),	&local_mem_type,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_LOCAL_MEM_SIZE,		sizeof(cl_ulong),	&local_mem_size,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_ERROR_CORRECTION_SUPPORT,	sizeof(cl_bool),	&error_correction,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_PROFILING_TIMER_RESOLUTION,	sizeof(size_t),		&profiling_timer_resolution,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_ENDIAN_LITTLE,		sizeof(cl_bool),	&endian_little,			NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_AVAILABLE,			sizeof(cl_bool),	&available,			NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_COMPILER_AVAILABLE,		sizeof(cl_bool),	&compiler_available,		NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_EXECUTION_CAPABILITIES,	sizeof(cl_device_exec_capabilities),	&execution_capabilities,	NULL));
			CL_CHECK_ERROR(clGetDeviceInfo(	device_id, CL_DEVICE_QUEUE_PROPERTIES,		sizeof(cl_command_queue_properties),	&queue_properties,		NULL));

			loadDeviceInfoString(device_id, CL_DEVICE_NAME, &name);
			loadDeviceInfoString(device_id, CL_DEVICE_VENDOR, &vendor);
			loadDeviceInfoString(device_id, CL_DRIVER_VERSION, &driver_version);
			loadDeviceInfoString(device_id, CL_DEVICE_PROFILE, &profile);
			loadDeviceInfoString(device_id, CL_DEVICE_VERSION, &version);
			loadDeviceInfoString(device_id, CL_DEVICE_EXTENSIONS, &extensions);
		}

		/**
		 * output previously loaded device information
		 */
		inline void printDeviceInfo(std::string sLinePrefix)
		{

			std::cout << sLinePrefix << "TYPE: ";

			if (device_type & CL_DEVICE_TYPE_CPU)	std::cout << "CPU ";
			if (device_type & CL_DEVICE_TYPE_GPU)	std::cout << "GPU ";
			if (device_type & CL_DEVICE_TYPE_ACCELERATOR)	std::cout << "ACCELERATOR ";
			if (device_type & CL_DEVICE_TYPE_DEFAULT)	std::cout << "DEFAULT ";
			std::cout << std::endl;

			std::cout << sLinePrefix << "VENDOR_ID: " << vendor_id << std::endl;
			std::cout << sLinePrefix << "MAX_COMPUTE_UNITS: " << max_compute_units << std::endl;
			std::cout << sLinePrefix << "MAX_WORK_ITEM_DIMENSIONS: " << max_work_item_dimensions << std::endl;

			for (cl_uint w = 0; w < max_work_item_dimensions; w++)
			{
				std::cout << sLinePrefix << "MAX_WORK_ITEM_SIZES[" << w << "]: " << max_work_item_sizes[w] << std::endl;
			}
			std::cout << sLinePrefix << "MAX_WORK_GROUP_SIZE: " << max_work_group_size << std::endl;
			std::cout << sLinePrefix << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_CHAR: " << preferred_vector_width_char << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_SHORT: " << preferred_vector_width_short << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_INT: " << preferred_vector_width_int << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_LONG: " << preferred_vector_width_long << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_FLOAT: " << preferred_vector_width_float << std::endl;
			std::cout << sLinePrefix << "PREFERRED_VECTOR_WIDTH_DOUBLE: " << preferred_vector_width_double << std::endl;
			std::cout << sLinePrefix << std::endl;
			std::cout << sLinePrefix << "MAX_CCLOCK_FREQUENCY: " << max_clock_frequency << std::endl;
			std::cout << sLinePrefix << "ADDRESS_BITS: " << address_bits << std::endl;
			std::cout << sLinePrefix << "MAX_MEM_ALLOC_SIZE: " << max_mem_alloc_size << std::endl;
			std::cout << sLinePrefix << std::endl;
			std::cout << sLinePrefix << "IMAGE_SUPPORT: " << image_support << std::endl;
			std::cout << sLinePrefix << "MAX_READ_IMAGE_ARGS: " << max_read_image_args << std::endl;
			std::cout << sLinePrefix << "MAX_WRITE_IMAGE_ARGS: " << max_write_image_args << std::endl;
			std::cout << sLinePrefix << "IMAGE2D_MAX_WIDTH: " << image2d_max_width << std::endl;
			std::cout << sLinePrefix << "IMAGE2D_MAX_HEIGHT: " << image2d_max_height << std::endl;
			std::cout << sLinePrefix << "IMAGE3D_MAX_WIDTH: " << image3d_max_width << std::endl;
			std::cout << sLinePrefix << "IMAGE3D_MAX_HEIGHT: " << image3d_max_height << std::endl;
			std::cout << sLinePrefix << "IMAGE3D_MAX_DEPTH: " << image3d_max_depth << std::endl;
			std::cout << sLinePrefix << std::endl;
			std::cout << sLinePrefix << "MAX_SAMPLERS: " << max_samplers << std::endl;
			std::cout << sLinePrefix << "MAX_PARAMETER_SIZE: " << max_parameter_size << std::endl;
			std::cout << sLinePrefix << "MEM_BASE_ADDR_ALIGN: " << mem_base_addr_align << std::endl;
			std::cout << sLinePrefix << "MIN_DATA_TYPE_ALIGN_SIZE: " << min_data_type_align_size << std::endl;
			std::cout << sLinePrefix << "SINGLE_FP_CONFIG: ";
			if (single_fp_config & CL_FP_DENORM)	std::cout << "FP_DENORM ";
			if (single_fp_config & CL_FP_INF_NAN)	std::cout << "FP_INF_NAN ";
			if (single_fp_config & CL_FP_ROUND_TO_NEAREST)	std::cout << "FP_ROUND_TO_NEAREST ";
			if (single_fp_config & CL_FP_ROUND_TO_ZERO)	std::cout << "FP_ROUND_TO_ZERO ";
			if (single_fp_config & CL_FP_ROUND_TO_INF)	std::cout << "FP_ROUND_TO_INF ";
			if (single_fp_config & CL_FP_FMA)	std::cout << "FP_FMA ";
			std::cout << std::endl;
			std::cout << sLinePrefix << std::endl;

			std::cout << sLinePrefix << "GLOBAL_MEM_CACHE_TYPE: ";
			switch(global_mem_cache_type)
			{
				case CL_NONE:			std::cout << "NONE";	break;
				case CL_READ_ONLY_CACHE:	std::cout << "CL_READ_ONLY_CACHE";	break;
				case CL_READ_WRITE_CACHE:	std::cout << "CL_READ_WRITE_CACHE";	break;
			}
			std::cout << std::endl;
			std::cout << sLinePrefix << "GLOBAL_MEM_CACHELINE_SIZE: " << global_mem_cacheline_size << std::endl;
			std::cout << sLinePrefix << "GLOBAL_MEM_CACHE_SIZE: " << global_mem_cache_size << std::endl;
			std::cout << sLinePrefix << "GLOBAL_MEM_SIZE: " << global_mem_size << " (" << (global_mem_size >> 20) << "MB)" << std::endl;
			std::cout << sLinePrefix << std::endl;
			std::cout << sLinePrefix << "MAX_CONSTANT_BUFFER_SIZE: " << max_constant_buffer_size << std::endl;
			std::cout << sLinePrefix << "MAX_CONSTANT_ARGS: " << max_constant_args << std::endl;

			std::cout << sLinePrefix << "LOCAL_MEM_TYPE: ";
			switch(local_mem_type)
			{
				case CL_LOCAL:	std::cout << "LOCAL";	break;
				case CL_GLOBAL:	std::cout << "GLOBAL";	break;
				default:	std::cout << "UNKNOWN";	break;
			}
			std::cout << std::endl;

			std::cout << sLinePrefix << "LOCAL_MEM_SIZE: " << local_mem_size << std::endl;
			std::cout << sLinePrefix << "ERROR_CORRECTION_SUPPORT: " << error_correction << std::endl;
			std::cout << sLinePrefix << "PROFILING_TIMER_RESOLUTION: " << profiling_timer_resolution << std::endl;
			std::cout << sLinePrefix << "ENDIAN_LITTLE: " << endian_little << std::endl;
			std::cout << sLinePrefix << "AVAILABLE: " << available << std::endl;
			std::cout << sLinePrefix << "COMPILER_AVAILABLE: " << compiler_available << std::endl;
			std::cout << sLinePrefix << "EXECUTION_CAPABILITIES: ";
			if (execution_capabilities & CL_EXEC_KERNEL)		std::cout << "EXEC_KERNEL ";
			if (execution_capabilities & CL_EXEC_NATIVE_KERNEL)	std::cout << "EXEC_NATIVE_KERNEL ";
			std::cout << std::endl;
			std::cout << sLinePrefix << "QUEUE_PROPERTIES: ";
			if (queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE)		std::cout << "QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE ";
			if (queue_properties & CL_QUEUE_PROFILING_ENABLE)	std::cout << "QUEUE_PROFILING_ENABLE";
			std::cout << std::endl;

			std::cout << sLinePrefix << std::endl;
			std::cout << sLinePrefix << "NAME: " << name << std::endl;
			std::cout << sLinePrefix << "VENDOR: " << vendor << std::endl;
			std::cout << sLinePrefix << "DRIVER_VERSION: " << driver_version << std::endl;
			std::cout << sLinePrefix << "PROFILE: " << profile << std::endl;
			std::cout << sLinePrefix << "VERSION: " << version << std::endl;
			std::cout << sLinePrefix << "EXTENSIONS: " << extensions << std::endl;
		}

private:
		inline void loadDeviceInfoString(	cl_device_id device_id,
						cl_device_info device_info,
						char **param_value)
		{
			size_t retval_size;
			CL_CHECK_ERROR(clGetDeviceInfo(device_id, device_info, 0, NULL, &retval_size));
			delete *param_value;
			*param_value = new char[retval_size];
			CL_CHECK_ERROR(clGetDeviceInfo(device_id, device_info, retval_size, *param_value, NULL));
		}

	};
};


#endif
