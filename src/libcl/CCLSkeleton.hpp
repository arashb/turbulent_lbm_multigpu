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


/*
 * CCLSkeleton.hpp
 *
 *  Created on: Jan 3, 2010
 *      Author: martin
 */

#ifndef CCLSKELETON_HPP_
#define CCLSKELETON_HPP_

#include "lib/CError.hpp"
#include "libcl/CCL.hpp"

/**
 * \brief OpenCL skeleton
 *
 * this class can be used to handle the OpenCL access variables/ids
 * conveniently to share them with other classes (e.g. simulation and gui)
 */
class CCLSkeleton
{
public:
    CError error;		///< error handler
    bool verbose;		///< verbose informations on/off

	CCL::CDevice cDevice;			///< OpenCL device for computations
	CCL::CPlatforms cPlatforms;		///< OpenCL platforms
	CCL::CDevices cDevices;			///< OpenCL devices
	CCL::CCommandQueue cCommandQueue;	///< OpenCL command queue for computations
	CCL::CContext cContext;			///< OpenCL context for computations

	/**
	 * initialize skeleton
	 */
	CCLSkeleton(bool p_verbose = false)
	{
		verbose = p_verbose;
	}


	/**
	 * initialize skeleton with existing skeleton
	 *
	 * only cDevice, cCommandQueue and cContext are copied
	 */
	CCLSkeleton(	const CCLSkeleton &cClSkeleton,
					bool p_verbose = false
	)
	{
		verbose = p_verbose;
		initCL(cClSkeleton);
	}

	/**
	 * initialize skeleton with existing skeleton
	 *
	 * the ids of cDevice, cCommandQueue and cContext are copied, not re-initialized
	 */
	void initCL(const CCLSkeleton &cClSkeleton)
	{
		cDevice.initWithExistingDevice(cClSkeleton.cDevice);
//		cPlatforms.init(cClSkeleton.cPlatforms);
//		cDevices.init(cClSkeleton.cDevices);
		cCommandQueue.initWithExistingCommandQueue(cClSkeleton.cCommandQueue);
		cContext.initWithExistingContext(cClSkeleton.cContext);
	}

	/**
	 * initialize OpenCL with given parameters
	 */
	void initCL(	const CCL::CContext &p_cContext,		///< existing OpenCL context
					const CCL::CDevice &p_cDevice,			///< existing OpenCL device
					const CCL::CCommandQueue &p_cCommandQueue	///< existing OpenCL command queue
			)
	{
		cDevice.initWithExistingDevice(p_cDevice);
//		cPlatforms.init(cClSkeleton.cPlatforms);
//		cDevices.init(cClSkeleton.cDevices);
		cCommandQueue.initWithExistingCommandQueue(p_cCommandQueue);
		cContext.initWithExistingContext(p_cContext);
	}

	/**
	 * initialize OpenCL skeleton
	 */
	void initCL(	int p_device_nr = -1,		///< device number (use -1 to list available devices)
					int p_platform_nr = -1		///< platform number (use -1 to list available platforms)
	)
	{
		// load platform information
		if (verbose)	std::cout << "loading platforms" << std::endl;
		cPlatforms.load();

		if (cPlatforms.platform_ids_count == 0)
		{
			error << "no platform found!" << std::endl;
			return;
		}

		if (p_platform_nr == -1)
		{
			for (int i = 0; i < (int)cPlatforms.platform_ids_count; i++)
			{
				CCL::CPlatform cPlatform(cPlatforms.platform_ids[i]);
				cPlatform.loadPlatformInfo();

				if (verbose)
				{
					std::cout << "Platform " << (i) << ":" << std::endl;
					std::cout << "        Name: " << cPlatform.name << std::endl;
					std::cout << "     Profile: " << cPlatform.profile << std::endl;
					std::cout << "     Version: " << cPlatform.version << std::endl;
					std::cout << "      Vendor: " << cPlatform.vendor << std::endl;
					std::cout << "  Extensions: " << cPlatform.extensions << std::endl;
					std::cout << std::endl;
				}

				/**
				 * autodetect platform
				 */
				if (p_platform_nr == -1)
					if (strcmp(cPlatform.profile, "FULL_PROFILE") == 0)
						p_platform_nr = i;

				if (verbose && (p_platform_nr == i))
				{
					std::cout << ">>> Using Platform " << (i) << " for simulation" << std::endl;
					std::cout << std::endl;
				}
			}
		}

		if (p_platform_nr == -1)
		{
			error << "no usable platform found... exiting" << std::endl;
			return;
		}

		if (p_platform_nr < 0 || p_platform_nr >= (int)cPlatforms.platform_ids_count)
		{
			error << "invalid platform number - use option \"-P -1\" to list all devices" << std::endl;
			return;
		}

		CCL::CPlatform cPlatform(cPlatforms.platform_ids[p_platform_nr]);

		// load devices belonging to platform
		if (verbose)	std::cout << "loading devices for platform " << p_platform_nr << std::endl;
		cDevices.load(cPlatform);

		if (cDevices.size() == 0)
		{
			error << "no devices found - aborting" << std::endl;
			return;
		}

		if (p_device_nr == -1)
		{
			// list available devices
			for (int i = 0; i < (int)cDevices.size(); i++)
			{
				CCL::CDeviceInfo cDeviceInfo(cDevices[i]);
				if (verbose)
				{
					std::cout << "Device " << (i) << ":" << std::endl;
					std::cout << "        Type: " << cDeviceInfo.getTypeString() << std::endl;
					std::cout << "        Name: " << cDeviceInfo.name << std::endl;
					std::cout << "     Profile: " << cDeviceInfo.profile << std::endl;
					std::cout << "     Version: " << cDeviceInfo.version << std::endl;
					std::cout << "      Vendor: " << cDeviceInfo.vendor << std::endl;
					std::cout << "  Extensions: " << cDeviceInfo.extensions << std::endl;
					std::cout << std::endl;
				}

				if (p_device_nr == -1)
					if (cDeviceInfo.device_type == CL_DEVICE_TYPE_GPU)
						p_device_nr = i;

				if (verbose && (p_device_nr == i))
				{
					std::cout << ">>> Using Device " << (i) << " for simulation" << std::endl;
					std::cout << std::endl;
				}
			}
		}

		if (p_device_nr < 0 || p_device_nr >= (int)cDevices.size())
		{
			error << "invalid device number - use option \"-D -1\" to list all devices" << std::endl;
			return;
		}

		cDevice.set(cDevices[p_device_nr]);
//		CCL::CDevice &cDevice = cDevices[p_device_nr];

		// load standard context for GPU devices
		if (verbose)	std::cout << "loading gpu context" << std::endl;
		cContext.load(cPlatform, cDevice);

		// initialize queue
		if (verbose)	std::cout << "creating command queue" << std::endl;
		cCommandQueue.init(cContext, cDevice);
	}



#ifdef C_GL_TEXTURE_HPP

	/**
	 * initialize OpenCL context from existing OpenGL context
	 */
	void initCLFromGLContext(	int p_device_nr = -1,		///< device number (use -1 to list available devices)
								int p_platform_nr = -1		///< platform number (use -1 to list available platforms)
	)
	{
		// load platform information
		if (verbose)	std::cout << "loading platforms" << std::endl;
		cPlatforms.load();

		if (cPlatforms.platform_ids_count == 0)
		{
			error << "no platform found!" << std::endl;
			return;
		}

		if (p_platform_nr == -1)
		{
			for (int i = 0; i < (int)cPlatforms.platform_ids_count; i++)
			{
				CCL::CPlatform cPlatform(cPlatforms.platform_ids[i]);
				cPlatform.loadPlatformInfo();

				if (verbose)
				{
					std::cout << "Platform " << (i) << ":" << std::endl;
					std::cout << "        Name: " << cPlatform.name << std::endl;
					std::cout << "     Profile: " << cPlatform.profile << std::endl;
					std::cout << "     Version: " << cPlatform.version << std::endl;
					std::cout << "      Vendor: " << cPlatform.vendor << std::endl;
					std::cout << "  Extensions: " << cPlatform.extensions << std::endl;
					std::cout << std::endl;
				}

				/**
				 * autodetect platform
				 */
				if (p_platform_nr == -1)
					if (strcmp(cPlatform.profile, "FULL_PROFILE") == 0)
						p_platform_nr = i;

				if (verbose && (p_platform_nr == i))
				{
					std::cout << ">>> Using Platform " << (i) << " for simulation" << std::endl;
					std::cout << std::endl;
				}
			}
		}

		if (p_platform_nr == -1)
		{
			error << "no usable platform found... exiting" << std::endl;
			return;
		}

		if (p_platform_nr < 0 || p_platform_nr >= (int)cPlatforms.platform_ids_count)
		{
			error << "invalid platform number - use option \"-P -1\" to list all devices" << std::endl;
			return;
		}

		CCL::CPlatform cPlatform(cPlatforms.platform_ids[p_platform_nr]);

		if (!cContext.loadCurrentGlContext(cPlatform))
			error << cContext.error.getString();

		cDevices.load(cContext);

		if (cDevices.size() == 0)
		{
			error << "no devices found - aborting" << std::endl;
			return;
		}

		if (p_device_nr == -1)
		{
			// list available devices
			for (int i = 0; i < (int)cDevices.size(); i++)
			{
				CCL::CDeviceInfo cDeviceInfo(cDevices[i]);
				if (verbose)
				{
					std::cout << "Device " << (i) << ":" << std::endl;
					std::cout << "        Type: " << cDeviceInfo.getTypeString() << std::endl;
					std::cout << "        Name: " << cDeviceInfo.name << std::endl;
					std::cout << "     Profile: " << cDeviceInfo.profile << std::endl;
					std::cout << "     Version: " << cDeviceInfo.version << std::endl;
					std::cout << "      Vendor: " << cDeviceInfo.vendor << std::endl;
					std::cout << "  Extensions: " << cDeviceInfo.extensions << std::endl;
					std::cout << std::endl;
				}

				if (p_device_nr == -1)
					if (cDeviceInfo.device_type == CL_DEVICE_TYPE_GPU)
						p_device_nr = i;

				if (verbose && (p_device_nr == i))
				{
					std::cout << ">>> Using Device " << (i) << " for simulation" << std::endl;
					std::cout << std::endl;
				}
			}
		}

		if (p_device_nr < 0 || p_device_nr >= (int)cDevices.size())
		{
			error << "invalid device number - use option \"-D -1\" to list all devices" << std::endl;
			return;
		}

		cDevice.set(cDevices[p_device_nr]);
//		CCL::CDevice &cDevice = cDevices[p_device_nr];

		// initialize queue
		if (verbose)	std::cout << "creating command queue" << std::endl;
		cCommandQueue.init(cContext, cDevice);
	}
#endif
};

#endif /* CCLSKELETON_HPP_ */
