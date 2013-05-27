#ifndef __CCONFIGURATION_HPP__
#define __CCONFIGURATION_HPP__

#include "tinyxml2.h"

/*
 * Class CConfiguration stores the necessary information for the simulation process.
 *
 */
template <typename T>
class CConfiguration
{
private:
	tinyxml2::XMLDocument doc;

public:
	 CVector<3,T> gravitation;		///< Specify the gravitation vector
	 T viscosity;
	 size_t computation_kernel_count;
	 int device_nr;
	 bool do_visualization;
	 T timestep;
	 int loops;
	 std::list<int> lbm_opencl_number_of_registers_list;
	 std::list<int> lbm_opencl_number_of_threads_list;
	 bool debug_mode;

	CConfiguration()
	{

	}

	CConfiguration( CVector<3,T> _gravitation,		///< Specify the gravitation vector
			T _viscosity,
			size_t _computation_kernel_count,
			int _device_nr,
			bool _do_visualization,
			T _timestep,
			int _loops
	) :
		gravitation(_gravitation),
		viscosity(_viscosity),
		computation_kernel_count(_computation_kernel_count),
		device_nr(_device_nr),
		do_visualization(_do_visualization),
		timestep(_timestep),
		loops(_loops)
		{

		}

	~CConfiguration() {

	}

	int load_file(std::string file_name)
	{
		doc.LoadFile( file_name.c_str() );
		return doc.ErrorID();
	}
};

#endif
