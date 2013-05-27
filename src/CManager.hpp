
#ifndef __CMANAGER_HPP__
#define __CMANAGER_HPP__


#include <vector>
#include <map>

#include "CDomain.hpp"
#include "CController.hpp"
#include "libmath/CVector.hpp"

/*
 * Class CManager is responsible for dividing and assigning the subdomains to different processors.
 *
 */
template <typename T>
class CManager
{
private:
	CDomain<T> _domain; 					///< The simulation domain.
	CVector<3,int> _subdomain_size; 			///< Each subdomain have the same size which is specified with this class member.
	CVector<3,T> _subdomain_length;			///< Each subdomain have the same lengthes which is specified with this class member.
	CVector<3,int> _subdomain_nums; 	///< number of subdomains in each direction.
	//std::map< int, CDomain<T>* > _subdomains_container;
	//std::map< int, CController<T>* > _controller_container; ///< A container for each subdomain simulation controller
	CController<T>* _lbm_controller;

public:

	CManager(CDomain<T> domain, CVector<3, int> subdomainNums):
		_domain(domain),
	//	_subdomain_nums(1,1,1),
		_lbm_controller(NULL)
	{
		this->setSubdomainNums(subdomainNums);
	}

	~CManager() {
	}

	CDomain<T> getDomain() const {
		return _domain;
	}

	void setDomain(CDomain<T> grid) {
		_domain = grid;
	}

	CVector<3, int> getSubdomainNums() const {
		return _subdomain_nums;
	}

	void setSubdomainNums(CVector<3, int> subdomainNums) {
		CVector<3,int> do_size = _domain.getSize();
		if ( (do_size[0] % subdomainNums[0] != 0) ||
			 (do_size[1] % subdomainNums[1] != 0) ||
			 (do_size[2] % subdomainNums[2] != 0)
			 ) {
			throw "Number of subdomains does not match with the grid size!";
		}

		CVector<3,int> tmpSD_size;
		tmpSD_size[0] = do_size[0] / subdomainNums[0];
		tmpSD_size[1] = do_size[1] / subdomainNums[1];
		tmpSD_size[2] = do_size[2] / subdomainNums[2];

		if ( (tmpSD_size[0] & 2 ) &&
			 (tmpSD_size[1] & 2 ) &&
			 (tmpSD_size[2] & 2 )
			 )
		{
			throw "Subdomain sizes should be power of 2!";
		}
		_subdomain_size = tmpSD_size;
		_subdomain_nums = subdomainNums;

		// subdomain lengths
		CVector<3,T> domain_length = _domain.getLength();
		_subdomain_length[0] = domain_length[0] / _subdomain_nums[0];
		_subdomain_length[1] = domain_length[1] / _subdomain_nums[1];
		_subdomain_length[2] = domain_length[2] / _subdomain_nums[2];

#if DEBUG
		std::cout << "NUMBER OF SUBDOMAINS: " << _subdo_subdomain_nums << std::endl;
		std::cout << "SUBDOMAIN_SIZE: " << _subdomain_size << std::endl;
		std::cout << "SUBDOMAIN_LENGTHS: " << _subdomain_length << std::endl;
#endif
	}


	void initSimulation(int my_rank) {
		// initialize the boundary condition
		int BC[3][2] = { /* x BC */FLAG_GHOST_LAYER,FLAG_GHOST_LAYER,
				/* y BC */FLAG_GHOST_LAYER,FLAG_GHOST_LAYER,
				/* z BC */FLAG_GHOST_LAYER,FLAG_GHOST_LAYER};

		int id = 0;
		// TODO: OPTIMIZATION: for huge number of subdomains three nested loops is slow.
		for( int nz = 0; nz < _subdomain_nums[2]; nz++) {
			for ( int ny = 0; ny < _subdomain_nums[1]; ny++ ) {
				for( int nx = 0; nx < _subdomain_nums[0]; nx++)
				{
					if ( id == my_rank)
					{
						// create the subdomains instances for the whole domain
						CVector<3,int> origin(nx*_subdomain_size[0],ny*_subdomain_size[1],nz*_subdomain_size[2]);
						CDomain<T> *subdomain = new CDomain<T>(id, _subdomain_size, origin, _subdomain_length);
						//_subdomains_container[id] = subdomain;

						// Setting the boundary conditions for the current Controller
						if ( nx == 0 )
							BC[0][0] = FLAG_OBSTACLE;
						if (nx == ( _subdomain_nums[0] - 1 ))
							BC[0][1] = FLAG_OBSTACLE;

						if ( ny == 0 )
							BC[1][0] = FLAG_OBSTACLE;
						if (ny == ( _subdomain_nums[1] - 1 ))
							BC[1][1] = FLAG_OBSTACLE;

						if ( nz == 0 )
							BC[2][0] = FLAG_OBSTACLE;
						if (nz == ( _subdomain_nums[2] - 1 ))
							BC[2][1] = FLAG_OBSTACLE;

						_lbm_controller = new CController<T>(id,*subdomain);
						_lbm_controller->setBC(BC);

						// Initializing the Controller's communication classes based on the already computed boundary conditions
						if (BC[0][0] == FLAG_GHOST_LAYER) {
							int comm_destination = id - 1;
							CVector<3,int> send_size(1,_subdomain_size[1],_subdomain_size[2]);
							CVector<3,int> recv_size(1,_subdomain_size[1],_subdomain_size[2]);
							CVector<3,int> send_origin(1,0,0);
							CVector<3,int> recv_origin(0,0,0);
							CVector<3,int> comm_direction(1,0,0);
							_lbm_controller->addCommunication(new CComm<T>(comm_destination,send_size,recv_size,send_origin,recv_origin,comm_direction));
						}
						if (BC[0][1] == FLAG_GHOST_LAYER) {
							int comm_destination = id + 1;
							CVector<3,int> send_size(1,_subdomain_size[1],_subdomain_size[2]);
							CVector<3,int> recv_size(1,_subdomain_size[1],_subdomain_size[2]);
							CVector<3,int> send_origin(_subdomain_size[0] - 2, 0, 0);
							CVector<3,int> recv_origin(_subdomain_size[0] - 1, 0, 0);
							CVector<3,int> comm_direction(-1,0,0);
							_lbm_controller->addCommunication(new CComm<T>(comm_destination,send_size,recv_size,send_origin,recv_origin,comm_direction));
						}
						if (BC[1][0] == FLAG_GHOST_LAYER) {
							int comm_destination = id - _subdomain_nums[0];
							CVector<3,int> send_size(_subdomain_size[0], 1, _subdomain_size[2]);
							CVector<3,int> recv_size(_subdomain_size[0], 1, _subdomain_size[2]);
							CVector<3,int> send_origin(0,1,0);
							CVector<3,int> recv_origin(0,0,0);
							CVector<3,int> comm_direction(0,1,0);
							_lbm_controller->addCommunication(new CComm<T>(comm_destination,send_size,recv_size,send_origin,recv_origin,comm_direction));
						}
						if (BC[1][1] == FLAG_GHOST_LAYER) {
							int comm_destination = id + _subdomain_nums[0];
							CVector<3,int> send_size(_subdomain_size[0], 1, _subdomain_size[2]);
							CVector<3,int> recv_size(_subdomain_size[0], 1, _subdomain_size[2]);
							CVector<3,int> send_origin(0, _subdomain_size[1] - 2, 0);
							CVector<3,int> recv_origin(0, _subdomain_size[1] - 1, 0);
							CVector<3,int> comm_direction(0,-1,0);
							_lbm_controller->addCommunication(new CComm<T>(comm_destination,send_size,recv_size,send_origin,recv_origin,comm_direction));
						}
						if (BC[2][0] == FLAG_GHOST_LAYER) {
							int comm_destination = id - _subdomain_nums[0]*_subdomain_nums[1];
							CVector<3,int> send_size(_subdomain_size[0], _subdomain_size[1], 1);
							CVector<3,int> recv_size(_subdomain_size[0], _subdomain_size[1], 1);
							CVector<3,int> send_origin(0,0,1);
							CVector<3,int> recv_origin(0,0,0);
							CVector<3,int> comm_direction(0,0,1);
							_lbm_controller->addCommunication(new CComm<T>(comm_destination,send_size,recv_size,send_origin,recv_origin,comm_direction));
						}
						if (BC[2][1] == FLAG_GHOST_LAYER) {
							int comm_destination = id + _subdomain_nums[0]*_subdomain_nums[1];
							CVector<3,int> send_size(_subdomain_size[0], _subdomain_size[1], 1);
							CVector<3,int> recv_size(_subdomain_size[0], _subdomain_size[1], 1);
							CVector<3,int> send_origin(0, 0, _subdomain_size[2] - 2);
							CVector<3,int> recv_origin(0, 0, _subdomain_size[2] - 1);
							CVector<3,int> comm_direction(0,0,-1);
							_lbm_controller->addCommunication(new CComm<T>(comm_destination,send_size,recv_size,send_origin,recv_origin,comm_direction));
						}
						return;
					}
					id++;
				}
			}
		}
	}

	// TODO: implement this function to use the configuration singleton
	void startSimulation() {
		if(!_lbm_controller)
			throw "CManager: Initialize the simulation before starting it!";

//		_lbm_controller->run(	debug,
//				gravitation,
//				viscosity,
//				computation_kernel_count,
//				device_nr,
//				gui,
//				pause,
//				timestep,
//				take_frame_screenshots,
//				steps,
//
//				lbm_opencl_number_of_threads_list,
//				lbm_opencl_number_of_registers_list
//		);

	}

	void startSimulation(
			bool debug, 				///< Set this variable to true to have a verbose output of simulation process.
			//CDomain<T> domain, 				///< Specify domain properties
			CVector<3,T> gravitation,		///< Specify the gravitation vector
			T viscosity,
			size_t computation_kernel_count,
			int device_nr,
			bool do_visualization,
			bool pause,
			T timestep,
			bool take_frame_screenshots,
			int steps,

			std::list<int> &lbm_opencl_number_of_threads_list,		///< List with number of threads for each successively created kernel
			std::list<int> &lbm_opencl_number_of_registers_list		///< List with number of registers for each thread threads for each successively created kernel
			)
	{
		if(!_lbm_controller)
			throw "CManager: Initialize the simulation before starting it!";
		_lbm_controller->run(	debug,
				gravitation,
				viscosity,
				computation_kernel_count,
				device_nr,
				do_visualization,
				pause,
				timestep,
				take_frame_screenshots,
				steps,

				lbm_opencl_number_of_threads_list,
				lbm_opencl_number_of_registers_list
		);
	}
};
#endif