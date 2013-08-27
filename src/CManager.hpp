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
template<typename T>
class CManager
{
private:
	CDomain<T> _domain; ///< The simulation domain.
	CVector<3, int> _subdomain_size; ///< Each subdomain have the same size which is specified with this class member.
	CVector<3, T> _subdomain_length; ///< Each subdomain have the same lengthes which is specified with this class member.
	CVector<3, int> _subdomain_nums; ///< number of subdomains in each direction.
	CController<T>* _lbm_controller;

public:

	CManager(CDomain<T> domain, CVector<3, int> subdomainNums) :
			_domain(domain), _lbm_controller(NULL) {
		this->setSubdomainNums(subdomainNums);
	}

	~CManager() {
		if (_lbm_controller)
			delete _lbm_controller;
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
		CVector<3, int> do_size = _domain.getSize();
		if ((do_size[0] % subdomainNums[0] != 0)
				|| (do_size[1] % subdomainNums[1] != 0)
				|| (do_size[2] % subdomainNums[2] != 0)) {
			throw "Number of subdomains does not match with the grid size!";
		}

		CVector<3, int> tmpSD_size;
		tmpSD_size[0] = do_size[0] / subdomainNums[0];
		tmpSD_size[1] = do_size[1] / subdomainNums[1];
		tmpSD_size[2] = do_size[2] / subdomainNums[2];

		_subdomain_size = tmpSD_size;
		_subdomain_nums = subdomainNums;

		// subdomain lengths
		CVector<3, T> domain_length = _domain.getLength();
		_subdomain_length[0] = domain_length[0] / _subdomain_nums[0];
		_subdomain_length[1] = domain_length[1] / _subdomain_nums[1];
		_subdomain_length[2] = domain_length[2] / _subdomain_nums[2];

#if DEBUG
		std::cout << "NUMBER OF SUBDOMAINS: " << _subdomain_nums << std::endl;
		std::cout << "SUBDOMAIN_SIZE: " << _subdomain_size << std::endl;
		std::cout << "SUBDOMAIN_LENGTHS: " << _subdomain_length << std::endl;
#endif
	}

	void initSimulation(int my_rank) {
		// initialize the boundary condition
		int BC[3][2] = { /* x BC */{ FLAG_GHOST_LAYER, FLAG_GHOST_LAYER },
		/* y BC */{ FLAG_GHOST_LAYER, FLAG_GHOST_LAYER },
		/* z BC */{ FLAG_GHOST_LAYER, FLAG_GHOST_LAYER } };
		int id = my_rank;
		if (id < 0)
			id = 0;

		int tmpid = id;
		int nx, ny, nz;
		nx = tmpid % _subdomain_nums[0];
		tmpid /= _subdomain_nums[0];
		ny = tmpid % _subdomain_nums[1];
		tmpid /= _subdomain_nums[1];
		nz = tmpid;
#if DEBUG
		std::cout << "ID: "<< id << " NX: " << nx << " NY: " << ny << " NZ: " << nz << std::endl;
#endif
		// create the subdomains instances for the whole domain
		CVector<3, int> origin(nx * _subdomain_size[0], ny * _subdomain_size[1],
				nz * _subdomain_size[2]);
		CDomain<T> *subdomain = new CDomain<T>(id, _subdomain_size, origin,
				_subdomain_length);

		// Setting the boundary conditions for the current Controller
		if (nx == 0)
			BC[0][0] = FLAG_OBSTACLE;
		if (nx == (_subdomain_nums[0] - 1))
			BC[0][1] = FLAG_OBSTACLE;

		if (ny == 0)
			BC[1][0] = FLAG_OBSTACLE;
		if (ny == (_subdomain_nums[1] - 1))
			BC[1][1] = FLAG_OBSTACLE;

		if (nz == 0)
			BC[2][0] = FLAG_OBSTACLE;
		if (nz == (_subdomain_nums[2] - 1))
			BC[2][1] = FLAG_OBSTACLE;

		_lbm_controller = new CController<T>(id, *subdomain, BC);

		// Initializing the Controller's communication classes based on the already computed boundary conditions
		if (BC[0][0] == FLAG_GHOST_LAYER) {
			int comm_destination = id - 1;
			CVector<3, int> send_size(1, _subdomain_size[1],
					_subdomain_size[2]);
			CVector<3, int> recv_size(1, _subdomain_size[1],
					_subdomain_size[2]);
			CVector<3, int> send_origin(1, 0, 0);
			CVector<3, int> recv_origin(0, 0, 0);
			CVector<3, int> comm_direction(1, 0, 0);
			_lbm_controller->addCommunication(
					new CComm<T>(comm_destination, send_size, recv_size,
							send_origin, recv_origin, comm_direction));
		}
		if (BC[0][1] == FLAG_GHOST_LAYER) {
			int comm_destination = id + 1;
			CVector<3, int> send_size(1, _subdomain_size[1],
					_subdomain_size[2]);
			CVector<3, int> recv_size(1, _subdomain_size[1],
					_subdomain_size[2]);
			CVector<3, int> send_origin(_subdomain_size[0] - 2, 0, 0);
			CVector<3, int> recv_origin(_subdomain_size[0] - 1, 0, 0);
			CVector<3, int> comm_direction(-1, 0, 0);
			_lbm_controller->addCommunication(
					new CComm<T>(comm_destination, send_size, recv_size,
							send_origin, recv_origin, comm_direction));
		}
		if (BC[1][0] == FLAG_GHOST_LAYER) {
			int comm_destination = id - _subdomain_nums[0];
			CVector<3, int> send_size(_subdomain_size[0], 1,
					_subdomain_size[2]);
			CVector<3, int> recv_size(_subdomain_size[0], 1,
					_subdomain_size[2]);
			CVector<3, int> send_origin(0, 1, 0);
			CVector<3, int> recv_origin(0, 0, 0);
			CVector<3, int> comm_direction(0, 1, 0);
			_lbm_controller->addCommunication(
					new CComm<T>(comm_destination, send_size, recv_size,
							send_origin, recv_origin, comm_direction));
		}
		if (BC[1][1] == FLAG_GHOST_LAYER) {
			int comm_destination = id + _subdomain_nums[0];
			CVector<3, int> send_size(_subdomain_size[0], 1,
					_subdomain_size[2]);
			CVector<3, int> recv_size(_subdomain_size[0], 1,
					_subdomain_size[2]);
			CVector<3, int> send_origin(0, _subdomain_size[1] - 2, 0);
			CVector<3, int> recv_origin(0, _subdomain_size[1] - 1, 0);
			CVector<3, int> comm_direction(0, -1, 0);
			_lbm_controller->addCommunication(
					new CComm<T>(comm_destination, send_size, recv_size,
							send_origin, recv_origin, comm_direction));
		}
		if (BC[2][0] == FLAG_GHOST_LAYER) {
			int comm_destination = id - _subdomain_nums[0] * _subdomain_nums[1];
			CVector<3, int> send_size(_subdomain_size[0], _subdomain_size[1],
					1);
			CVector<3, int> recv_size(_subdomain_size[0], _subdomain_size[1],
					1);
			CVector<3, int> send_origin(0, 0, 1);
			CVector<3, int> recv_origin(0, 0, 0);
			CVector<3, int> comm_direction(0, 0, 1);
			_lbm_controller->addCommunication(
					new CComm<T>(comm_destination, send_size, recv_size,
							send_origin, recv_origin, comm_direction));
		}
		if (BC[2][1] == FLAG_GHOST_LAYER) {
			int comm_destination = id + _subdomain_nums[0] * _subdomain_nums[1];
			CVector<3, int> send_size(_subdomain_size[0], _subdomain_size[1],
					1);
			CVector<3, int> recv_size(_subdomain_size[0], _subdomain_size[1],
					1);
			CVector<3, int> send_origin(0, 0, _subdomain_size[2] - 2);
			CVector<3, int> recv_origin(0, 0, _subdomain_size[2] - 1);
			CVector<3, int> comm_direction(0, 0, -1);
			_lbm_controller->addCommunication(
					new CComm<T>(comm_destination, send_size, recv_size,
							send_origin, recv_origin, comm_direction));
		}
		if (ny == _subdomain_nums[1] - 1) {
			_lbm_controller->setGeometry();
		}

	}

	void startSimulation() {
		if (!_lbm_controller)
			throw "CManager: Initialize the simulation before starting it!";
		_lbm_controller->run();

	}

	CController<T>* getController() const {
		return _lbm_controller;
	}

	void setController(CController<T>* lbmController) {
		_lbm_controller = lbmController;
	}

	CVector<3, int> getSubdomainSize() const {
		return _subdomain_size;
	}
};
#endif
