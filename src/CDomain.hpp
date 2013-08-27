#ifndef __CDOMAIN_HPP__
#define __CDOMAIN_HPP__

#include <stdlib.h>
#include <iostream>
#include "libmath/CVector.hpp"

// simulation type
typedef float T;
//typedef double T;

/*
 * Class CDomain contains data related to a subdomain of simulation grid
 *
 */

template<typename T>
class CDomain {
	int _UID; ///< Unique id of the sub/domain
	CVector<3, int> _size; ///< Size of the domain
	CVector<3, int> _origin_cell; ///< Origin of the domain points if it is part of a bigger domain
	CVector<3, T> _length; ///< Length of the domain in each direction

public:
	CDomain(int UID, CVector<3, int> size, CVector<3, int> origin_cell,
			CVector<3, T> length) :
			_UID(UID), _size(size), _origin_cell(origin_cell), _length(length) {

	}
	CDomain(int UID, CVector<3, int> size) :
			_UID(UID), _size(size) {
		_origin_cell = CVector<3, int>(0, 0, 0);
		_length = CVector<3, T>(0.05, 0.05, 0.05);
	}

	~CDomain() {

	}

public:
	CVector<3, int> getOrigin() const {
		return _origin_cell;
	}

	CVector<3, int> getSize() const {
		return _size;
	}

	int getUid() const {
		return _UID;
	}

	CVector<3, T> getLength() const {
		return _length;
	}

};

#endif
