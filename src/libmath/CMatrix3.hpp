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


#ifndef __CMATRIX_HH
	#error "dont include CMatrix3.hpp directly - use CMatrix.hpp instead!"
#endif

#ifndef TMATRIX3_H
#define TMATRIX3_H
/**
 * CHANGELOG:
 *
 * 2008-10-28: martin schreiber
 * 	created
 *
 */

#include <limits>
#include <iostream>
#include "CMath.hpp"
#include "CVector.hpp"

template <typename T>
class CMatrix4;

/**
 * \brief	3x3 Matrix Handler
 */
template <typename T>
class CMatrix3
{
public:
	CVector<3,T> matrix[3];	///< matrix array

	/**
	 * constructor: load identity matrix
	 */
	inline CMatrix3()
	{
		loadIdentity();
	}

	/**
	 * set values to 0
	 */
	inline void setZero()
	{
		matrix[0][0] = (T)0;
		matrix[0][1] = (T)0;
		matrix[0][2] = (T)0;

		matrix[1][0] = (T)0;
		matrix[1][1] = (T)0;
		matrix[1][2] = (T)0;

		matrix[2][0] = (T)0;
		matrix[2][1] = (T)0;
		matrix[2][2] = (T)0;
	}


	/**
	 * load identity matrix
	 */
	inline void loadIdentity()
	{
		matrix[0][0] = (T)1;
		matrix[0][1] = (T)0;
		matrix[0][2] = (T)0;

		matrix[1][0] = (T)0;
		matrix[1][1] = (T)1;
		matrix[1][2] = (T)0;

		matrix[2][0] = (T)0;
		matrix[2][1] = (T)0;
		matrix[2][2] = (T)1;
	}


	/**
	 * copy data from matrix m
	 * \param m	source matrix
	 */
	inline CMatrix3<T>& operator=(const CMatrix3<T> &m)
	{
		matrix[0][0] = m.matrix[0][0];
		matrix[0][1] = m.matrix[0][1];
		matrix[0][2] = m.matrix[0][2];

		matrix[1][0] = m.matrix[1][0];
		matrix[1][1] = m.matrix[1][1];
		matrix[1][2] = m.matrix[1][2];

		matrix[2][0] = m.matrix[2][0];
		matrix[2][1] = m.matrix[2][1];
		matrix[2][2] = m.matrix[2][2];

		return *this;
	};

	/**
	 * copy data from matrix m
	 * \param m	source matrix
	 */
	inline CMatrix3<T>& operator=(const CMatrix4<T> &m)
	{
		matrix[0][0] = m.matrix[0][0];
		matrix[0][1] = m.matrix[0][1];
		matrix[0][2] = m.matrix[0][2];

		matrix[1][0] = m.matrix[1][0];
		matrix[1][1] = m.matrix[1][1];
		matrix[1][2] = m.matrix[1][2];

		matrix[2][0] = m.matrix[2][0];
		matrix[2][1] = m.matrix[2][1];
		matrix[2][2] = m.matrix[2][2];
		return *this;
	};

	/**
	 * copy data from matrix m
	 * \param m	source matrix
	 */
	inline CMatrix3(const CMatrix4<T> &m)
	{
		matrix[0][0] = m.matrix[0][0];
		matrix[0][1] = m.matrix[0][1];
		matrix[0][2] = m.matrix[0][2];

		matrix[1][0] = m.matrix[1][0];
		matrix[1][1] = m.matrix[1][1];
		matrix[1][2] = m.matrix[1][2];

		matrix[2][0] = m.matrix[2][0];
		matrix[2][1] = m.matrix[2][1];
		matrix[2][2] = m.matrix[2][2];
	};

	/**
	 * matrix vector product
	 *
	 * use 1 as 3rd component of input vector v
	 *
	 * avoid using parameters by references because this cannot be handled by
	 * subsequenced orders like a+(o*i)
	 */
	inline CVector<2,T>	operator*(CVector<2,T> v)
	{
		return CVector<2,T>(	matrix[0][0]*v[0] + matrix[0][1]*v[1] + matrix[0][2],
								matrix[1][0]*v[0] + matrix[1][1]*v[1] + matrix[1][2]
				);
	}

	/**
	 * matrix vector product
	 */
	inline CVector<3,T>	operator*(CVector<3,T> &v)
	{
		return CVector<3,T>(
					matrix[0][0]*v[0] + matrix[0][1]*v[1] + matrix[0][2]*v[2],
					matrix[1][0]*v[0] + matrix[1][1]*v[1] + matrix[1][2]*v[2],
					matrix[2][0]*v[0] + matrix[2][1]*v[1] + matrix[2][2]*v[2]
				);
	}

	/**
	 * multiply this matrix with m2 and return result
	 * \return this->matrix * m2
	 */
	inline CMatrix3<T>	operator*(CMatrix3<T> &m2)
	{
		CMatrix3 m;

		m[0][0] = matrix[0][0]*m2[0][0] + matrix[0][1]*m2[1][0] + matrix[0][2]*m2[2][0];
		m[0][1] = matrix[0][0]*m2[0][1] + matrix[0][1]*m2[1][1] + matrix[0][2]*m2[2][1];
		m[0][2] = matrix[0][0]*m2[0][2] + matrix[0][1]*m2[1][2] + matrix[0][2]*m2[2][2];

		m[1][0] = matrix[1][0]*m2[0][0] + matrix[1][1]*m2[1][0] + matrix[1][2]*m2[2][0];
		m[1][1] = matrix[1][0]*m2[0][1] + matrix[1][1]*m2[1][1] + matrix[1][2]*m2[2][1];
		m[1][2] = matrix[1][0]*m2[0][2] + matrix[1][1]*m2[1][2] + matrix[1][2]*m2[2][2];

		m[2][0] = matrix[2][0]*m2[0][0] + matrix[2][1]*m2[1][0] + matrix[2][2]*m2[2][0];
		m[2][1] = matrix[2][0]*m2[0][1] + matrix[2][1]*m2[1][1] + matrix[2][2]*m2[2][1];
		m[2][2] = matrix[2][0]*m2[0][2] + matrix[2][1]*m2[1][2] + matrix[2][2]*m2[2][2];

		return m;
	}

private:
	///< switch 3 colums and vector elements
	inline void solve_switch_row(CMatrix3<T> &m, CVector<3,T> &b, int ra, int rb)
	{
		T a;
		// switch 1st and 3nd row
		a = m[ra][0];	m[ra][0] = m[rb][0];	m[rb][0] = a;
		a = m[ra][1];	m[ra][0] = m[rb][1];	m[rb][1] = a;
		a = m[ra][2];	m[ra][0] = m[rb][2];	m[rb][2] = a;
		a = b[ra];	b[ra] = b[rb];		b[rb] = a;
	}

	///< switch just 2 columns and vector elements
	inline void solve_switch_row2(CMatrix3<T> &m, CVector<3,T> &b, int ra, int rb)
	{
		T a;
		// switch 1st and 3nd row
		a = m[ra][1];	m[ra][0] = m[rb][1];	m[rb][1] = a;
		a = m[ra][2];	m[ra][0] = m[rb][2];	m[rb][2] = a;
		a = b[ra];	b[ra] = b[rb];		b[rb] = a;
	}

public:
	/**
	 * create a rotation matrix
	 *
	 * taken from opengl manpage:
	 *
	 *	x^2(1-c)+c	xy(1-c)-zs	xz(1-c)+ys
	 *	yx(1-c)+zs	y^2(1-c)+c	yz(1-c)-xs
	 *	xz(1-c)-ys	yz(1-c)+xs	z^2(1-c)+c
	 */
	void genRotation(	T angle, 		// rotation angle in radians
						const CVector<3,T> axis	// axis of rotation
			)
	{
		T c = cos(angle);
		T s = sin(angle);
		T cm = T(1)-c;

		CVector<3,T> naxis = axis;
		naxis.normalize();

		matrix[0][0] = naxis[0]*naxis[0]*cm + c;
		matrix[0][1] = naxis[0]*naxis[1]*cm - naxis[2]*s;
		matrix[0][2] = naxis[0]*naxis[2]*cm + naxis[1]*s;

		matrix[1][0] = naxis[1]*naxis[0]*cm + naxis[2]*s;
		matrix[1][1] = naxis[1]*naxis[1]*cm + c;
		matrix[1][2] = naxis[1]*naxis[2]*cm - naxis[0]*s;

		matrix[2][0] = naxis[2]*naxis[0]*cm - naxis[1]*s;
		matrix[2][1] = naxis[2]*naxis[1]*cm + naxis[0]*s;
		matrix[2][2] = naxis[2]*naxis[2]*cm + c;
	}

	/**
	 * linear array access to matrix components
	 */
	inline CVector<3,T>& operator[](const int i)
	{
		return matrix[i];
	}
};



template <class T>
inline
std::ostream&
operator<<(std::ostream &co, CMatrix3<T> &m)
{
	return co	<< "[" << m[0][0] << ", " << m[0][1] << ", " << m[0][2] << "]" << std::endl
				<< "[" << m[1][0] << ", " << m[1][1] << ", " << m[1][2] << "]" << std::endl
				<< "[" << m[2][0] << ", " << m[2][1] << ", " << m[2][2] << "]" << std::endl	;
}



#endif
