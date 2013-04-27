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
	#error "dont include CMatrix4.hpp directly - use CMatrix.hpp instead!"
#endif

#ifndef TMATRIX4_HH
#define TMATRIX4_HH

/**
 * CHANGELOG:
 *
 * 2008-10-28: martin schreiber
 * 	created
 *
 * 2009-12-15: martin schreiber
 * 	avoiding initialization with zero values with default constructor
 * 	values are set to 0, if some integer parameter is given to the constructor
 *
 * 2009-12-22: martin schreiber
 *  modifications for matrix inversion / transpose
 */

#include <limits>
#include <iostream>
#include "CMath.hpp"
#include "CVector.hpp"

/**
 * \brief	4x4 matrix handler
 */
template <typename T>
class CMatrix4
{
public:
	T matrix[4][4];		///< matrix array

	/**
	 * constructor: load identity matrix
	 */
	inline CMatrix4()
	{
	}

	/**
	 * initialize only with zero, if parameter is set
	 */
	inline CMatrix4(int i)
	{
		setZero();
	}

	/**
	 * initialize matrix with given scalar values
	 */
	inline CMatrix4(	T m00, T m01, T m02, T m03,
						T m10, T m11, T m12, T m13,
						T m20, T m21, T m22, T m23,
						T m30, T m31, T m32, T m33
			)
	{
		matrix[0][0] = m00;
		matrix[0][1] = m01;
		matrix[0][2] = m02;
		matrix[0][3] = m03;

		matrix[1][0] = m10;
		matrix[1][1] = m11;
		matrix[1][2] = m12;
		matrix[1][3] = m13;

		matrix[2][0] = m20;
		matrix[2][1] = m21;
		matrix[2][2] = m22;
		matrix[2][3] = m23;

		matrix[3][0] = m30;
		matrix[3][1] = m31;
		matrix[3][2] = m32;
		matrix[3][3] = m33;
	}

	/**
	 * set values to 0
	 */
	inline void setZero()
	{
		matrix[0][0] = (T)0;
		matrix[0][1] = (T)0;
		matrix[0][2] = (T)0;
		matrix[0][3] = (T)0;
		matrix[1][0] = (T)0;
		matrix[1][1] = (T)0;
		matrix[1][2] = (T)0;
		matrix[1][3] = (T)0;
		matrix[2][0] = (T)0;
		matrix[2][1] = (T)0;
		matrix[2][2] = (T)0;
		matrix[2][3] = (T)0;
		matrix[3][0] = (T)0;
		matrix[3][1] = (T)0;
		matrix[3][2] = (T)0;
		matrix[3][3] = (T)0;
	}

	/**
	 * linear array access to matrix components
	 */
	inline T* operator[](const int i)
	{
		return &(matrix[i][0]);
	}

	/**
	 * copy data from matrix m
	 * \param m	source matrix
	 */
	inline CMatrix4<T>& operator=(const CMatrix4<T> &m)
	{
		matrix[0][0] = m.matrix[0][0];
		matrix[0][1] = m.matrix[0][1];
		matrix[0][2] = m.matrix[0][2];
		matrix[0][3] = m.matrix[0][3];

		matrix[1][0] = m.matrix[1][0];
		matrix[1][1] = m.matrix[1][1];
		matrix[1][2] = m.matrix[1][2];
		matrix[1][3] = m.matrix[1][3];

		matrix[2][0] = m.matrix[2][0];
		matrix[2][1] = m.matrix[2][1];
		matrix[2][2] = m.matrix[2][2];
		matrix[2][3] = m.matrix[2][3];

		matrix[3][0] = m.matrix[3][0];
		matrix[3][1] = m.matrix[3][1];
		matrix[3][2] = m.matrix[3][2];
		matrix[3][3] = m.matrix[3][3];

		return *this;
	};

	/**
	 * load identity matrix
	 */
	inline void loadIdentity()
	{
		matrix[0][0] = (T)1;
		matrix[0][1] = (T)0;
		matrix[0][2] = (T)0;
		matrix[0][3] = (T)0;

		matrix[1][0] = (T)0;
		matrix[1][1] = (T)1;
		matrix[1][2] = (T)0;
		matrix[1][3] = (T)0;

		matrix[2][0] = (T)0;
		matrix[2][1] = (T)0;
		matrix[2][2] = (T)1;
		matrix[2][3] = (T)0;

		matrix[3][0] = (T)0;
		matrix[3][1] = (T)0;
		matrix[3][2] = (T)0;
		matrix[3][3] = (T)1;
	}

	/**
	 * matrix vector product
	 * use 1 as 4th component of input vector v
	 * avoid using parameters by references because this cannot be handled by
	 * subsequenced  orders like a+(o*i)
	 */
	inline CVector<3,T>	operator*(CVector<3,T> v)
	{
		return CVector<3,T>(	matrix[0][0]*v[0] + matrix[0][1]*v[1] + matrix[0][2]*v[2] + matrix[0][3],
								matrix[1][0]*v[0] + matrix[1][1]*v[1] + matrix[1][2]*v[2] + matrix[1][3],
								matrix[2][0]*v[0] + matrix[2][1]*v[1] + matrix[2][2]*v[2] + matrix[2][3]
				);
	}

	/**
	 * matrix vector product
	 */
	inline CVector<4,T>	operator*(CVector<4,T> v)
	{
		return CVector<4,T>(	matrix[0][0]*v[0] + matrix[0][1]*v[1] + matrix[0][2]*v[2] + matrix[0][3]*v[3],
								matrix[1][0]*v[0] + matrix[1][1]*v[1] + matrix[1][2]*v[2] + matrix[1][3]*v[3],
								matrix[2][0]*v[0] + matrix[2][1]*v[1] + matrix[2][2]*v[2] + matrix[2][3]*v[3],
								matrix[3][0]*v[0] + matrix[3][1]*v[1] + matrix[3][2]*v[2] + matrix[3][3]*v[3]
			);
	}


	/**
	 * matrix - scalar product
	 */
	inline CMatrix4<T>&	operator*=(T v)
	{
		matrix[0][0] *= v;
		matrix[0][1] *= v;
		matrix[0][2] *= v;
		matrix[0][3] *= v;
		matrix[1][0] *= v;
		matrix[1][1] *= v;
		matrix[1][2] *= v;
		matrix[1][3] *= v;
		matrix[2][0] *= v;
		matrix[2][1] *= v;
		matrix[2][2] *= v;
		matrix[2][3] *= v;
		matrix[3][0] *= v;
		matrix[3][1] *= v;
		matrix[3][2] *= v;
		matrix[3][3] *= v;

		return *this;
	}

	/**
	 * matrix - scalar division
	 */
	inline CMatrix4<T>&	operator/=(const T v)
	{
		matrix[0][0] /= v;
		matrix[0][1] /= v;
		matrix[0][2] /= v;
		matrix[0][3] /= v;
		matrix[1][0] /= v;
		matrix[1][1] /= v;
		matrix[1][2] /= v;
		matrix[1][3] /= v;
		matrix[2][0] /= v;
		matrix[2][1] /= v;
		matrix[2][2] /= v;
		matrix[2][3] /= v;
		matrix[3][0] /= v;
		matrix[3][1] /= v;
		matrix[3][2] /= v;
		matrix[3][3] /= v;

		return *this;
	}

	/**
	 * multiply this matrix with m2 and return result
	 * \return this->matrix * m2
	 */
	inline CMatrix4<T>	operator*(CMatrix4<T> &m2)
	{
		CMatrix4 m;

		m[0][0] = matrix[0][0]*m2[0][0] + matrix[0][1]*m2[1][0] + matrix[0][2]*m2[2][0] + matrix[0][3]*m2[3][0];
		m[0][1] = matrix[0][0]*m2[0][1] + matrix[0][1]*m2[1][1] + matrix[0][2]*m2[2][1] + matrix[0][3]*m2[3][1];
		m[0][2] = matrix[0][0]*m2[0][2] + matrix[0][1]*m2[1][2] + matrix[0][2]*m2[2][2] + matrix[0][3]*m2[3][2];
		m[0][3] = matrix[0][0]*m2[0][3] + matrix[0][1]*m2[1][3] + matrix[0][2]*m2[2][3] + matrix[0][3]*m2[3][3];

		m[1][0] = matrix[1][0]*m2[0][0] + matrix[1][1]*m2[1][0] + matrix[1][2]*m2[2][0] + matrix[1][3]*m2[3][0];
		m[1][1] = matrix[1][0]*m2[0][1] + matrix[1][1]*m2[1][1] + matrix[1][2]*m2[2][1] + matrix[1][3]*m2[3][1];
		m[1][2] = matrix[1][0]*m2[0][2] + matrix[1][1]*m2[1][2] + matrix[1][2]*m2[2][2] + matrix[1][3]*m2[3][2];
		m[1][3] = matrix[1][0]*m2[0][3] + matrix[1][1]*m2[1][3] + matrix[1][2]*m2[2][3] + matrix[1][3]*m2[3][3];

		m[2][0] = matrix[2][0]*m2[0][0] + matrix[2][1]*m2[1][0] + matrix[2][2]*m2[2][0] + matrix[2][3]*m2[3][0];
		m[2][1] = matrix[2][0]*m2[0][1] + matrix[2][1]*m2[1][1] + matrix[2][2]*m2[2][1] + matrix[2][3]*m2[3][1];
		m[2][2] = matrix[2][0]*m2[0][2] + matrix[2][1]*m2[1][2] + matrix[2][2]*m2[2][2] + matrix[2][3]*m2[3][2];
		m[2][3] = matrix[2][0]*m2[0][3] + matrix[2][1]*m2[1][3] + matrix[2][2]*m2[2][3] + matrix[2][3]*m2[3][3];

		m[3][0] = matrix[3][0]*m2[0][0] + matrix[3][1]*m2[1][0] + matrix[3][2]*m2[2][0] + matrix[3][3]*m2[3][0];
		m[3][1] = matrix[3][0]*m2[0][1] + matrix[3][1]*m2[1][1] + matrix[3][2]*m2[2][1] + matrix[3][3]*m2[3][1];
		m[3][2] = matrix[3][0]*m2[0][2] + matrix[3][1]*m2[1][2] + matrix[3][2]*m2[2][2] + matrix[3][3]*m2[3][2];
		m[3][3] = matrix[3][0]*m2[0][3] + matrix[3][1]*m2[1][3] + matrix[3][2]*m2[2][3] + matrix[3][3]*m2[3][3];

		return m;
	}

	/**
	 * multiply this matrix with m2 and store result to dst without modification of m2 and this->matrix
	 * \return this->matrix * m2
	 */
	inline void multiplyAndStore(CMatrix3<T> &m2, CMatrix4<T> &dst)
	{
		dst[0][0] = matrix[0][0]*m2[0][0] + matrix[0][1]*m2[1][0] + matrix[0][2]*m2[2][0] + matrix[0][3]*m2[3][0];
		dst[0][1] = matrix[0][0]*m2[0][1] + matrix[0][1]*m2[1][1] + matrix[0][2]*m2[2][1] + matrix[0][3]*m2[3][1];
		dst[0][2] = matrix[0][0]*m2[0][2] + matrix[0][1]*m2[1][2] + matrix[0][2]*m2[2][2] + matrix[0][3]*m2[3][2];
		dst[0][3] = matrix[0][0]*m2[0][3] + matrix[0][1]*m2[1][3] + matrix[0][2]*m2[2][3] + matrix[0][3]*m2[3][3];

		dst[1][0] = matrix[1][0]*m2[0][0] + matrix[1][1]*m2[1][0] + matrix[1][2]*m2[2][0] + matrix[1][3]*m2[3][0];
		dst[1][1] = matrix[1][0]*m2[0][1] + matrix[1][1]*m2[1][1] + matrix[1][2]*m2[2][1] + matrix[1][3]*m2[3][1];
		dst[1][2] = matrix[1][0]*m2[0][2] + matrix[1][1]*m2[1][2] + matrix[1][2]*m2[2][2] + matrix[1][3]*m2[3][2];
		dst[1][3] = matrix[1][0]*m2[0][3] + matrix[1][1]*m2[1][3] + matrix[1][2]*m2[2][3] + matrix[1][3]*m2[3][3];

		dst[2][0] = matrix[2][0]*m2[0][0] + matrix[2][1]*m2[1][0] + matrix[2][2]*m2[2][0] + matrix[2][3]*m2[3][0];
		dst[2][1] = matrix[2][0]*m2[0][1] + matrix[2][1]*m2[1][1] + matrix[2][2]*m2[2][1] + matrix[2][3]*m2[3][1];
		dst[2][2] = matrix[2][0]*m2[0][2] + matrix[2][1]*m2[1][2] + matrix[2][2]*m2[2][2] + matrix[2][3]*m2[3][2];
		dst[2][3] = matrix[2][0]*m2[0][3] + matrix[2][1]*m2[1][3] + matrix[2][2]*m2[2][3] + matrix[2][3]*m2[3][3];

		dst[3][0] = matrix[3][0]*m2[0][0] + matrix[3][1]*m2[1][0] + matrix[3][2]*m2[2][0] + matrix[3][3]*m2[3][0];
		dst[3][1] = matrix[3][0]*m2[0][1] + matrix[3][1]*m2[1][1] + matrix[3][2]*m2[2][1] + matrix[3][3]*m2[3][1];
		dst[3][2] = matrix[3][0]*m2[0][2] + matrix[3][1]*m2[1][2] + matrix[3][2]*m2[2][2] + matrix[3][3]*m2[3][2];
		dst[3][3] = matrix[3][0]*m2[0][3] + matrix[3][1]*m2[1][3] + matrix[3][2]*m2[2][3] + matrix[3][3]*m2[3][3];
	}


	/**
	 * compute determinant
	 * taken from crystalspace engine
	 * http://www.crystalspace3d.org/docs/online/api/matrix4_8h-source.html
	 */
	inline float getDeterminant() const
	{
#define m matrix
		return	  m[0][3] * m[1][2] * m[2][1] * m[3][0] - m[0][2] * m[1][3] * m[2][1] * m[3][0] - m[0][3] * m[1][1] * m[2][2] * m[3][0] + m[0][1] * m[1][3] * m[2][2] * m[3][0]
				+ m[0][2] * m[1][1] * m[2][3] * m[3][0] - m[0][1] * m[1][2] * m[2][3] * m[3][0] - m[0][3] * m[1][2] * m[2][0] * m[3][1] + m[0][2] * m[1][3] * m[2][0] * m[3][1]
				+ m[0][3] * m[1][0] * m[2][2] * m[3][1] - m[0][0] * m[1][3] * m[2][2] * m[3][1] - m[0][2] * m[1][0] * m[2][3] * m[3][1] + m[0][0] * m[1][2] * m[2][3] * m[3][1]
				+ m[0][3] * m[1][1] * m[2][0] * m[3][2] - m[0][1] * m[1][3] * m[2][0] * m[3][2] - m[0][3] * m[1][0] * m[2][1] * m[3][2] + m[0][0] * m[1][3] * m[2][1] * m[3][2]
				+ m[0][1] * m[1][0] * m[2][3] * m[3][2] - m[0][0] * m[1][1] * m[2][3] * m[3][2] - m[0][2] * m[1][1] * m[2][0] * m[3][3] + m[0][1] * m[1][2] * m[2][0] * m[3][3]
				+ m[0][2] * m[1][0] * m[2][1] * m[3][3] - m[0][0] * m[1][2] * m[2][1] * m[3][3] - m[0][1] * m[1][0] * m[2][2] * m[3][3] + m[0][0] * m[1][1] * m[2][2] * m[3][3];
#undef m
	}

	/**
	 * return inverse of matrix
	 * taken from crystalspace engine
	 * http://www.crystalspace3d.org/docs/online/api/matrix4_8h-source.html
	 *
	 * the inversion makes use of cramers rule
	 */
	CMatrix4<T> getInverse() const
	{
	#define m matrix
		CMatrix4<T> sol (
			m[1][2]*m[2][3]*m[3][1] - m[1][3]*m[2][2]*m[3][1] + m[1][3]*m[2][1]*m[3][2] - m[1][1]*m[2][3]*m[3][2] - m[1][2]*m[2][1]*m[3][3] + m[1][1]*m[2][2]*m[3][3],
			m[0][3]*m[2][2]*m[3][1] - m[0][2]*m[2][3]*m[3][1] - m[0][3]*m[2][1]*m[3][2] + m[0][1]*m[2][3]*m[3][2] + m[0][2]*m[2][1]*m[3][3] - m[0][1]*m[2][2]*m[3][3],
			m[0][2]*m[1][3]*m[3][1] - m[0][3]*m[1][2]*m[3][1] + m[0][3]*m[1][1]*m[3][2] - m[0][1]*m[1][3]*m[3][2] - m[0][2]*m[1][1]*m[3][3] + m[0][1]*m[1][2]*m[3][3],
			m[0][3]*m[1][2]*m[2][1] - m[0][2]*m[1][3]*m[2][1] - m[0][3]*m[1][1]*m[2][2] + m[0][1]*m[1][3]*m[2][2] + m[0][2]*m[1][1]*m[2][3] - m[0][1]*m[1][2]*m[2][3],

			m[1][3]*m[2][2]*m[3][0] - m[1][2]*m[2][3]*m[3][0] - m[1][3]*m[2][0]*m[3][2] + m[1][0]*m[2][3]*m[3][2] + m[1][2]*m[2][0]*m[3][3] - m[1][0]*m[2][2]*m[3][3],
			m[0][2]*m[2][3]*m[3][0] - m[0][3]*m[2][2]*m[3][0] + m[0][3]*m[2][0]*m[3][2] - m[0][0]*m[2][3]*m[3][2] - m[0][2]*m[2][0]*m[3][3] + m[0][0]*m[2][2]*m[3][3],
			m[0][3]*m[1][2]*m[3][0] - m[0][2]*m[1][3]*m[3][0] - m[0][3]*m[1][0]*m[3][2] + m[0][0]*m[1][3]*m[3][2] + m[0][2]*m[1][0]*m[3][3] - m[0][0]*m[1][2]*m[3][3],
			m[0][2]*m[1][3]*m[2][0] - m[0][3]*m[1][2]*m[2][0] + m[0][3]*m[1][0]*m[2][2] - m[0][0]*m[1][3]*m[2][2] - m[0][2]*m[1][0]*m[2][3] + m[0][0]*m[1][2]*m[2][3],

			m[1][1]*m[2][3]*m[3][0] - m[1][3]*m[2][1]*m[3][0] + m[1][3]*m[2][0]*m[3][1] - m[1][0]*m[2][3]*m[3][1] - m[1][1]*m[2][0]*m[3][3] + m[1][0]*m[2][1]*m[3][3],
			m[0][3]*m[2][1]*m[3][0] - m[0][1]*m[2][3]*m[3][0] - m[0][3]*m[2][0]*m[3][1] + m[0][0]*m[2][3]*m[3][1] + m[0][1]*m[2][0]*m[3][3] - m[0][0]*m[2][1]*m[3][3],
			m[0][1]*m[1][3]*m[3][0] - m[0][3]*m[1][1]*m[3][0] + m[0][3]*m[1][0]*m[3][1] - m[0][0]*m[1][3]*m[3][1] - m[0][1]*m[1][0]*m[3][3] + m[0][0]*m[1][1]*m[3][3],
			m[0][3]*m[1][1]*m[2][0] - m[0][1]*m[1][3]*m[2][0] - m[0][3]*m[1][0]*m[2][1] + m[0][0]*m[1][3]*m[2][1] + m[0][1]*m[1][0]*m[2][3] - m[0][0]*m[1][1]*m[2][3],

			m[1][2]*m[2][1]*m[3][0] - m[1][1]*m[2][2]*m[3][0] - m[1][2]*m[2][0]*m[3][1] + m[1][0]*m[2][2]*m[3][1] + m[1][1]*m[2][0]*m[3][2] - m[1][0]*m[2][1]*m[3][2],
			m[0][1]*m[2][2]*m[3][0] - m[0][2]*m[2][1]*m[3][0] + m[0][2]*m[2][0]*m[3][1] - m[0][0]*m[2][2]*m[3][1] - m[0][1]*m[2][0]*m[3][2] + m[0][0]*m[2][1]*m[3][2],
			m[0][2]*m[1][1]*m[3][0] - m[0][1]*m[1][2]*m[3][0] - m[0][2]*m[1][0]*m[3][1] + m[0][0]*m[1][2]*m[3][1] + m[0][1]*m[1][0]*m[3][2] - m[0][0]*m[1][1]*m[3][2],
			m[0][1]*m[1][2]*m[2][0] - m[0][2]*m[1][1]*m[2][0] + m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2]
		);

		sol /= getDeterminant();
		return sol;
	#undef m
	}

	/**
	 * return transpose of matrix
	 */
	CMatrix4<T> getTranspose() const
	{
		return CMatrix4<T> (
			matrix[0][0], matrix[1][0], matrix[2][0], matrix[3][0],
			matrix[0][1], matrix[1][1], matrix[2][1], matrix[3][1],
			matrix[0][2], matrix[1][2], matrix[2][2], matrix[3][2],
			matrix[0][3], matrix[1][3], matrix[2][3], matrix[3][3]
		);
	}

	/**
	 * return the inverse transpose of the matrix
	 */
	CMatrix4<T> getInverseTranspose() const
	{
		CMatrix4<T> m = getInverse();
		return m.getTranspose();
	}

    /**
     * generate a rotation matrix
     *
     * taken from opengl manpage:
     *
     *      x^2(1-c)+c      xy(1-c)-zs      xz(1-c)+ys
     *      yx(1-c)+zs      y^2(1-c)+c      yz(1-c)-xs
     *      xz(1-c)-ys      yz(1-c)+xs      z^2(1-c)+c
     */
	void genRotation(   T angle,			// rotation angle in radians
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
		matrix[0][3] = 0;

		matrix[1][0] = naxis[1]*naxis[0]*cm + naxis[2]*s;
		matrix[1][1] = naxis[1]*naxis[1]*cm + c;
		matrix[1][2] = naxis[1]*naxis[2]*cm - naxis[0]*s;
		matrix[1][3] = 0;

		matrix[2][0] = naxis[2]*naxis[0]*cm - naxis[1]*s;
		matrix[2][1] = naxis[2]*naxis[1]*cm + naxis[0]*s;
		matrix[2][2] = naxis[2]*naxis[2]*cm + c;
		matrix[2][3] = 0;

		matrix[3][0] = 0;
		matrix[3][1] = 0;
		matrix[3][2] = 0;
		matrix[3][3] = 1;
    }
};

typedef CMatrix4<float>	mat4f;


template <class T>
inline
std::ostream&
operator<<(std::ostream &co, CMatrix4<T> &m)
{
	return co	<< "[" << m[0][0] << ", " << m[0][1] << ", " << m[0][2] << ", " << m[0][3] << "]" << std::endl
			<< "[" << m[1][0] << ", " << m[1][1] << ", " << m[1][2] << ", " << m[1][3] << "]" << std::endl
			<< "[" << m[2][0] << ", " << m[2][1] << ", " << m[2][2] << ", " << m[2][3] << "]" << std::endl
			<< "[" << m[3][0] << ", " << m[3][1] << ", " << m[3][2] << ", " << m[3][3] << "]" << std::endl	;
}



#endif
