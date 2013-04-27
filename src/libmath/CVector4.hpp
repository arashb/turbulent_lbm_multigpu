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


#ifndef __CVECTOR_HH
	#error "dont include CVector4.hpp directly!"
#endif

#ifndef __CVECTOR4_HH
#define __CVECTOR4_HH
/**
 * CHANGELOG:
 *
 * 2008-02-27: martin schreiber
 * 	dist2();
 *
 * 2010-01-18: martin schreiber
 *  updates for doxygen
 */

#include <iostream>
#include "CMath.hpp"

/**
 * \brief	4D Vector handler
 */
template <typename T>
class CVector<4,T>
{
public:
	T data[4];		///< vector data

	/*******************
	 * CONSTRUCTURS
	 *******************/
	inline CVector()
	{
		setZero();
	}

	/**
	 * initialize vector with (x0, x1, x2, x3)
	 */
	inline CVector(const T x0, const T x1, const T x2, const T x3)
	{
		data[0] = x0;
		data[1] = x1;
		data[2] = x2;
		data[3] = x3;
	}

	/**
	 * initialize all vector components with the scalar value 'x'
	 */
	inline CVector(const T x)
	{
		data[0] = x;
		data[1] = x;
		data[2] = x;
		data[3] = x;
	}

	/**
	 * initialize integer vector components with the array 'v' converting to the type (T)
	 */
	inline CVector(const T v[4])
	{
		data[0] = v[0];
		data[1] = v[1];
		data[2] = v[2];
		data[3] = v[3];
	}


	/**
	 * set all components of vector to 0
	 */
	inline void setZero()
	{
		data[0] = T(0);
		data[1] = T(0);
		data[2] = T(0);	
		data[3] = T(0);	
	}

	/**
	 * return length of vector
	 */
	inline T length()
	{
		return CMath<T>::sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2] + data[3]*data[3]);
	}

	/**
	 * return (length*length) of vector
	 */
	inline T length2()
	{
		return data[0]*data[0] + data[1]*data[1] + data[2]*data[2] + data[3]*data[3];
	}



	/*******************
	 * OPERATORS
	 *******************/
	/// assign values of a[4] to this vector and return reference to this vector
	inline CVector<4,T>&	operator=(const T a[4])	{	data[0] = a[0]; data[1] = a[1]; data[2] = a[2];	data[3] = a[3]; return *this;	};
	/// assign values of vector a to this vector and return reference to this vector
	inline CVector<4,T>&	operator=(CVector<4,T> const & a)	{	data[0] = a.data[0]; data[1] = a.data[1]; data[2] = a.data[2];	data[3] = a.data[3]; return *this;	};
	/// assign values of vector a to this vector, set the 4th component to 1 and return reference to this vector
	inline CVector<4,T>&	operator=(CVector<3,T> const & a)	{	data[0] = a.data[0]; data[1] = a.data[1]; data[2] = a.data[2];	data[3] = 1; return *this;	};


// T
	/// return new vector (this+a)
	inline CVector<4,T>	operator+(const T a)	{	return CVector<4,T>(data[0]+a, data[1]+a, data[2]+a, data[3]+a);	}
	/// return new vector (this-a)
	inline CVector<4,T>	operator-(const T a)	{	return CVector<4,T>(data[0]-a, data[1]-a, data[2]-a, data[3]-a);	}
	/// return new vector with component wise (this*a)
	inline CVector<4,T>	operator*(const T a)	{	return CVector<4,T>(data[0]*a, data[1]*a, data[2]*a, data[3]*a);	}
	/// return new vector with component wise (this/a)
	inline CVector<4,T>	operator/(const T a)	{	return CVector<4,T>(data[0]/a, data[1]/a, data[2]/a, data[3]/a);	}
	/// add a to this vector and return reference to this vector
	inline CVector<4,T>& operator+=(const T a)	{	data[0] += a; data[1] += a; data[2] += a; data[3] += a;	return *this;	}
	/// subtract a from this vector and return reference to this vector
	inline CVector<4,T>& operator-=(const T a)	{	data[0] -= a; data[1] -= a; data[2] -= a; data[3] -= a;	return *this;	}
	/// multiply each component of this vector with scalar a and return reference to this vector
	inline CVector<4,T>& operator*=(const T a)	{	data[0] *= a; data[1] *= a; data[2] *= a; data[3] *= a;	return *this;	}
	/// divide each component of this vector by scalar a and return reference to this vector
	inline CVector<4,T>& operator/=(const T a)	{	data[0] /= a; data[1] /= a; data[2] /= a; data[3] /= a;	return *this;	}

	/// return new vector with sum of this vector and v
	inline CVector<4,T>	operator+(const CVector<4,T> &v)	{	return CVector<4,T>(data[0]+v.data[0], data[1]+v.data[1], data[2]+v.data[2], data[3]+v.data[3]);	}
	/// return new vector with subtraction of vector v from this vector
	inline CVector<4,T>	operator-(const CVector<4,T> &v)	{	return CVector<4,T>(data[0]-v.data[0], data[1]-v.data[1], data[2]-v.data[2], data[3]-v.data[3]);	}
	/// return new vector with values of this vector multiplied component wise with vector v
	inline CVector<4,T>	operator*(const CVector<4,T> &v)	{	return CVector<4,T>(data[0]*v.data[0], data[1]*v.data[1], data[2]*v.data[2], data[3]*v.data[3]);	}
	/// return new vector with values of this vector divided component wise by components of vector v
	inline CVector<4,T>	operator/(const CVector<4,T> &v)	{	return CVector<4,T>(data[0]/v.data[0], data[1]/v.data[1], data[2]/v.data[2], data[3]/v.data[3]);	}

	/// return this vector after adding v
	inline CVector<4,T>&	operator+=(const CVector<4,T> &v)	{	data[0] += v.data[0]; data[1] += v.data[1]; data[2] += v.data[2]; data[3] += v.data[3];	return *this;	}
	/// return this vector after subtracting v
	inline CVector<4,T>&	operator-=(const CVector<4,T> &v)	{	data[0] -= v.data[0]; data[1] -= v.data[1]; data[2] -= v.data[2]; data[3] -= v.data[3];	return *this;	}

	/// return true, if each component of the vector is equal to the corresponding component of vector v
	inline bool	operator==(const CVector<4,T> &v)	{	return bool(data[0] == v.data[0] && data[1] == v.data[1] && data[2] == v.data[2] && data[3] == v.data[3]);	}
	/// return true, if at lease component of the vector is not equal to the corresponding component of vector v
	inline bool	operator!=(const CVector<4,T> &v)	{	return bool(data[0] != v.data[0] || data[1] != v.data[1] || data[2] != v.data[2] || data[3] != v.data[3]);	}


	/**
	 * access element i
	 */
	inline T& operator[](const int i)
	{
		return data[i];
	}

	/**
	 * \brief	compare set for sort operation
	 */
	struct compareSet
	{
		/**
		 * compare set operator
		 */
		inline bool operator()(CVector<4,T> *v1, CVector<4,T> *v2)
		{
			if ((*v1)[0] != (*v2)[0])
				return (*v1)[0] < (*v2)[0];
			if ((*v1)[1] != (*v2)[1])
				return (*v1)[1] < (*v2)[1];
			if ((*v1)[2] != (*v2)[2])
				return (*v1)[2] < (*v2)[2];
			if ((*v1)[3] != (*v2)[3])
				return (*v1)[3] < (*v2)[3];
			return false;
		}

		/**
		 * compare set operator
		 */
		inline bool operator()(const CVector<4,T> &v1, const CVector<4,T> &v2)
		{
			if (v1.data[0] != v2.data[0])
				return v1.data[0] < v2.data[0];
			if (v1.data[1] != v2.data[1])
				return v1.data[1] < v2.data[1];
			if (v1.data[2] != v2.data[2])
				return v1.data[2] < v2.data[2];
			if (v1.data[3] != v2.data[3])
				return v1.data[3] < v2.data[3];
			return false;
		}
	};
};


typedef CVector<4,double>	vec4d;
typedef CVector<4,float>	vec4f;
typedef CVector<4,int>	vec4i;


template <class T>
inline
std::ostream&
operator<<(std::ostream &co, CVector<4,T> &v)
{
	return co << "[" << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << "]";
}


#endif
