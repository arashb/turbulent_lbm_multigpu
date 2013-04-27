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
 * a small math lib with some useful functions and a c++ abstraction layer
 */
#ifndef CMATH_HPP
#define CMATH_HPP

#include <iostream>
#include <limits>

extern "C" {
	#include <math.h>
	#include <stdlib.h>
}

#ifndef fabsf
#define	m_fabsf(x)	(::fabs(x))
#else
#define	m_fabsf(x)	(::fabsf(x))
#endif

#ifndef floorf
#define	m_floorf(x)	(::floor(x))
#else
#define	m_floorf(x)	(::floorf(x))
#endif

#ifndef ceilf
#define	m_ceilf(x)	(::ceil(x))
#else
#define	m_ceilf(x)	(::ceilf(x))
#endif

#ifndef sqrtf
#define	m_sqrtf(x)	(::sqrt(x))
#else
#define	m_sqrtf(x)	(::sqrtf(x))
#endif

#ifndef sqrtl
#define	m_sqrtl(x)	((unsigned long)::sqrt((double)x))
#else
#define	m_sqrtl(x)	((unsigned long)::sqrtl((unsigned long)x))
#endif


#ifndef roundf
#define m_roundf(x)	((float)(int)(x+0.5f))
#else
#define m_roundf(x)	(::roundf(x))
#endif

#ifndef round
#define m_round(x)	((double)(int)(x+0.5))
#else
#define m_round(x)	(round(x))
#endif


/**
 * \brief	math handler to use same function names for different types
 */
template <typename T>
class CMath
{
public:

	inline CMath()	{
	}

	/// assign value
	inline CMath(T a)	{	value = a;	}

// CMath numbers
	/// add operator
	inline CMath	operator+(const CMath &a)	{	return CMath(value + a.value);	}
	/// sub operator
	inline CMath	operator-(const CMath &a)	{	return CMath(value - a.value);	}

	/// multiply operator
	inline CMath	operator*(const CMath a)	{	return CMath(value * a.value);	}

	/// div operator
	inline CMath	operator/(const CMath &a)	{	return CMath(value / a.value);	}

	/// add and return this operator
	inline CMath&	operator+=(const CMath &a)	{	value += a.value;		return *this;	}
	/// sub and return this operator
	inline CMath&	operator-=(const CMath &a)	{	value -= a.value;		return *this;	}
	/// multiply and return this operator
	inline CMath&	operator*=(const CMath &a)	{	value *= a.value;		return *this;	}
	/// div and return this operator
	inline CMath&	operator/=(const CMath &a)	{	value /= a.value;		return *this;	}

// T numbers
	/// add operator
	inline CMath	operator+(const T a)		{	return CMath(value + a);	}
	/// sub operator
	inline CMath	operator-(const T a)		{	return CMath(value - a);	}
	/// multiply operator
	inline CMath	operator*(const T a)		{	return CMath(value * a);	}
	/// div operator
	inline CMath	operator/(const T a)		{	return CMath(value / a);	}
	/// add and return this operator
	inline CMath&	operator+=(const T a)		{	value += a;		return *this;	}
	/// sub and return this operator
	inline CMath&	operator-=(const T a)		{	value -= a;		return *this;	}
	/// mul and return this operator
	inline CMath&	operator*=(const T a)		{	value *= a;		return *this;	}
	/// div and return this operator
	inline CMath&	operator/=(const T a)		{	value /= a;		return *this;	}

	/// return PI
	static inline T PI();

	/// return absolute value
	static inline T abs(T a);
	/// return power of value
	static inline T pow(T base, T exp);
	/// return floor value
	static inline T floor(T a);
	/// return ceil value
	static inline T ceil(T a);

	/// return ceiled value in binary system
	static inline T ceil2(T a);
	/// return sqrt value
	static inline T sqrt(T a);
	/// return rounded value
	static inline T round(T a);

	/// return digits of a in binary system
	static inline T digits2(T a);

	/// return exp value
	static inline T exp(T a);

	/// return sin value
	static inline T sin(T a);
	/// return cos value
	static inline T cos(T a);
	/// return tan value
	static inline T tan(T a);

	/// return max of both values
	static inline T max(T a, T b)		{	return (a < b ? b : a);	}
	/// return min of both values
	static inline T min(T a, T b)		{	return (a > b ? b : a);	}

	/// return value of type T for string s
	static inline T aton(const char *s);

	/// return maximum available finite number
	static inline T max();
	/// return minimum available finite number (positive number closest to 0)
	static inline T min();
	/// return the value for infinity
	static inline T inf();

	/**
	 * greatest common divisor
	 * http://en.wikipedia.org/wiki/Euclidean_algorithm
	 */
	static inline T gcd(T a, T b)
	{
		if (a == 0)	return b;
		while (b != 0)
		{
			if (a > b)
				a -= b;
			else
				b -= a;
		}
		return a;
	}

	T value;	///< the value itself!
};


/** return PI */
template <>	inline float	CMath<float>::PI()
{
	return M_PI;
}

/** return PI of type double */
template <>	inline double	CMath<double>::PI()				{	return M_PI;		}

/** return absolute value of a */
template <>	inline float	CMath<float>::abs(float a)		{	return m_fabsf(a);	}
/** return absolute value of a */
template <>	inline double	CMath<double>::abs(double a)	{	return ::fabs(a);	}
/** return absolute value of a */
template <>	inline int		CMath<int>::abs(int a)			{	return ::abs(a);	}

/** return the power of 'base' to 'exp' */
template <>	inline float	CMath<float>::pow(float base, float exp)	{	return ::powf(base, exp);	}
/** return the power of 'base' to 'exp' */
template <>	inline double	CMath<double>::pow(double base, double exp)	{	return ::pow(base, exp);	}
/** return the power of 'base' to 'exp' */
//template <>	inline int		CMath<int>::pow(int base, int exp)			{	return ::powl(base, exp);	}

/** return floored value of a */
template <>	inline float	CMath<float>::floor(float a)	{	return m_floorf(a);	}
/** return floored value of a */
template <>	inline double	CMath<double>::floor(double a)	{	return ::floor(a);	}
/** return floored value of a */
template <>	inline int		CMath<int>::floor(int a)		{	return a;			}

/** return ceiled value of a */
template <>	inline float	CMath<float>::ceil(float a)		{	return m_ceilf(a);	}
/** return ceiled value of a */
template <>	inline double	CMath<double>::ceil(double a)	{	return ::ceil(a);	}
/** return ceiled value of a */
template <>	inline int		CMath<int>::ceil(int a)			{	return a;			}

/** return maximum finite value of a */
template <>	inline float	CMath<float>::max()				{	return std::numeric_limits<float>::max();		}
/** return maximum finite value of a */
template <>	inline double	CMath<double>::max()			{	return std::numeric_limits<double>::max();		}
/** return maximum finite value of a */
template <>	inline int		CMath<int>::max()				{	return std::numeric_limits<int>::max();			}

/** return minimum available finite number (positive number closest to 0) */
template <>	inline float	CMath<float>::min()				{	return std::numeric_limits<float>::min();		}
/** return minimum available finite number (positive number closest to 0) */
template <>	inline double	CMath<double>::min()			{	return std::numeric_limits<double>::min();		}
/** return minimum available finite number (positive number closest to 0) */
template <>	inline int		CMath<int>::min()				{	return std::numeric_limits<int>::min();			}

/** return minimum available finite number (positive number closest to 0) */
template <>	inline float	CMath<float>::inf()				{	return std::numeric_limits<float>::infinity();		}
/** return minimum available finite number (positive number closest to 0) */
template <>	inline double	CMath<double>::inf()			{	return std::numeric_limits<double>::infinity();		}


/**
 * return ceiled value of a in binary system
 *
 * safety check only valid for unsigned data
 */
template <typename T>	inline T		CMath<T>::ceil2(T a)
		{
			if (a > ((T)1<<(sizeof(T)*8-2)))	return 0;

			T r = 1;
			while (r < a)
			{
				r <<= 1;
			}
			return r;
		}

/**
 * return digits of a in binary system
 *
 * a must not be larger than 2^30!!!
 */
template <>	inline int		CMath<int>::digits2(int a)
		{
			if (a > 0x40000000)	return 0;

			if (a == 0)		return 0;

			int r = 1;
			int c = 1;

			while (r < a)
			{
				r <<= 1;
				c++;
			}
			return c;
		}


/** return square root of a */
template <>	inline float	CMath<float>::sqrt(float a)		{	return m_sqrtf(a);	}
/** return square root of a */
template <>	inline double	CMath<double>::sqrt(double a)	{	return ::sqrt(a);	}
/** return square root of a */
template <>	inline int		CMath<int>::sqrt(int a)			{	return m_sqrtl((unsigned long)a);	}


/** return rounded value of a */
template <>	inline float	CMath<float>::round(float a)	{	return m_roundf(a);	}
/** return rounded value of a */
template <>	inline double	CMath<double>::round(double a)	{	return m_round(a);	}
/** return rounded value of a */
template <>	inline int		CMath<int>::round(int a)		{	return a;			}

/** convert string to type float */
template <>	inline float	CMath<float>::aton(const char *s)	{	return ::atof(s);	}
/** convert string to type double */
template <>	inline double	CMath<double>::aton(const char *s)	{	return ::atof(s);	}
/** convert string to type int */
template <>	inline int		CMath<int>::aton(const char *s)		{	return ::atoi(s);	}

/** sinus of s */
template <>	inline float	CMath<float>::sin(float s)	{	return ::sinf(s);	}
/** cosinus of s */
template <>	inline float	CMath<float>::cos(float s)	{	return ::cosf(s);	}
/** tangens of s */
template <>	inline float	CMath<float>::tan(float s)	{	return ::tanf(s);	}

/** exponential of s */
template <>	inline float	CMath<float>::exp(float s)	{	return ::expf(s);	}


template <typename T>
std::ostream&
operator<<(std::ostream& os, const CMath<T> &lm)
{
	return os << lm.value;
}


#endif
