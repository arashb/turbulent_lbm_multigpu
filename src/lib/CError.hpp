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


#ifndef CERROR_HPP
#define CERROR_HPP

#include <string>
#include <iostream>
#include <sstream>

/**
 * \brief	error message class
 *
 * this class helps the programmer to output error (or information) messages.
 *
 * the error messages have to be read by the calling class using getString().
 *
 * helper macros exist to avoid duplicate code.
 *
 * examples:
 *
 * 		// pipe error message to the 'error' class
 * 		error << "This is an Error Message" << std::endl;
 *
 * 		// pipe className.error message to error class and return in case of an error in className.error
 * 		CError_AppendReturn(className);
 */
class CError
{
protected:
	std::ostringstream errorStream;		///< error output stream
	bool errorFlag;						///< true, if error message was written to 'errorStream'

public:
	/**
	 * endl - use CError::endl for line termination
	 *
	 * \param __os	output stream
	 */
	static std::ostream& endl(std::ostream& __os)	///< pointer to endl (use CError::endl). Same as std::endl
	{
		return __os;
	}

	/**
	 * return error flag
	 */
	bool operator()()
	{
		return errorFlag;
	}

	CError()
	{
		errorFlag = false;
	}

	/**
	 * append string s
	 *
	 * \param s	string to append
	 */
	CError& operator<<(const std::string& s)
	{
		errorStream << s;
		errorFlag = true;

		return *this;
	}

	/**
	 * append string pointed by s
	 *
	 * \param s	string to append
	 */
	CError& operator<<(const char *s)
	{
		errorStream << s;
		errorFlag = true;

		return *this;
	}

	/**
	 * append number of type size_t
	 *
	 * \param i	number to append
	 */
#if __LP64__
	CError& operator<<(const size_t i)
	{
		errorStream << i;
		errorFlag = true;

		return *this;
	}
#endif

#if !WIN32
	/**
	 * append number of type unsigned int
	 *
	 * \param i	number to append
	 */
	CError& operator<<(const unsigned int i)
	{
		errorStream << i;
		errorFlag = true;

		return *this;
	}
#endif

	/**
	 * append integer number
	 *
	 * \param i	number to append
	 */
	CError& operator<<(const int i)
	{
		errorStream << i;
		errorFlag = true;

		return *this;
	}

	/**
	 * append float number
	 *
	 * \param f	float number to append
	 */
	CError& operator<<(const float f)
	{
		errorStream << f;
		errorFlag = true;

		return *this;
	}


	/**
	 * append double number to error string
	 *
	 * \param f	double number to append
	 */
	CError& operator<<(const double f)
	{
		errorStream << f;
		errorFlag = true;

		return *this;
	}

	/**
	 * catch std::endl
	 *
	 * \param pf	outputstream
	 */
	CError& operator<<(std::ostream &(*pf)(std::ostream &))
	{
		if (pf == &CError::endl)
			errorStream << std::endl;
		return *this;
	}

	/**
	 * return the error string and set error flag to false
	 */
	std::string getString()
	{
		std::string errorMessage = errorStream.str();
		errorStream.str("");

		errorFlag = false;
		return errorMessage;
	}
};


typedef CError CMessage;	///< create CMessages with same functionality as CError class to handle simple messages

/**
 * little helpers for convenience
 */

#define CError_Location 	__FILE__ << ": " << __FUNCTION__ << " (" << __LINE__ <<"):"
#define CError_AppendCode(class_name)	error << "Error in class '" << #class_name "' created in "<< CError_Location << CError::endl
#define CError_AppendCodeThis(class_name)	this->error << "Error in class '" << #class_name "' created in "<< CError_Location << CError::endl

#define CError_AppendReturn(class_name)					\
		if (class_name.error())							\
		{												\
			CError_AppendCode(class_name) << class_name.error.getString();				\
			return;										\
		}

#define CError_AppendReturnThis(class_name)					\
		if (class_name.error())							\
		{												\
			CError_AppendCodeThis(class_name) << class_name.error.getString();				\
			return;										\
		}

// error reporting with prefix string
#define CError_PrefixAppendReturn(prefix, class_name)					\
		if (class_name.error())							\
		{												\
			CError_AppendCode(class_name) << prefix << class_name.error.getString();				\
			return;										\
		}

#define CError_PtrAppendReturn(class_name)				\
		if (class_name->error())						\
		{												\
			CError_AppendCode(class_name) << class_name->error.getString();				\
			return;										\
		}


#define CError_Append(class_name)						\
		if (class_name.error())							\
		{												\
			CError_AppendCode(class_name) << class_name.error.getString();				\
		}


#define CError_PtrAppend(class_name)					\
		if (class_name->error())						\
		{												\
			CError_AppendCode(class_name) << class_name->error.getString();				\
		}

#endif
