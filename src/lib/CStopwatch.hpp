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


#ifndef CSTOPWATCH_HPP
#define CSTOPWATCH_HPP

#include <sys/time.h>


/**
 * \brief start, stop, continue and restart virtual stopwatch
 */
class CStopwatch
{
private:
	struct timeval timevalue_start;	///< time value of last start
	struct timeval timevalue_stop;	///< time value of last stop

public:
	double time;		///< stopped time

	/**
	 * reset time with 0
	 */
	void reset()
	{
		time = 0.0f;
	}

	/**
	 * initialize time with 0
	 */
	CStopwatch()
	{
		reset();
	}

	/**
	 * start counter
	 */
	inline void start()
	{
		gettimeofday(&timevalue_start, NULL);
	}

	/**
	 * stop counter and add the delta value to the time counter
	 */
	inline void stop()
	{
		gettimeofday(&timevalue_stop, NULL);

		time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
		suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

		time += (double)dsec + (double)dsusec/1000000.0f;
	}

	/**
	 * return the time in seconds
	 */
	inline double getTime()
	{
		gettimeofday(&timevalue_stop, NULL);

		time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
		suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

		return time + (double)dsec + (double)dsusec/1000000.0f;
	}

	/**
	 * return the time in seconds
	 */
	inline double operator()()
	{
		return time;
	}
};

#endif
