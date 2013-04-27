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


class CStopwatch
{
private:
	struct timeval timevalue_start;
	struct timeval timevalue_stop;

public:
	double time;

	void reset()
	{
		time = 0.0f;
	}

	CStopwatch()
	{
		reset();
	}

	inline void start()
	{
		gettimeofday(&timevalue_start, NULL);
	}

	inline void stop()
	{
		gettimeofday(&timevalue_stop, NULL);

		time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
		suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

		time += (double)dsec + (double)dsusec/1000000.0f;
	}

	inline double getTime()
	{
		gettimeofday(&timevalue_stop, NULL);

		time_t dsec = timevalue_stop.tv_sec - timevalue_start.tv_sec;
		suseconds_t dsusec = timevalue_stop.tv_usec - timevalue_start.tv_usec;

		return time + (double)dsec + (double)dsusec/1000000.0f;
	}
};

#endif
