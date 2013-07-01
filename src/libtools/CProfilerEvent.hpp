
#ifndef CPROFILEREVENT_HPP
#define CPROFILEREVENT_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "../libcl/CCLErrors.hpp"
#include "CStopwatch.hpp"

typedef enum {
	EVENT_TYPE_UNKNOWN = 0,
	EVENT_TYPE_DEVICE_KERNEL,
	EVENT_TYPE_HOST_FUNCTION,
} EVENT_TYPE;

class CProfilerEvent
{
private:
  static int event_counter;
  int _uuid;
  EVENT_TYPE _type;
  std::string _event_id;
  cl_ulong _start; // start time in milliseconds
  cl_ulong _end;  // end time in milliseconds
  cl_float  _duration; // duration time in milliseconds

  void get_profile_info(cl_event* _clevent) {

    CL_CHECK_ERROR(clGetEventProfilingInfo(*_clevent, CL_PROFILING_COMMAND_START,
					   sizeof(_start), &_start, NULL));

    CL_CHECK_ERROR(clGetEventProfilingInfo(*_clevent, CL_PROFILING_COMMAND_END,
					   sizeof(_end), &_end, NULL));

    _duration = (_end - _start)/1000000.0;
  }

public:
  CProfilerEvent(std::string event_id, cl_event* event)
  {
    if (event_id.empty())
      throw "CProfilerEvent: ID of the event is unknown!";
    if (!event)
      throw "CProfilerEvent: event is not valid!";

    event_counter++;
    _type = EVENT_TYPE_DEVICE_KERNEL;
    _event_id = event_id;
    _uuid = event_counter;
    get_profile_info(event);
  }

  ~CProfilerEvent()
  {
    event_counter--;
  }

  int getUuid() {
    return _uuid;
  }

  std::string getEventId() {
    return _event_id;
  }
  
  cl_ulong getEventStartTime() {
    return _start;
  }

  cl_ulong getEventEndTime() {
    return _end;
  }

  cl_ulong getEventDuration() {
    return _duration;
  }
  
  void printEvent(std::ofstream& prof_file) {
    	prof_file << "[EVENT" << _uuid << "]" << std::endl;
	prof_file << "TYPE : " << ((_type == EVENT_TYPE_DEVICE_KERNEL)? "DEVICE_KERNEL":"HOST_FUNCTION") << std::endl;;
	prof_file << "NAME : " << _event_id << std::endl;
	prof_file << "# start/end in nanoseconds" << std::endl;
	prof_file << "START : " << _start << std::endl;
	prof_file << "END : " << _end << std::endl;
	prof_file << "# duration in milliseconds" << std::endl;
	prof_file << "DURATION : " << _duration << std::endl;
	prof_file << std::endl;
  }
};

#endif
