
#ifndef CPROFILER_HPP
#define CPROFILER_HPP

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <fstream>

#include "CProfilerEvent.hpp"

typedef std::multimap<int, CProfilerEvent*> event_map;
typedef std::multimap<int, CProfilerEvent*>::value_type event_map_pair;
typedef std::multimap<int, CProfilerEvent*>::iterator event_map_ptr;

class CProfiler
{
private:
  event_map _event_container;
  
public:
  CProfiler()
  {

  }

  ~CProfiler()
  {
    event_map_ptr it = _event_container.begin();
    for( ;it != _event_container.end(); it++){
      delete (*it).second;
    }
  }

  void addProfilerEvent(CProfilerEvent* profEvent) {
    _event_container.insert(event_map_pair(profEvent->getUuid(), profEvent));
  }

  void saveEvents(std::string file_name) {
    std::ofstream prof_file (file_name.c_str(), std::ios::out | std::ios::app );
    if (prof_file.is_open()) {
      event_map_ptr it = _event_container.begin();
      for( ;it != _event_container.end(); it++) {
	(*it).second->printEvent(prof_file);
      }
    } else std::cout << "Unable to open file: " << file_name << std::endl;
  }

};

#endif
