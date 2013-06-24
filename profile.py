#!/usr/bin/python
import os, subprocess, sys, ConfigParser, json

keys = ['NAME', 'START', 'END', 'DURATION']
SECTION_BASE = 'EVENT'
INI_FILE_DIR = "./profileOutput/"

METADATA_SECTION = 'METADATA'

class ProfilerEvent:
    def __init__(self, uuid, event_type, name, start, end, duration):
        self.uuid = uuid
        self.event_type = event_type
        self.name = name
        self.start = start
        self.end = end
        self.duration = duration

    def __repr__(self):
        return "ProfilerEvent()"
    
    def __str__(self):
        event_str =  "ProfilerEvent:\n"
        event_str += "\tuuid : " + str(self.uuid) + "\n"
        event_str += "\ttype : " + self.event_type + "\n"
        event_str += "\tname : " + self.name + "\n"
        event_str += "\ttime : [" + str(self.start) + ", " + str(self.end) + "]\n"
        event_str += "\tduration : " + str(self.duration) + "\n"
        return event_str

    def overlap(self, event):
        if self.start < event.end and self.end > event.start:
            return True
        return False

def get_current_proc_id_events(config):
    current_events = []
    for event_counter in range( 1, len(config.sections()) ): # first element is a meta section
        event_section_name = "EVENT"+str(event_counter) 
        name = config.get( event_section_name, 'NAME')
        start = config.getint( event_section_name, 'START' )
        end = config.getint( event_section_name, 'END')
        duration = config.getfloat( event_section_name, 'DURATION' )
        event_type = config.get( event_section_name, 'TYPE' )
        # todo: implement other types of events (e.g. host functions)
        event = ProfilerEvent(event_counter, event_type, name, start, end, duration)
        current_events.append(event)
    return current_events

def get_overlapping_events(events_list):
    overlapping_events = []
    for cur_event_idx in range(0, len(events_list)):
        for next_event_idx in range(cur_event_idx+1, len(events_list)):
            current_event = events_list[cur_event_idx]
            next_event = events_list[next_event_idx]
            if current_event.overlap(next_event):
                overlapping_events.append((current_event, next_event))
    return overlapping_events

def analyse(filenames):
    print "start analysing " , len(filenames), " files." 
    print "files: " , str(filenames)
    config = ConfigParser.ConfigParser()
    for filename in filenames:
        print "\nanalysing file: ", filename
        config.read(filename)
        total_num_proc = config.getint( METADATA_SECTION, 'TOTAL_NUM_PROC' )
        current_proc_id = config.getint( METADATA_SECTION, 'CURRENT_PROC_ID' )
        print "TOTAL NUMBER OF PROCESSES: ", total_num_proc
        print "CURRENT PROCCESSOR ID: ", current_proc_id
        current_proc_events = get_current_proc_id_events(config)
        print "# EVENTS: ", str(len(current_proc_events))
        overlapping_events = get_overlapping_events(current_proc_events)
        print "# OEVERLAPPING EVENTS FOUND: ", str(len(overlapping_events))

def visualize(results):
    import pprint
    print "profiling results pretty print: "
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(profiling_res)

if __name__ == "__main__":
    import glob
    filenames = glob.glob(INI_FILE_DIR+"*.ini")
    ask = True
    do_benchmark = True
    if len(filenames) != 0:
        while ask:
            print "The profile output directory is not empty. What should I do with older files?"
            input_variable = raw_input("(a)ppend, a(r)chive, (d)elete, (u)se?")
            if input_variable is "a":
                ask = False
            elif input_variable is "r":
                ask = False
                try:
                    import tarfile
                except ImportError:
                    print "not archiving module available."
                import datetime
                now = datetime.datetime.now()
                tar = tarfile.open(INI_FILE_DIR+"archive" + now.strftime("%Y%m%d_%H%M")+ ".tar.gz", "w:gz")
                for f in filenames:
                    print "archiving file: ", f
                    tar.add(f)
                tar.close()
                for f in filenames:
                    os.remove(f)
            elif input_variable is "d":
                ask = False
                for f in filenames:
                    print "removing file: ", f
                    os.remove(f)
            elif input_variable is "u":
                ask = False
                do_benchmark = False
            else:
                print "wrong input."
    filenames = glob.glob(INI_FILE_DIR+"*.ini")
    analyse(filenames)
