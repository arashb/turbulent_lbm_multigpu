import os, subprocess, sys, ConfigParser

keys = ['SECONDS', 'FPS', 'MLUPS', 'BANDWIDTH']
SECTION_BASE = 'EXP'
INI_FILE_DIR = "./profileOutput/"

def get_average_value( config, key, num_exp ):
    sum = 0.0
    for exp in range(1,num_exp+1):
        sum += float( config.get( SECTION_BASE+str(exp), key ) )
    return sum / num_exp
            
def analyse(filenames):
    print "start analysing " , len(filenames), " files." 
    print "files: " , str(filenames)
    config = ConfigParser.ConfigParser()
    # gres: contains the analyzation results of all profiling cases
    gres = dict();
    for filename in filenames:
        print "analysing file: ", filename
        config.read(filename)
        num_exp = len(config.sections())
        # pres: contains the average values of the metrics for profiling case
        pres = dict()
        # get the key of current profiling case
        proc_key = config.get( SECTION_BASE+str(1), 'NP' )
        # compute the average value of profiling metrics for the current profiling case
        for key in keys:
             pres[key]= str(get_average_value( config, key, num_exp))
        gres[proc_key] = pres
    return gres

def execute(command):
    print "executing command: ", command
    os.system(command)

def profile(max_num_proc = 1, num_exp = 1):
    SRUN = "srun -n "
    LBM_COMMAND = "./build/lbm_opencl_dc_mpicxx_release"
    DOMAIN_LENGTH = 0.1;
    for num_proc in range(1,max_num_proc+1):
        x_size = 32*num_proc;
        command_str = SRUN + str(num_proc)+ " " + LBM_COMMAND +  " -x " + str(x_size) +  " -X " +  str(num_proc) + " -l 100"
        command_str +=  " -n " + str(DOMAIN_LENGTH*num_proc) + " -m " + str(DOMAIN_LENGTH) + " -p " + str(DOMAIN_LENGTH)
        ini_filename = INI_FILE_DIR + "profile_" + str(num_proc)+ ".ini"
        for exp_counter in range(1,num_exp+1):
            f = open(ini_filename, 'a')
            f.write("["+SECTION_BASE+str(exp_counter)+"]\n")
            f.write("NP : " + str(num_proc))
            f.write("\n")
            f.close()
            execute(command_str)

def visualize(profiling_res):
    try:
        // TODO: test plotting of the files
        import numpy as np
        import matplotlib.pyplot as plt
        ngpus = profiling_res.keys()
        for key in keys:
            fig = plt.figure()
            values = [res[key] for proc, res in profiling_res.iteritems()]
            plt.plot(ngpus, values,'g^')
            fig.savefig('plot_'+ key+ '.png')

    except ImportError:
        print "Could not import numpy/matplotlib module"
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(profiling_res)
        for key in keys:
            values = [res[key] for proc, res in profiling_res.iteritems()]
            print key, ":", values

if __name__ == "__main__":
    import glob
    filenames = glob.glob(INI_FILE_DIR+"*.ini")
    ask = True
    do_profile = True
    while ask:
        if len(filenames) != 0:
            print "The profile output directory is not empty. what should I do with older files?"
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
                do_profile = False
            else:
                print "wrong input."
    max_num_proc = int(sys.argv[1])
    num_exp = int(sys.argv[2])
    if do_profile:
        profile(max_num_proc, num_exp)
    filenames = glob.glob(INI_FILE_DIR+"*.ini")
    res = analyse(filenames)
    visualize(res)
