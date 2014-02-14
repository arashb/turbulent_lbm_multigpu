#!/usr/bin/python
import os, subprocess, sys, ConfigParser, json, math

keys = ['CUBE_X', 'CUBE_Y', 'CUBE_Z', 'SECONDS', 'FPS']#, 'MLUPS', 'BANDWIDTH']
SECTION_BASE = 'EXP'
INI_FILE_DIR = "./output/benchmark/"
MPI_COMMAND = "mpirun"
LBM_COMMAND = "./build/lbm_opencl_dc_mpicxx_release"

def get_average_value( config, key, num_exp ):
    sum = 0.0
    for exp in range(1,num_exp+1):
        sum += config.getfloat( SECTION_BASE+str(exp), key ) 
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
        proc_key = config.getint( SECTION_BASE+str(1), 'NP' )
        # compute the average value of profiling metrics for the current profiling case
        for key in keys:
             pres[key]= get_average_value( config, key, num_exp)
        gres[proc_key] = pres
    return gres

def execute(command):
    print "executing command: ", command
    os.system(command)

def weak_scaling_benchmark_2d(max_num, num_exp, loops = 100, grid_size_increase_step = 1024, base_domain_length = 0.1 ):
    DOMAIN_LENGTH = base_domain_length
    for num_increase in range(1,max_num+1):
        num_proc = num_increase*num_increase
        # computing the grid size for current benchmark
        x_size = grid_size_increase_step*num_increase
        y_size = x_size
        # generating the string of execution command  
        command_str = MPI_COMMAND + " -n " + str(num_proc)+ " -N 1 " + LBM_COMMAND 
        command_str += " -x " + str(x_size) +  " -X " +  str(num_increase) 
        command_str += " -y " + str(y_size) +  " -Y " +  str(num_increase) 
        command_str += " -l " + str(loops)
        command_str +=  " -n " + str(DOMAIN_LENGTH*num_increase) + " -m " + str(DOMAIN_LENGTH*num_increase) + " -p " + str(DOMAIN_LENGTH)
        ini_filename = INI_FILE_DIR + "benchmark_" + str(num_proc)+ ".ini"
        # performing the experiments
        for exp_counter in range(1,num_exp+1):
            f = open(ini_filename, 'a')
            f.write("["+SECTION_BASE+str(exp_counter)+"]\n")
            f.write("NP : " + str(num_proc))
            f.write("\n")
            f.close()
            execute(command_str)

def weak_scaling_benchmark_1d(max_num, num_exp, loops = 100, grid_size_increase_step = 1024, base_domain_length = 0.1 ):
    DOMAIN_LENGTH = base_domain_length
    for num_increase in range(1,max_num+1):
        num_proc = num_increase
        # computing the grid size for current benchmark
        x_size = grid_size_increase_step*num_increase
        y_size = 1024
        # generating the string of execution command  
        command_str = MPI_COMMAND + " -n " + str(num_proc)+ " -N 1 " + LBM_COMMAND 
        command_str += " -x " + str(x_size) +  " -X " +  str(num_increase) 
        command_str += " -y " + str(y_size) #+  " -Y " +  str(num_increase) 
        command_str += " -l " + str(loops)
        command_str +=  " -n " + str(DOMAIN_LENGTH*num_increase) + " -m " + str(DOMAIN_LENGTH) + " -p " + str(DOMAIN_LENGTH)
        ini_filename = INI_FILE_DIR + "benchmark_" + str(num_proc)+ ".ini"
        # performing the experiments
        for exp_counter in range(1,num_exp+1):
            f = open(ini_filename, 'a')
            f.write("["+SECTION_BASE+str(exp_counter)+"]\n")
            f.write("NP : " + str(num_proc))
            f.write("\n")
            f.close()
            execute(command_str)

def strong_scaling_benchmark_2d(max_num, num_exp, loops = 100, grid_size = 128, base_domain_length = 0.1 ):
    DOMAIN_LENGTH = base_domain_length
    for num_increase in range(1,max_num+1):
        num_proc = num_increase*num_increase
        # computing the grid size for current benchmark
        x_size = grid_size
        y_size = grid_size
        # generating the string of execution command  
        command_str = MPI_COMMAND + " -n " + str(num_proc)+ " " + LBM_COMMAND 
        command_str += " -x " + str(x_size) +  " -X " +  str(num_increase) 
        command_str += " -y " + str(y_size) +  " -Y " +  str(num_increase) 
        command_str += " -l " + str(loops)
        command_str +=  " -n " + str(DOMAIN_LENGTH) + " -m " + str(DOMAIN_LENGTH) + " -p " + str(DOMAIN_LENGTH)
        ini_filename = INI_FILE_DIR + "benchmark_" + str(num_proc)+ ".ini"
        # performing the experiments
        for exp_counter in range(1,num_exp+1):
            f = open(ini_filename, 'a')
            f.write("["+SECTION_BASE+str(exp_counter)+"]\n")
            f.write("NP : " + str(num_proc))
            f.write("\n")
            f.close()
            execute(command_str)

def strong_scaling_benchmark_1d(max_num, num_exp, loops = 100, grid_size = 1024, base_domain_length = 0.1 ):
    DOMAIN_LENGTH = base_domain_length
    for num_increase in range(0,max_num+1):
        num_proc = int(math.pow(2,num_increase))#num_increase*num_increase
        # computing the grid size for current benchmark
        x_size = grid_size
        y_size = grid_size
        # generating the string of execution command  
        command_str = MPI_COMMAND + " -n " + str(num_proc)+ " " + LBM_COMMAND 
        command_str += " -x " + str(x_size) +  " -X " +  str(num_proc) 
        command_str += " -y " + str(y_size)# +  " -Y " +  str(num_increase) 
        command_str += " -l " + str(loops)
        command_str +=  " -n " + str(DOMAIN_LENGTH) + " -m " + str(DOMAIN_LENGTH) + " -p " + str(DOMAIN_LENGTH)
        ini_filename = INI_FILE_DIR + "benchmark_" + str(num_proc)+ ".ini"
        # performing the experiments
        for exp_counter in range(1,num_exp+1):
            f = open(ini_filename, 'a')
            f.write("["+SECTION_BASE+str(exp_counter)+"]\n")
            f.write("NP : " + str(num_proc))
            f.write("\n")
            f.close()
            execute(command_str)

def benchmark(benchmark_strategy, max_num = 1, num_exp = 1):
    benchmark_strategy(max_num, num_exp)

def visualize(profiling_res):
    import pprint
    print "profiling results pretty print: "
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(profiling_res)
    try:
        import numpy as np
        import matplotlib.pyplot as plt
        ngpus = sorted(profiling_res.keys())
        runtime_values = []
        print "visualizing keys: ", ngpus
        for key in keys:
            values = [float(res[key]) for proc, res in sorted(profiling_res.iteritems())]
            if key is "SECONDS":
                runtime_values = values
            fig = plt.figure()
            plt.title(key + ' scaling')
            plt.xlabel('# GPUs')
            plt.ylabel(key)
            plt.grid(True)
            plt.plot(ngpus, values, marker='^', linestyle='--', color='g' )
            x1,x2,y1,y2 = plt.axis()
            plt.axis((min(ngpus)-1,max(ngpus)+1,y1,y2))
            filename = INI_FILE_DIR+'plot_'+ key+ '.png'
            fig.savefig(filename)
            print "saved graph: ", filename
        # computing speedup
        speedup_values = []
        for val in runtime_values:
            speedup_values.append(runtime_values[0]/val)
        fig = plt.figure()
        # TODO: implement diffent speedup for weak/strong scaling
        plt.title('Speedup Scaling')
        plt.xlabel('# GPUs')
        plt.ylabel("speedup")
        plt.grid(True)
        #values = [res[key] for proc, res in sorted(profiling_res.iteritems())]
        plt.plot(ngpus, speedup_values, marker='^', linestyle='--', color='g' )
        x1,x2,y1,y2 = plt.axis()
        plt.axis((min(ngpus)-1,max(ngpus)+1,y1,y2))
        filename = INI_FILE_DIR+'plot_speedup.png'
        fig.savefig(filename)
        print "saved graph: ", filename

    except ImportError:
        print "Could not import numpy/matplotlib module"
        for key in keys:
            values = [float(res[key]) for proc, res in sorted(profiling_res.iteritems())]
            print key, ":", values
        print 'MLUPS: ' , [float(res['MLUPS']) for proc, res in sorted(profiling_res.iteritems())]

if __name__ == "__main__":
    import glob
    filenames = glob.glob(INI_FILE_DIR+"*.ini")
    ask = True
    do_benchmark = True
    if len(filenames) != 0:
        while ask:
            print "The benchmark output directory is not empty. What should I do with older files?"
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
    if not os.path.exists(INI_FILE_DIR):
        os.makedirs(INI_FILE_DIR)
    if do_benchmark:
        max_num = int(sys.argv[1])
        num_exp = int(sys.argv[2])
        MPI_COMMAND = sys.argv[3]
        LBM_COMMAND = sys.argv[4]
        BENCHMARK_STRATEGY = str(sys.argv[5])
        print "benchmark strategy: ", BENCHMARK_STRATEGY
        if BENCHMARK_STRATEGY == "weak":
            #benchmark(weak_scaling_benchmark_2d, max_num, num_exp)
            benchmark(weak_scaling_benchmark_1d, max_num, num_exp)
        elif BENCHMARK_STRATEGY == "strong":
            benchmark(strong_scaling_benchmark_1d, max_num, num_exp)
        else:
            print "Unknown benchmarking strategy!"
            sys.exit(0)
    filenames = glob.glob(INI_FILE_DIR+"*.ini")
    res = analyse(filenames)
    print "RES"
    for p in range(1,len(res)+1):
        res[p]['MLUPS'] = res[p]['CUBE_X'] * res[p]['CUBE_Y'] * res[p]['CUBE_Z'] * res[p]['FPS'] * 0.000001
    # saving the results in json format
    analysation_file_name = INI_FILE_DIR + "results.txt"
    with open(analysation_file_name,'w') as outfile:
        json.dump(res,outfile, indent = 4)
        print "saved analysation results in file:", analysation_file_name
    visualize(res)
