'''
Created on 20 févr. 2020

@author: Abdelkader Ouali
update: Arnold Hien (09/03/2020)

'''

from include import *
import multiprocessing as mp

# number parallel launches
npr = 1

# program parameters 
fthresholds = {
#                "BMS1.dat" : ["0.0015", "0.0014", "0.0012"],
#                "connect.dat" : ["0.3", "0.18", "0.15"],
#                "hepatitis.dat" : ["0.3", "0.2", "0.1"],
#                "mushroom.dat" : ["0.05", "0.01", "0.008", "0.005"],
#                "pumsb.dat" : [ "0.4", "0.3", "0.2"],
#                "retail.dat" : ["0.05", "0.01"],
#                "T10I4D100K.dat" : [ "0.05", "0.01", "0.005" ],
#                "chess.dat" : ["0.2", "0.15", "0.1"],
#                "heart-cleveland.dat" : [ "0.2", "0.1", "0.08", "0.06"],
#                "kr-vs-kp.dat" : ["0.4", "0.3", "0.2"],
#                "splice1.dat" : [ "0.1", "0.05", "0.02", "0.01" ],
#                "T40I10D100K.dat" : ["0.1", "0.08", "0.05", "0.01"]
#                "ijcai16.dat" : ["0.1"],
                "these22.dat" : ["0.1"],
               }

def nbrTrans(absolute_input_dataset_file):
    with open(absolute_input_dataset_file, "r") as infile:
        nbr_trans = 0
        for line in infile:
            if line.rstrip() and line[0] != '@':
                nbr_trans += 1
        return nbr_trans


def krimp_convert(queue, project_bin_file, config):
    '''
    config is a list of (filename, freq_threshold, absolute_freq_threshold, jmax, historyAggregator, branchStrategy, number_threads and timeout))
    '''
    try:
        # preparing directories
        absolute_bin_file = os.path.abspath(os.path.join(project_bin_dir, project_bin_file))
        
        if (not os.path.exists(absolute_bin_file)):
            print("Error: needed files for launching are missing!")
            sys.exit(1)
        
        convertText = ""
        
        # print(absolute_bin_file)
        with open(absolute_bin_file, 'r') as binfile:
            for line in binfile:
                line = line.replace("\n", "")
                if line.startswith("dbName = "):
                    convertText = convertText + "dbName = " + config[0].split('.')[0] + "\n"
                    # convertText = convertText + "-closed-" + config[2] + "d" + "\n"
                else:
                    convertText += line + "\n"
            
        # print("\n\n#########################\n#########################\n", compressText, "\n\n#########################\n#########################\n")
        with open(absolute_bin_file, 'w') as binfile:
            binfile.write(convertText)
        
        print("COMMAND = \'../bin/krimp convertdb ", config[0], "\'")
        # subprocess.run(["../bin/krimp compress"], shell=True, check=True, stdout=launch_log_file, stderr=launch_log_file)
        subprocess.run(["../bin/krimp convertdb"], shell=True, check=True)
    finally:
        queue.put(mp.current_process().name)    


if __name__ == '__main__':
    
    # print date and time
    timenow = datetime.datetime.now()
    print(timenow.strftime("%d/%m/%Y %H:%M:%S"), "\nnpr =", npr, "\n", fthresholds, "\n")
    # checking required needed directories and files 
    if (not os.path.exists(project_data_dir)):
        print("Error: the data directory was not found")
        sys.exit(1)
    # if (not os.path.exists(project_res_dir)):
    #    os.makedirs(project_res_dir, exist_ok=True)
    # prepare the different configurations
    configs = []
    for filename in os.listdir(project_data_dir):
        if filename.endswith(".dat"):
            if filename in fthresholds:
                configs.append((filename,))
#                for threshold in fthresholds[filename]:
#                    absolute_dataset_file = os.path.abspath(os.path.join(project_data_dir, filename))
#                    absolute_threshold = int(math.ceil(float(threshold) * nbrTrans(absolute_dataset_file)))
#                    configs.append((filename, threshold, str(absolute_threshold)))
    # launch a chunk of processes all at once
    curr = 0
    procs = dict()
    queue = mp.Queue()
    convert_file = "convertdb.conf"
    for i in range(0, npr):
        if curr < len(configs):
            # launch
            print(configs[i][0], configs[curr])
            proc = Process(target=krimp_convert, args=(queue, convert_file, configs[curr],))
            proc.start()
            procs[proc.name] = proc
            curr += 1
    # using the queue, launch a new process whenever an old process finishes its workload
    while procs:
        name = queue.get()
        proc = procs[name]
        print(proc) 
        proc.join()
        del procs[name]
        if curr < len(configs):
            print(configs[i][0], configs[curr])
            proc = Process(target=krimp_convert, args=(queue, convert_file, configs[curr],))
            proc.start()
            procs[proc.name] = proc
            curr += 1
    print("finished")
    # subprocess.run(["/bin/chmod 770 -R" + project_res_dir + "; echo \"Exit status: $?\" "], shell=True, check=True)
