# import modules
import melp
import multiprocessing as mp
import subprocess
from glob import glob

#global variables
input_files = []
output_files = []

arg_txt = False
arg_npz = False
arg_binned = False
arg_angle = 0


#defining the functions with multithreading support
def mt_rate(i, npz = arg_npz, binned = arg_binned):
    calls = melp.TileHitRate(input_files[i], output_files[i], output_files[i]+"edep")
    print("read file ", i+1)
    print("started thread ", i+1)
    calls.tileHitRateHID()  
    if npz == True:
        calls.saveNpz()
    elif binned == True:
        calls.saveBinned()


def mt_angle(i):
    txt = arg_txt
    npz = arg_npz
    binned = arg_binned
    angle = arg_angle
    
    calls = melp.TileHitAngle(input_files[i], output_files[i])
    print("read file ", i+1)
    print("started thread ", i+1)
    calls.hitAngleTID(angle)
    if txt == True:
        calls.saveTxt()
    elif npz == True:
        calls.saveNpz()
    elif binned == True:
        calls.saveBinned()
    

def mt_angleHelix(i, txt = arg_txt, npz = arg_npz, binned = arg_binned):
    calls = melp.TileHitAngle(input_files[i], output_files[i])
    print("read file ", i+1)
    print("started thread ", i+1)
    calls.hitAngleHelix()
    if txt == True:
        calls.saveTxt()
    elif npz == True:
        calls.saveNpz()
    elif binned == True:
        calls.saveBinned()

    
# defining which function to run
def run_mt(function_str, src, args): #set func to the one you want to use from above and pass args
    # set used threads
    unused_threads = 2 #set the number of threads you don't want to use
    print("-----------------------")
    print("Available threads = ",mp.cpu_count())
    print("Used threads = ",mp.cpu_count() - unused_threads)
    print("-----------------------")

    #redefine global variables
    arg_txt = args[0]
    arg_npz = args[1]
    arg_binned = args[2]
    arg_angle = args[3]

    # get list of input files
    input_files = glob(src) #set to directory with .root files you want to analyse

    # generate list of output files 
    output_files = []
    for j in range(len(input_files)):
        output_files.append("out"+str(j+1))

    print(input_files)
    print(output_files)


    # map multiprocessing pool
    pool = mp.Pool(mp.cpu_count() - unused_threads) 
    if function_str == "mt_rate":
        pool.map(mt_rate, [i for i in range(len(input_files))])

    elif function_str == "mt_angle":
        pool.map(mt_angle, [i for i in range(len(input_files))])

    elif function_str == "mt_angleHelix": 
        pool.map(mt_angleHelix, [i for i in range(len(input_files))])




#if __name__ == "__main__":
#  run_mt()


