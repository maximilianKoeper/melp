# import modules
import melp
import multiprocessing as mp
import subprocess
from glob import glob
from functools import partial
import numpy as np

#global variables
input_files = []
output_files = []


#defining the functions with multithreading support
def mt_rate(npz, binned, i):
    calls = melp.TileHitRate(input_files[i], output_files[i], output_files[i]+"edep")
    print("read file ", i+1)
    print("started thread ", i+1)
    calls.tileHitRateHID()
    if npz == True:
        calls.saveNpz()
    if binned == True:
        calls.saveBinned()


def mt_angle(txt, npz, binned, angle, i):
    calls = melp.TileHitAngle(input_files[0][i], output_files[i])
    print("read file ", i+1)
    print("started thread ", i+1)
    calls.hitAngleTID(angle=angle)
    if txt == True:
        calls.saveTxt()
    if npz == True:
        calls.saveNpz()
    if binned == True:
        calls.saveBinned()


def mt_angleHelix(txt, npz, binned, i):
    calls = melp.TileHitAngle(input_files[0][i], output_files[i])
    print("read file ", i+1)
    print("started thread ", i+1)
    calls.hitAngleHelix()
    if txt == True:
        calls.saveTxt()
    if npz == True:
        calls.saveNpz()
    if binned == True:
        calls.saveBinned()


# defining which function to run
def run_mt(function_str, src, args): #set func to the one you want to use from above and pass args
    # set used threads
    unused_threads = 2 #set the number of threads you don't want to use
    print("-----------------------")
    print("Available threads = ",mp.cpu_count())
    print("Used threads = ",mp.cpu_count() - unused_threads)
    print("-----------------------")


    # get list of input files
    input_files.append(glob(src)) #set to directory with .root files you want to analyse

    # generate list of output files
    for j in range(len(input_files[0])):
        output_files.append("out"+str(j+1))


    # map multiprocessing pool
    pool = mp.Pool(mp.cpu_count() - unused_threads)
    if function_str == "mt_rate":
        func = partial(mt_rate, *args)
        pool.map(func, [i for i in range(len(input_files[0]))])

    elif function_str == "mt_angle":
        func = partial(mt_angle, *args)
        pool.map(func, [i for i in range(len(input_files[0]))])

    elif function_str == "mt_angleHelix":
        func = partial(mt_angleHelix, *args)
        pool.map(func, [i for i in range(len(input_files[0]))])




#if __name__ == "__main__":
#  run_mt()
