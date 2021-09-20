# import modules
import melp
import multiprocessing as mp
import subprocess

# set used threads
unused_threads = 2 #set the number of threads you don't want to use
print("Available threads = ",mp.cpu_count())
print("Used threads = ",mp.cpu_count() - unused_threads)

# list of input files
input_files =  ["sorted1.root", 
                "sorted2.root",
                "sorted3.root",
                "sorted4.root",
                "sorted5.root",
                "sorted6.root",
                "sorted7.root",
                "sorted8.root",
                "sorted9.root",
                "sorted10.root"]

# list of output files 
output_files = ["out1", "out2", "out3", "out4", "out5", "out6", "out7", "out8", "out9", "out10"]


#defining the functions with multithreading support
def multithreading_angle(i, angle = "norm"):
    calls = melp.TileHitAngle(input_files[i], output_files[i])
    print("read file ", i+1)
    print("started thread ", i+1)
    calls.hitAngleTID(angle)
    calls.saveTxt()

def multithreading_rate(i):
    calls = melp.input_files[i], output_files[i], output_files[i]+"edep")
    print("read file ", i+1)
    print("started thread ", i+1)
    calls.tileHitRateHID()  
    calls.saveNpz()

def multithreading_angleHelix(i):
    calls = melp.TileHitAngle(input_files[i], output_files[i])
    print("read file ", i+1)
    print("started thread ", i+1)
    calls.hitAngleHelix()
    calls.saveBinned()

    
# defining which function to run
def main(func = multithreading_angleHelix):
    pool = mp.Pool(mp.cpu_count() - unused_threads) 
    pool.map(func, [i for i in range(len(input_files))])

if __name__ == "__main__":
  main()


