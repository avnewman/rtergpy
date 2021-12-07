#!/usr/bin/env python3 
# 
from joblib import Parallel, delayed
import numpy as np
import multiprocessing as mp
from obspy import UTCDateTime


def random_square(seed):
    np.random.seed(seed)
    random_num = np.random.randint(0, 10)
    return random_num**2

n_cpu = mp.cpu_count()
print(n_cpu, "CPUs available")
#n_cpu=8
#pool = mp.Pool(processes=n_cpu)
#results = [pool.map(random_square, range(int(1e3)))]

count=int(5e6)
ncores=6
then=UTCDateTime()
results = Parallel(n_jobs=ncores, verbose=1) (delayed(random_square)(i) for i in range(count))
now=UTCDateTime()
print("Parlallel Process took", now-then, "seconds with", ncores)

then=UTCDateTime()
#results = [random_square(i) for i in range(count)]
now=UTCDateTime()
print("Serial Process took", now-then, "seconds")

# sure in this case parallel is a bit faster, but user results vary.  Depending on size, I get better times with 4,6, or 8 cores.
# I get worse times with 16 cores, but I suspect taht's becaue the mac is double-counting actual cores.