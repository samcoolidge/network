# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 15:42:06 2016

@author: samuelcoolidge
"""

from MPI_equilibrate import *
import sys
import time 
from create_network import *

nodes = sys.argv[1]
reps = int(sys.argv[2])/cpu_count() - 1

dense_list = [.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0]

time_list = []

for density in dense_list:
    create(int(nodes),float(density))
    adj_Matrix, k = parse('network_with_' + str(nodes) + '_nodes_and_density_' + str(density) + ".dat")
    t0 = time.time()
    reliability(adj_Matrix, k, reps)
    t1 = time.time()
    t =  t1 - t0
    time_list.append([(density, t1-t0)])
    
f = open("speed_density_test" + nodes + ".dat", "w") 

for time in time_list:
    f.write(str(time) + '\n')
f.close()