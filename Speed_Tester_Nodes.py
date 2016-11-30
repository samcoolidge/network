# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 11:05:13 2016

@author: samuelcoolidge
"""

from Par_Equilibrate import *
import sys
import time 

density = sys.argv[1]
reps = int(sys.argv[2])/cpu_count() - 1

from create_network import *



node_list = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

time_list = []

for nodes in node_list:
    create(nodes,float(density))
    adj_Matrix, k = parse('network_with_' + str(nodes) + '_nodes_and_density_' + str(density) + ".dat")
    t0 = time.time()
    reliability(adj_Matrix, k, reps)
    t1 = time.time()
    t =  t1 - t0
    time_list.append([(nodes, t1-t0)])
    
f = open("speed_density_test" + str(density) + ".dat", "w") 

for time in time_list:
    f.write(str(time) + '\n')
f.close()