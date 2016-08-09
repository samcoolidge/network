# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 14:33:13 2016

@author: samuelcoolidge
"""
from Functions import *
from MC_step_cython import *
import numpy as np
import operator as op
import time
import sys


file_name = sys.argv[1]
reps = int(sys.argv[2])


# Find number of nodes
maxim = 0;    
for line in open(file_name, 'r'):
    edge = map(int, line.split())
    if max(edge) > maxim:
        maxim = max(edge)

k = maxim; 
 

# initialize adjacency matrix 
adj_Matrix = set();

# initialize link reliabilities
reliability_list_seq = np.zeros((k,k))

#convert file to adjacency matrix   
for line in open(file_name, 'r'):
    edge = map(int, line.split())
    adj_Matrix.add((edge[0]-1,edge[1]-1))
    adj_Matrix.add((edge[1]-1,edge[0]-1))
    
for i in range(k):
    if (i,i) in adj_Matrix:
        print "Error: Network Contains Self-Linking Node"
        exit()


#create numpy adjacency matrix
adj_Matrix_np = np.zeros(shape = (k,k), dtype = int)
for (i,j) in adj_Matrix:
    adj_Matrix_np[i][j] = 1

#get decorrelation step
decor_step = get_decorrelation_step(adj_Matrix,adj_Matrix_np,k)

print "decorellation step" , decor_step

t_10 = time.time()

#get equilibrated starting partition
start_parts, g2g, group_size_vector,node_to_group, group_of_node, nonempty_groups, H = get_sample(adj_Matrix,adj_Matrix_np,decor_step, k)  

t_11 = time.time()

print "time for starting partition/equilibration" , t_11-t_10

t20 = time.time()

#get the decorralated partitions
sample,reliability_list = add_sample(start_parts,adj_Matrix,adj_Matrix_np, k, reps, decor_step , reliability_list_seq ,g2g, group_size_vector,node_to_group, group_of_node, nonempty_groups,H)

t21 = time.time()

print "time for decorraleted partitions" , t21 - t20

print "reliability" , reliability_list

missing_list = []

spurious_list = []

for i in range(k):
    for j in range(i):
        if (i,j) in adj_Matrix:
            spurious_list.append([(i, j) , reliability_list[i,j]])
        else :
            missing_list.append([(i, j) , reliability_list[i,j]])

spurious_list.sort(key = op.itemgetter(1)) 

missing_list.sort(key = op.itemgetter(1), reverse = True) 

f = open("spurious.dat", "w")
for pair in spurious_list:
    f.write(str(pair[0]) + ' ' + str(pair[1]) + "\n")
f.close()              

f = open("missing.dat", "w")
for pair in missing_list:
    f.write(str(pair[0]) + ' ' + str(pair[1]) + "\n")
f.close()  

