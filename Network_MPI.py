# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 22:29:02 2016

@author: samuelcoolidge
"""

from Functions import *
from MC_step_cython import *
# from multiprocessing import cpu_count
import numpy as np
import operator as op
import time
import sys
# from joblib import Parallel,delayed
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def parse(file_name):
    adj_Matrix = set()
    maxim = 0
    for line in open(file_name, 'r'):
        edge = map(int, line.split())
        if edge[0] == edge[1]:
            continue
        adj_Matrix.add((edge[0]-1,edge[1]-1))
        adj_Matrix.add((edge[1]-1,edge[0]-1))
        if max(edge) > maxim:
            maxim = max(edge)
    return adj_Matrix, maxim
    
def main(file_name, reps):
    # partition_number = cpu_count() #how many new partitions should we make
    adj_Matrix, k = parse(file_name)
    return reliability(adj_Matrix, k, reps)

   
def reliability(adj_Matrix, k, reps):
    
    
    if rank==0:
        
        time1 = time.time()
        #initialize reliability list
        reliability_list_seq = np.zeros((k,k))
           
        for i in range(k):
            if (i,i) in adj_Matrix:
                print "Error: Network Contains Self-Linking Node"
                exit()
        
        #create numpy adjacency matrix
        adj_Matrix_np = np.zeros(shape = (k,k), dtype = int)
        for (i,j) in adj_Matrix:
            adj_Matrix_np[i][j] = 1
        
        
        decor_step = get_decorrelation_step(adj_Matrix,adj_Matrix_np,k)
        
        print "decorellation step" , decor_step
        
    
    else :
        adj_Matrix_np, decor_step, reliability_list_seq = [None for _ in range(3)]
    
    t_10 = time.time()  
    
    #broadcast data to all nodes
    adj_Matrix_np = comm.bcast(adj_Matrix_np, root = 0)
    decor_step = comm.bcast(decor_step, root = 0)
    reliability_list_seq = comm.bcast(reliability_list_seq, root = 0)
        
    #get equilibrated starting partitions
        
    start_parts, g2g, group_size_vector,node_to_group, group_of_node, nonempty_groups, H = get_sample(adj_Matrix, adj_Matrix_np, decor_step,k)
        
    t_11 = time.time()
        
    print "time for starting partitions", t_11-t_10, rank
    
    #list of variables : parameters - adj_Matrix, k, reps
    #other inputs to get_sample - adj_Matrix_np, decor_step
    #outputs of get_sample - start_parts, g2g, group_size_vector,node_to_group, group_of_node, nonempty_groups, H
    #misc - reliability_list_seq
    
    #print "REPS" , reps
    
    t20 = time.time()
    sample_dirty, reliability_dirty = add_sample(
        start_parts, adj_Matrix, adj_Matrix_np, k, reps, decor_step, reliability_list_seq, 
        g2g, group_size_vector,node_to_group, group_of_node, nonempty_groups,H 
    )
    t21 = time.time()
    
    #print "reliability_dirty" , reliability_dirty 
    
    #sample_dirty, reliability_dirty = add_sample(start_parts,adj_Matrix,adj_Matrix_np, k, reps, decor_step, reliability_list_seq ,
    #    g2g, group_size_vector,node_to_group, group_of_node, nonempty_groups,H)
    
    print "time for decor parts" , t21 - t20
    

    sample_dirty = comm.gather(sample_dirty, root=0)
    reliability_dirty = comm.gather(reliability_dirty, root=0)
    
    
    
    #print "shape sample" ,np.shape(sample_dirty)
    #print "shape reliability" ,np.shape(reliability_dirty)
    
    if rank == 0:
        reliability_list = np.sum(reliability_dirty, axis=0) / len(reliability_dirty)
        
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
        
        time2 = time.time()
        
        print "total run time" , time2-time1
        
        return missing_list, spurious_list
    
if __name__ == '__main__':
    time0 = time.time()
    file_name = sys.argv[1]
    reps = int(sys.argv[2]) / size - 1
    # print "CPU_COUNT" , cpu_count()
    main(file_name, reps)
