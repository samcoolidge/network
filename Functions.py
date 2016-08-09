# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 14:31:41 2016

@author: samuelcoolidge

"""
from MC_step_cython import *
import math
import numpy as np
from copy import copy,deepcopy

def compute_entropy(part, size):
    h = 0
    for i in xrange(size):
        i_size = np.sum(part[i])
        if i_size == 0:
            continue
        h += i_size * math.log(float(i_size) / size)
    return h


def compute_joint_entropy(part_1, part_2, size):
    h12 = 0
    for i in xrange(size):
        i_size = np.sum(part_1[i])
        if i_size == 0:
            continue
        for j in xrange(size):
            j_size = np.sum(part_2[j])
            if j_size == 0:
                continue
            s12 = np.sum(part_1[i] * part_2[j])
            if s12 > 0:
                # print(s12, size, i_size, j_size, float((s12 * size)) / (i_size * j_size))
                h12 += s12 * math.log(
                    float((s12 * size)) / (i_size * j_size))
    return h12

def mutual_information(part_1, part_2, size):
    h1 = compute_entropy(part_1, size)
    h2 = compute_entropy(part_2, size)
    h12 = compute_joint_entropy(part_1, part_2, size)
    return -2.0 * h12 / (h1 + h2)
            
def calculate_decay(k,x1,y1,x2,y2):
    
    def func(p,*data):
        a, b = p
        x1, y1, x2, y2 = data
        out = [y1 - a - (1. - a) * math.exp(-x1 / b)]
        out.append(y2 - a - (1. - a) * math.exp(-x2 / b))
        return out
    try:
        root = fsolve(func,[y1,math.sqrt(k)],args = (x1,y1,x2,y2))
    except OverflowError:
        return 'failed'

    return root[1]
    
def get_decorrelation_step(adj_Matrix,adj_Matrix_np,k):
    
    partition = np.zeros(shape=(k,k),dtype = np.int64)
    
    for x in xrange(k):
            partition[x,x]=1
            
    g2g = deepcopy(adj_Matrix_np)
                        
    node_to_group = deepcopy(adj_Matrix_np)
        
    group_of_node = np.arange(k, dtype = np.int64)
        
    nonempty_groups = set(xrange(k))  
    
    group_size_vector = group_size_vector = np.ones(k,dtype=np.int64)
    
    
    x2 = k/5
    x1 = x2 /4
    iteration_number = 10
    decay = []
    rep = 0
    H = [calculate_Hvalue(group_size_vector,nonempty_groups,g2g, k)]
    while rep < iteration_number:
        print "decor iteration" , rep
        part_ref = deepcopy(partition)
        for step in range(x2):
           factor_MC_step(partition,adj_Matrix_np,k,1,None,g2g,group_size_vector,node_to_group,group_of_node,nonempty_groups,H)
           if step == x1 :
               y1 = mutual_information(part_ref,partition,k)
        y2 = mutual_information(part_ref,partition,k)
        
        z = calculate_decay(k,x1,y1,x2,y2)
        if z == 'failed' or z < 0:
            continue
        else:
            decay.append(z)
            rep += 1
    

    result = sum(decay)
    norm = 0
    for rep in range(iteration_number):
        if decay[rep] > k: #abs(decay[rep] - mean_decay) / std_decay > 2:
            result -= decay[rep]
        else:
            norm +=1
    
            
    return result/norm
    
def get_sample(adj_Matrix,adj_Matrix_np,decor_step,k): 

    partition = np.eye(k, dtype = np.int)
    
    group_size_vector = np.ones(k,dtype=np.int)
    
    g2g = deepcopy(adj_Matrix_np)

    node_to_group = deepcopy(adj_Matrix_np)
    
    group_of_node = np.arange(k)
    
    nonempty_groups = set(xrange(k))     
    
    H_values = np.zeros(shape=(20,))
    
    H = [calculate_Hvalue(group_size_vector,nonempty_groups,g2g,k)]
    print "initial h value" , H[0]
    Hmean0 = 1e10
    Hstd0 = 1e-10
    equilibrated = 0
    while equilibrated < 5:
        for rep in xrange(20): 
            factor_MC_step(partition,adj_Matrix_np,k,decor_step, None, g2g,group_size_vector,node_to_group, group_of_node, nonempty_groups,H)
            H_values[rep] = H[0]
            print "H value while equilibrating" , H[0]
        Hmean1 = np.mean(H_values)
        Hstd1 = np.std(H_values)
        print "H mean" ,Hmean1 , "H std" , Hstd1
        if (Hmean0 - Hstd0/math.sqrt(20)) - (Hmean1 + Hstd1 / math.sqrt(20)) < 1e-10:
            equilibrated += 1
            print "equilibrated" + str(equilibrated) + "/5"
        else:
            print "not equilibrated"
            equilibrated = 0
            Hmean0 = Hmean1
            Hstd0 = Hstd1
    return partition, g2g, group_size_vector, node_to_group, group_of_node, nonempty_groups, H;
    
def add_sample(part, adj_Matrix,adj_Matrix_np, k, r, decor_step, reliability_list,g2g,group_size_vector,node_to_group, group_of_node,nonempty_groups,H):
    l = [part]
    reliability_copy = deepcopy(reliability_list)
    for rep in xrange(r):
       part_copy = deepcopy(part)
       factor_MC_step(part_copy,adj_Matrix_np,k,decor_step, reliability_copy, g2g,group_size_vector,node_to_group, group_of_node, nonempty_groups,H)
       l.append(part_copy)
       part = part_copy
    reliability_copy = reliability_copy / float(r)
    return l, reliability_copy