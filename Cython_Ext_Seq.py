# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:37:18 2016

@author: samuelcoolidge
"""

# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 13:31:37 2016

@author: samuelcoolidge
"""

# -*- coding: utf-8 -*-

"""
Created on Fri May 27 12:12:37 2016

@author: samuelcoolidge
"""

# import pandas as pd

from multiprocessing import Process, Queue, cpu_count
q = Queue()
import numpy as np
import line_profiler
from  MC_step_cython43 import *


from scipy.optimize import fsolve
#from scipy import sparse

import collections;
import math;
import operator as op
import random

# import matplotlib.pyplot as plt
from copy import copy, deepcopy
from decimal import *

import cPickle as pickle
import json
import time
from joblib import Parallel,delayed
#import sympy as sp

def log_ncr(n, r, d={}):
    if (n, r) in d:
        return d[(n, r)]
    #print n , r 
    result = math.lgamma(n + 1) - math.lgamma(r + 1) - math.lgamma(n - r + 1)
    d[(n, r)] = result
    return result

    
def factorial(n, d={}):
    if n in d:
        return d[n]
    result = math.factorial(n)
    d[n] = result
    return result


def ncr(n, r, d={}):
    if (n, r) in d:
        return d[(n, r)]
    if r == 0:
        return 1
    if n < r or n < 0 or r < 0:
        print n, r
        raise Exception()
    result = factorial(n) / factorial(r) / factorial(n - r)
    d[(n, r)] = result
    return result

'''
def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom


import math
import time
import random
# I tried "import operator" and the speed difference was insignificant.
import operator as op

def factorial(n, d={}):
    if n in d:
        return d[n]
    result = math.factorial(n)
    d[n] = result
    return result


def ncr(n, r, d={}):
    if (n, r) in d:
        return d[(n, r)]
    result = factorial(n) / (factorial(r) * factorial(n - r))
    d[(n, r)] = result
    return result

def old_ncr(n, r):
    return math.factorial(n) / (math.factorial(r) * math.factorial(n - r))

def current_ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

def time_ncr(f, rands):
    t0 = time.time()
    for item in rands:
        f(item[0] + item[1], item[0])
    t1 = time.time()
    return t1 - t0

def get_rands(r, n):
    return [(random.randint(0, r), random.randint(0, r)) for _ in range(n)]

# Two tests. The first number is hopefully lower than the second and third.

time_ncr(ncr, get_rands(5, 10 ** 4))
time_ncr(old_ncr, get_rands(5, 10 ** 4))
time_ncr(current_ncr, get_rands(5, 10 ** 4))

time_ncr(ncr, get_rands(100, 10 ** 4))
time_ncr(old_ncr, get_rands(100, 10 ** 4))
time_ncr(current_ncr, get_rands(100, 10 ** 4))

"""
results

>>> time_ncr(ncr, get_rands(5, 10 ** 4))
0.005514860153198242
>>> time_ncr(old_ncr, get_rands(5, 10 ** 4))
0.007944822311401367
>>> time_ncr(current_ncr, get_rands(5, 10 ** 4))
0.01625990867614746
>>> 
>>> time_ncr(ncr, get_rands(100, 10 ** 4))
0.020054101943969727
>>> time_ncr(old_ncr, get_rands(100, 10 ** 4))
0.16937613487243652
>>> time_ncr(current_ncr, get_rands(100, 10 ** 4))
0.09639692306518555
"""
'''
    

def calculate_Rvalue(a,b,len_a,len_b):
    if a!=b:
        r=len_a*len_b;
    else:
        r=(len_a * (len_b - 1))/2;
    return r


def calculate_Lvalue_orig(list_a,list_b,adj_Matrix):
    l=0;
    for i in list_a:
        for j in list_b:
            if (i, j) in adj_Matrix:
                l += 1
    
    
    
    if list_a == list_b:
        return int(l/2);
    else:
        return int(l);    
        

def calculate_Lvalue(a,b,adj_Matrix_np,grp_Matrix,k):
    
    #t10 = time.time()
    l = 0
    
    #t1 = time.time()
    list_a = np.flatnonzero(grp_Matrix[a])
    for i in list_a:
        dot = np.dot(grp_Matrix[b],adj_Matrix_np[i])
        l += dot
                
    if a == b:
        return int(l/2);
    else:
        return int(l); 

#profile
'''
def calculate_Hvalue(adj_Matrix_np,grp_Matrix,group_size_vector,k):       
    
    h=0;
    for b in range(k):
        for a in range(b+1):
           # if len(grp_Matrix[a]) != 0 and len(grp_Matrix[b]) !=0:
           len_a = group_size_vector[a]
           len_b = group_size_vector[b]
           l=calculate_Lvalue(a,b,adj_Matrix_np,grp_Matrix,k)
           r=calculate_Rvalue(a,b,len_a,len_b)
           h = h + math.log(r+1) + math.log(ncr(r,l));
            
            
    return h;

    
def calculate_Hvalue(adj_Matrix_np,grp_Matrix,group_size_vector,nonempty_groups,k):       
    s = []
    h = 0
    for b in nonempty_groups:
        s.append(b)
        for a in s:
           # if len(grp_Matrix[a]) != 0 and len(grp_Matrix[b]) !=0:
           len_a = group_size_vector[a]
           len_b = group_size_vector[b]
           l=calculate_Lvalue(a,b,adj_Matrix_np,grp_Matrix,k)
           r=calculate_Rvalue(a,b,len_a,len_b)
           h = h + math.log(r+1) + log_ncr(r,l);
            
            
    return h;
'''    
    
def decision(probability):
    return random.random() < probability;
    
    
def max(l):
    if l[0]>l[1]:
        return l[0]
    else:
        return l[1]
    
def find_group(partition,node):
    for x in range(len(partition)):
        if partition[x].count(node) == 1:
            return x
            break
        else:
            continue   


#dont use this
def link_reliability(node_1,node_2,sample,adj_Matrix,k):
   global time_find
   global time_rel
   global time_l
   t10 = time.time()
   s = 0         
   for part in sample:
       t20 = time.time()
       l1 = part[find_group(part, node_1)]
       l2 = part[find_group(part, node_2)]
       t21 = time.time() 
       time_find += t21 - t20
       t30 = time.time()
       l = calculate_Lvalue(l1, l2, adj_Matrix_np,k)
       t31 = time.time()
       time_l += t31 - t30
       r = calculate_Rvalue(l1, l2)
       s = s + (l+1)/float((r+2))
       
   l_r = s/float(len(sample))
   t11 = time.time()
   time_rel += t11 - t10
   return l_r

def remove_links(adj_Matrix,p,k):
    adj_Matrix_F = deepcopy(adj_Matrix)
    link_list = []
    for pair in adj_Matrix:
        link_list.append(pair)
        
    remove = int(math.floor(p*len(adj_Matrix)))


    for x in range(remove):

        link = random.choice(link_list)
        link_list.remove(link)
        if (link[0], link[1]) in adj_Matrix_F:
            adj_Matrix_F.remove((link[0], link[1]))
        if (link[1], link[0]) in adj_Matrix_F:
            adj_Matrix_F.remove((link[1], link[0]))
    return adj_Matrix_F
        
    
             
def add_links(adj_Matrix,p,k):
    adj_Matrix_F = deepcopy(adj_Matrix)
    link_list = []
    
    for i in range(k):
        for j in range(i):
            if adj_Matrix[i][j] == 0:
                link_list.append([i,j])
    add = int(math.floor(p*len(link_list)))
    for x in range(add):
        link = random.choice(link_list)
        link_list.remove(link)
        adj_Matrix_F.add((link[0], link[1]))
        adj_Matrix_F.add((link[1], link[0]))
    
    return adj_Matrix_F
    
def modify(adj_Matrix,p,k):
    adj_Matrix_M = remove_links(adj_Matrix,p,k)
    removed = len(adj_Matrix) - len(adj_Matrix_M)
    count = 0
    while True:
        i = random.randint(0,k-1)
        j = random.randint(0,k-1)
        if (i, j) not in adj_Matrix_M and (i, j) not in adj_Matrix and i!=j:
            adj_Matrix_M.add((i, j))
            adj_Matrix_M.add((j, i))
            count = count + 1
        if count == removed:
            break
    return adj_Matrix_M
    
            

def link_checker(adj_Matrix,adj_Matrix_modified,k):
    spurious_count = 0
    missing_count = 0
    links_missing = []
    links_spurious = []
    for i in range(k):
        for j in range(i):
            if (i, j) in adj_Matrix and (i, j) not in adj_Matrix_modified:
                missing_count = missing_count + 1
                links_missing.append([[i,j]])
            if (i, j) not in adj_Matrix and (i, j) in adj_Matrix_modified:
                spurious_count = spurious_count + 1
                links_spurious.append([i,j])
    count = [missing_count,spurious_count,links_missing,links_spurious]
    return [missing_count,spurious_count]
      
    
def test_rel_spurious(adj_Matrix_added,adj_Matrix_T,rel,k):
    Positives_T = []
    Positives_F = []
    
    #print rel, rel.shape

    for a in range(k):
        for b in range(a):
            link_rel = rel[a,b]
            if (a, b) in adj_Matrix_added:
               if (a, b) in adj_Matrix_T:
                   Positives_T.append(link_rel)
               else:
                   Positives_F.append(link_rel)  
    prob_list_less = []
    prob_list_equal = []
    prob_list_greater = []
    
    for i in range(len(Positives_F)):
        num_greater = 0
        num_equal = 0
        num_less = 0
        l = len(Positives_T)
        for j in range(l):
            if Positives_F[i]< Positives_T[j]:
               num_less = num_less+1
            elif Positives_F[i] == Positives_T[j] :
                 num_equal = num_equal+1
            elif Positives_F[i]>Positives_T[j]:
                 num_greater = num_greater+1
        prob_greater = num_greater/float(l)
        prob_list_greater.append(prob_greater)
        prob_equal = num_equal/float(l)
        prob_list_equal.append(prob_equal)
        prob_less = num_less/float(l)
        prob_list_less.append(prob_less)   
        mean_list = [np.mean(prob_list_less),np.mean(prob_list_equal),np.mean(prob_list_greater)]
         
    return mean_list
    
def test_rel_missing(adj_Matrix_removed,adj_Matrix_T,rel,k):
    Negatives_T = []
    Negatives_F = []

    for a in range(k):
        for b in range(a):
            link_rel = rel[a,b]
            if (a, b) not in adj_Matrix_removed:
                if (a, b) in adj_Matrix_T:
                    Negatives_F.append(link_rel)
                else:
                    Negatives_T.append(link_rel)
    
    prob_list_greater = []
    prob_list_equal = []
    prob_list_less = []
    
    for i in range(len(Negatives_F)):
        num_greater = 0
        num_equal = 0
        num_less = 0
        l = len(Negatives_T)
        for j in range(l):
            if Negatives_F[i]>Negatives_T[j]:
               num_greater = num_greater+1
            elif Negatives_F[i] == Negatives_T[j]:
                 num_equal = num_equal+1
            elif Negatives_F[i]<Negatives_T[j]:
                 num_less = num_less+1
        prob_greater = num_greater/float(l)
        prob_list_greater.append(prob_greater)
        prob_equal = num_equal/float(l)
        prob_list_equal.append(prob_equal)
        prob_less = num_less/float(l)
        prob_list_less.append(prob_less)   
        mean_list = [np.mean(prob_list_greater),np.mean(prob_list_equal),np.mean(prob_list_less)] 
        
    return mean_list

   
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
       #profile = line_profiler.LineProfiler(factor_MC_step)
       #profile.runcall(factor_MC_step, part_copy,adj_Matrix_np,k,decor_step, reliability_copy, g2g,group_size_vector,node_to_group, group_of_node, nonempty_groups,H)
       #profile.print_stats()
       l.append(part_copy)
       part = part_copy
    reliability_copy = reliability_copy / float(r)
    return l, reliability_copy

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

                   
    #mean_decay = np.mean(decay)
    #std_decay = np.std(decay)
    

    result = sum(decay)
    norm = 0
    for rep in range(iteration_number):
        if decay[rep] > k: #abs(decay[rep] - mean_decay) / std_decay > 2:
            result -= decay[rep]
        else:
            norm +=1
    
            
    return result/norm        

#@profile        
def main():
    file_name = 'blogs0.2_1224'
    #file_name = 'eu2_net'
    reps = 1
    make_new_mods = True #do we want to make a new modified adjacency matrix
    p = .2 #what percent of mods should we make
    make_new_samples = True #do we want to make new partitions
    partition_number = cpu_count() #how many new partitions should we make
    t0 = time.time()
    listHValues = []
    
    maxim = 0;    
    for line in open(file_name + ".dat", 'r'):
        #item = line.rstrip()
        edge = map(int, line.split())
        if max(edge) > maxim:
            maxim = max(edge)
        #print y

    k = maxim;   
    
    adj_Matrix = set();
    
    #reliability_list_all = np.zeros((cpu_count(), k, k))
    
    reliability_list_seq = np.zeros((k,k))
    
    #print "starting rel list" , reliability_list

       
    for line in open(file_name + ".dat", 'r'):
        #item = line.rstrip()
        edge = map(int, line.split())
        adj_Matrix.add((edge[0]-1,edge[1]-1))
        adj_Matrix.add((edge[1]-1,edge[0]-1))
        

    #in this section of code you can decide your modified matrix to test link reliability
    
    '''    
    #just use these lines for ecoli
    print len(adj_Matrix)
    adj_copy = deepcopy(adj_Matrix)  
    for link in adj_Matrix:
        if link[0] > 7086 or link[1]>7086:
            adj_copy.remove(link)
    adj_Matrix = adj_copy
    k = 1099
    #use above lines for ecoli
    ''' 
           
    adj_Matrix_T = deepcopy(adj_Matrix)
    

    p_string = str(p)
    matrix_file = "adj_Matrix_" + file_name + "_modify_" + p_string + ".p"
    if make_new_mods:  
        adj_Matrix = modify(adj_Matrix_T, p, k)
        output = open(matrix_file,'wb')
        pickle.dump(adj_Matrix,output)
    else:
        adj_Matrix = pickle.load(open(matrix_file,'rb')) 
       
    adj_Matrix_O = deepcopy(adj_Matrix)
    

    
    adj_Matrix_np = np.zeros(shape = (k,k), dtype = int)
    for (i,j) in adj_Matrix:
        adj_Matrix_np[i][j] = 1
    
    t_10 = time.time()
    
    print "got to here"
    
    decor_step = get_decorrelation_step(adj_Matrix,adj_Matrix_np,k)
    print "decorellation step" , decor_step
    
    #start_parts = Parallel(n_jobs=4)(delayed(sqrt)(i ** 2) for i in range(10))
    #get_sample(adj_Matrix,adj_Matrix_np,1,k)
    #start_parts = Parallel(n_jobs=cpu_count())(delayed(get_sample)(adj_Matrix,adj_Matrix_np, 1, k) for i in range(cpu_count()))
    
    start_parts, g2g, group_size_vector,node_to_group, group_of_node, nonempty_groups, H = get_sample(adj_Matrix,adj_Matrix_np,decor_step, k)  

    t_11 = time.time()
    
    print "time for starting partitions" , t_11-t_10
    
    t20 = time.time()
    sample,reliability_list = add_sample(start_parts,adj_Matrix,adj_Matrix_np, k, reps, decor_step , reliability_list_seq ,g2g, group_size_vector,node_to_group, group_of_node, nonempty_groups,H)
    t21 = time.time()
    
    print "time for decor parts" , t21 - t20
    
    
    #print start_parts[0]
    
    print "reliability" , reliability_list
    
    
    time100 = time.time()
    test_prob = test_rel_spurious(adj_Matrix,adj_Matrix_T,reliability_list,k)
    print "probabiliy that a false positive has (lower, equal, greater) reliability than a true positive"
    print test_prob
    print
        
    test_prob = test_rel_missing(adj_Matrix, adj_Matrix_T,reliability_list,k)
    print "probability that a false negative has  (greater, equal, lower) reliability than a true negative"
    print test_prob
    print
    time101 = time.time()

    
    
if __name__ == "__main__":
    main()

