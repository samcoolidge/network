
import numpy as np
cimport numpy as np
import math
import random
cimport cython


import line_profiler

from libc.math cimport log, lgamma

def log_c(int n):
    return log(n)


def log_ncr_c(int n, int r, dict d={}):
    cdef double result
    if (n, r) in d:
        return d[(n, r)]
    #print n , r 
    result = lgamma(n + 1) - lgamma(r + 1) - lgamma(n - r + 1)
    d[(n, r)] = result
    return result
    
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
    
def decision(probability):
    return random.random() < probability;
    
def calculate_Rvalue(int a,int b,int len_a,int len_b):
    cdef int r     
    if a!=b:
        r=len_a*len_b;
    else:
        r=(len_a * (len_b - 1))/2;
    return r

def calculate_Lvalue(int a,int b,np.ndarray[np.int64_t, ndim=2] adj_Matrix_np,np.ndarray[np.int64_t, ndim=2]grp_Matrix,int k):
    cdef int l, dot 
    cdef np.ndarray[np.int64_t, ndim = 1] list_a 
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


def calculate_Hvalue(np.ndarray[np.int64_t, ndim=1] group_size_vector,set nonempty_groups, np.ndarray[np.int64_t, ndim=2] g2g, int k):       
    cdef list s 
    cdef float h
    cdef int l, r
    
    s = []
    h = 0
    for b in nonempty_groups:
        s.append(b)
        for a in s:
           # if len(grp_Matrix[a]) != 0 and len(grp_Matrix[b]) !=0:
           len_a = group_size_vector[a]
           len_b = group_size_vector[b]
           l=g2g[a,b]
           r=calculate_Rvalue(a,b,len_a,len_b)
           h += log_c(r+1) + log_ncr_c(r,l)
            
    return h;



def factor_MC_step(np.ndarray[np.int64_t, ndim=2] partition,np.ndarray[np.int64_t, ndim=2]adj_Matrix_np,int k,int factor,np.ndarray[np.float64_t, ndim=2] reliability_list,np.ndarray[np.int64_t, ndim=2] g2g, np.ndarray[np.int64_t, ndim=1] group_size_vector,np.ndarray[np.int64_t, ndim=2] node_to_group, np.ndarray[np.int64_t, ndim=1] group_of_node, set nonempty_groups, list H):
    
    cdef int z, x, y, a, b, g, lx, ly, links,links_x, links_xtoz, links_y, links_ytoz, links_xtoy, i, j, l, r, m
    cdef np.ndarray[np.int64_t, ndim = 1] links_grp_to_z, adj_z
    # cdef float current_H_value, swap_contribution, swap_contribution_1 , swap_contribution_2, swap_contribution_3, swap_contribution_4, swap_contribution_5, prob
    cdef float swap_contribution, swap_contribution_1 , swap_contribution_2, swap_contribution_3, swap_contribution_4, swap_contribution_5, prob
    
    # current_H_value = H[0]
    
    for steps in xrange(int(math.floor(factor*k))):
            z = random.randint(0, k - 1)
            x = group_of_node[z]
            y = random.randint(0, k - 2)
            if y >= x:
                y += 1
            
            a = group_size_vector[x]
            b = group_size_vector[y]
            swap_contribution_1 = 0
            swap_contribution_2 = 0
            
            nonempty_groups.remove(x)
            nonempty_groups.discard(y)

            for group in nonempty_groups:
                g = group_size_vector[group]

                                       
                swap_contribution_1 += log_c(g*(a-1)+1)+log_c(g*(b+1)+1)-log_c(g*a+1)-log_c(g*b+1)

                lx = g2g[group,x]
                ly = g2g[group,y]


                links = node_to_group[group, z]
                
                #if lx < links:
                #    print g2g[group], group_size_vector[group], np.sum(g2g[group]), node_to_group[group]
                 
                #print g*(a-1),int(lx - links),lx, links
                swap_contribution_2 += log_ncr_c(g*(a-1),int(lx - links)) + log_ncr_c(g*(b+1),int(ly + links))- log_ncr_c(g*a,lx)- log_ncr_c(g*b,ly)  
                #swap_contribution_2 += log_ncr(g*(a-1),int(lx - links)) + log_ncr(g*(b+1),int(ly + links))- log_ncr(g*a,lx)- log_ncr(g*b,ly)
                
                
            
            links_x = g2g[x,x]
            links_xtoz = np.dot(partition[x],adj_Matrix_np[z])
            
            swap_contribution_3 = log_c((((a-1)*(a-2))/2)+1)-log_c((a*(a-1)/2)+1)+ log_ncr_c((((a-1)*(a-2))/2),int(links_x - links_xtoz)) - log_ncr_c((a*(a-1))/2,links_x)
            
            links_y = g2g[y,y]
            links_ytoz = np.dot(partition[y],adj_Matrix_np[z])
        
            swap_contribution_4 = math.log((((b+1)*b)/2)+1)-math.log((b*(b-1)/2)+1)+ log_ncr_c(((b+1)*b/2),int(links_y + links_ytoz)) - log_ncr_c((b*(b-1))/2,links_y)
            
            links_xtoy = g2g[x,y]
            
            swap_contribution_5 = log_c((a-1)*(b+1)+1)-log_c(a*b+1)+log_ncr_c(((a-1)*(b+1)),int(links_xtoy - links_ytoz + links_xtoz))-log_ncr_c(a*b,links_xtoy)
            swap_contribution = swap_contribution_1 + swap_contribution_2+ swap_contribution_3 + swap_contribution_4 +swap_contribution_5
            
            
            if (swap_contribution < 0) or (decision(math.exp(-swap_contribution))==True):
                partition[y,z]=1;
                partition[x,z]=0;
                H[0] += swap_contribution
                

                links_grp_to_z = node_to_group[:,z]
                g2g[x] -= links_grp_to_z 
                g2g[:,x] = g2g[x]
                g2g[y] += links_grp_to_z 
                g2g[:,y] = g2g[y]
                        
                adj_z = adj_Matrix_np[z]
                node_to_group[x] -= adj_z
                node_to_group[y] += adj_z
                        
                group_size_vector[x] -= 1
                group_size_vector[y] += 1
                
                group_of_node[z] = y
                
            if group_size_vector[x] != 0:
                nonempty_groups.add(x)
            if group_size_vector[y] != 0:
                nonempty_groups.add(y)
           
    
    # print "# accpeted moves" ,accepted_moves
    if reliability_list is None:
        return partition

    for i in xrange(k):
        if group_size_vector[i] == 0:
            continue
        for j in xrange(k):
            if group_size_vector[j] == 0:
                 continue
            l = g2g[i,j]
            if i==j:
                l /= 2
                r = (group_size_vector[i] * (group_size_vector[j] - 1)) / 2
            else:
                r = group_size_vector[i] * group_size_vector[j]
            
            s = (l+1) / float(r+2)
            

            for k in np.flatnonzero(partition[i]):
                for m in np.flatnonzero(partition[j]):
                    if k > m:
                        reliability_list[k][m] += s
           
    return partition
    
