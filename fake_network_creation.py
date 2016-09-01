# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 13:06:19 2016

@author: dan
"""

# Note: this file will add about 15 test networks of sizes up to 30MB
# (although most are much smaller) to the folder you run it in.

import random

def random_group_size(average):
    return int(max([random.gauss(average, average / 3), 0]))

def create_group_sizes(n, number_of_groups, density):
    total = 0
    groups = []
    average = n / number_of_groups * density
    for i in range(0, number_of_groups):
        group_size = random_group_size(average)
        if group_size == 0:
            continue
        total += group_size
        groups.append(group_size)
    return groups


def break_group(group_size):
    part_size = int(random.gauss(group_size / 2, group_size / 10))
    if part_size < 0:
        part_size = 0
    if part_size > group_size:
        part_size = group_size
    return (part_size, group_size - part_size, group_size)


def get_group_splits(n, number_of_groups, density):
    return map(break_group, create_group_sizes(n, number_of_groups, density))

def add_group(s, nodes, i):
    group = random.sample(nodes, i[2])
    part_1 = set(random.sample(group, i[0]))
    part_2 = set(group) - part_1
    for i in part_1:
        for j in part_2:
            s.add((i, j))
            s.add((j, i))

def as_adj_matrix(group_splits, n):
    s = set()
    nodes = range(1, n + 1)
    for i in group_splits:
        add_group(s, nodes, i)
    return s

'''
# Use this function if you want more memorable but longer file names.
def make_file_name(n, number_of_groups, density):
    first = ['dolphin', 'hamster', 'karate', 'neural', 'celegans', 'drug']
    second = ['network', 'colony', 'friendships', 'airports',
              'connections', 'interactions']
    return '_'.join([random.choice(first), random.choice(second),
                     str(n), str(number_of_groups), str(density)])
'''

def make_file_name(n, number_of_groups, density):
    return '_'.join(['network', str(n), str(number_of_groups), str(density)])

def adj_matrix(n, number_of_groups, density):
    group_splits = get_group_splits(n, number_of_groups, density)
    return as_adj_matrix(group_splits, n)

file_extension = '.dat'

def write_to_file(file_name, s):
    with open(file_name + file_extension, 'w') as f:
        for i in s:
            f.write(' '.join(map(str, i)) + '\n')

def main(n, number_of_groups, density):
    my_adj_matrix = adj_matrix(n, number_of_groups, density)
    file_name = make_file_name(n, number_of_groups, density)
    write_to_file(file_name, my_adj_matrix)

def make_tests():
    i = 20
    while i < 30000:
        main(int(i), int(i ** .5), 1)
        i *= 1.5

make_tests()