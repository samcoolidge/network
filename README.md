# Network

A parallel version of seeslab/rgraph.

##Installation

**Mac**

Clone Github

1) git clone http://github.com/samcoolidge/network 


Install Python Packages

1) sudo easy_install pip

2) pip install scipy

3) pip install joblib


**Linux**

Clone Github

1) git clone http://github.com/samcoolidge/network


Install Python Packages

1) sudo apt-get install python-numpy python-scipy

2) apt-get install python-pip

3) pip install joblib


Install Cython

1) download tarball from https://pypi.python.org/pypi/Cython/

2) tar –xzvf Cython-0.24.1.tar.gz

3) cd cython-0.24.1

4) python setup.py install

5) pip install cython

6) python setup.py build_ext —inplace

**HPC**

If you want to run on multiple nodes of a High Performance Computing Cluster, please install Open MPI and mpi4py

##Usage

1) Network_Seq: evaluates the reliability of links using sequential code 

2) Network_Parallel: evaluates the reliability of links using parallel code

3) Network_MPI: evaluates the reliability of links using parallel code and is compatible with a HPC Cluster running on multiple nodes. 

For Network_MPI run the program using the following command

mpiexec -np numberofnodes python Network_MPI.py filename samplesize

Replace numberofnodes, filename and samplesize with the desired values

*Arguments*

The program takes two arguments

arg1 = filename

The file must be a list of positive integer node pairs representing links. Each node must be separated by spaces and each link must be separated by a line. Self-linking nodes are not allowed i.e. a link between node 1 and node 1. The format is shown below

N1 N2

N3 N4

arg2 = sample size (10,000 is recommended)


*Output*

The program produces two output files, Missing.dat and Spurious.dat 

Spurious.dat lists all the reliabilities of links that are present in the network. An element of the list is written as (i,j), r where i and j are nodes and r is the reliability of the corresponding link. The list is sorted from low to high where a low score indicates that the link is likely spurious.

Missing.dat lists all the reliabilities of links that are not present in the network. An element of the list is written as (i,j), r where i and j are nodes and r is the reliability of the corresponding link. The list is sorted from high to low where a high score indicates that the link is likely missing.

##License

The program is licensed under the MIT licence.

## Speed

The program is written in Python and Cython and, for large enough datasets (about 1000 nodes and above), is only 2-3 times slower than the non-parallel code it is based on. This means that even running it on a typical computer, you will generally achieve a speed benefit over the original code.



