# network

##Installation

**Mac**

Clone Github

1) git clone http://github.com/samcoolidge/network 


Install Python Packages

1) sudo easy_install pip

2) pip install scipy

3) pip install joblib


*Linux*

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


##Usage

1) Network_Seq: evaluates the reliability of links using sequential code 2) Network_Parallel: evaluates the reliability of links using parallel code
*Arguments*
The program takes two arguments
arg1 = filenameThe file must be a list of positive integer node pairs representing links. Each node must be separated by spaces and each link must be separated by a line. Self-linking nodes are not allowed i.e. a link between node 1 and node 1. The format is shown belowN1  N2N3  N4arg2 = sample size (10,000 is recommended)
*Output*The program produces two output files, Missing.dat and Spurious.dat Spurious.dat lists all the reliabilities of links that are present in the network. An element of the list is written as (i,j), r where i and j are nodes and r is the reliability of the corresponding link. The list is sorted from low to high where a low score indicates that the link is likely spurious.Missing.dat lists all the reliabilities of links that are not present in the network. An element of the list is written as (i,j), r where i and j are nodes and r is the reliability of the corresponding link. The list is sorted from high to low where a high score indicates that the link is likely missing.
