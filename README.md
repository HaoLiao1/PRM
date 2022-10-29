# Title

PRM-IMM algorithm from our SIGMOD 2023 paper: Popularity Ratio Maximization: Surpassing Competitors through Influence Propagation.
For some reasons, part of the basic graph algorithm code in this project cannot be open sourced. 
So we compiled an executable file and provided the corresponding python script to run this program. You can run this program on windows 10.

## Install

We implement PRM-IMM algorithms in Visual C++, compiled in Visual Studio 2019, and run our tests on a computer with 3.6GHz Intel(R) Core(TM) i9-9900K CPU (16 cores), 128G memory, and Windows 10 professional (64 bits).

## Usage
Input
// read from file, there are also a example data in dm_real.txt.
	// Example input:          // Nodes starts from 1, and <edge_size>*2 lines follows
	//	15233 58891 ---------> #nodes #edges(undirected)
	//	1 2 2.380952e-002 2.500000e-001
	//	2 1 2.500000e-001 2.380952e-002
	//	3 4 7.692308e-002 5.000000e-001
	//	4 3 5.000000e-001 7.692308e-002
	//	...
	//	15232 15233 3.127006e-001 3.127006e-001
	//	15233 15232 3.127006e-001 3.127006e-001

### max_inf_origin.exe
The file "max_inf_origin.exe" is the main executable file for PRM-IMM(OINS) algorithm. It contains the PRM-IMM algorithm.

-t seeds_file <num_iter=10000> <seed_set_size = 50> <output_file=GC_spread.txt> <nthreads=1> <mode=0>: test influence spread with seeds \n
example: max_inf.exe -t < dm_real.txt

-rr : reverse influence maximization algorithms.
	   -rr5 <eps=0.1> <ell=1.0>	<k = 50> <mode = 1> <round = 10> <dp0 = 400> <dn0 = 10> <a = 10> (PRM-IMM) 
example: max_inf.exe -rr3o 0.1 1.0 10 1 10 400 10 50 < dm_real.txt

### max_inf6.exe
The file "max_inf6.exe" is the main executable file for PRM-IMM(NIOS) algorithm. It contains the PRM-IMM algorithm.

-t seeds_file <num_iter=10000> <seed_set_size = 50> <output_file=GC_spread.txt> <nthreads=1> <mode=0>: test influence spread with seeds \n
example: max_inf.exe -t < dm_real.txt

-rr : reverse influence maximization algorithms.
	   -rr3o <eps=0.1> <ell=1.0>	<k = 50> <mode = 1> <round = 10> <dp0 = 400> <dn0 = 10> <a = 10> (PRM-IMM) 
example: max_inf.exe -rr3o 0.1 1.0 10 1 10 400 10 50 < dm_real.txt


## Contributing


