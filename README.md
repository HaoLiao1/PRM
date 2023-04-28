# Title

PRM-IMM algorithm from our SIGMOD 2023 paper: Popularity Ratio Maximization: Surpassing Competitors through Influence Propagation.

## Environment

We implement PRM-IMM algorithms in Visual C++, compiled in Visual Studio 2019, and run our tests on a computer with 3.6GHz Intel(R) Core(TM) i9-9900K CPU (16 cores), 128G memory, and Windows 10 professional (64 bits).
This code can also be compiled and run on Linux. But the omp setting is hard to use on linux. The OpenMP is an API specification for parallel programming. So Windows is recommended.

## Code Structure
We have two folders for the PRM-OINS and NIOS setting. You need to put the code from these two folders into two projects and compile.
```
├── PRM_OINS
│   ├── code
│   │   ├── mi_command_line.cpp // process command line parameter
│   │	├── reverse_general_cascade.cpp // generate rr set
│   │	├── rr_infl.cpp // main algorithm
│   │	├── event_timer.cpp // record running time
│   │	└── ....			
│   └── data
│       └── dm_wc.txt	
├── PRM_NIOS
│   ├── code
│   │   ├── mi_command_line.cpp // process command line parameter
│   │	├── reverse_general_cascade.cpp // generate rr set
│   │	├── rr_infl.cpp // main algorithm
│   │	├── event_timer.cpp // record running time
│   │	├── greedy.cpp // greedy algorithm
│   │	├── simulate.cpp // simulate the PA-IC NIOS setting
│   │	└── ....			
│   └── data
│       └── dm_wc.txt	
└── readme.md
```
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

### PRM_OINS folder
The file "PRM_OINS.exe" is the main executable file for PRM-IMM(OINS) algorithm. It contains the PRM-IMM algorithm.

	-h: print the help
	-rr5 <eps=0.1> <ell=1.0> <k = 50> <mode = 1> <round = 10> <dp0 = 400> <dn0 = 10> <a = 50> (PRM-IMM).

example: PRM_OINS.exe -rr3o 0.1 1.0 10 1 10 400 10 50 < dm_real.txt

### PRM_NIOS folder
The file "PRM_NIOS.exe" is the main executable file for PRM-IMM(NIOS) algorithm. It contains the PRM-IMM algorithm.

	-g : greedy algorithm for PRM NIOS and OINS setting.
	-tp simulate the process of PA-IC in NIOS setting and evaluate the result of different algorithm.
	-rr5 <eps=0.1> <ell=1.0> <k = 50> <mode = 1> <round = 10> <dp0 = 400> <dn0 = 10> <a = 50> (PRM-IMM).

example: PRM_NIOS.exe -rr5o 0.1 1 10 1 10 400 10 50 < dm_real.txt > out.txt

### MonteCarlo-Test.py
The simulation of the process of PA-IC in OINS is easier than NIOS. So we write an python script to simulate.

### PRM__appendix_proof.pdf

This file contains the proof for all the lemmas in the paper. And we will put full version on arxiv.


