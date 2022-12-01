#include "mi_command_line.h"
#include "rr_infl.h"

using namespace std;



std::string MICommandLine::Help()
{
	std::string help =
		"[Help]\n"
		MI_SOFTWARE_META "\n"
		"\n"
		"-h: print the help \n"
		"-g : greedy algorithm for IM and PRM\n"
		"-tp simulate the process of PA-IC in NIOS setting and evaluate the result of different algorithm. \n"
		"-t seeds_file <num_iter=10000> <seed_set_size = 50> <output_file=GC_spread.txt> <nthreads=1> <mode=0>: test influence spread with seeds \n"
		"-rr : reverse influence maximization algorithms.\n"
		"	   -rr5 <eps=0.1> <ell=1.0>	<k = 50> <mode = 1> <round = 10> <dp0 = 400> <dn0 = 10> <a = 10> (PRM-IMM) \n"
		"\n"
		"example: PRM_OINS.exe -rr5o 0.1 1 10 1 10 400 10 50 < dm_real.txt > out.txt \n"
	;

	// read from file
	// Example input:          // Nodes starts from 1, and <edge_size>*2 lines follows
	//	15233 58891 ---------> #nodes #edges(undirected)
	//	1 2 2.380952e-002 2.500000e-001
	//	2 1 2.500000e-001 2.380952e-002
	//	3 4 7.692308e-002 5.000000e-001
	//	4 3 5.000000e-001 7.692308e-002
	//	бн
	//	15232 15233 3.127006e-001 3.127006e-001
	//	15233 15232 3.127006e-001 3.127006e-001
	//
	//  printf("example: python gen_graph.py <parm1> <parm2> ... | max_influence -p 20 2000 \n"); // read from other program
	
	return help;
}

int MICommandLine::Main(int argc, char * argv[])
{
	// create an arguments vector to include all strings
	std::vector<std::string> argVec;
	for (int i = 0; i < argc; i++) {
		std::string param = argv[i];
		argVec.push_back(param);
	}

	return Main(argc, argVec);
}

int MICommandLine::Main(int argc, std::vector<std::string>& argv)
{
	srand((unsigned)time(NULL));
	
	if (argc <= 1) {
		std::cout << Help() << std::endl;
		return 0;
	}

	std::string arg1 = argv[1]; // the switch string like "-abc"
	std::string s;

	// for -h, print help
	s = "-h";
	if (s.compare(arg1) == 0) {
		std::cout << Help() << std::endl;
		return 0;
	}

	// the following contains switches for algorithms
	system("mkdir tmp");
	system("cd tmp");

	// create empty _running_.log to indicate running
	system("del /Q _finished_.log");
	system("echo. 2> _running_.log");

	s = "-r";
	if (s.compare(arg1) == 0){
		BuildRanking(argc, argv);
	}

	s = "-tp";
	if (s.compare(arg1) == 0) {
		TestPRM2(argc, argv);
	}

	s = "-t";
	if (s.compare(arg1) == 0){
		TestSeeds(argc, argv);
	}

	s = "-st";
	if (s.compare(arg1) == 0) {
		GraphStat(argc, argv);
	}

	s = "-g";
	if (s.compare(arg1) == 0) {
		GreedyAlg(argc, argv);
	}

	s = "-rr";
	if (s.compare(arg1.substr(0, 3)) == 0) {
		RRAlg(argc, argv);
	}

	// delete _running_.log to indicate finish
	system("del /Q _running_.log");
	system("echo. 2> _finished_.log");

	return 0;
}

void MICommandLine::TestSeeds(int argc, std::vector<std::string>& argv)
{
	//EventTimer timer;
	//timer.SetTimeEvent("start");
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);
	string seeds_file = "Test.txt";
	if (argc >= 3) {
		seeds_file = argv[2];
	}
	int num_iter = NUM_ITER;
	if (argc >= 4) {
		num_iter = std::stoi(argv[3]);
	}
	int seed_set_size = SET_SIZE;
	if (argc >= 5) {
		seed_set_size = std::stoi(argv[4]);
	}
	string outfile = "GC_spread.txt";
	if (argc >= 6) {
		outfile = argv[5];
	}
	int nthreads = 10;
	if (argc >= 7) {
		nthreads = std::stoi(argv[6]);
	}

	int mode = 0;
	if (argc >= 8) {
		mode = std::stoi(argv[7]);
	}

	cascade.nthreads = (nthreads <= 1) ? 1 : nthreads;

	SeedIO io;
	std::vector<int> seeds = io.Read(seeds_file, gf);
	int size = min(int(seeds.size()), seed_set_size);

	EventTimer timer;
	timer.SetTimeEvent("start");

	if (mode == 0) {
		Simu(seeds, outfile, size, num_iter).toSimulate(cascade);
	}
	else {
		Simu(seeds, outfile, size, num_iter).toSimulateOnceFile(cascade);
	}

	timer.SetTimeEvent("end");
	char timefilename[] = "time_test.txt";
	FILE *out;
	fopen_s(&out, timefilename, "w");
	fprintf(out, "%g\n", timer.TimeSpan("start", "end"));
	fclose(out);
}

void MICommandLine::TestPRM2(int argc, std::vector<std::string>& argv){
	int simu_time = 500;
	//EventTimer timer;
	//timer.SetTimeEvent("start");
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);
	string seeds_list_file = "Test_seeds.txt";
	if (argc >= 3) {
		seeds_list_file = argv[2];
	}
	double dp = 400;
	if (argc >= 4) {
		dp = std::stoi(argv[3]);
	}
	double dn = 10;
	if (argc >= 5) {
		dn = std::stoi(argv[4]);
	}
	double a = 50;
	if (argc >= 6) {
		a = std::stoi(argv[5]);
	}
	string outfile = "GC_spread.txt";

	int nthreads = 10;
	int num_iter = 1;

	cascade.nthreads = (nthreads <= 1) ? 1 : nthreads;
	std::ifstream inFile(seeds_list_file, ios_base::in);
	//inFile.open(seeds_list_file, ios_base::in);
	string seeds_file;
	vector<string> seeds_file_list;
	while (std::getline(inFile, seeds_file)){
		seeds_file_list.push_back(seeds_file);
	}

	vector<bool> active_old;
	vector<bool> active_new;
	EventTimer timer;
	double all_ratio = 0;
	for (int i = 0; i < simu_time; i++) {
		double new_dp = dp;
		double new_dn = dn;
		int marginal_index = 0;
		for (int i = 0; i < seeds_file_list.size(); i++) {
			SeedIO io;
			string seeds_file_name = seeds_file_list[i];
			std::vector<int> seeds = io.Read(seeds_file_name, gf);
			int size = int(seeds.size());

			timer.SetTimeEvent("start");
			int marginal = 0;
			active_new = Simu(seeds, outfile, size, num_iter).toSimulatePRM2(cascade);
			if (marginal_index == 0) {
				for (int acti = 0; acti < active_new.size();acti++) {
					if (active_new[acti]) {
						marginal++;
						//active_old[i] = true;
					}
				}
				active_old = active_new;
				marginal_index = 1;
			}
			else {
				for (int it = 0; it < active_new.size(); it++) {
					if (active_new[it] && !active_old[it]) {
						marginal++;
						active_old[it] = true;
					}
				}

			}
			double rate = new_dn / (new_dp + new_dn);
			new_dp += (1 - rate) * a;
			new_dn += rate * a + marginal;

		}
		all_ratio += new_dn / new_dp;
	}
	cout << "R:" << all_ratio / simu_time << endl;
	timer.SetTimeEvent("end");


	char timefilename[] = "time_test.txt";
	FILE* out;
	fopen_s(&out, timefilename, "w");
	fprintf(out, "%g\n", timer.TimeSpan("start", "end"));
	fclose(out);
}

void MICommandLine::BuildRanking(int argc, std::vector<std::string>& argv)
{
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);
	Greedy alg;
	alg.BuildRanking(gf, 100, cascade);
}


void MICommandLine::GraphStat(int argc, std::vector<std::string>& argv)
{
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GraphStatistics sta;
	sta.Stats(gf, std::cout);
}


void MICommandLine::GreedyAlg(int argc, std::vector<std::string>& argv)
{
	// GreedyGC (improved by CELF)
	int d_p, d_n, a, time, topk, mode=0;
	if (argc >= 3) topk = std::stoi(argv[2]);
	if (argc >= 4) time = std::stoi(argv[3]);
	if (argc >= 5) d_p = std::stoi(argv[4]);
	if (argc >= 6) d_n = std::stoi(argv[5]);
	if (argc >= 7) a = std::stoi(argv[6]);
	if (argc >= 8) mode = std::stoi(argv[7]);

	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	GeneralCascade cascade;
	cascade.Build(gf);

	EventTimer timer;
	timer.SetTimeEvent("start");
	Greedy alg;
	if (mode == 0) {
		alg.Build(gf, min(topk, gf.GetN()), cascade, d_p, d_n, a, time);
	}
	else {
		alg._Build(gf, min(topk, gf.GetN()), cascade, d_p, d_n, a, time);
	}
	
	timer.SetTimeEvent("end");

	/*
	FILE* timetmpfile;
	fopen_s(&timetmpfile, "time_greedy_gc.txt", "w");
	fprintf(timetmpfile, "%g\n", timer.TimeSpan("start", "end"));
	fclose(timetmpfile);
	system("copy greedy.txt greedy_gc.txt");
	//Simu(alg.GetSeedList(), "GC_Greedy.txt").toSimulate(cascade);
	system("del /Q tmp\\*");
	
	//system("pause");

	
	// GreedyGC_SPM (improved by CELF)
	QueryPerformanceCounter(&start);
	Greedy::Build(SET_SIZE,SPM_gc::Run);
	QueryPerformanceCounter(&Eend);
	timer = (double)(Eend.QuadPart - start.QuadPart) / freq.QuadPart;
	fopen_s(&timetmpfile, "time_greedy_gc_spm.txt", "w");
	fprintf(timetmpfile,"%g\n", timer);
	fclose(timetmpfile);
	system("copy greedy.txt greedy_gc_spm.txt");
	system("del /Q tmp\\*");
	toSimulate("GC_SPM.txt", Greedy::GetNode, GeneralCascade::Run);

	// GreedyGC_SP1M (improved by CELF)
	QueryPerformanceCounter(&start);
	Greedy::Build(SET_SIZE,SP1M_gc::Run);
	QueryPerformanceCounter(&Eend);
	timer = (double)(Eend.QuadPart - start.QuadPart) / freq.QuadPart;
	fopen_s(&timetmpfile,"time_greedy_gc_sp1m.txt", "w");
	fprintf(timetmpfile,"%g\n", timer);
	fclose(timetmpfile);
	system("copy greedy.txt greedy_gc_sp1m.txt");
	system("del /Q tmp\\*");
	toSimulate("GC_SP1M.txt", Greedy::GetNode, GeneralCascade::Run);*/
}



void MICommandLine::RRAlg(int argc, std::vector<std::string>& argv)
{
	string arg1(argv[1]);

	bool isSODA14 = false, isSIGMOD14 = false, isSIGMOD15 = false, isWIMM = false, isMultiIMM = false;
	bool isCIMM = false;
	int num_iter = 1000000;
	double eps = 0.1;
	double ell = 1.0;
	int mode = 1;
	int time = 10;
	double delta = 1.0;
	int topk = SET_SIZE;
	int nout = 50;
	bool isConcurrent = false;
	bool isSingleInf = false;
	bool isShapley = false;
	double d_n = 0, d_p = 0, a = 0;
	// switches:
	// -rr  -rro
	// -rr1  -rr1o
	// -rr2  -rr2o
	// -rr3  -rr3o 
	// -rrs  -rrsn -rrso
	/// -rr4  -rr4o
	if (arg1.compare("-rr") == 0) {
		isSODA14 = isSIGMOD14 = isSIGMOD15 = isCIMM = true; // test algorithms
	}
	else if (arg1.compare("-rro") == 0) {
		isSODA14 = isSIGMOD14 = isSIGMOD15 = isCIMM = true; // test algorithms with parallel optimization
		isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rr1") == 0) {
		isSODA14 = true;
		if (argc >= 3) num_iter = std::stoi(argv[2]);
		if (argc >= 4) topk = std::stoi(argv[3]);
		if (arg1.compare("-rr1o") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rr2") == 0) {
		isSIGMOD14 = true;
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (argc >= 5) topk = std::stoi(argv[4]);
		if (arg1.compare("-rr2o") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rr3") == 0) {
		isSIGMOD15 = true;
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (argc >= 5) topk = std::stoi(argv[4]);
		if (argc >= 6) mode = std::stoi(argv[5]);
		if (arg1.compare("-rr3o") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rrs") == 0) {
		isShapley = true;
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (argc >= 5) topk = std::stoi(argv[4]);
		else if (arg1.compare("-rrsn") == 0){
			isSingleInf = true;
		}
		if (arg1.compare("-rrso") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rr4") == 0){
		isCIMM = true;
		//parameters
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (argc >= 5) delta = std::stod(argv[4]);
	}
	else if (arg1.substr(0, 4).compare("-rr5") == 0) {
		isWIMM = true;
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (argc >= 5) topk = std::stoi(argv[4]);
		if (argc >= 6) mode = std::stoi(argv[5]);
		if (argc >= 6) time = std::stoi(argv[6]);
		if (argc >= 7) d_p = std::stoi(argv[7]);
		if (argc >= 8) d_n = std::stoi(argv[8]);
		if (argc >= 9) a = std::stoi(argv[9]);
		if (arg1.compare("-rr5o") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rr6") == 0) {
		isMultiIMM = true;
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (argc >= 5) topk = std::stoi(argv[4]);
		if (argc >= 6) mode = std::stoi(argv[5]);
		if (argc >= 6) time = std::stoi(argv[6]);
		if (argc >= 7) d_p = std::stoi(argv[7]);
		if (argc >= 8) d_n = std::stoi(argv[8]);
		if (argc >= 9) a = std::stoi(argv[9]);
		if (arg1.compare("-rr6o") == 0)
			isConcurrent = true;
	}
	else if (arg1.substr(0, 4).compare("-rr7") == 0) {
		isWIMM = true;
		isMultiIMM = true;
		if (argc >= 3) eps = std::stod(argv[2]);
		if (argc >= 4) ell = std::stod(argv[3]);
		if (argc >= 5) topk = std::stoi(argv[4]);
		if (argc >= 6) mode = std::stoi(argv[5]);
		if (argc >= 6) time = std::stoi(argv[6]);
		if (arg1.compare("-rr7o") == 0)
			isConcurrent = true;
	}
	GraphFactory fact;
	Graph gf = fact.Build(std::cin);
	ReverseGCascade cascade;
	cascade.Build(gf);
	

	if (isSODA14) {
		int maxK = min(topk, gf.GetN());
		std::cout << "=== Algorithm 1: SODA'14 ===" << endl;
		std::cout << "#seeds = " << maxK << endl;
		RRInfl infl;
		infl.isConcurrent = isConcurrent;
		infl.Build(gf, maxK, cascade, num_iter);
		char rrinfl_simu_file[] = "GC_rr_infl.txt";
		// toSimulate(rrinfl_simu_file, RRInfl::GetNode, GeneralCascade::Run);
	}

	if (isSIGMOD14) {
		int maxK = min(topk, gf.GetN());
		std::cout << "=== Algorithm 2: TimPlus, SIGMOD'14 ===" << endl;
		std::cout << "#seeds = " << maxK << endl;
		std::cout << "eps = " << eps << endl;
		std::cout << "ell = " << ell << endl;
		TimPlus infl;
		infl.isConcurrent = isConcurrent;
		infl.Build(gf, maxK, cascade, eps, ell);
		//char rrinfl_simu_file[] = "GC_rr_timplus_infl.txt";
		// toSimulate(rrinfl_simu_file, TimPlus::GetNode, GeneralCascade::Run);
	}

	if (isSIGMOD15) {
		int maxK = min(topk, gf.GetN());
		std::cout << "=== Algorithm 3: IMM, SIGMOD'15 ===" << endl;
		std::cout << "#seeds = " << maxK << endl;
		std::cout << "eps = " << eps << endl;
		std::cout << "ell = " << ell << endl;
		std::cout << "mode = " << mode << endl;
		std::cout << "isConcurrent = " << isConcurrent << endl;
		IMM infl;
		infl.isConcurrent = isConcurrent;
		//infl.Build(gf, maxK, cascade, eps, ell, mode);
		// char rrinfl_simu_file[] = "GC_rr_imm_infl.txt";
		// toSimulate(rrinfl_simu_file, IMM::GetNode, GeneralCascade::Run);
	}

	if (isWIMM) {
		int maxK = min(topk, gf.GetN());
		std::cout << "=== Algorithm 4: WIMM ===" << endl;
		std::cout << "#seeds = " << maxK << endl;
		std::cout << "eps = " << eps << endl;
		std::cout << "ell = " << ell << endl;
		std::cout << "mode = " << mode << endl;
		std::cout << "isConcurrent = " << isConcurrent << endl;
		std::cout << "dp0 = " << d_p << endl;
		std::cout << "dn0 = " << d_n << endl;
		std::cout << "a = " << a << endl;
		IMM infl;
		infl.kp_0 = d_p;
		infl.kb_0= d_n;
		infl.m_0 = a;
		infl.isConcurrent = isConcurrent;
		infl.Build(gf, maxK, time, cascade, eps, ell, mode);
		// char rrinfl_simu_file[] = "GC_rr_imm_infl.txt";
		// toSimulate(rrinfl_simu_file, IMM::GetNode, GeneralCascade::Run);
	}

	if (isMultiIMM) {
		int maxK = min(topk, gf.GetN());
		std::cout << "=== Algorithm 5: MultiIMM ===" << endl;
		std::cout << "#seeds = " << maxK << endl;
		std::cout << "eps = " << eps << endl;
		std::cout << "ell = " << ell << endl;
		std::cout << "mode = " << mode << endl;
		std::cout << "isConcurrent = " << isConcurrent << endl;
		std::cout << "dp0 = " << d_p << endl;
		std::cout << "dn0 = " << d_n << endl;
		std::cout << "a = " << a << endl;
		IMM infl;
		infl.kp_0 = d_p;
		infl.kb_0 = d_n;
		infl.m_0 = a;
		infl.isConcurrent = isConcurrent;
		infl._Build(gf, maxK, time, cascade, eps, ell, mode);
		// char rrinfl_simu_file[] = "GC_rr_imm_infl.txt";
		// toSimulate(rrinfl_simu_file, IMM::GetNode, GeneralCascade::Run);
	}


}





