#ifndef mi_command_line_h__
#define mi_command_line_h__

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <set>
#include <fstream>
#include <sstream>
#include <random>

#include "common.h"
#include "graph.h"
#include "graph_stat.h"
#include "greedy.h"
#include "independ_cascade.h"



#include "event_timer.h"
#include "simulate.h"


#include "rr_infl.h"



/// Command line for a set of max_influence algorithms
class MICommandLine
{
public:
	int Main(int argc, char* argv[]);
	int Main(int argc, std::vector<std::string>& argv);
	std::string Help();
	void BuildRanking(int argc, std::vector<std::string>& argv);
	void TestSeeds(int argc, std::vector<std::string>& argv);
	void TestPRM2(int argc, std::vector<std::string>& argv);
	void GraphStat(int argc, std::vector<std::string>& argv);
	void GreedyAlg(int argc, std::vector<std::string>& argv);
	void RRAlg(int argc, std::vector<std::string>& argv);
};


#endif // mi_command_line_h__
