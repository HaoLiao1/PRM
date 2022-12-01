#ifndef GREEDY_H
#define GREEDY_H

#include <vector>
#include <set>
#include <string>
#include "algo_base.h"
#include "common.h"
#include "graph.h"
#include "cascade.h"
#include "general_cascade.h"

/// Greedy algorithm with lazy-forward optimization
class Greedy
	: public AlgoBase
{
protected:
	int nthreads = 16;
	std::string file;
	void InitializeConcurrent()
	{
		if (IsConcurrent())
		{
#ifdef MI_USE_OMP
			/////////////////////////////////////////////
			// run concurrently
			const double DYNAMIC_RATIO = 0.25;
			omp_set_num_threads(nthreads);
			int dynamicThreads = (int)(nthreads * DYNAMIC_RATIO);
			omp_set_dynamic(dynamicThreads);

			std::cout << "== Turn on omp optimization: == " << std::endl;
			std::cout << "#Max Threads = " << omp_get_max_threads() << "\t#Dynamic Threads = " << omp_get_dynamic() << std::endl;
#else
			std::cout << "== omp is not supported or enabled == " << std::endl;
#endif
		}
	}

	inline bool IsConcurrent()
	{
		return (nthreads > 1);
	}

public:
	Greedy();

	void Build(IGraph& gf, int k, ICascade& cascade);
	void Build(IGraph& gf, int k, ICascade& cascade, int dp, int dn, int a, int t);
	void _Build(IGraph& gf, int k, GeneralCascade& cascade, int dp, int dn, int a, int t);
	void BuildRanking(IGraph& gf, int k, ICascade& cascade);
	void BuildFromFile(IGraph& gf, const char* filename);
	void WriteToFileWithTime(const std::string& filename, IGraph& gf, std::vector<std::pair<int, int>> listWithTime);
	void Write(const std::string& filename, const std::vector<std::pair<int, int>>& seeds, const std::vector<double>& infl, IGraph& gf);
	void Write(std::ostream& out, const std::vector<std::pair<int, int>>& seeds, const std::vector<double>& infl, IGraph& gf);
	double Weight_iter(int dn, int dp, int a, int time);
};

#endif

