///
/// rr_infl.h

#ifndef __rr_infl_h__
#define __rr_infl_h__


#include <set>
#include <vector>
#include "graph.h"
#include "common.h"
#include "reverse_general_cascade.h"
#include "algo_base.h"
#include "general_cascade.h"


/// Base class for Reverse Influence Maximization 
class RRInflBase
	: public AlgoBase
{
public:
	typedef Graph graph_type;
	typedef ReverseGCascade cascade_type;

	RRInflBase() : m(0),
		isConcurrent(false) {
	}

	// for concurrent optimization: using omp
	bool isConcurrent; // turn on to use openmp

protected:
	int m;

	std::vector< RRVec > table;
	std::vector<int> targets;
	// degree of hyper-edges v, where e(u, v) in the hyper graph
	// source id --> degrees
	std::vector<int> degrees;
	std::vector< std::vector<int> > degreeRRIndices;
	std::set<int> sourceSet; // all the source node ids

	void InitializeConcurrent();

	void _AddRRSimulation(size_t num_iter,
		cascade_type& cascade,
		std::vector< RRVec >& refTable,
		std::vector<int>& refTargets);
	void _AddRRSimulation(size_t num_iter,
		cascade_type& cascade,
		std::vector< RRVec >& refTable,
		std::vector<int>& refTargets,
		std::vector<int>& refEdgeVisited);

	double _RunGreedy(int seed_size,
		std::vector<int>& outSeeds,
		std::vector<double>& outMarginalCounts);

	void _RebuildRRIndices();
	double _EstimateInfl(const std::vector<int>& seeds, std::vector<double>& out_cumuInfl);
	void _SetResults(const std::vector<int>& seeds, const std::vector<double>& cumu_spread);

};


/// Reverse Influence Maximization
class RRInfl :
	public RRInflBase
{
protected:
	std::string file;
	std::string time_file;

public:
	RRInfl() : RRInflBase(),
		file("rr_infl.txt"),
		time_file("time_rr_infl.txt")
	{}


public:
	virtual void Build(graph_type& gf, int k, cascade_type& cascade, size_t num_iter = 1000000); // [1]
	virtual void BuildInError(graph_type& gf, int k, cascade_type& cascade, double epsilon = 0.1); // [1] 0 < epsilon < 1

protected:
	void _Build(graph_type& gf, int k, cascade_type& cascade, size_t num_iter = 0); // internal
	// methods for finding the solution
	double DefaultRounds(int n, int m, double epsilon = 0.2); // [1]

};

/// TimPlus Algorithm
class TimPlus :
	public RRInflBase
{
protected:
	std::string file;
	std::string time_file;

public:
	TimPlus() : RRInflBase(),
		file("rr_timplus_infl.txt"),
		time_file("time_rr_timplus_infl.txt") {
	}

public:
	virtual void Build(graph_type& gf, int k, cascade_type& cascade, double eps = 0.1, double ell = 1.0); // [2]

protected:
	double StepThreshold(int n, double lb, double ell = 1.0);
	double RThreshold_0(double eps, double opt, double ell = 1.0);
	double LogNChooseK(int n, int k);
	double RThreshold(double eps, double opt, int k, double ell = 1.0); // Lemma 3 in [2]
	double EpsPrime(double eps, int k, double ell = 1.0); // Last equation in Section 4.1 [2]
};

/// IMM algorithm
class IMM :
	public TimPlus
{
public:
	IMM() : TimPlus() {
		file = "rr_imm_infl.txt";
		time_file = "time_rr_imm_infl.txt";
	}

public:
	/// override Build
	void Build(graph_type& gf, int k, cascade_type& cascade, double eps = 0.1, double ell = 1.0, int mode = 0); // [3]

	double LambdaPrime(double epsprime, int k, double ell, int n); // Equation(9) in [3]
	double LambdaStar(double eps, int k, double ell, int n); // Equation(6) in [3]
};



/// PRM-IMM algorithm
/// Implementation of PRM (see: rr_inf.h)
class PRM_IMM :
	public TimPlus
{
public:
	PRM_IMM() : TimPlus() {
		file = "rr_imm_infl.txt";
		time_file = "time_rr_imm_infl.txt";
	}

public:
	float kp_0 = 990;
	float kb_0 = 10;
	float m_0 = 50;
	/// override Build
	void _Build(graph_type& gf, int k, int time, cascade_type& cascade, double eps = 0.1, double ell = 1.0, int mode = 0); // [3]
	void Build(graph_type& gf, int k, int time, cascade_type& cascade, double eps = 0.1, double ell = 1.0, int mode = 0); // [3]

	double LambdaPrime(double epsprime, int k, double ell, int n);
	double LambdaStar(double eps, int k, double ell, int n);
	double LambdaPrimeOrigin(double epsprime, int k, double ell, int n);
	double LambdaStarOrigin(double eps, int k, double ell, int n);
	double _RunGreedy1(int seed_size,
		std::vector< std::pair< int, int > >& outSeeds,
		std::vector<double>& outEstSpread);
	void _AddRRSimulation1(size_t num_iter,
		cascade_type& cascade,
		std::vector< std::pair< RRVec, int > >& refTable,
		std::vector<int>& refTargets,
		int k);

	void _RebuildRRIndicesWithTime();
	void _RebuildRRIndicesWithReuse();
	double _RunGreedyTest(int seed_size,
		std::vector< std::pair< int, int > >& outSeeds,
		std::vector<double>& outEstSpread, float& ratio);

	double FindTopK(int seed_size,
		std::vector< std::pair< int, int > >& outSeeds,
		std::vector<double>& outEstSpread);
	double uniformChoose(int seed_size,
		std::vector< std::pair< int, int > >& outSeeds,
		std::vector<double>& outEstSpread);
	double DecreasingChoose(int seed_size,
		std::vector< std::pair< int, int > >& outSeeds,
		std::vector<double>& outEstSpread);
	double RandomChoose(int seed_size,
		std::vector< std::pair< int, int > >& outSeeds,
		std::vector<double>& outEstSpread);
	double reuseRunGreedy(int seed_size,
		std::vector< std::pair< int, int > >& outSeeds,
		std::vector<double>& outEstSpread);
	double Weight_iter(int weight_mode, int time);
	void _SetResults1(const std::vector<std::pair<int, int>>& seeds, const std::vector<double>& cumu_spread);
	virtual void WriteToFileWithTime(const std::string& filename, IGraph& gf);
	void Write(const std::string& filename,
		const std::vector<std::pair<int, int>>& seeds,
		const std::vector<double>& infl, IGraph& gf);
	void Write(std::ostream& out,
		const std::vector<std::pair<int, int>>& seeds,
		const std::vector<double>& infl, IGraph& gf);
protected:
	int max_time = 0;

	int weight_mode = 1;
	/// seedsWithTime
	std::vector<std::pair<int, int>> listWithTime;

	std::vector< std::vector<double> > degreesWithTime; //c_t[v]
	std::vector< std::pair< RRVec, int > > tableWithTime; // set of RR-set with label
	std::vector< std::vector< std::vector<int> > > degreeRRIndicesWithTime; //RR_t[v]
	std::set<std::pair<int, int>> sourceSetWithTime;
	std::vector<int> RR_number;
};



#endif ///:~ __rr_infl_h__

