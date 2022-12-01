
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <functional>
#include <algorithm>
#include <cassert>

#include "rr_infl.h"
#include "reverse_general_cascade.h"
#include "graph.h"
#include "event_timer.h"
#include "common.h"

using namespace std;



// define a comparator for counts
struct CountComparator
{
public:
	vector<int>& counts;
	CountComparator(vector<int>& c) :counts(c){}
public:
	bool operator () (int a, int b) {
		return (counts[a] < counts[b]);
	}
};

struct dCountComparator
{
public:
	vector<double>& counts;
	dCountComparator(vector<double>& c) :counts(c){}
public:
	bool operator () (double a, double b) {
		return (counts[a] < counts[b]);
	}
};

struct PairdCountComparator
{
public:
	vector<pair<pair<double, int>, int>>& counts;
	PairdCountComparator(vector<pair<pair<double, int>, int>>& c) :counts(c) {}
public:
	bool operator () (double a, double b) {
		
		return (counts[a].first.first <= counts[b].first.first);
	}
};



void RRInflBase::InitializeConcurrent()
{
	if (isConcurrent) 
	{
#ifdef MI_USE_OMP
		/////////////////////////////////////////////
		// run concurrently
		const double DYNAMIC_RATIO = 1.0;
		int maxThreads = omp_get_max_threads();
		omp_set_num_threads(maxThreads);
		int dynamicThreads = (int)(maxThreads * DYNAMIC_RATIO);
		omp_set_dynamic(dynamicThreads);
		
		cout << "== Turn on omp optimization: == " << endl;
		cout << "#Max Threads = " << omp_get_max_threads() << "\t#Dynamic Threads = " << omp_get_dynamic() << endl;
#else
        cout << "== omp is not supported or enabled == " << endl;
#endif
    }
}


void RRInflBase::_AddRRSimulation(size_t num_iter, 
							cascade_type& cascade, 
							std::vector< RRVec >& refTable,
							std::vector<int>& refTargets) 
{
	vector<int> edgeVisited; // discard
	_AddRRSimulation(num_iter, cascade, refTable, refTargets, edgeVisited);
}

void RRInflBase::_AddRRSimulation(size_t num_iter, 
							  cascade_type& cascade, 
							  std::vector< RRVec >& refTable,
							  std::vector<int>& refTargets,
							  std::vector<int>& refEdgeVisited)
{
#ifdef MI_USE_OMP
	if (!isConcurrent) {
#endif
		// run single thread

		for (size_t iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited;
			cascade.ReversePropagate(1, id, refTable, edgeVisited);

			refTargets.push_back(id);
			refEdgeVisited.push_back(edgeVisited);
		}

#ifdef MI_USE_OMP
	} else {
		// run concurrently
		// #pragma omp parallel for private(iter, id, edge_visited_count, tmpTable) shared(refTable, refTargets, refEdgeVisited)

		#pragma omp parallel for ordered
		for (int iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited;
			vector<RRVec> tmpTable;
			cascade.ReversePropagate(1, id, tmpTable, edgeVisited);
			#pragma omp critical
			{
				refTable.push_back(tmpTable[0]);
				refTargets.push_back(id);
				refEdgeVisited.push_back(edgeVisited);
			}
		}
	}
#endif

}

// Apply Greedy to solve Max Cover
double RRInflBase::_RunGreedy(int seed_size,
	vector<int>& outSeeds,
	vector<double>& outEstSpread)
{
	outSeeds.clear();
	outEstSpread.clear();

	// set enables for table
	vector<bool> enables;
	enables.resize(table.size(), true);

	set<int> candidates(sourceSet);
	CountComparator comp(degrees);

	double spread = 0;
	for (int iter = 0; iter < seed_size; ++iter) {
		set<int>::const_iterator maxPt = max_element(candidates.begin(), candidates.end(), comp);
		int maxSource = *maxPt;
		// cout << "" << maxSource << "\t";
		assert(degrees[maxSource] >= 0);

		// selected one node
		outSeeds.push_back(maxSource);

		// estimate spread
		spread = spread + ((double)n * degrees[maxSource] / table.size());

		// if (iter==0)
		// {
		// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
		// }

		outEstSpread.push_back(spread);

		// clear values
		candidates.erase(maxPt);
		degrees[maxSource] = -1;

		// deduct the counts from the rest nodes
		const vector<int>& idxList = degreeRRIndices[maxSource];
		if (!isConcurrent) {
			for (int idx : idxList) {
				if (enables[idx]) {
					const RRVec& RRset = table[idx];
					for (int rr : RRset) {
						if (rr == maxSource) continue;
						degrees[rr]--; // deduct
					}
					enables[idx] = false;
				}
			}
		}
		else {
			// run concurrently
#pragma omp parallel for
			for (int idxIter = 0; idxIter < idxList.size(); ++idxIter) {
				int idx = idxList[idxIter];
				if (enables[idx]) {
					const RRVec& RRset = table[idx];
					for (int rr : RRset) {
						if (rr == maxSource) continue;

#pragma omp atomic
						degrees[rr]--; // deduct
					}
					enables[idx] = false;
				}
			}
		}
	}

	assert(outSeeds.size() == seed_size);
	assert(outEstSpread.size() == seed_size);

	return spread;
}

void RRInflBase::_RebuildRRIndices()
{
	degrees.clear();
	degrees.resize(n, 0);
	degreeRRIndices.clear();

	for (int i = 0; i < n; ++i) {
		degreeRRIndices.push_back(vector<int>());
	}

	// to count hyper edges:
	for (size_t i = 0; i < table.size(); ++i) {
		const RRVec& RR = table[i];
		for (int source : RR) {
			degrees[source]++;
			degreeRRIndices[source].push_back(i); // add index of table
		}
	}

	// add to sourceSet where node's degree > 0
	sourceSet.clear();
	for (size_t i = 0; i < degrees.size(); ++i) {
		if (degrees[i] >= 0) {
			sourceSet.insert(i);
		}
	}
}

void IMM::_AddRRSimulation1(size_t num_iter,
	cascade_type& cascade,
	std::vector< std::vector< RRVec > >& refTable,
	std::vector<int>& refTargets,
	int k)
{
#ifdef MI_USE_OMP
	if (!isConcurrent) {
#endif
		// run single thread

		for (size_t iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited; //好像没用
			std::vector< RRVec > tmpRefTable;
			for (size_t i = 0; i < k; ++i) {
				cascade.ReversePropagate(1, id, tmpRefTable, edgeVisited);
			}
			refTable.push_back(tmpRefTable);
			//cascade.ReversePropagate(1, id, refTable, edgeVisited, k - 1); //refTable 就是tablewithtime，也就是返回的反向可达集pair的集合。

			refTargets.push_back(id);
			// refEdgeVisited.push_back(edgeVisited);
		}

#ifdef MI_USE_OMP
	}
	else {
		// run concurrently
		// #pragma omp parallel for private(iter, id, edge_visited_count, tmpTable) shared(refTable, refTargets, refEdgeVisited)

#pragma omp parallel for ordered
		for (int iter = 0; iter < num_iter; ++iter) {
			int id = cascade.GenRandomNode();
			int edgeVisited;
			std::vector< RRVec > tmpRefTable;
			for (size_t i = 0; i < k; ++i) {
				cascade.ReversePropagate(1, id, tmpRefTable, edgeVisited);
			}
#pragma omp critical
			{
				refTable.push_back(tmpRefTable);
				refTargets.push_back(id);
			}
		}
	}
#endif

}


void IMM::_RebuildRRIndicesWithTime()
{
	//RR_number.clear();
	//RR_number.resize(top, 0);   //保存不同时间的反向可达集的数量
	degreesWithTime.clear();
	vector<double> temp(n);
	degreesWithTime.resize(top, temp); //k's value  分别保存不同时间，不同节点的影响力扩展度。
	degreeRRIndicesWithTime.clear();  //分别保存不同时间，不同节点cover的反向可达集
	vector< vector<int> > temp1;
	temp1.clear();
	for (int i = 0; i < n; ++i) {
		temp1.push_back(vector<int>());
	}
	degreeRRIndicesWithTime.resize(top, temp1);


	// to count hyper edges:
	for (size_t i = 0; i < tableWithTime.size(); ++i) {
		const vector< RRVec >& RR = tableWithTime[i];
		//RR_number[RR.second]++;
#pragma omp parallel for ordered
		for (int T = 0; T < RR.size(); T++) {
			//double weight = 1.0 / (k_0 + m_0 * RR.second);
			//double weight = 1.0;
			//double weight = top - RR.second ;
#pragma omp critical
			for (int source : RR[T]) {
				double weight = Weight_iter(weight_mode, T + 1);
				degreesWithTime[T][source] += weight;
				degreeRRIndicesWithTime[T][source].push_back(i); // add index of table
			}
			
		}
	}

	// add to sourceSet where node's degree > 0
	sourceSetWithTime.clear(); //
	sourceSet.clear();
	for (size_t i = 0; i < degreesWithTime.size(); ++i) {
//#pragma omp parallel for ordered
		for (int j = 0; j < degreesWithTime[i].size(); ++j) {
			if (degreesWithTime[i][j] >= 0) {
				sourceSet.insert(j);
				pair<int, int> NodeAndTime;
				NodeAndTime.first = j; //node
				NodeAndTime.second = i; //time
//#pragma omp critical
				sourceSetWithTime.insert(NodeAndTime);
			}
		}
	}
}

void IMM::_RebuildRRIndicesWithReuse()
{
	degreesWithTime.clear();
	vector<double> temp(n);
	degreesWithTime.resize(max_time, temp); //k's value  分别保存不同时间，不同节点的影响力扩展度。
	degreeRRIndicesWithTime.clear();  //分别保存不同时间，不同节点cover的反向可达集
	vector< vector<int> > temp1;
	temp1.clear();
	for (int i = 0; i < n; ++i) {
		temp1.push_back(vector<int>());
	}
	degreeRRIndicesWithTime.resize(max_time, temp1);


	// to count hyper edges:
	for (size_t i = 0; i < table.size(); ++i) {
		const RRVec& RR = table[i];
		for (int source : RR) {
			//double weight = 1.0 / (k_0 + m_0 * RR.second);
			//double weight = 1.0;
			//double weight = top - RR.second ;
			//double weight = Weight_iter(weight_mode, 0 + 1);
			degreesWithTime[0][source] += 1;
			degreeRRIndicesWithTime[0][source].push_back(i); // add index of table
		}
	}

	for (int j = 1; j < max_time; j++)
	{
		degreesWithTime[j] = degreesWithTime[0];
		degreeRRIndicesWithTime[j] = degreeRRIndicesWithTime[0];

	}


	// add to sourceSet where node's degree > 0
	sourceSetWithTime.clear(); //
	sourceSet.clear();
	for (size_t i = 0; i < degreesWithTime.size(); ++i) {
		for (size_t j = 0; j < degreesWithTime[i].size(); ++j) {
			if (degreesWithTime[i][j] >= 0) {
				sourceSet.insert(j);
				pair<int, int> NodeAndTime;
				NodeAndTime.first = j; //node
				NodeAndTime.second = i; //time
				sourceSetWithTime.insert(NodeAndTime);
			}
		}
	}
}



double IMM::_RunGreedy1(int seed_size,
	vector< pair< int, int > >& outSeeds,
	vector<double>& outEstSpread)
{
	outSeeds.clear();
	outEstSpread.clear();

	// set enables for table
	vector<bool> enables;
	enables.resize(tableWithTime.size(), true);
	vector<int> cover_round;
	cover_round.resize(tableWithTime.size(), 0);

	set<int> candidates(sourceSet);
	vector<set<int>> candidatesWithTime(top, candidates);
	// dCountComparator comp(degreesWithTime[0]);  //edited
	vector<dCountComparator> camp;
	for (int i = 0; i < top; i++) {
		dCountComparator comp(degreesWithTime[i]);
		camp.push_back(comp);
	}

	double spread = 0;
	for (int iter = 0; iter < seed_size; ++iter) {
		vector<pair<pair<double, int>, int>> winner; //记录每个时间上最大的节点
		for (int i = 0; i < top; i++) {
			set<int>::const_iterator maxPtIn = max_element(candidatesWithTime[i].begin(), candidatesWithTime[i].end(), camp[i]);
			pair<pair<double, int>, int> maxSource;
			maxSource.first.second = *maxPtIn;
			maxSource.first.first = degreesWithTime[i][*maxPtIn];
			maxSource.second = i;
			winner.push_back(maxSource);
		}
		PairdCountComparator comp_end(winner);
		set<int> Timecandidates;
		for (int i = 0; i < top; i++)
		{
			Timecandidates.insert(i);
		}
		set<int>::const_iterator maxPt = max_element(Timecandidates.begin(), Timecandidates.end(), comp_end);
		pair<int, int> maxSourceWithTime;
		maxSourceWithTime.first = winner[*maxPt].first.second; // node id
		maxSourceWithTime.second = winner[*maxPt].second; // time
		// cout << "" << maxSource << "\t";
		assert(degreesWithTime[winner[*maxPt].second][maxSourceWithTime.first] >= 0);

		// selected one node
		outSeeds.push_back(maxSourceWithTime);

		// estimate spread
		spread = spread + ((double)n * degreesWithTime[maxSourceWithTime.second][maxSourceWithTime.first] / tableWithTime.size());

		// if (iter==0)
		// {
		// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
		// }

		outEstSpread.push_back(spread);

		// clear values
		candidatesWithTime[maxSourceWithTime.second].erase(maxSourceWithTime.first);   //把选中的节点从预备节点中移除
		//degreesWithTime[maxSourceWithTime.second][maxSourceWithTime.first] = -1;   似乎是多此一举

		// deduct the counts from the rest nodes
		const vector<int>& idxList = degreeRRIndicesWithTime[maxSourceWithTime.second][maxSourceWithTime.first];
		if (isConcurrent) {
			int count = 0;
			for (int idx : idxList) {
				if (cover_round[idx]==0) {
					cover_round[idx] = maxSourceWithTime.second;
					const vector<RRVec>& RRset = tableWithTime[idx];
#pragma omp parallel for ordered
					for (int T = 0; T < RRset.size();T++) {
						for (int node : RRset[T]) {
							if (node == maxSourceWithTime.first && T == maxSourceWithTime.second) continue;
#pragma omp critical
							if (T < maxSourceWithTime.second) {
								count++;
								degreesWithTime[T][node] -= Weight_iter(weight_mode, maxSourceWithTime.second + 1);
							}
							else {
								degreesWithTime[T][node] -= Weight_iter(weight_mode, T + 1); // deduct
							}
						}
						
					}
					
				}
				else if (cover_round[idx] > maxSourceWithTime.second) {
					int old_round = cover_round[idx];
					cover_round[idx] = maxSourceWithTime.second;
					const vector<RRVec>& RRset = tableWithTime[idx];
#pragma omp parallel for ordered
					for (int T = 0; T < RRset.size(); T++) {
						for (int node : RRset[T]) {
							if (node == maxSourceWithTime.first && T == maxSourceWithTime.second) continue;
#pragma omp critical
							if (T < maxSourceWithTime.second) {
								count++;
								degreesWithTime[T][node] -= Weight_iter(weight_mode, maxSourceWithTime.second + 1) - Weight_iter(weight_mode, old_round + 1);
							}
							else if(old_round > T){
								degreesWithTime[T][node] -= Weight_iter(weight_mode, T + 1) - Weight_iter(weight_mode, old_round + 1);
							}
						}

					}
				}
				if(count == 0) enables[idx] = false;
			}
		}
		else {
			int count = 0;
			for (int idx : idxList) {
				if (cover_round[idx] == 0) {
					cover_round[idx] = maxSourceWithTime.second;
					const vector<RRVec>& RRset = tableWithTime[idx];
#pragma omp parallel for ordered
					for (int T = 0; T < RRset.size(); T++) {
						for (int node : RRset[T]) {
							if (node == maxSourceWithTime.first && T == maxSourceWithTime.second) continue;
#pragma omp critical
							if (T < maxSourceWithTime.second) {
								count++;
								degreesWithTime[T][node] -= Weight_iter(weight_mode, maxSourceWithTime.second + 1);
							}
							else {
								degreesWithTime[T][node] -= Weight_iter(weight_mode, T + 1); // deduct
							}
						}

					}

				}
				else if (cover_round[idx] > maxSourceWithTime.second) {
					int old_round = cover_round[idx];
					cover_round[idx] = maxSourceWithTime.second;
					const vector<RRVec>& RRset = tableWithTime[idx];
#pragma omp parallel for ordered
					for (int T = 0; T < RRset.size(); T++) {
						for (int node : RRset[T]) {
							if (node == maxSourceWithTime.first && T == maxSourceWithTime.second) continue;
#pragma omp critical
							if (T < maxSourceWithTime.second) {
								count++;
								degreesWithTime[T][node] -= Weight_iter(weight_mode, maxSourceWithTime.second + 1) - Weight_iter(weight_mode, old_round + 1);
							}
							else if (old_round > T) {
								degreesWithTime[T][node] -= Weight_iter(weight_mode, T + 1) - Weight_iter(weight_mode, old_round + 1);
							}
						}

					}
				}
				if (count == 0) enables[idx] = false;
			}
			/*
			// run concurrently
#pragma omp parallel for
			for (int idxIter = 0; idxIter < idxList.size(); ++idxIter) {
				int idx = idxList[idxIter];
				if (enables[idx]) {
					const pair<RRVec, int>& RRset = tableWithTime[idx];
					for (int rr : RRset.first) {
						if (rr == maxSourceWithTime.first) continue;

#pragma omp atomic
						degreesWithTime[maxSourceWithTime.second][rr] -= Weight_iter(weight_mode, maxSourceWithTime.second + 1); // deduct
					}
					enables[idx] = false;
				}
			}
			*/
		
		}
		
	}

	assert(outSeeds.size() == seed_size);
	assert(outEstSpread.size() == seed_size);

	return spread;
}

double IMM::_RunGreedyTest(int seed_size,
	vector< pair< int, int > >& outSeeds,
	vector<double>& outEstSpread,
	float& ratio)
{
	outSeeds.clear();
	outEstSpread.clear();

	// set enables for table
	vector<bool> enables;
	enables.resize(tableWithTime.size(), true);

	set<int> candidates(sourceSet);
	// dCountComparator comp(degreesWithTime[0]);  //edited
	vector<dCountComparator> camp;
	for (int i = 0; i < top; i++) {
		dCountComparator comp(degreesWithTime[i]);
		camp.push_back(comp);
	}


	double spread1 = 0;
	double spread2 = 0;
	for (int iter = 0; iter < seed_size; ++iter) {
		vector<pair<pair<double, int>, int>> winner; //记录每个时间上最大的节点
		for (int i = 0; i < top; i++) {
			set<int>::const_iterator maxPtIn = max_element(candidates.begin(), candidates.end(), camp[i]);
			pair<pair<double, int>, int> maxSource;
			maxSource.first.second = *maxPtIn;
			maxSource.first.first = degreesWithTime[i][*maxPtIn];
			maxSource.second = i;
			winner.push_back(maxSource);
		}
		PairdCountComparator comp_end(winner);
		set<int> Timecandidates;
		for (int i = 0; i < top; i++)
		{
			Timecandidates.insert(i);
		}
		set<int>::const_iterator maxPt = max_element(Timecandidates.begin(), Timecandidates.end(), comp_end);
		pair<int, int> maxSourceWithTime;
		maxSourceWithTime.first = winner[*maxPt].first.second; //node id
		maxSourceWithTime.second = winner[*maxPt].second; //time
		// cout << "" << maxSource << "\t";
		assert(degreesWithTime[winner[*maxPt].second][maxSourceWithTime.first] >= 0);

		// selected one node
		outSeeds.push_back(maxSourceWithTime);

		// estimate spread
		spread1 = spread1 + ((double)n * degreesWithTime[maxSourceWithTime.second][maxSourceWithTime.first] / RR_number[maxSourceWithTime.second]);
		spread2 = spread2 + ((double)n * degreesWithTime[maxSourceWithTime.second][maxSourceWithTime.first] / RR_number[maxSourceWithTime.second]);
		// if (iter==0)
		// {
		// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
		// }

		outEstSpread.push_back(spread1);

		// clear values
		candidates.erase(maxSourceWithTime.first);   //把选中的节点从预备节点中移除
		//degreesWithTime[maxSourceWithTime.second][maxSourceWithTime.first] = -1;   似乎是多此一举

		// deduct the counts from the rest nodes
		/*
		const vector<int>& idxList = degreeRRIndicesWithTime[maxSourceWithTime.second][maxSourceWithTime.first];
		if (!isConcurrent) {
			for (int idx : idxList) {
				if (enables[idx]) {
					const pair<RRVec, int>& RRset = tableWithTime[idx];
					for (int rr : RRset.first) {
						if (rr == maxSourceWithTime.first) continue;
						degreesWithTime[maxSourceWithTime.second][rr] -= Weight_iter(weight_mode, maxSourceWithTime.second + 1); // deduct
					}
					enables[idx] = false;
				}
			}
		}
		else {
			// run concurrently
#pragma omp parallel for
			for (int idxIter = 0; idxIter < idxList.size(); ++idxIter) {
				int idx = idxList[idxIter];
				if (enables[idx]) {
					const pair<RRVec, int>& RRset = tableWithTime[idx];
					for (int rr : RRset.first) {
						if (rr == maxSourceWithTime.first) continue;

#pragma omp atomic
						degreesWithTime[maxSourceWithTime.second][rr] -= Weight_iter(weight_mode, maxSourceWithTime.second + 1); // deduct
					}
					enables[idx] = false;
				}
			}
		}
		*/
	}
	for (int iter = 0; iter < seed_size; ++iter) {
		vector<pair<pair<double, int>, int>> winner; //记录每个时间上最大的节点
		for (int i = 0; i < top; i++) {
			set<int>::const_iterator maxPtIn = max_element(candidates.begin(), candidates.end(), camp[i]);
			pair<pair<double, int>, int> maxSource;
			maxSource.first.second = *maxPtIn;
			maxSource.first.first = degreesWithTime[i][*maxPtIn];
			maxSource.second = i;
			winner.push_back(maxSource);
		}
		PairdCountComparator comp_end(winner);
		set<int> Timecandidates;
		for (int i = 0; i < top; i++)
		{
			Timecandidates.insert(i);
		}
		set<int>::const_iterator maxPt = max_element(Timecandidates.begin(), Timecandidates.end(), comp_end);
		pair<int, int> maxSourceWithTime;
		maxSourceWithTime.first = winner[*maxPt].first.second; //node id
		maxSourceWithTime.second = winner[*maxPt].second; //time
		// cout << "" << maxSource << "\t";
		assert(degreesWithTime[winner[*maxPt].second][maxSourceWithTime.first] >= 0);

		// selected one node
		//outSeeds.push_back(maxSourceWithTime);

		// estimate spread
		spread2 = spread2 + ((double)n * degreesWithTime[maxSourceWithTime.second][maxSourceWithTime.first] / RR_number[maxSourceWithTime.second]);
		// if (iter==0)
		// {
		// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
		// }

		//outEstSpread.push_back(spread1);

		// clear values
		candidates.erase(maxSourceWithTime.first);   //把选中的节点从预备节点中移除
		//degreesWithTime[maxSourceWithTime.second][maxSourceWithTime.first] = -1;   似乎是多此一举


	}
	assert(outSeeds.size() == seed_size);
	assert(outEstSpread.size() == seed_size);
	ratio = spread1 / spread2;
	return spread1;
}

// random choose
double IMM::RandomChoose(int seed_size,
	vector< pair< int, int > >& outSeeds,
	vector<double>& outEstSpread)
{
	outSeeds.clear();
	outEstSpread.clear();

	// set enables for table
	vector<bool> enables;
	enables.resize(table.size(), true);

	set<int> candidates(sourceSet);
	CountComparator comp(degrees);

	double spread = 0;
	vector<int> nodes;
	nodes.resize(n);
	for (int i = 0; i < n; i++)
		nodes[i] = i;
	std::random_device rd;
	std::mt19937 g(rd());

	std::shuffle(nodes.begin(), nodes.end(), g);

	MIRandom random;

	// selected one node
	for (int i = 0; i < top; i++)
	{
		pair< int, int > node;
		node.first = nodes[i];
		int time = random.RandInt(1, top);
		node.second = time;
		outSeeds.push_back(node);
		// estimate spread
		int iter = random.RandInt(1, top);
		spread = spread + Weight_iter(weight_mode, time + 1) * ((double)n * degrees[nodes[i]] / table.size());

		outEstSpread.push_back(spread);
	}
	// clear values
	assert(outSeeds.size() == seed_size);
	assert(outEstSpread.size() == seed_size);

	return spread;
}

// Find Top K
double IMM::FindTopK(int seed_size,
	vector< pair< int, int > >& outSeeds,
	vector<double>& outEstSpread)
{
	outSeeds.clear();
	outEstSpread.clear();


	// set enables for table
	vector<bool> enables;
	enables.resize(table.size(), true);

	set<int> candidates(sourceSet);
	CountComparator comp(degrees);

	double spread = 0;
	for (int iter = 0; iter < seed_size; ++iter) {
		set<int>::const_iterator maxPt = max_element(candidates.begin(), candidates.end(), comp);
		int maxSource = *maxPt;
		// cout << "" << maxSource << "\t";
		assert(degrees[maxSource] >= 0);

		// selected one node
		pair< int, int > node;
		node.first = maxSource;
		node.second = iter;
		outSeeds.push_back(node);

		// estimate spread
		spread = spread + (1 / (kb_0 + kp_0 + m_0 * (iter + 1))) * ((double)n * degrees[maxSource] / table.size());

		// if (iter==0)
		// {
		// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
		// }

		outEstSpread.push_back(spread);

		// clear values
		candidates.erase(maxPt);
		degrees[maxSource] = -1;

		// deduct the counts from the rest nodes

	}

	assert(outSeeds.size() == seed_size);
	assert(outEstSpread.size() == seed_size);

	return spread;
}

// uniform choose
double IMM::uniformChoose(int seed_size,
	vector< pair< int, int > >& outSeeds,
	vector<double>& outEstSpread)
{
	outSeeds.clear();
	outEstSpread.clear();
	int less = top / max_time;
	int more_number = top % max_time;
	int more = less + 1;



	set<int> candidates(sourceSet);
	//CountComparator comp(degrees);

	double spread = 0;
	for (int iter = 0; iter < max_time; ++iter) {
		dCountComparator comp(degreesWithTime[iter]);
		if (iter < more_number)
		{
			// set enables for table
			vector<bool> enables;
			enables.resize(table.size(), true);
			for (int j = 0; j < more; j++)
			{

				set<int>::const_iterator maxPt = max_element(candidates.begin(), candidates.end(), comp);
				int maxSource = *maxPt;
				// cout << "" << maxSource << "\t";
				assert(degreesWithTime[iter][maxSource] >= 0);

				// selected one node
				pair< int, int > node;
				node.first = maxSource;
				node.second = iter;
				outSeeds.push_back(node);

				// estimate spread
				spread = spread + Weight_iter(weight_mode, iter + 1) * ((double)n * degreesWithTime[iter][maxSource] / table.size());

				// if (iter==0)
				// {
				// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
				// }

				outEstSpread.push_back(spread);

				// clear values
				candidates.erase(maxSource);
				degreesWithTime[iter][maxSource] = -1;

				// deduct the counts from the rest nodes
				const vector<int>& idxList = degreeRRIndicesWithTime[iter][maxSource];
				if (!isConcurrent) {
					for (int idx : idxList) {
						if (enables[idx]) {
							const RRVec& RRset = table[idx];
							for (int rr : RRset) {
								if (rr == maxSource) continue;
								degreesWithTime[iter][rr] -= 1; // deduct
							}
							enables[idx] = false;
						}
					}
				}
				else {
					// run concurrently
#pragma omp parallel for
					for (int idxIter = 0; idxIter < idxList.size(); ++idxIter) {
						int idx = idxList[idxIter];
						if (enables[idx]) {
							const RRVec& RRset = table[idx];
							for (int rr : RRset) {
								if (rr == maxSource) continue;

#pragma omp atomic
								degreesWithTime[iter][rr] -= 1; // deduct
							}
							enables[idx] = false;
						}
					}
				}
			}
		}
		else
		{
			// set enables for table
			vector<bool> enables;
			enables.resize(table.size(), true);
			for (int j = 0; j < less; j++)
			{

				set<int>::const_iterator maxPt = max_element(candidates.begin(), candidates.end(), comp);
				int maxSource = *maxPt;
				// cout << "" << maxSource << "\t";
				assert(degreesWithTime[iter][maxSource] >= 0);

				// selected one node
				pair< int, int > node;
				node.first = maxSource;
				node.second = iter;
				outSeeds.push_back(node);

				// estimate spread
				spread = spread + Weight_iter(weight_mode, iter + 1) * ((double)n * degreesWithTime[iter][maxSource] / table.size());

				// if (iter==0)
				// {
				// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
				// }

				outEstSpread.push_back(spread);

				// clear values
				candidates.erase(maxSource);
				degreesWithTime[iter][maxSource] = -1;

				// deduct the counts from the rest nodes
				const vector<int>& idxList = degreeRRIndicesWithTime[iter][maxSource];
				if (!isConcurrent) {
					for (int idx : idxList) {
						if (enables[idx]) {
							const RRVec& RRset = table[idx];
							for (int rr : RRset) {
								if (rr == maxSource) continue;
								degreesWithTime[iter][rr] -= 1; // deduct
							}
							enables[idx] = false;
						}
					}
				}
				else {
					// run concurrently
#pragma omp parallel for
					for (int idxIter = 0; idxIter < idxList.size(); ++idxIter) {
						int idx = idxList[idxIter];
						if (enables[idx]) {
							const RRVec& RRset = table[idx];
							for (int rr : RRset) {
								if (rr == maxSource) continue;

#pragma omp atomic
								degreesWithTime[iter][rr] -= 1; // deduct
							}
							enables[idx] = false;
						}
					}
				}
			}
		}

	}

	assert(outSeeds.size() == seed_size);
	assert(outEstSpread.size() == seed_size);

	return spread;
}
double IMM::reuseRunGreedy(int seed_size,
	vector< pair< int, int > >& outSeeds,
	vector<double>& outEstSpread)
{
	outSeeds.clear();
	outEstSpread.clear();

	set<int> candidates(sourceSet);
	double spread = 0;
	int seed_size_temp = 0;
	seed_size_temp = seed_size;
	std::random_device rd;
	std::mt19937 g(rd());
	MIRandom random;
	std::vector<int> seize_size;
	seize_size.resize(max_time, 0);
	for (int i = 0; i < seed_size; i++)
	{
		seize_size[random.RandInt(0, max_time - 1)]++;
	}
	for (int iter = 0; iter < max_time; iter++)
	{
		dCountComparator comp(degreesWithTime[iter]);
		// set enables for table
		vector<bool> enables;
		enables.resize(table.size(), true);
		for (int j = 0; j < seize_size[iter]; j++)
		{
			set<int>::const_iterator maxPt = max_element(candidates.begin(), candidates.end(), comp);
			int maxSource = *maxPt;
			// cout << "" << maxSource << "\t";
			assert(degreesWithTime[iter][maxSource] >= 0);

			// selected one node
			pair< int, int > node;
			node.first = maxSource;
			node.second = iter;
			outSeeds.push_back(node);

			// estimate spread
			spread = spread + Weight_iter(weight_mode, iter + 1) * ((double)n * degreesWithTime[iter][maxSource] / table.size());

			// if (iter==0)
			// {
			// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
			// }

			outEstSpread.push_back(spread);

			// clear values
			candidates.erase(maxPt);
			degreesWithTime[iter][maxSource] = -1;

			// deduct the counts from the rest nodes
			const vector<int>& idxList = degreeRRIndicesWithTime[iter][maxSource];
			if (!isConcurrent) {
				for (int idx : idxList) {
					if (enables[idx]) {
						const RRVec& RRset = table[idx];
						for (int rr : RRset) {
							if (rr == maxSource) continue;
							degreesWithTime[iter][rr] -= 1; // deduct
						}
						enables[idx] = false;
					}
				}
			}
			else {
				// run concurrently
#pragma omp parallel for
				for (int idxIter = 0; idxIter < idxList.size(); ++idxIter) {
					int idx = idxList[idxIter];
					if (enables[idx]) {
						const RRVec& RRset = table[idx];
						for (int rr : RRset) {
							if (rr == maxSource) continue;

#pragma omp atomic
							degreesWithTime[iter][rr] -= 1; // deduct
						}
						enables[idx] = false;
					}
				}
			}
		}
	}


	assert(outSeeds.size() == seed_size);
	assert(outEstSpread.size() == seed_size);

	return spread;
}

double IMM::Weight_iter(int weight_mode, int time)
{
	if (weight_mode == 1)
	{
		return 1.0 / (kb_0 + kp_0 + m_0 * float(time));

	}
	else if (weight_mode == 2)
	{
		return 1.0 / (200 + 50 * float(time));
	}
	else if (weight_mode == 3)
	{
		return 1.0 / (50 + 50 * float(time));
	}
}

void IMM::WriteToFileWithTime(const std::string& filename, IGraph& gf)
{
	
	Write(filename, listWithTime, d, gf);
}

void IMM::Write(const std::string& filename,
	const std::vector<std::pair<int, int>>& seeds,
	const std::vector<double>& infl, IGraph& gf)
{

	std::ofstream ft(filename.c_str());
	Write(ft, seeds, infl, gf);
}

void IMM::Write(std::ostream& out,
	const std::vector<std::pair<int, int>>& seeds,
	const std::vector<double>& infl, IGraph& gf)
{
	size_t num = seeds.size();
	out << num << std::endl;
	int MaxTime = 0;
	for (size_t i = 0; i < num; i++) {
		if (seeds[i].second > MaxTime)
		{
			MaxTime = seeds[i].second;
		}
	}
	out << "MaxTime: " << MaxTime << std::endl;

	for (size_t i = 0; i < num; i++) {
		out << "Node: " << gf.MapIndexToNodeName(seeds[i].first) << "\t" << "Time: " << seeds[i].second << "\t" << infl[i] << std::endl;
	}
}



vector<int> vectors_intersection(vector<int> v1, vector<int> v2){
	vector<int> v;
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v));//求交集   
	return v;
}






void RRInflBase::_SetResults(const vector<int>& seeds, const vector<double>& cumu_spread)
{
	for (int i = 0; i < top; ++i) {
		list[i] = seeds[i];
		// d[i] = cumu_spread[i];  // use cumulative spread
		d[i] = (i > 0) ? (cumu_spread[i] - cumu_spread[i - 1]) : cumu_spread[i]; // use marginal spread
	}
}
void IMM::_SetResults1(const vector<pair<int, int>>& seeds, const vector<double>& cumu_spread)
{
	for (int i = 0; i < listWithTime.size(); i++)
	{
		listWithTime[i] = seeds[i];
		listWithTime[i].second += 1;
		d[i] = (i>0)? (cumu_spread[i] - cumu_spread[i - 1]) : cumu_spread[i];
	}
}

// Same thing has been implemented as a part of _Greedy
double RRInflBase::_EstimateInfl(const vector<int>& seeds, vector<double>& out_cumu_spread)
{
	set<int> covered;
	double spd = 0;

	vector<bool> enables;
	enables.resize(table.size(), true);
	for (size_t i = 0; i < seeds.size(); ++i) {
		int sd = seeds[i];
		for (int idx : degreeRRIndices[sd]) {
			if (enables[idx]) {
				covered.insert(idx);
				enables[idx] = false;
			}
		}
		spd = (double)(n * covered.size()) / table.size();
		out_cumu_spread.push_back(spd);
	}
	return spd;
}



////////////////////////////////////////////////////////////////////
// RRInfl: paper [1]
double RRInfl::DefaultRounds(int n, int m, double epsilon)
{
	return max(144.0 * (n + m) / pow(epsilon, 3) * log( max(n,1) ), 1.0); // to make it positive
}

void RRInfl::BuildInError(graph_type& gf, int k, cascade_type& cascade, double epsilon/*=0.1*/)
{
	n = gf.GetN();
	m = gf.GetM();
	size_t num_iter = (size_t)(ceil(DefaultRounds(n, m, epsilon)));
	_Build(gf, k, cascade, num_iter);
}

void RRInfl::Build(graph_type& gf, int k, cascade_type& cascade, size_t num_iter/*=1000000*/)
{
	n = gf.GetN();
	m = gf.GetM();
	_Build(gf, k, cascade, num_iter);
}
	

void RRInfl::_Build(graph_type& gf, int k, cascade_type& cascade, size_t num_iter)
{
	InitializeConcurrent();

	n = gf.GetN();
	m = gf.GetM();

	top = k;
	d.resize(top, 0.0);
	list.resize(top, 0);

	cascade.Build(gf);

	cout << "#round = " << num_iter << endl;

	// table: target_id  RRset
	// id1  RR: rr1, rr2 ...
	// id2  RR: rr1, rr2 ...
	// ...
	table.clear();
	targets.clear();

	EventTimer pctimer;
	pctimer.SetTimeEvent("start");

	// Step 1:
	_AddRRSimulation(num_iter, cascade, table, targets);
	assert(targets.size() == num_iter);
	assert(table.size() == num_iter);
	pctimer.SetTimeEvent("step1");

	// Step 2:
	vector<int> seeds;
	vector<double> est_spread;
	_RebuildRRIndices();
	pctimer.SetTimeEvent("step2");

	// Step 3:
	double spd = _RunGreedy(k, seeds, est_spread);
	_SetResults(seeds, est_spread);
	pctimer.SetTimeEvent("end");

	cout << "  final (estimated) spread = " << spd << "\t round = " << num_iter << endl;
	
	// Write results to file:
	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i=0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);
	//fclose(out);
	WriteToFile(file, gf);

	FILE *timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile,"%g\n", pctimer.TimeSpan("start", "end"));
	fprintf(timetmpfile,"Gen graph: %g\n", pctimer.TimeSpan("start", "step1"));
	fprintf(timetmpfile,"Build RR: %g\n", pctimer.TimeSpan("step1", "step2"));
	fprintf(timetmpfile,"Greedy: %g\n", pctimer.TimeSpan("step2", "end"));
	fclose(timetmpfile);
}




///////////////////////////////////////////////////////////////////////////
/// TimPlus: paper 2
void TimPlus::Build(graph_type& gf, int k, cascade_type& cascade, double eps/*=0.1*/, double ell/*=1.0*/)
{
	InitializeConcurrent();

	n = gf.GetN();
	m = gf.GetM();
	top = k;
	d.resize(top, 0.0);
	list.resize(top, 0);

	ell = ell + log(2) / log(n); // Original IMM has failure probability 2/n^ell, so use this transformation to 
	// make the failure probability 1/n^ell

	cascade.Build(gf);
	table.clear();
	targets.clear();

	EventTimer pctimer;
	EventTimer pctimer2;
	pctimer.SetTimeEvent("start");

	// Step1: 
	double est_opt1 = 0.0;
	size_t round1 = 0;
	{
		cout << "Step 1: estimate EPT" << endl;
		pctimer2.SetTimeEvent("est_start");
		int maxIth = (int) ceil(log(n) / log(2) - 1);
		for (int ith=1; ith < maxIth; ++ith) {
			double lb = pow((double)0.5, ith);

			int loop = (int)(StepThreshold(n, lb, ell) + 1);

			vector<int> edgeVisited;
			_AddRRSimulation(loop, cascade, table, targets, edgeVisited);
			round1 += loop;

			double sumKappa = 0.0;
			for (int visit : edgeVisited) {
				double width = (double)visit / m;
				double kappa = 1.0 - pow(1.0 - width, k);
				sumKappa += kappa;
			}
			if (sumKappa / loop > lb) {
				// Algorithm 2 in [2] uses  (n * sumKappa) / (2.0 * loop), but their code does not have 2.0
				est_opt1 = (n * sumKappa) / (loop); 
				break;
			}
		}
		est_opt1 = max(est_opt1, (double)k);
		pctimer2.SetTimeEvent("est_end");
	}
	cout << "  est_opt1 = " << est_opt1 << "\t round1 = " << round1 << endl;
	pctimer.SetTimeEvent("step1");
	double time_per_round = pctimer2.TimeSpan("est_start", "est_end") / ((double) round1);
	cout << "  (Time per round = " << time_per_round << ")" << endl;

	// Step2: Use greedy to estimate opt lower bound
	double eps_step2, spd2, est_opt2;
	size_t round2;
	{
		eps_step2 = EpsPrime(eps, k, ell);
		vector<int> seeds2;
		vector<double> est_spread2;
		double theta = RThreshold_0(eps, est_opt1, ell);
		round2 = (size_t)max(theta+1.0, 1.0);
		cout << "Step 2: estimate opt by greedy. round = " << round2 << endl
			 <<  "  (Estimate time = " << time_per_round * round2 << ")"  << endl;
		_AddRRSimulation(round2, cascade, table, targets);
		_RebuildRRIndices();
		spd2 = _RunGreedy(k, seeds2, est_spread2);
		est_opt2 = spd2 / (1+eps_step2);
		est_opt2 = max(est_opt2, est_opt1);
	}
	cout << "  est_opt2 = " << est_opt2 << "\t round2 = " << round2 << endl;
	pctimer.SetTimeEvent("step2");

	// Step3: Set final runs
	vector<int> seeds3;
	vector<double> est_spread3;
	double spd3;
	size_t round3;
	{
		double theta = RThreshold(eps, est_opt2, k, ell);
		round3 = (size_t)max(theta+1.0, 1.0);
		cout << "Step 3: final greedy. round = " << round3 << endl
			 <<  "  (Estimate time = " << time_per_round * round3 << ")" << endl;
		pctimer.SetTimeEvent("step3-1");
		_AddRRSimulation(round3, cascade, table, targets);
		pctimer.SetTimeEvent("step3-2");
		_RebuildRRIndices();
		pctimer.SetTimeEvent("step3-3");
		spd3 = _RunGreedy(top, seeds3, est_spread3);
	}
	cout << "  final spread = " << spd3 << "\t round3 = " << round3 << endl;
	_SetResults(seeds3, est_spread3);
	pctimer.SetTimeEvent("end");

	// Write results to file:
	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i=0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);
	//fclose(out);
	WriteToFile(file, gf);

	FILE *timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile,"%g\n", pctimer.TimeSpan("start", "end"));
	fprintf(timetmpfile,"Step1: %g\n", pctimer.TimeSpan("start", "step1"));
	fprintf(timetmpfile,"Step2: %g\n", pctimer.TimeSpan("step1", "step2"));
	fprintf(timetmpfile,"Step3: %g\n", pctimer.TimeSpan("step2", "end"));
	fprintf(timetmpfile,"   Step3-1: %g\n", pctimer.TimeSpan("step3-1", "step3-2"));
	fprintf(timetmpfile,"   Step3-2: %g\n", pctimer.TimeSpan("step3-2", "step3-3"));
	fprintf(timetmpfile,"   Step3-3: %g\n", pctimer.TimeSpan("step3-3", "end"));
	fclose(timetmpfile);
}

double TimPlus::LogNChooseK(int n, int k){
	double ans = 0.0;
	assert(n >= k && k >= 0);
	for(int i = 0; i < k; i++){
		ans += log(n-i);
	}
	for(int i = 1; i <= k; i++){
		ans -= log(i);
	}
	return ans;
}

double TimPlus::RThreshold_0(double eps, double opt, double ell/*=1.0*/)
{
	double lambda = (double)(8 + 2*eps) * n * (ell*log(n) + log(2)) / (eps*eps);
	double theta = lambda / (opt*4.0);
	return theta;
}

double TimPlus::RThreshold(double eps, double opt, int k, double ell/*=1.0*/)
{
	double lambda = (double)(8 + 2*eps) * n * (ell*log(n) + LogNChooseK(n,k) + log(2)) / (eps*eps);
	double theta = lambda / opt;
	return theta;
}


double TimPlus::EpsPrime(double eps, int k, double ell/*=1.0*/) 
{
	return 5.0 * pow(eps*eps*ell/((double)ell+k), 1.0/3.0);
}

double TimPlus::StepThreshold(int n, double lb, double ell/*=1.0*/)
{
	double logValue = max(log(n) / log(2.0), 1.0);
	double loop = ((double)6 * ell * log(n) + 6 * log(logValue)) / lb;
	return loop;
}


///////////////////////////////////////////////////////////////////////////
/// IMM: paper 3
/// mode = 0, the original IMM algorithm
///      = 1, the first workaround fixing the IMM issue reported in arXiv:1808.09363
///      = 2, the second workaround fixing the IMM issue reported in arXiv:1808.09363
void IMM::_Build(graph_type& gf, int k, int time, cascade_type& cascade, double eps /* = 0.1 */, double ell /* = 1.0 */, int mode /* = 0 */)
{
	InitializeConcurrent();

	EventTimer pctimer;
	vector<EventTimer> steptimers;
	pctimer.SetTimeEvent("start");

	n = gf.GetN();
	m = gf.GetM();
	max_time = time;
	top = k;
	d.resize(top, 0.0);
	list.resize(top, 0);
	pair<int, int> listWT;
	listWithTime.resize(top, listWT);

	cascade.Build(gf);

	table.clear();
	targets.clear();

	double epsprime = eps * sqrt(2.0); // eps'
	double LB = 1.0;   // lower bound
	double maxRounds = max(max(log2((double)n), 1.0) - 1.0, 1.0);


	ell = ell + log(2) / log(n); // Original IMM has failure probability 2/n^ell, so use this transformation to 
	// make the failure probability 1/n^ell

	// first a few passes of sampling
	double spread = 0.0, spread1 = 0.0, spread2 = 0.0, spread3 = 0.0;
	for (int r = 1; r < maxRounds; r++) {
		EventTimer stept;
		stept.SetTimeEvent("step_start");

		cout << "  Step" << r << ":" << endl;
		double x = max(double(n) / pow(2, r), 1.0);
		double theta = LambdaPrimeOrigin(epsprime, 10, ell, n) / x;

		// select nodes
		vector<int> seeds;
		vector<double> est_spread;
		LARGE_INT64 nNewSamples = LARGE_INT64(theta - table.size() + 1);
		if (table.size() < theta) {
			// generate samples
			_AddRRSimulation(nNewSamples, cascade, table, targets);
			_RebuildRRIndices();
			spread = _RunGreedy(10, seeds, est_spread);
			//_SetResults(seeds, est_spread);
		}
		cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;

		stept.SetTimeEvent("step_end");
		steptimers.push_back(stept);

		// check influence and set lower bound
		if (spread >= ((1.0 + epsprime) * x)) {
			LB = spread / (1.0 + epsprime);
			break;
		}
	}

	// final pass of sampling
	{
		EventTimer stept;
		stept.SetTimeEvent("step_start");
		cout << "  Estimated Lower bound: " << LB << endl;
		cout << "  Final Step:" << endl;
		double theta = LambdaStarOrigin(eps, 10, ell, n) / LB;

		vector<int> seeds;
		vector<double> est_spread;

		vector< pair< int, int > > seeds1, seeds2, seeds3, seeds4;
		vector<double> est_spread1, est_spread2, est_spread3, est_spread4;

	
		if (mode == 1) {
			LARGE_INT64 nNewSamples = LARGE_INT64(theta + 1);
			cout << "  IMM Workaround 1 --- Regenerating RR sets, # RR sets = " << nNewSamples << endl;
			table.clear();
			targets.clear();
			// generate samples
			_AddRRSimulation(nNewSamples, cascade, table, targets);
			_RebuildRRIndices();

			//random choose
			double random_spread_accmulate=0;
			for (int i = 0; i < 100; i++)
			{
				spread2 = RandomChoose(top, seeds2, est_spread2);
				random_spread_accmulate += spread2;
			}
			spread2 = random_spread_accmulate / 100.0;
			_SetResults1(seeds2, est_spread2);
			WriteToFileWithTime("rr_random_infl.txt", gf);

			//find top k
			if (k != max_time)
			{
				_RebuildRRIndicesWithReuse();
				spread1 = uniformChoose(top, seeds1, est_spread1);
				_SetResults1(seeds1, est_spread1);
				WriteToFileWithTime("rr_unniform_infl.txt", gf);
			}
			else
			{
				spread1 = FindTopK(top, seeds1, est_spread1);
				_SetResults1(seeds1, est_spread1);
				WriteToFileWithTime("rr_unniform_infl.txt", gf);
			}
			//random seize
			random_spread_accmulate = 0;
			for (int i = 0; i < 10; i++)
			{
				_RebuildRRIndicesWithReuse();
				spread3 = reuseRunGreedy(top, seeds3, est_spread3);
				random_spread_accmulate += spread3;
			}
			spread3 = random_spread_accmulate / 10.0;
			_SetResults1(seeds3, est_spread3);
			WriteToFileWithTime("rr_reuse_infl.txt", gf);
			//all in first
			_RebuildRRIndices();
			spread = _RunGreedy(top, seeds, est_spread);
			pair<int, int> seed;
			seeds4.resize(top, seed);
			for (int i = 0; i < seeds.size(); i++)
			{
				seeds4[i].first = seeds[i];
				seeds4[i].second = 0;
				est_spread4.push_back(Weight_iter(weight_mode, 1)*est_spread[i]);
			}
			_SetResults1(seeds4, est_spread4);
			WriteToFileWithTime("rr_first_infl.txt", gf);
			_SetResults(seeds, est_spread);

		}
		else {
			LARGE_INT64 nNewSamples = LARGE_INT64(theta - table.size() + 1);
			nNewSamples = ((nNewSamples > 0) ? nNewSamples : 0);
			cout << "  Original IMM without regenerating RR sets, # new RR sets needed = " << nNewSamples << endl;
			if (table.size() < theta) {
				// generate samples
				_AddRRSimulation(nNewSamples, cascade, table, targets);
				_RebuildRRIndices();
				spread = _RunGreedy(top, seeds, est_spread);
				_SetResults(seeds, est_spread);
			}
		}


		cout << " spread0 = " << (spread*Weight_iter(weight_mode, 1)+1.0)*(kb_0/kp_0+1.0)-1.0 << endl;  //first
		cout << " spread1 = " << (spread1+1) * (kb_0 / kp_0 + 1.0) - 1.0 << endl;  //uniform
		cout << " spread2 = " << (spread2+1) * (kb_0 / kp_0 + 1.0) - 1.0 << endl;
		cout << " spread3 = " << (spread3+1) * (kb_0 / kp_0 + 1.0) - 1.0 << endl;  //random round
		stept.SetTimeEvent("step_end");
		steptimers.push_back(stept);
	}
	pctimer.SetTimeEvent("end");

	// Write results to file:
	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i = 0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);
	//fclose(out);

	//WriteToFile(file, gf);

	FILE* timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%g\n", pctimer.TimeSpan("start", "end"));
	for (size_t t = 0; t < steptimers.size(); t++) {
		if (t != steptimers.size() - 1) {
			fprintf(timetmpfile, "  Step%d: %g\n", t + 1, steptimers[t].TimeSpan("step_start", "step_end"));
		}
		else {
			fprintf(timetmpfile, "  Final Step: %g\n", steptimers[t].TimeSpan("step_start", "step_end"));
		}
	}
	fclose(timetmpfile);
}

void IMM::Build(graph_type& gf, int k, int time, cascade_type& cascade, double eps /* = 0.1 */, double ell /* = 1.0 */, int mode /* = 0 */)
{
	InitializeConcurrent();
	float ratio = 0.0;
	EventTimer pctimer;
	vector<EventTimer> steptimers;
	pctimer.SetTimeEvent("start");

	n = gf.GetN();
	m = gf.GetM();
	
	top = time;
	d.resize(k, 0.0);
	list.resize(k, 0);
	pair<int, int> listWT;
	listWithTime.resize(k, listWT);
	cascade.Build(gf);

	table.clear();
	tableWithTime.clear();
	targets.clear();

	double sum_weight = 0.0;
	for (int i = 0; i < time; i++)
	{
		sum_weight += Weight_iter(weight_mode, i+1);
	}
	double epsprime = eps * sqrt(2.0); // eps'
	double LB = 1.0 ;   // lower bound
	double maxRounds = max(max(log2((double)n), 1.0) - 1.0, 1.0); //x_i  edited
	

	ell = ell + log(2) / log(n); // Original IMM has failure probability 2/n^ell, so use this transformation to 
	// make the failure probability 1/n^ell

	// first a few passes of sampling
	double spread = 0.0;
	for (int r = 1; r < maxRounds; r++) {
		EventTimer stept;
		stept.SetTimeEvent("step_start");

		std::cout << "  Step" << r << ":" <<endl;
		double x = max(double(n) / pow(2, r), 1.0 );
		//double x = max(double(n) * sum_weight / pow(2, r), 1.0 * sum_weight);
		//double theta = Weight_iter(weight_mode, 1) * LambdaPrime(epsprime, k, ell, n, time) / x;// / 10;// *sum_weight;
		double theta = LambdaPrime(epsprime, k, ell, n, time) / x;

		// select nodes
		vector<int> seeds;
		vector< pair< int, int > > seedsWithTime;
		vector<double> est_spread;
		LARGE_INT64 nNewSamples = LARGE_INT64(theta - tableWithTime.size() + 1);
		if (tableWithTime.size() < theta) {
			// generate samples
			_AddRRSimulation1(nNewSamples, cascade, tableWithTime, targets, time);
			_RebuildRRIndicesWithTime();
			spread = _RunGreedy1(k, seedsWithTime, est_spread);
			_SetResults1(seedsWithTime, est_spread);
		}
		//spread = (spread + 1.0) * (1.0+kb_0/kp_0) - 1.0;
		std::cout << " spread = " << spread << "\t round = " << ((nNewSamples > 0) ? nNewSamples : 0) << endl;

		stept.SetTimeEvent("step_end");
		steptimers.push_back(stept);
		
		// check influence and set lower bound
		spread = spread;// / sum_weight;

		if (spread/Weight_iter(weight_mode,1) >= ((1.0 + epsprime) * x)) {
			LB =  spread / (1.0 + epsprime);
			break;
		}
		/*if (spread >= x) {
			LB = spread;
			break;
		}*/
	}

	// final pass of sampling
	{
		EventTimer stept;
		stept.SetTimeEvent("step_start");
		cout << "  Estimated Lower bound: " << LB  << endl;
		cout << "  Final Step:" << endl;
		double theta = Weight_iter(weight_mode, 1)*LambdaStar(eps, k, ell, n, time) / LB;// / 10;// *sum_weight;

		vector<int> seeds;
		vector< pair< int, int > > seedsWithTime;
		vector<double> est_spread; 



		if (mode == 1) {
			LARGE_INT64 nNewSamples = LARGE_INT64(theta + 1);
			cout << "  IMM Workaround 1 --- Regenerating RR sets, # RR sets = " << nNewSamples  << endl;
			tableWithTime.clear();
			targets.clear();
			// generate samples
			_AddRRSimulation1(nNewSamples, cascade, tableWithTime, targets, time);
			_RebuildRRIndicesWithTime();
			//spread = _RunGreedyTest(k, seedsWithTime, est_spread, ratio);
			spread = _RunGreedy1(k, seedsWithTime, est_spread);
			_SetResults1(seedsWithTime, est_spread);

		}
		else {
			LARGE_INT64 nNewSamples = LARGE_INT64(theta - table.size() + 1);
			nNewSamples = ((nNewSamples > 0) ? nNewSamples : 0);
			cout << "  Original IMM without regenerating RR sets, # new RR sets needed = " << nNewSamples << endl;
			if (table.size() < theta) {
				// generate samples
				_AddRRSimulation(nNewSamples, cascade, table, targets);
				_RebuildRRIndicesWithTime();
				//spread = _RunGreedyTest(k, seedsWithTime, est_spread, ratio);
				spread = _RunGreedy1(k, seedsWithTime, est_spread);
				_SetResults1(seedsWithTime, est_spread);
			}
		}

		cout << " spread(final ratio) = " << (spread+1.0)*(1.0+kb_0/kp_0)-1.0 << endl;

		cout << " Approximate ratio = " << ratio << endl;

		stept.SetTimeEvent("step_end");
		steptimers.push_back(stept);
	}
	pctimer.SetTimeEvent("end");

	// Write results to file:
	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i = 0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);
	//fclose(out);

	WriteToFileWithTime(file, gf);

	FILE *timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%g\n", pctimer.TimeSpan("start", "end"));
	for (size_t t = 0; t < steptimers.size(); t++) {
		if (t != steptimers.size() - 1) {
			fprintf(timetmpfile, "  Step%d: %g\n", t + 1, steptimers[t].TimeSpan("step_start", "step_end"));
		}
		else {
			fprintf(timetmpfile, "  Final Step: %g\n", steptimers[t].TimeSpan("step_start", "step_end"));
		}
	}
	fclose(timetmpfile);
}

double IMM::LambdaPrime(double epsprime, int k, double ell, int n, int time)  
{
	static double SMALL_DOUBLE = 1e-16;
	double cst = (2.0 + 2.0 / 3.0 * epsprime) / max(epsprime * epsprime, SMALL_DOUBLE);
	double part2 = time * LogNChooseK(n, k);
	part2 += ell * log(max((double)n, 1.0));  
	double sum_weight = 0.0;
	part2 += log(max(log2(max((double)n, 1.0 )), 1.0 ));
	part2 += log(4); //edited
	// calc \lambda'
	double lambda = cst * part2 * n;  //edited
	return lambda;
}

double IMM::LambdaStar(double eps, int k, double ell, int n, int time)
{
	static double SMALL_DOUBLE = 1e-16;
	static double APPRO_RATIO = (1.0 - 1.0 / 2);
	double logsum = ell * log(n) + log(4);
	// calc \alpha and \beta
	double alpha = sqrt(max(logsum, SMALL_DOUBLE));
	double beta = sqrt(APPRO_RATIO * (LogNChooseK(n * time, k) + alpha * alpha));
	// calc \lambda*
	double lambda = 2.0 * n / max(pow(eps, 2), SMALL_DOUBLE);
	lambda = lambda * pow(APPRO_RATIO * alpha + beta, 2);
	return lambda;
}

double IMM::LambdaPrimeOrigin(double epsprime, int k, double ell, int n)
{
	static double SMALL_DOUBLE = 1e-16;
	double cst = (2.0 + 2.0 / 3.0 * epsprime) / max(epsprime * epsprime, SMALL_DOUBLE);
	double part2 = LogNChooseK(n, k);
	part2 += ell * log(max((double)n, 1.0));
	part2 += log(max(log2(max((double)n, 1.0)), 1.0));
	// calc \lambda'
	double lambda = cst * part2 * n;
	return lambda;
}

double IMM::LambdaStarOrigin(double eps, int k, double ell, int n)
{
	static double SMALL_DOUBLE = 1e-16;
	static double APPRO_RATIO = (1.0 - 1.0 / exp(1));
	double logsum = ell * log(n) + log(2);
	// calc \alpha and \beta
	double alpha = sqrt(max(logsum, SMALL_DOUBLE));
	double beta = sqrt(APPRO_RATIO * (LogNChooseK(n, k) + logsum));
	// calc \lambda*
	double lambda = 2.0 * n / max(pow(eps, 2), SMALL_DOUBLE);
	lambda = lambda * pow(APPRO_RATIO * alpha + beta, 2);
	return lambda;
}










