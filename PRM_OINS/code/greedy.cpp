#include <cstdio>
#include <string>
#include <cstdlib>
#include "greedy.h"
#include "graph.h"
#include "common.h"
#include "cascade.h"
#include <cstring>
#include <set>
#include "event_timer.h"

using namespace std;



struct dCountComparator
{
public:
	vector<double>& counts;
	dCountComparator(vector<double>& c) :counts(c) {}
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

Greedy::Greedy() 
{
	file = "greedy.txt";
}

void Greedy::Build(IGraph& gf, int num, ICascade& cascade)
{
	n = gf.GetN();
	top = num;
	d.resize(num, 0);
	list.resize(num, 0);

	bool *used= new bool[n];
	memset(used, 0, sizeof(bool)*n);
	int* set = new int[num];  //seed set

	double old = 0.0;

	double *improve = new double[n];
	int *lastupdate = new int[n];
	int *heap = new int[n];
	for (int i=0; i<n; i++)
	{
		heap[i] = i;  //heap all seed set sorted by marginal
		lastupdate[i] = -1;
		improve[i] = (double)(n+1);//initialize largest
	}

	for (int i=0; i<top; i++)
	{
		int ccc = 0;
		//printf("%d\n",i);
		while (lastupdate[heap[0]] != i)
		{
			//printf("%d %d %d\n", i, heap[0], ccc);
			ccc++;
			lastupdate[heap[0]] = i;
			set[i] = heap[0];
			//printf("GreedyGC_SPM %d %d\n",heap[0],improve[heap[0]]);
			improve[heap[0]] = cascade.Run(NUM_ITER, i+1, set) - old;
			
			char tmpfilename[200];
			sprintf_s(tmpfilename, "tmp/%02d%05d.txt", i, heap[0]);
			//FILE *tmpfile;
			//fopen_s(&tmpfile, tmpfilename,"w");
			//fprintf(tmpfile, "%g\n", improve[heap[0]]); 
			//fclose(tmpfile);

			int x = 0;
			while (x*2+2<=n-i)
			{
				int newx=x*2+1;
				if ((newx+1<n-i) && (improve[heap[newx]]<improve[heap[newx+1]]))
					newx++;
				if (improve[heap[x]]<improve[heap[newx]])
				{
					int t=heap[x];
					heap[x] = heap[newx];
					heap[newx] = t;
					x = newx;
				}
				else
					break;
			}
		}

		used[heap[0]] = true;
		set[i] = heap[0];
		list[i] = heap[0];
		d[i] = improve[heap[0]];
		old+=d[i];

		//char bakname[200];
		//sprintf(bakname, "greedychoice%02d.txt", i+1);
		//FILE *bak = fopen(bakname, "w");
		//fprintf(bak, "%6d\t%d\t%g\n", i+1, heap[0], improve[heap[0]]);
		//fclose(bak);
		//printf("%d\t%g\n", i+1, improve[131]);

		heap[0] = heap[n-i-1];
		int x = 0;
		while (x*2+2<=n-i)//bug should-1
		{
			int newx=x*2+1;
			if ((newx+1<n-i) && (improve[heap[newx]]<improve[heap[newx+1]]))	//bug should-1
				newx++;
			if (improve[heap[x]]<improve[heap[newx]])
			{
				int t=heap[x];
				heap[x] = heap[newx];
				heap[newx] = t;
				x = newx;
			}
			else
				break;
		}
	}

	//FILE *out;
	//fopen_s(&out, file.c_str(), "w");
	//fprintf(out, "%d\n", top);
	//for (int i=0; i<top; i++)
	//	fprintf(out, "%d\t%g\n", list[i], d[i]);	//the nodes we want!
	//fclose(out);
	WriteToFile(file, gf);

	SAFE_DELETE_ARRAY(set);
	SAFE_DELETE_ARRAY(heap);
	SAFE_DELETE_ARRAY(lastupdate);
	SAFE_DELETE_ARRAY(improve);
	SAFE_DELETE_ARRAY(used);
}

void Greedy::Build(IGraph& gf, int k, ICascade& cascade, int dp, int dn, int a, int t){
	n = gf.GetN();
	top = k;
	d.resize(k, 0);
	list.resize(k, 0);

	pair<int, int> listWT;
	std::vector<std::pair<int, int>> listWithTime;
	listWithTime.resize(k, listWT);

	vector< pair< int, int > > outSeeds;
	vector<double> outEstSpread;

	vector<int*> setWithTime = vector<int*>(t);
	vector<int> seedNumber = vector<int>(top, 0);
	for (int i = 0; i < t; i++) {
		int* seedSet = new int[top];
		setWithTime[i] = seedSet;
	}

	vector<double> oldWithTime = vector<double>(t);
	for (int i = 0; i < t; i++) {
		double old = 0.0;
		oldWithTime[i] = old;
	}

	vector<vector<double>> imp = vector<vector<double>>(t, vector<double>(n));

	/*算法开始*/
	EventTimer pctimer;
	pctimer.SetTimeEvent("start");

	std::set<int> sourceSet;
	for (int i = 0; i < n; i++) {
		sourceSet.insert(i);
	}

	set<int> candidates(sourceSet);
	vector<std::set<int>> candidatesWithTime(top, candidates);

	double spread = 0;
	for (int iter = 0; iter < top; ++iter) {

		for (int i = 0; i < t; i++) {
			for (int j = 0; j < n; j++) {
				setWithTime[i][seedNumber[i]] = j;
				double inf = cascade.Run(500, seedNumber[i] + 1, setWithTime[i]);
				imp[i][j] = Weight_iter(dn, dp, a, i+1) * inf - oldWithTime[i];
			}
		}

		vector<dCountComparator> camp;
		for (int i = 0; i < top; i++) {
			dCountComparator comp(imp[i]);
			camp.push_back(comp);
		}

		vector<pair<pair<double, int>, int>> winner; //记录每个时间上最大的节点
		for (int i = 0; i < top; i++) {
			set<int>::const_iterator maxPtIn = max_element(candidatesWithTime[i].begin(), candidatesWithTime[i].end(), camp[i]);
			pair<pair<double, int>, int> maxSource;
			maxSource.first.second = *maxPtIn;
			maxSource.first.first = imp[i][*maxPtIn];
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
		assert(imp[winner[*maxPt].second][maxSourceWithTime.first] >= 0);

		// selected one node
		listWithTime[iter].first = maxSourceWithTime.first;
		listWithTime[iter].second = maxSourceWithTime.second + 1;

		// 记录用来计算边际效益
		setWithTime[maxSourceWithTime.second][seedNumber[maxSourceWithTime.second]] = maxSourceWithTime.first;
		seedNumber[maxSourceWithTime.second] += 1;

		// estimate spread
		spread = spread + imp[maxSourceWithTime.second][maxSourceWithTime.first];

		// if (iter==0)
		// {
		// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
		// }

		d[iter] = imp[maxSourceWithTime.second][maxSourceWithTime.first];
		oldWithTime[maxSourceWithTime.second] += imp[maxSourceWithTime.second][maxSourceWithTime.first];

		// 删除所有轮的节点
		for (int i = 0; i < t; i++) {
			candidatesWithTime[i].erase(maxSourceWithTime.first);
		}
		//candidatesWithTime[maxSourceWithTime.second].erase(maxSourceWithTime.first);   //只删除这一轮的节点
		//degreesWithTime[maxSourceWithTime.second][maxSourceWithTime.first] = -1;   似乎是多此一举

	}
	pctimer.SetTimeEvent("end");
	/*结束*/
	file = "rr_imm_infl.txt";
	WriteToFileWithTime(file, gf, listWithTime);

	std::string time_file = "time_rr_imm_infl.txt";
	FILE* timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%g\n", pctimer.TimeSpan("start", "end"));
	fclose(timetmpfile);
}

void Greedy::_Build(IGraph& gf, int k, GeneralCascade& cascade, int dp, int dn, int a, int t) {
	n = gf.GetN();
	top = k;
	d.resize(k, 0);
	list.resize(k, 0);

	pair<int, int> listWT;
	std::vector<std::pair<int, int>> listWithTime;
	listWithTime.resize(k, listWT);

	vector< pair< int, int > > outSeeds;
	vector<double> outEstSpread;

	vector<int*> setWithTime = vector<int*>(t);
	vector<int> seedNumber = vector<int>(top, 0);
	for (int i = 0; i < t; i++) {
		int* seedSet = new int[top];
		setWithTime[i] = seedSet;
	}

	vector<double> oldWithTime = vector<double>(t);
	for (int i = 0; i < t; i++) {
		double old = 0.0;
		oldWithTime[i] = old;
	}

	vector<vector<double>> imp = vector<vector<double>>(t, vector<double>(n));

	/*算法开始*/
	EventTimer pctimer;
	pctimer.SetTimeEvent("start");

	std::set<int> sourceSet;
	for (int i = 0; i < n; i++) {
		sourceSet.insert(i);
	}

	set<int> candidates(sourceSet);
	vector<std::set<int>> candidatesWithTime(top, candidates);

	double spread = 0;
	for (int iter = 0; iter < top; ++iter) {

		for (int i = 0; i < t; i++) {
			for (int j = 0; j < n; j++) {
				//模拟500次
				int round = 500;
				vector<int> all_count = vector<int>(round, 0);
#pragma omp parallel for
				for (int r = 0; r < round; r++) {
					// 计算在没有此种子前影响到的节点集合
					std::vector<bool> old_active = vector<bool>(n, false);
					for (int k = 0; k <= i; k++) {
						std::vector<bool> active;
						cascade.Run(1, seedNumber[k], setWithTime[k], active);
						for (int single_node = 0; single_node < active.size(); single_node++) {
							if (active[single_node]) {
								old_active[single_node] = true;
							}
						}
					}
					// 计算加上此种子后影响到的节点集合
					setWithTime[i][seedNumber[i]] = j;
					std::vector<bool> new_active = vector<bool>(n, false);
					int active_number = 0;
					for (int k = 0; k <= i; k++) {
						std::vector<bool> active;
						cascade.Run(1, seedNumber[k] + 1, setWithTime[k], active);
						for (int single_node = 0; single_node < active.size(); single_node++) {
							if (active[single_node]) {
								new_active[single_node] = true;
							}
						}
					}
					// 计算多出来的节点
					for (int single_node = 0; single_node < new_active.size(); single_node++) {
						if (new_active[single_node] && !old_active[single_node]) {
							all_count[r]++;
						}
					}
				}
				//求平均值
				double inf = 0;
				for (int r = 0; r < round; r++) {
					inf += all_count[r];
				}
				inf /= (double)round;
				imp[i][j] = Weight_iter(dn, dp, a, i + 1) * inf;
			}
		}

		vector<dCountComparator> camp;
		for (int i = 0; i < top; i++) {
			dCountComparator comp(imp[i]);
			camp.push_back(comp);
		}

		vector<pair<pair<double, int>, int>> winner; //记录每个时间上最大的节点
		for (int i = 0; i < top; i++) {
			set<int>::const_iterator maxPtIn = max_element(candidatesWithTime[i].begin(), candidatesWithTime[i].end(), camp[i]);
			pair<pair<double, int>, int> maxSource;
			maxSource.first.second = *maxPtIn;
			maxSource.first.first = imp[i][*maxPtIn];
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
		assert(imp[winner[*maxPt].second][maxSourceWithTime.first] >= 0);

		// selected one node
		listWithTime[iter].first = maxSourceWithTime.first;
		listWithTime[iter].second = maxSourceWithTime.second + 1;

		// 记录用来计算边际效益
		setWithTime[maxSourceWithTime.second][seedNumber[maxSourceWithTime.second]] = maxSourceWithTime.first;
		seedNumber[maxSourceWithTime.second] += 1;

		// estimate spread
		spread = spread + imp[maxSourceWithTime.second][maxSourceWithTime.first];

		// if (iter==0)
		// {
		// 	  cout << "hhh" << maxSource << (double)n * degrees[28] / table.size() << "hhh" << endl;
		// }

		d[iter] = imp[maxSourceWithTime.second][maxSourceWithTime.first];
		oldWithTime[maxSourceWithTime.second] += imp[maxSourceWithTime.second][maxSourceWithTime.first];

		// 删除所有轮的节点
		/*
		for (int i = 0; i < t; i++) {
			candidatesWithTime[i].erase(maxSourceWithTime.first);
		}
		*/
		candidatesWithTime[maxSourceWithTime.second].erase(maxSourceWithTime.first);   //只删除这一轮的节点
		//degreesWithTime[maxSourceWithTime.second][maxSourceWithTime.first] = -1;   似乎是多此一举

	}
	pctimer.SetTimeEvent("end");
	/*结束*/
	file = "rr_imm_infl.txt";
	WriteToFileWithTime(file, gf, listWithTime);

	std::string time_file = "time_rr_imm_infl.txt";
	FILE* timetmpfile;
	fopen_s(&timetmpfile, time_file.c_str(), "w");
	fprintf(timetmpfile, "%g\n", pctimer.TimeSpan("start", "end"));
	fclose(timetmpfile);
}


void Greedy::BuildRanking(IGraph& gf, int num, ICascade& cascade)
{
	this->n = gf.GetN();
	this->top = num;
	d.resize(num, 0);
	list.resize(num, 0);

	/*int* set = new int[num];
	double* inf = new double[top];*/
	//for(int i=0;i < top;i++)
	//{
	//	set[i] =0;
	//	inf[i] =0;
	//}

	vector<double> improve(n, 0.0);
	int tmp[1];
	for (int i=0; i<n; i++)
	{
		tmp[0] = i;
		improve[i] = cascade.Run(NUM_ITER, 1, tmp);
		if (improve[i] > list[top - 1])
		{
			d[top - 1] = improve[i];
			list[top - 1] = i;
			int j = top - 2;
			while(j>=0)
			{
				if (improve[i] > d[j])
				{
					int int_tmp = list[j];
					double inf_tmp = d[j];
					list[j] = i;
					d[j] = improve[i];
					list[j + 1] = int_tmp;
					d[j + 1] = inf_tmp;
				}
				else
					break;
				j--;
			}
		}
	}

	WriteToFile(file, gf);
	//SAFE_DELETE_ARRAY(set);
	//SAFE_DELETE_ARRAY(improve);
	//SAFE_DELETE_ARRAY(inf);
	//delete[] used;
}


void Greedy::BuildFromFile(IGraph& gf, const char* name)
{
	ReadFromFile(name, gf);
}

void Greedy::WriteToFileWithTime(const std::string& filename, IGraph& gf, std::vector<std::pair<int, int>> listWithTime)
{

	Write(filename, listWithTime, d, gf);
}

void Greedy::Write(const std::string& filename,
	const std::vector<std::pair<int, int>>& seeds,
	const std::vector<double>& infl, IGraph& gf)
{

	std::ofstream ft(filename.c_str());
	Write(ft, seeds, infl, gf);
}

void Greedy::Write(std::ostream& out,
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

double Greedy::Weight_iter(int dn, int dp, int a, int time)
{
	return 1.0 / (dn + dp + a * float(time));
}