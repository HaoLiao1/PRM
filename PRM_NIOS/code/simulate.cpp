#include "simulate.h"



Simu::Simu() 
{
	file = "";
	simuSize = 0;
	simuIterNum = NUM_ITER;
}

Simu::Simu(std::vector<int>& seeds, 
			std::string filename /*= ""*/, 
			int simulateSize /*= SET_SIZE*/, 
			int simuIterNum /*= NUM_ITER*/)
{
	this->file = filename;
	this->seeds = seeds;
	this->simuSize = ((int)seeds.size() < simulateSize) ? (int)seeds.size() : simulateSize;
	this->simuIterNum = simuIterNum;
}

Simu::Simu(int* seeds, 
			std::string filename /*= ""*/, 
			int simulateSize /*= SET_SIZE*/, 
			int simuIterNum /*= NUM_ITER*/)
{
	this->file = filename;
	for (int i = 0; i < simulateSize; i++) {
		this->seeds.push_back(seeds[i]);
	}
	this->simuSize = simulateSize;
	this->simuIterNum = simuIterNum;
}


Simu::Simu(int(*GetNode)(int i), 
			std::string filename /*= ""*/, 
			int simulateSize /*= SET_SIZE*/, 
			int simuIterNum /*= NUM_ITER*/)
{
	this->file = filename;
	for (int i = 0; i < simulateSize; i++) {
		seeds.push_back(GetNode(i));
	}
	this->simuSize = simulateSize;
	this->simuIterNum = simuIterNum;
}

double Simu::toSimulate(ICascade& cascade)
{
	std::ofstream out;
	if (isWriteFile()) {
		out.open(file.c_str());
	}

	int* set = new int[simuSize];
	double spread = 0.0;
	for (int t = 0; t < simuSize; t++)
	{
		set[t] = seeds[t];
		spread = cascade.Run(simuIterNum, t + 1, set);

		printf("%02d \t %10g\n", t + 1, spread);
		if (isWriteFile()) {
			out << (t + 1) << "\t" << spread << std::endl;
		}
	}

	if (isWriteFile()) {
		out.close();
	}
	SAFE_DELETE_ARRAY(set);
	return spread;
}


double Simu::toSimulateOnce(ICascade& cacade)
{
	int* set = new int[simuSize];
	for (int t = 0; t < simuSize; t++)
	{
		set[t] = seeds[t];
	}
	double spread = cacade.Run(simuIterNum, simuSize, set);
	SAFE_DELETE_ARRAY(set);
	return spread;
}

double Simu::toSimulateOnceFile(ICascade& cacade)
{
	std::ofstream out;
	if (isWriteFile()) {
		out.open(file.c_str());
	}

	int* set = new int[simuSize];
	for (int t = 0; t < simuSize; t++)
	{
		set[t] = seeds[t];
	}
	double spread = cacade.Run(simuIterNum, simuSize, set);

	printf("seed set size = %d, spread = %10g\n", simuSize, spread);
	if (isWriteFile()) {
		out << simuSize << "\t" << spread << std::endl;
	}

	SAFE_DELETE_ARRAY(set);
	return spread;
}

std::vector<bool> Simu::toSimulatePRM2(GeneralCascade& cascade)
{
	std::ofstream out;
	if (isWriteFile()) {
		out.open(file.c_str());
	}

	int* set = new int[simuSize];
	for (int t = 0; t < simuSize; t++)
	{
		set[t] = seeds[t];
	}
	std::vector<bool> active;
	double spread = cascade.Run(simuIterNum, simuSize, set, active);



	printf("seed set size = %d, spread = %10g\n", simuSize, spread);
	if (isWriteFile()) {
		out << simuSize << "\t" << spread << std::endl;
	}

	SAFE_DELETE_ARRAY(set);
	return active;
}
