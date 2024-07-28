#include "option.h"
#include "graph.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>
#include <omp.h>

using namespace std;

// compute eps-delta approx of INF_F using generalized SRA
float estimateFakeInfluence(Graph &g, int n, double epsilon_prime, double delta_prime) {
	int b = n - g.getNumFakeSeeds();
	double eps = epsilon_prime * (1 - (epsilon_prime * b)/((2 + 2.0/3.0 * epsilon_prime) * log(2.0 / delta_prime) * b));
	long long int gamma = (1 + epsilon_prime) * (2 + 2.0/3.0 * eps) * log(2.0 / delta_prime) * (1.0 / (eps*eps)) * b;

	long long int counter = 0;
	long long int global_running_total = 0;
    int c = 5000;

    omp_set_num_threads(14);
    #pragma omp parallel
    {
        vector<bool> visit(n,false);
        vector<int> visit_index(n,0);
        vector<int> dist(n,1000);
        priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > pq;
        long long int running_total;

        while (running_total < gamma) {
            running_total = 0;
            for (int i = 0; i < c; i++) {
                running_total += g.generateFakeInfluenceSample(visit, visit_index, dist, pq);
            }
            #pragma omp critical
            {
                global_running_total += running_total;
                counter += c;
            }
        }
	}

	return (float)global_running_total / counter;
}

// compute eps-delta approx of INF_F using generalized SRA
float estimateFakeInfluenceParallel(Graph &g, int n) {
	long long int counter = 0;
	long long int global_running_total = 0;
	int c = 100;

	#pragma omp parallel
	{
        vector<bool> visit(n,false);
        vector<int> visit_index(n,0);
        vector<double> threshold(n,2.0);
        vector<int> dist(n,100000);
        priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > pq;

		long long int running_total;
		while (counter < 10000) {
			running_total = 0;
//            #pragma omp parallel for
			for (int i = 0; i < c; i++) {
                long long int oneRes = g.generateFakeInfluenceSample(visit, visit_index,  dist, pq);
                running_total += oneRes;
//                if(oneRes > 5)
//                    cout << oneRes << endl;
			}
			#pragma omp critical
			{
				global_running_total += running_total;
				counter += c;
			}
		}
	}
    cout << global_running_total << " " << counter << endl;
	return (float)global_running_total / counter;
}

int main(int argc, char ** argv)
{
	srand(time(NULL));
	
	OptionParser op(argc, argv);
	if (!op.validCheck()){
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}

	char * inFile = op.getPara("-i");
	if (inFile == NULL){
		inFile = (char*)"network";
	}

	char * outFile = op.getPara("-o");
	if (outFile == NULL){
		outFile = (char*)"fake.inf";
	}

	char * fakeSeedsFile = op.getPara("-fakeseeds");
	if (fakeSeedsFile == NULL){
		fakeSeedsFile = (char*)"fake.seeds";
	}

	char * tmp = op.getPara("-epsilon");
	float epsilon = 0.1;
	if (tmp != NULL){
		epsilon = atof(tmp);
	}

	float ew = -1.0;
	tmp = op.getPara("-ew");
	if (tmp != NULL){
		ew = atof(tmp);
	}
	bool fixed = (ew < 0.0) ? false : true;

	Graph g(8,32);
	g.readGraph(inFile, fixed, ew);
	g.readFakeSeeds(fakeSeedsFile);

	int n = g.getNumNodes();

	float delta = 1.0/n;
	tmp = op.getPara("-delta");
    if (tmp != NULL){
    	delta = atof(tmp);
    }

    cout << "\n*******************" << endl;
	cout << "\tSTART" << endl;
	cout << "*******************\n" << endl;


    double start = omp_get_wtime();
    float inf_f = estimateFakeInfluenceParallel(g, n);
    cout << "Time to estimate INF_F: " << omp_get_wtime()-start << "s" << endl;
    cout << "INF_F = " << inf_f << endl;

    cout << "fake seed set: ";
    const vector<int> &fs = g.getFakeSeeds();
	for (unsigned int s = 0; s < g.getNumFakeSeeds(); s++) {
		cout << fs[s] << " ";
	}
	cout << endl;

	ofstream out(outFile, ios::app);
	out << inf_f << endl;
	out.close();

	cout << "\n*******************" << endl;
	cout << "\tALL DONE" << endl;
	cout << "*******************" << endl;

	//
	// ALL DONE
	//
}