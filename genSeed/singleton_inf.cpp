#include "option.h"
#include "graph.h"
#include <algorithm>
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

typedef pair<float,int> fi;
typedef vector<fi> vfi;

// estimate influence of singleton seeds using 10K monte carlo simulations
void estimateInfluence(Graph &g, int n, int r, vfi &inf) {
    omp_set_num_threads(16);
	#pragma omp parallel
	{
		long long int running_total;
		float influence;
        int cnt = 0;
		#pragma omp for
		for (int i = 0; i < n; i++) {
//            vector<bool> visit(n,false);
//            vector<int> visit_index(n,0);
//            vector<double> threshold(n,2.0);
//            vector<int> dist(n,100000);
			running_total = 0;
			for (int j = 0; j < r; j++) {
                vector<bool> visit(n,false);
                vector<int> visit_index(n,0);
                vector<double> threshold(n,2.0);
                vector<int> dist(n,100000);
				auto res = g.generateInfluenceSample(visit, visit_index, threshold, dist, i);
//                if(j == 0)
//                    cout << res << endl;
                if(res < 100000000){
                    running_total += res * r;
                    break;
                }
                else{
                    running_total += res;
//                    cout << res << endl;
                }

//                cout << running_total << endl;
			}
			influence = (float)running_total / r;
//            #pragma omp critical
//            {
//                cnt++;
//                if( cnt % 100 == 0)
//                    cout << cnt << endl;
//            }
			inf[i] = fi(influence, i);
		}
	}

	sort(inf.begin(), inf.end());
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
		outFile = (char*)"singleton.inf";
	}

	float ew = -1.0;
	char * tmp = op.getPara("-ew");
	if (tmp != NULL){
		ew = atof(tmp);
	}

    int t = 8;
    tmp = op.getPara("-t");
    if (tmp != NULL){
        t = atoi(tmp);
    }

    int T = 32;
    tmp = op.getPara("-T");
    if (tmp != NULL){
        T = atoi(tmp);
    }

    int r = 1000;
    tmp = op.getPara("-r");
    if (tmp != NULL){
        r = atoi(tmp);
    }

	bool fixed = (ew < 0.0) ? false : true;

	Graph g(t,T);
	g.readGraph(inFile, fixed, ew);
    g.init_visit_thresh_hold();

	int n = g.getNumNodes();

	vfi inf(n);

    double start = omp_get_wtime();
    estimateInfluence(g, n, r, inf);
    cout << "Time to estimate singleton influences: " << (omp_get_wtime()-start) << "s" << endl;

    ofstream of(outFile);
	for (int i = n-1; i >= 0; i--){
		of << inf[i].second << "\t" << inf[i].first << endl;
	}
	of.close();

	//
	// ALL DONE
	//
}