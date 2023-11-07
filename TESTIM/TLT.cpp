#include "option.h"
#include "hypergraph.hpp"
#include "graph.h"
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cstring>
#include <omp.h>

using namespace std;


//// -i ../../datasets/DBLP/network -fakeseeds ../../datasets/DBLP/fake_20.seeds -k 10

int main(int argc, char ** argv)
{
	srand(time(NULL));
    bool time_generate_flag = false;
	
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
		outFile = (char*)"results.txt";
	}

	char * fakeSeedsFile = op.getPara("-fakeseeds");
	if (fakeSeedsFile == NULL){
		fakeSeedsFile = (char*)"fake.seeds";
	}



	char * tmp = op.getPara("-epsilon");
	float epsilon = 0.3;
	if (tmp != NULL){
		epsilon = atof(tmp);
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

	int k = 1;
	tmp = op.getPara("-k");
	if (tmp != NULL){
		k = atoi(tmp);
	}

    int lamda= 1;
    tmp = op.getPara("-lamda");
    if (tmp != NULL){
        lamda = atoi(tmp);
        time_generate_flag = false;
    }

    float p = 0.3;
    tmp = op.getPara("-p");
    if (tmp != NULL){
        p = atof(tmp);
        time_generate_flag = true;
    }

	float ew = -1.0;
	tmp = op.getPara("-ew");
	if (tmp != NULL){
		ew = atof(tmp);
	}
	bool fixed = (ew < 0.0) ? false : true;

	cout << "\n*******************" << endl;
	cout << "\tSTART" << endl;
	cout << "*******************\n" << endl;

	Graph g(t,T);
	g.readGraph(inFile, fixed, ew);
	g.readFakeSeeds(fakeSeedsFile);
    g.init_visit_thresh_hold();
	int n = g.getNumNodes();


	float delta = 1.0/n;
	tmp = op.getPara("-delta");
    if (tmp != NULL){
    	delta = atof(tmp);
    }


    double start_total = omp_get_wtime();


    const vi &fs = g.getFakeSeeds();
    cout << "fake seed set ["<< fs.size()<<"]: ";
	for (unsigned int s = 0; s < g.getNumFakeSeeds(); s++) {
		cout << fs[s] << " ";
	}
	cout << endl;
    /*
     * test code begin
     */
//    vector<bool> isSeedNeighbor(n,false);
//    int isFakeNeighborCount = 0;
//    for(unsigned int s = 0; s < g.getNumFakeSeeds(); s++) {
//        auto fkn = g.getOutNeighbours(fs[s]);
//        for(unsigned int f = 0; f < g.getOutDegree(fs[s]); f++){
//            isSeedNeighbor[fkn[f]] = true;
//        }
//        cout<<"The degree of "<<fs[s]<<"is"<<g.getOutDegree(fs[s]);
//    }
    /*
     * test code end
     */
    const int try_times = 1;
    double final_inf = 0.0, final_save = 0.0, final_reduction = 0.0;
    long long cur_samples;
    vector<int> RR_block_set;
//    double primary_bound = (8 * e * e + e * 4 + 2 * e * e * epsilon) / (epsilon * epsilon * e * e + e * 4 * epsilon + 4) * n * (log(6 / delta)+ n * log(2)) / epsilon / epsilon;
    double primary_bound = (8  + 2 * epsilon)  * n * (log(6 / delta)+ n * log(2)) / epsilon / epsilon;
    priority_queue<ii,vector<ii>,less<ii>> pq;
    map<int,int>cnt_map;
    for (int try_time = 0; try_time < try_times; try_time++) {
        HyperGraph coll_one(n);
        coll_one.hyperGT.clear();
        coll_one.hyperG.clear();
        coll_one.hyperGB.clear();
        for(int i=0; i<n; i++)
            coll_one.hyperG.push_back(vi());
        for(int i=0; i<n; i++)
            coll_one.hyperGB.push_back(vi());
        cur_samples = g.getNumNodes();
        long long samples_done = 0;
        unsigned int count_cnt=0;
        vector<int> is_RB_set_last;
        vector<int> set_dist_last;
        int usefulTR = 0, usefulBR = 0;
        do{
            cout << "================================  round " << count_cnt+1 << "  ================================" << endl;
            count_cnt++;

            //// begin generate rr sets and select blocking set
            auto bd_return = addHyperedgeParallel(g, coll_one, cur_samples, samples_done, RR_block_set, count_cnt, delta, epsilon, is_RB_set_last, set_dist_last, usefulTR, usefulBR, k, time_generate_flag, lamda, p);
            cout << "[R1] current misinformation INF (after blocking): " << bd_return.second - coll_one.EstPro1 << endl;

            //// bd_return.first = [true]: current result satisfy constrains
            if(bd_return.first){
                cout << "[R1] misinformation INF: "<< bd_return.second << endl;
                cout << "[R1] saved nodes: "<< coll_one.EstPro1 << endl;
                final_inf += bd_return.second;
                final_save += coll_one.EstPro1 ;
                cout << "Reduction ratio = " << coll_one.EstPro1 / bd_return.second << endl;
                final_reduction += coll_one.EstPro1 / bd_return.second;
                cout << "Running time: " << (omp_get_wtime()-start_total) / 60  << " min is gone\n";
                break;
            }
            //// bd_return.first = [false]: double the rr sets
            samples_done = cur_samples;
            cur_samples *= 2;
            if(usefulBR  > (long long)primary_bound * (1 + coll_one.eps1) / coll_one.EstPro2){
                cout << "[Exit] because rr set sample reach the maximal number\n";
            }
        }while(usefulBR  < (long long)primary_bound * (1 + coll_one.eps1) / coll_one.EstPro2);
        for(auto bs : RR_block_set){
            if(cnt_map.find(bs)!=cnt_map.end()){
                cnt_map[bs]++;
            }
            else
                cnt_map[bs] = 1;
        }

    }
    RR_block_set.clear();
    for(auto item : cnt_map){
        pq.push(make_pair(item.second,item.first));
    }
    for(int i = 0; i < k; i++){
        RR_block_set.emplace_back(pq.top().second);
        pq.pop();
    }
    double time_total = omp_get_wtime()-start_total;
    cout << "lamda = " << lamda << endl;

    ofstream out(outFile);
    out <<  cur_samples << "\t" << "RR set Number" << endl;
    out << "The final inf = " << final_inf / try_times << endl;
    out << "Finally saved " << final_save / try_times << "nodes" << endl;
    out << "Final inf(after block ) = " << (final_inf / try_times) - (final_save / try_times) << endl;
    out << "Reduction ratio = " << final_reduction / try_times << endl;
    out << "budget:" << k << '\t';
    for(int write_node = 0; write_node < k; write_node++){
        out << RR_block_set[write_node] << ' ';
    }
    out << '\n';
    out << "The average cost time is " << time_total / try_times << " s" << endl;
    out.close();


	cout << "\nTotal Time: " << time_total << "s" << endl;
	cout << "\n*******************" << endl;
	cout << "\tALL DONE" << endl;
	cout << "*******************" << endl;

	//
	// ALL DONE
	//
}