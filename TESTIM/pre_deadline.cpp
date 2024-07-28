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


/*
* read fake seeds
*/
float readFakeInfluence(const char* filename)
{
    float fi;
    ifstream in(filename);
    in >> fi;
    in.close();
    return fi;
}


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
    float epsilon = 0.2;
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

    int k = 20;
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
    bool fixed = (ew < 0.0) ? false : true;

    cout << "\n*******************" << endl;
    cout << "\tSTART" << endl;
    cout << "*******************\n" << endl;

    Graph g(t,T);
    g.readGraph(inFile, fixed, ew);
    g.readFakeSeeds(fakeSeedsFile);
    g.init_visit_thresh_hold();
    int n = g.getNumNodes();

    float delta = 1.0 / n;
    tmp = op.getPara("-delta");
    if (tmp != NULL){
        delta = atof(tmp);
    }

    float precision = 1-1/exp(1);
    tmp = op.getPara("-precision");
    if (tmp != NULL){
        precision = atof(tmp);
    }



    double start_total = omp_get_wtime();
    double start = omp_get_wtime();


    cout << "fake seed set: ";
    const vi &fs = g.getFakeSeeds();
    for (unsigned int s = 0; s < g.getNumFakeSeeds(); s++) {
        cout << fs[s] << " ";
    }
    cout << endl;

    const int try_times = 2;
    start = omp_get_wtime();
    double e = exp(1);
//    double primary_bound = (8 * e * e + e * 4 + 2 * e * e * epsilon) / (epsilon * epsilon * e * e + e * 4 * epsilon + 4) * n * (log(6 / delta)+ n * log(2)) / epsilon / epsilon;
    double primary_bound = (8  + 2 * epsilon)  * n * (log(6 / delta)+ n * log(2)) / epsilon / epsilon;
    float ratio;
    int Final_MIN_left = 1e5;
    float final_ratio;
    for (int try_time = 0; try_time < try_times; try_time++) {
        int MIN_left = 18;
        HyperGraph coll_one(n);
        coll_one.hyperGT.clear();
        vector<int> RR_block_set;
        long long cur_samples = g.getNumNodes();
        long long samples_done = 0;
        unsigned int count_cnt=0;
        vector<int> is_RB_set_last;
        vector<int> set_dist_last;
        int usefulTR = 0, usefulBR = 0;
        float inf_number_record = 0.0;
        bool change_again = false;
        vector<bool>mark_RB_set, mark_dist_last;
        float last_ratio = 0.0;
        while(MIN_left = MIN_left + 3){
            coll_one.hyperG.clear();
            coll_one.hyperGB.clear();
            for(int i=0; i<n; i++)
                coll_one.hyperG.push_back(vi());
            for(int i=0; i<n; i++)
                coll_one.hyperGB.push_back(vi());
            g.setDeadline(MIN_left);
            if(samples_done != 0)
                change_again = true;
            do{
                if(!change_again)
                    count_cnt++;
                else{
                    change_again = false;
                    usefulTR = cur_samples;
                    usefulBR = 0;
                }
                auto bd_return = addHyperedgeParallel(g, coll_one, cur_samples, samples_done, RR_block_set, count_cnt, delta, epsilon, is_RB_set_last, set_dist_last, usefulTR, usefulBR, k, time_generate_flag, lamda, p);
                samples_done = cur_samples;
                if(bd_return.first){
                    inf_number_record = bd_return.second - coll_one.EstPro1;
                    break;
                }
                cur_samples *= 2;
                if(cur_samples  > (long long)primary_bound * (1 + coll_one.eps1) / coll_one.EstPro2){
                    inf_number_record = bd_return.second - coll_one.EstPro1;
                    cout << "It is because of sample too small\n";
                }
            }while(cur_samples  < (long long)primary_bound * (1 + coll_one.eps1) / coll_one.EstPro2);
            ratio = inf_number_record * 1.0 / n;
            cout << ratio << endl;
            if(ratio - last_ratio >= 0.005)
                last_ratio = ratio;
            else{
                if(MIN_left < Final_MIN_left){
                    Final_MIN_left = MIN_left;
                    final_ratio = last_ratio;
                }
                break;
            }
        }
        cout << "Now try the " << try_time << " times\n" ;
        cout<<"The max deadline is " << MIN_left << "\n" ;
        cout<<"The max ratio is " << last_ratio << "\n" ;
    }
    double time_total = omp_get_wtime()-start_total;
    ofstream out(outFile);
    out << Final_MIN_left << endl;
    out << final_ratio << endl;
    out << "The average time cost is " << time_total / try_times << "s";
    out.close();
    cout << "\nTotal Time: " << time_total << "s" << endl;
    cout << "\n*******************" << endl;
    cout << "\tALL DONE" << endl;
    cout << "*******************" << endl;

    //
    // ALL DONE
    //
}


