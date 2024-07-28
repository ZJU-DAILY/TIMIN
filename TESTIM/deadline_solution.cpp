#include "option.h"
#include "hypergraph.hpp"
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

    char * preFile = op.getPara("-preFile");
    if (preFile == NULL){
        preFile = (char*)"preFile";
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

    int k = 10;
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

    float alpha = 0.25;
    tmp = op.getPara("-alpha");
    if (tmp != NULL){
        alpha = atof(tmp);
    }

    bool nub = false;
    tmp = op.getPara("-nub");
    if (tmp != NULL){
        nub = atoi(tmp);
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

    float precision = 1-1/exp(1);
    tmp = op.getPara("-precision");
    if (tmp != NULL){
        precision = atof(tmp);
    }



    double start_total = omp_get_wtime();
    double start = omp_get_wtime();
    int max_time_right;
    ifstream preIn(preFile);
    preIn >> max_time_right;
    preIn.close();

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
    float final_ratio;
    int final_left = max_time_right + 1;
    for (int try_time = 0; try_time < try_times; try_time++) {
        int MIN_left = 1, MAX_right = max_time_right;
//        cout << "max_time_right= " << max_time_right << endl;
        float ratio;
        HyperGraph coll_one(n);
        coll_one.hyperGT.clear();
        vector<int> RR_block_set;
        long long cur_samples = g.getNumNodes();
        long long samples_done = 0;
        unsigned int count_cnt=0;
        vector<int> is_RB_set_last;
        vector<int> set_dist_last;
        int usefulTR = 0, usefulBR = 0;
        int mid;
        float inf_number_record = 0.0;
        bool change_again = false;
        vector<bool>mark_RB_set, mark_dist_last;
        while(MIN_left < MAX_right){
            coll_one.hyperG.clear();
            coll_one.hyperGB.clear();
            for(int i=0; i<n; i++)
                coll_one.hyperG.push_back(vi());
            for(int i=0; i<n; i++)
                coll_one.hyperGB.push_back(vi());
            mid = MIN_left + (MAX_right-MIN_left+1) / 2;
//            cout << "mid= " << mid << endl;
            g.setDeadline(mid);
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
            ratio = inf_number_record / n;
//            cout << ratio << endl;
            if(ratio > alpha){
                MAX_right = mid - 1;
            }
            else
                MIN_left =mid;
        }
        cout<<"The max deadline is " << MIN_left << "\n" ;
        while(ratio > alpha){
            coll_one.hyperG.clear();
            coll_one.hyperGB.clear();
            for(int i=0; i<n; i++)
                coll_one.hyperG.push_back(vi());
            for(int i=0; i<n; i++)
                coll_one.hyperGB.push_back(vi());
            MIN_left = MIN_left - 1;
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
            ratio = inf_number_record / n;
        }
        if(MIN_left < final_left && ratio < alpha){
            final_left = MIN_left;
            final_ratio = ratio;
        }
        cout<<"The " << try_time + 1 << " times max deadline is " << MIN_left <<"\n";
        cout<<"The " << try_time + 1 << "  times last ratio is " << ratio <<"\n";
    }
    double time_total = omp_get_wtime() - start_total;
    ofstream out(outFile);
    out << "The: max deadline is " << final_left << '\n';
    out << "The last ratio is " << final_ratio << '\n';
    out << "The average time is " << time_total / try_times << "s\n";
    cout << "\nAverage Time: " << time_total / try_times << "s" << endl;
    cout << "\n*******************" << endl;
    cout << "\tALL DONE" << endl;
    cout << "*******************" << endl;

    //
    // ALL DONE
    //
}

