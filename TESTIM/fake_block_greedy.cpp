#include "option.h"
#include "graph.h"
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <fstream>
#include <cstring>
#include <omp.h>

using namespace std;

// generate a single forward monte carlo estimate of the influence of F/B
float blockGreedyParallel(Graph &g, int n, const vector<int> &block) {
    vector<bool> visit(n,false);
    int b = n - g.getNumFakeSeeds();
    const vector<int> &fs = g.getFakeSeeds();
    long long int counter = 0;
    long long int global_running_total = 0;
    int c = 100;
    for(int i=0 ; i < g.getNumFakeSeeds(); i++){
        visit[fs[i]] = true;
    }
    #pragma omp parallel
    {
        vector<bool> visit(n,false);
        vector<int> visit_index(n,0);
        vector<double> threshold(n,2.0);
        vector<int> dist(n,100000);
        priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > pq;

        long long int running_total;
        while (counter < 1000) {
            running_total = 0;
            #pragma omp parallel for
            for (int i = 0; i < c; i++) {
//                long long int oneRes = g.generateBlockSampleGreedy(visit, visit_index, threshold, dist, block, pq);
                long long int oneRes = g.liveEdgeSampleGreedy(visit, visit_index, dist, block, pq);
                running_total += oneRes;
//                if(oneRes > 1000)
//                    cout << oneRes << endl;
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

//greedily choose budget node to block using monte carlo
void GreedyChooseBlock(Graph &g, int n, vector<int> &block, int budget, vector<double> &feedback, vector<bool>not_inf){
    int i, j;
    const vector<int> &fs = g.getFakeSeeds();
    vector<bool> notBlock(n,false);
    int fakeNum = (int)g.getNumFakeSeeds();
    for (i = 0; i < fakeNum; i++)
        notBlock[fs[i]] = true;
    for (i = 0; i < budget; i++){
        int x = -1;
        double MIN = 1e6;
        for (j = 0; j < n; j++){
            if (notBlock[j] == false && not_inf[j]==false){
                block.emplace_back(j);
                float res = blockGreedyParallel(g, n ,block);
                if (res < MIN){
                    MIN = res;
                    x = j;
                }
                block.pop_back();
            }
            if(j%10==0)
                cout<<j<<" "<<MIN<<endl;
        }
        cout << "now block " << x << " inf = " << MIN << endl;
        block.emplace_back(x);
        feedback.emplace_back(MIN);
        notBlock[x] = true;
    }
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
        outFile = (char*)"block.txt";
    }

    int budget = 1;
    char * tmp = op.getPara("-k");
    if (tmp != NULL){
        budget = atoi(tmp);
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

    char * fakeSeedsFile = op.getPara("-fakeseeds");
    if (fakeSeedsFile == NULL){
        fakeSeedsFile = (char*)"fake.seeds";
    }

    char * notInfFile = op.getPara("-notInf");
    if (notInfFile == NULL){
        notInfFile = (char*)"./WCModel/notInf.txt";
    }

    float ew = -1.0;
    tmp = op.getPara("-ew");
    if (tmp != NULL){
        ew = atof(tmp);
    }
    bool fixed = (ew < 0.0) ? false : true;

    Graph g(t,T);
    g.readGraph(inFile, fixed, ew);
    g.readFakeSeeds(fakeSeedsFile);
    int n = g.getNumNodes();
    vector<bool> not_inf(n,false);
    ifstream not_inf_file(notInfFile);
    unsigned int not_inf_size;
    int not_inf_node;
    not_inf_file >> not_inf_size;
    for (unsigned int i = 0; i < not_inf_size; i++){
        not_inf_file >> not_inf_node;
        not_inf[not_inf_node] = true;
    }
    vector<int> block;
    vector<double> feedback;
    cout << "\n*******************" << endl;
    cout << "\tSTART" << endl;
    cout << "*******************\n" << endl;

    double start = omp_get_wtime();
    GreedyChooseBlock(g, n ,block, budget, feedback, not_inf);
    double time_mark = omp_get_wtime()-start;
    cout << "Time to Block Node: " << time_mark << "s" << endl;

    cout << "fake seed set: ";
    const vector<int> &fs = g.getFakeSeeds();
    for (unsigned int s = 0; s < g.getNumFakeSeeds(); s++) {
        cout << fs[s] << " ";
    }
    cout << endl;

    ofstream out(outFile, ios::app);
    out << "block budget = " << budget << endl;
    out << "The program costs " << time_mark << "s to finish" << endl;
    for (int o = 0; o< budget; o++){
        out << block[o] << "  " << feedback[o] << endl;
    }
    out.close();

    cout << "\n*******************" << endl;
    cout << "\tALL DONE" << endl;
    cout << "*******************" << endl;

    //
    // ALL DONE
    //
}