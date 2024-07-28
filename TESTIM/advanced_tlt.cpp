#include "option.h"
#include "hypergraph.hpp"
#include "graph.h"
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cstring>
#include <omp.h>

using namespace std;

void bfs(const Graph &g, vector<vector<int>> &neighbor){
    vector<bool> visit(g.getNumNodes(),false);
    auto fs = g.getFakeSeeds();
    unsigned int fake_num = g.getNumFakeSeeds();
//    int cnt = 0;
    queue<pair<int,int>>q;
    for(unsigned int i=0; i < fake_num; i++){
        q.push(make_pair(fs[i],0));
    }
//    cout <<1 << endl;
    while(!q.empty()){
        auto f = q.front();
        q.pop();
        if(visit[f.first])
            continue;
        visit[f.first] = true;
        neighbor[f.second].emplace_back(f.first);
        auto child = g.getOutNeighbours(f.first);
        auto child_num = g.getOutDegree(f.first);
        for(unsigned int j = 0; j < child_num; j++){
            q.push(make_pair( child[j] , f.second + 1 ));
        }
//        cout << cnt++ << endl;
    }
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

//    char * RRFile = op.getPara("-ir");
//    if (RRFile == NULL){
//        RRFile = (char*)"rrSetNumber.txt";
//    }

    char * outFile = op.getPara("-o");
    if (outFile == NULL){
        outFile = (char*)"results.txt";
    }

    char * fakeSeedsFile = op.getPara("-fakeseeds");
    if (fakeSeedsFile == NULL){
        fakeSeedsFile = (char*)"fake.seeds";
    }



    char * tmp ;


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

    int rr = 1;
    tmp = op.getPara("-rr");
    if (tmp != NULL){
        rr = atoi(tmp);
    }

    float ew = -1.0;
    tmp = op.getPara("-ew");
    if (tmp != NULL){
        ew = atof(tmp);
    }
    bool fixed = (ew < 0.0) ? false : true;   //determine whether the edge weight is a fixed number or not

    cout << "\n*******************" << endl;
    cout << "\tSTART" << endl;
    cout << "*******************\n" << endl;

    Graph g(t,T);
    g.readGraph(inFile, fixed, ew);
    g.readFakeSeeds(fakeSeedsFile);
    int n = g.getNumNodes();
//    cout << n <<endl;

    float delta = 1.0/n;
    tmp = op.getPara("-delta");
    if (tmp != NULL){
        delta = atof(tmp);
    }


    double interval;
    double time_pre = 0.0;
    double start_total = omp_get_wtime();
    double start = omp_get_wtime();

    interval = omp_get_wtime()-start;
    time_pre += interval;

    const vi &fs = g.getFakeSeeds();
//    cout << 1 << endl;
    cout << "fake seed set ["<< fs.size()<<"]: ";
    for (unsigned int s = 0; s < g.getNumFakeSeeds(); s++) {
        cout << fs[s] << " ";
    }
    cout << endl;

    const int try_times = 1;
    start = omp_get_wtime();
    vector<vector<int>> neighbor;
    for(int i=0; i < 40; i++){
        neighbor.emplace_back(vector<int>());
    }

    cout << "=====================begin BFS =================\n";
    bfs(g, neighbor);
    cout << "=====================end BFS =================\n";

    int max_height = 0;
    for(int i=0; i < 40; i++){
        if(neighbor[i].size() == 0){
            max_height = i - 1;
            break;
        }
    }
    cout << "max height = " << max_height << endl;
    priority_queue<ii,vector<ii>,less<ii>> pq;
    map<int,int>cnt_map;
    vector<int> block_set;
    double final_save = 0.0;
    vector<bool> is_larger(max_height+1,false);
    vector<int>replace_num(max_height+1,0);
    vector<double>saved_vector(max_height+1, 0.0);
    float final_inf = 0.0, final_reduction = 0.0;
//    ifstream in2(RRFile);
//    long long initial_rr_set;
//    in2 >> initial_rr_set;
//    in2.close();
    for(int try_time = 0; try_time < try_times; try_time++){
        block_set.clear();
        std::fill(is_larger.begin(), is_larger.end(), false);
        std::fill(replace_num.begin(), replace_num.end(), 0);
        std::fill(saved_vector.begin(), saved_vector.end(), 0.0);
        HyperGraph coll_one(n);
        coll_one.hyperGT.clear();
        coll_one.hyperG.clear();
        coll_one.hyperGB.clear();
        for(int i=0; i<n; i++)
            coll_one.hyperG.push_back(vi());
        for(int i=0; i<n; i++)
            coll_one.hyperGB.push_back(vi());
        //pow(2,8) * 1e5
        long long rr_set_num = rr * n;
        cout << rr_set_num << endl;
        vector<int> block_set_coverage;
        double INF ;
        auto res = simpleBuild(g, coll_one, rr_set_num, INF, time_generate_flag, lamda, p);
        cout << "INF = " << INF << endl;
        final_inf += INF;
        int useful_BR = res.second;
        float max_saved;
        for(int i = 1; i <= max_height; i++){
            auto ret = max_cover(g, coll_one, block_set, k, neighbor[i], block_set_coverage, INF, useful_BR, i, replace_num);
            max_saved = ret.second ;
            saved_vector[i] = ret.second ;
            is_larger[i] = ret.first;
        }
        cout << "saved " << max_saved << endl;
        final_save += max_saved;
        final_reduction += (max_saved / INF);
        for(auto bs : block_set){
            if(cnt_map.find(bs)!=cnt_map.end()){
                cnt_map[bs]++;
            }
            else
                cnt_map[bs] = 1;
        }
    }
    block_set.clear();
    for(auto item : cnt_map){
        pq.push(make_pair(item.second,item.first));
    }
    for(int i = 0; i < k; i++){
        block_set.emplace_back(pq.top().second);
        pq.pop();
    }

//    for(unsigned i = 1; i < 20; i++){
//        if(neighbor[i].size()){
//            cout << "This is depth " << i << ", the size is " << neighbor[i].size() << endl;
//        }
//    }
    double time_total = omp_get_wtime()-start_total;
    ofstream out(outFile);
    out << "The rr set number = " << rr * n << "\n";
    out << "The final inf  = " << final_inf / try_times << "\n";
    out << "Finally saved " << final_save / try_times << " nodes\n";
    out << "Final inf(after block) = " << (final_inf / try_times) - (final_save / try_times) << "\n";
    out << "The reduction ratio = " << final_reduction / try_times << "\n";
    out << "The block budget is " << k << '\t';
    for(int write_node = 0; write_node < k; write_node++){
        out << block_set[write_node] << ' ';
    }
    out << '\n';
    out << "The average cost time is " << time_total / try_times << " s\n";
    out << "The change vector is " ;
    for(int i = 1; i <= max_height; i++)
        out << is_larger[i]  << ' ';
    out << '\n';
    out << "The replace number vector is " ;
    for(int i = 1; i <= max_height; i++)
        out << replace_num[i]  << ' ';
    out << '\n';
    out << "The save vector is ";
    for(int i = 1; i <= max_height; i++)
        out << saved_vector[i]  << ' ';
    out << '\n';
    out.close();
    cout << "\nTotal Time: " << time_total << "s" << endl;
    cout << "\n*******************" << endl;
    cout << "\tALL DONE" << endl;
    cout << "*******************" << endl;

    //
    // ALL DONE
    //
}