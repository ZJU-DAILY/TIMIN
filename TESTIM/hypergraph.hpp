#include "graph.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <functional>
#include <fstream>
#include <cstring>
#include <random>
#include <climits>
#include <omp.h>

using namespace std;

typedef pair<int,int> ii;
typedef pair<float,float> ff;
typedef pair<float,int> fi;
typedef vector<int> vi;
typedef vector<bool> vb;
typedef vector<float> vf;
typedef vector<ii> vii;
typedef vector<vi> vvi;
typedef vector<unordered_set<int>> vus;
typedef vector<vii> vvii;
typedef pair<bool,bool> bb;
typedef pair<bool,double>bd;

struct CompareBySecond {
    bool operator()(pair<int, int> a, pair<int, int> b)
    {
        return a.second < b.second;
    }
};
/*
 *
 */
Time_list delay_get(int len, bool flag, int lamda, float p){
    static std::default_random_engine e(time(0));
    static std::poisson_distribution<UI> distribution(lamda);
    static std::geometric_distribution<uint32_t> distribution_geo(p);
    Time_list a;
    for(size_t i=0;i<len;i++){
        UI number;
        if(!flag)
            number = distribution(e)+1;
        else
            number = distribution_geo(e)+1;
//        cout << number << " ";
        a.emplace_back(number);
    }
//    cout << endl;
    return a;
}

Time_list geometric_delay(size_t len, double p){
    static std::default_random_engine e(time(0));
    std::geometric_distribution<uint32_t> distribution(p);
    Time_list a;
    for(size_t i=0;i<len;i++){
        uint32_t number = distribution(e)+1;
        a.emplace_back(number);
    }
    return a;
}

double CalSeedSetInf(Graph &g, HyperGraph &hg, int usefulTR);

bool BuildBlockSet(Graph &g, HyperGraph &hg, int usefulBR, vi &RR_block_set, unsigned int count_cnt, double delta, double epsilon, int budget, int totalINF, bool fflag, int lamda, float p);

//// generating R1 set function
ii BuildHypergraphNode(Graph &g, HyperGraph &hg, int uStart, int hyperiiid, bool addHyperEdge ,vb &visit, vi &visit_mark, bool fflag, int lamda, float p){
    UI UI_MAX = 4294967295U;
    int cnt = 0;
    if(addHyperEdge) ////always true
    {
        if(!g.isFakeSeed(uStart))
            hg.hyperGT[hyperiiid].push_back(uStart);
        else
            return make_pair(1,0);
    }

    int n_visit_mark =0;
    int ddl = g.getDeadline();
    int max_edge_delay = g.getMaxEdgeDelay();
    priority_queue<ii, vii, greater<ii> > pq;
    pq.push(make_pair(0,uStart));
//    ASSERT(n_visit_mark < g.getNumNodes());
    visit_mark[n_visit_mark++]=uStart;
    visit[uStart]=true;
    int last_dist = 0;
    while(!pq.empty()) {
        int expand = pq.top().second;
        // total time delay
        int dist = pq.top().first;
        //update last dist
        last_dist = dist;
        // check if expand was fake seed
        if(g.isFakeSeed(expand))
            break;
        pq.pop();
        if(dist > ddl)
            continue;
        const vi &neigh = g.getInNeighbours(expand);
        const vf &w = g.getInWeights(expand);
        Time_list expand_time = delay_get(g.getInDegree(expand), fflag, lamda, p);
        if(neigh.size()==0)
            continue;
        double randDouble = sfmt_genrand_real1(&g.sfmt_seed) ;
//        double randDouble = sfmt_genrand_uint32(&g.sfmt_seed) / (double)UI_MAX ;
        for(int i=0; i<(int)neigh.size(); i++){
            if(expand_time[i] <= max_edge_delay)
                randDouble -= w[i] * hg.f_t[expand_time[i]];
            if(randDouble>0)
                continue;
            int v=neigh[i];
            if(visit[v]){
                continue;
//                break;
            }
            if(!visit[v])
            {
                visit_mark[n_visit_mark++]=v;
                visit[v]=true;
            }
            pq.push(make_pair(expand_time[i]+dist,v));
            if(addHyperEdge)
            {
                if(g.isFakeSeed(v))
                    cnt = 1;
                hg.hyperGT[hyperiiid].push_back(v);
            }
            break;
        }
    }
    for(int i=0; i<n_visit_mark; i++)
        visit[visit_mark[i]]=false;
//    cout << "n_visit_mark=" << n_visit_mark << endl;
    return make_pair(cnt,last_dist);
}


/*
* generate hyperedges in parallel following TLT model
*/
bd addHyperedgeParallel(Graph &g, HyperGraph &hg, long long num, long long already_done, vi &RR_block_set, unsigned int count_cnt, double delta, double epsilon, vi &is_RB_set, vi &set_dist, int &usefulTR,int &usefulBR, int budget, bool time_generate_flag, int lamda, float p)
{
	int numNodes = g.getNumNodes();
    vb visit(numNodes,false);
    vi visit_mark(numNodes,0);
    is_RB_set.resize(num);
    set_dist.resize(num);
    int ddl = g.getDeadline();
    while((long long)hg.hyperGT.size() <= num)
        hg.hyperGT.push_back(vi());

    //// generating  TR set  (R1) and picking RB set
    for(int i = already_done; i< num; i++){
        auto ress = BuildHypergraphNode(g, hg, sfmt_genrand_uint32(&g.sfmt_seed) % numNodes, i, true, visit, visit_mark, time_generate_flag, lamda, p);
        is_RB_set[i] = ress.first;
        set_dist[i] = ress.second;
    }
    usefulTR += (num - already_done);

    for(int i = (num==already_done?0:already_done); i < num; i++){
        if (set_dist[i] <= ddl){
            for(int t:hg.hyperGT[i])
            {
                hg.hyperG[t].push_back(i);
            }
            if (is_RB_set[i]){
                for(int t:hg.hyperGT[i])
                    hg.hyperGB[t].push_back(i);
                usefulBR++;
            }
        }
//        else
//            usefulTR--;
    }

    cout<<"total sampled T-RR set number: " << usefulTR <<endl;
    cout<<"total covered T-RR set (RB set) number: " << usefulBR <<endl; //// the number of RB-set
    // CalSeedSetInf(g, hg, usefulTR);
    double totalINF =  usefulBR * 1.0 / usefulTR * numNodes ; //// mis_inf
    cout<< "total misinformation INF: " <<totalINF<<endl;

//// -------------------------- select blocking set on R1 and evaluating effect on R2 -------------------------------------------

    if(BuildBlockSet(g, hg, usefulBR, RR_block_set, count_cnt, delta, epsilon, budget, totalINF, time_generate_flag, lamda, p)){
        cout << "It is because the constrain is satisfied\n";
        return bd(true,totalINF);
    }
    else
        return bd(false,totalINF);
//        cout << "total_num="<<num<<endl<<"testcnt="<<testcnt<<endl;
//    }

}




//// calculate the influence of the seed set use T-RR set (no use)
double CalSeedSetInf(Graph &g, HyperGraph &hg, int usefulTR)
{
    int n = g.getNumNodes();
    vector<int>coverage(n, 0);

    for (int i = 0; i < n; i++)
    {
        coverage[i] = (int)hg.hyperG[i].size();
    }
    const vi &fs = g.getFakeSeeds();

    long long influence = 0;
    long long numEdge = hg.hyperGT.size();

    // check if an edge is removed
    vector<bool> edgeMark(numEdge, false);
    // check if an node is remained in the heap
    vector<bool> nodeMark(n + 1, true);

    for (unsigned int i = 0; i < fs.size(); i++)
    {
        influence += coverage[fs[i]];
        nodeMark[fs[i]] = false;
        vector<int>e = hg.hyperG[fs[i]];
        for (unsigned int j = 0; j < e.size(); ++j){
            if (edgeMark[e[j]])continue;

            vector<int>nList = hg.hyperGT[e[j]];
            for (unsigned int l = 0; l < nList.size(); ++l){
                if (nodeMark[nList[l]])
                    coverage[nList[l]]--;
            }
            edgeMark[e[j]] = true;
        }
    }
//    cout << "Influence = " << influence << endl;
    return 1.0*influence / usefulTR * n;
}


//// generating R2 set function
bb BuildHypergraphNodeSingle(vector<bool>&tag, Graph &g, HyperGraph &hg, vb &visit, vi &visit_mark, bool fflag, int lamda, float p)
{
    UI UI_MAX = 4294967295U;
    int deadline = g.getDeadline();
    int max_edge_delay = g.getMaxEdgeDelay();
    int n = g.getNumNodes();
    auto uStart = sfmt_genrand_uint32(&g.sfmt_seed) % n;
//    if(tag[uStart] == false)
//        return true;
    if(g.isFakeSeed(uStart))
        return bb(true,true); //// bb.first (flag && flag2): covered RB set  bb.second (isRBset): RB set
    unsigned int n_visit_mark = 0;
    visit_mark[n_visit_mark++] = uStart;
    visit[uStart] = true;
    priority_queue<ii, vii, greater<ii> > pq;
    pq.push(make_pair(0,uStart));

    bool isBlock = false;
    bool iscoveredRB = false;
    bool isRBset = false;
    if(tag[uStart] == false)
        isBlock = true;

    while(!pq.empty()) {
        if(iscoveredRB || (iscoveredRB == false && isRBset == true))
            break;
        int expand = pq.top().second;
        // total time delay
        int dist = pq.top().first;
        if(dist > deadline){
            break;
        }
        pq.pop();
        const vi &neigh = g.getInNeighbours(expand);
        const vf &w = g.getInWeights(expand);
        Time_list expand_time = delay_get(g.getInDegree(expand), fflag, lamda, p);
        if(neigh.size()==0)
            continue;
        double randDouble = sfmt_genrand_real1(&g.sfmt_seed);
//        double randDouble= sfmt_genrand_uint32(&g.sfmt_seed) / (double)UI_MAX ;
        for(int i=0; i<(int)neigh.size(); i++){
            if(expand_time[i] <= max_edge_delay)
                randDouble -= w[i] * decayFunc(expand_time[i]);
            if(randDouble>0)
                continue;
            int v=neigh[i];
            if(visit[v]){
                continue;
//                break;
            }
            if(!visit[v])
            {
                visit_mark[n_visit_mark++]=v;
                visit[v]=true;

                if(tag[v]== false && expand_time[i] + dist <= deadline)//// finding blocking nodes
                    isBlock = true;

                if(isBlock == true)
                    if(g.isFakeSeed(v) && expand_time[i] + dist <= deadline) //// meet fake seeds after finding blocking nodes
                        iscoveredRB=true;

                if(g.isFakeSeed(v) && expand_time[i] + dist <= deadline)
                    isRBset = true;


            }
            pq.push(make_pair(expand_time[i]+dist,v));
            break;
        }
    }
    for (unsigned int i = 0; i < n_visit_mark; i++)visit[visit_mark[i]] = false;

    //hyperGT.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark));

    return bb(iscoveredRB, isRBset);
}
//build block set
bool BuildBlockSet(Graph &g, HyperGraph &hg, int usefulBR, vi &RR_block_set, unsigned int count_cnt, double delta, double epsilon, int budget, int totalINF, bool fflag, int lamda, float p)
{

    int n = g.getNumNodes();
    priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
    vector<int>coverage(n, 0);
    for (int i = 0; i < n; i++)
    {
        if(g.isFakeSeed(i))
            continue;
        pair<int, int>tep(make_pair(i, (int)hg.hyperGB[i].size()));
        heap.push(tep);
        coverage[i] = (int)hg.hyperGB[i].size();
    }
    int maxInd;

    long long coveredRB = 0;
    long long numEdge = hg.hyperGT.size();

    // check if an edge is removed
    vector<bool> edgeMark(numEdge, false);
    // check if an node is remained in the heap
    vector<bool> nodeMark(n + 1, true);

    RR_block_set.clear();
    while ((int)RR_block_set.size() < budget)
    {
        pair<int, int>ele = heap.top();
        heap.pop();
        if (ele.second > coverage[ele.first])
        {
            ele.second = coverage[ele.first];
            heap.push(ele);
            continue;
        }

        maxInd = ele.first;
        vector<int>e = hg.hyperGB[maxInd];  //the edge influence
        coveredRB += coverage[maxInd]; //// total number of covered RB set
//        cout << "coverage i = " << coverage[maxInd] <<endl;
        RR_block_set.push_back(maxInd);
        nodeMark[maxInd] = false;  //// nodeMark=false: maxInd is a blocking node
//        cout << 1.0 * coverage[maxInd] / usefulBR * totalINF << endl;     //saved nodes by each block node
        for (unsigned int j = 0; j < e.size(); ++j){
            if (edgeMark[e[j]])continue;

            vector<int>nList = hg.hyperGT[e[j]];
            for (unsigned int l = 0; l < nList.size(); ++l){
                if (nodeMark[nList[l]])
                    coverage[nList[l]]--;
            }
            edgeMark[e[j]] = true;
        }
    }
    vb visit(n,false);
    vi visit_mark(n,0);

    cout << "[R1] num of covered RB: " << coveredRB << "\t num of total RB: " << usefulBR << endl;

    hg.EstPro1 = 1.0 * coveredRB / usefulBR * totalINF;  //// the expected number of saved nodes on R1


////-------------------------------------- estimate on R2 ----------------------------------------------
    unsigned int i = 0;
    hg.EstPro2 = 0.0;
    int num_rbset = 0;
    double valid_sum = 0.0;
    // pow(2,count_cnt-1)
    while(i++ < numEdge) {
        auto rrr = BuildHypergraphNodeSingle(nodeMark, g, hg, visit, visit_mark, fflag, lamda, p);
        valid_sum += rrr.first;
        num_rbset += rrr.second;
    }
    cout << "[R2] num of covered RB: " << valid_sum << "\t num of total RB: " << num_rbset << endl;

    hg.EstPro2 = valid_sum / num_rbset * totalINF;

    cout << "R1 saved " << hg.EstPro1 << " nodes\n";
    cout << "R2 saved " << hg.EstPro2  << " nodes\n";

////-------------------------------------- evaluate current result ----------------------------------------------

//    double x = hg.EstPro2 * numEdge / log(5. * count_cnt * count_cnt / delta) / n;
    double x1 = hg.EstPro2 * coveredRB / log(5. * count_cnt * count_cnt / delta) / n;
    double x2 = hg.EstPro2 * valid_sum / log(5. * count_cnt * count_cnt / delta) / n;
    hg.eps1=4./(sqrt(1.+8.*x1)-3);
    //ASSERT(eps1>0);
    hg.eps2=sqrt((2.*hg.eps1+2.)/x2);
    double t=hg.EstPro1 / hg.EstPro2;


    double cons=(1-1/exp(1))/(1-1/exp(1)-epsilon) * (1-hg.eps2) / (1+hg.eps1);
    cout << "eps1: " << hg.eps1<< "\t eps2: " << hg.eps2<< "\t t: " << t<< "\t cons: " << cons<<endl;

//&& ((hg.eps1 + hg.eps2) <= epsilon)
    if(t <= (1-1/exp(1))/(1-1/exp(1)-epsilon) * (1-hg.eps2) / (1+hg.eps1)){
        cout << "[satisfy constrains] current misinformation INF on R1: " << totalINF - hg.EstPro1 << endl;
        return true;
    }
    return false;
//    return 1.0 * influence / usefulBR * n;    //the number of nodes saved
}




double Simple_BuildBlockSet(Graph &g, HyperGraph &hg, int usefulBR, vi &RR_block_set, int budget, int totalINF)
{

    int n = g.getNumNodes();
    priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
    vector<int>coverage(n, 0);
    for (int i = 0; i < n; i++)
    {
        if(g.isFakeSeed(i))
            continue;
        pair<int, int>tep(make_pair(i, (int)hg.hyperGB[i].size()));
        heap.push(tep);
        coverage[i] = (int)hg.hyperGB[i].size();
    }

    int maxInd;

    long long influence = 0;
    long long numEdge = hg.hyperGT.size();

    // check if an edge is removed
    vector<bool> edgeMark(numEdge, false);
    // check if an node is remained in the heap
    vector<bool> nodeMark(n + 1, true);

    RR_block_set.clear();
    while ((int)RR_block_set.size() < budget)
    {
        pair<int, int>ele = heap.top();
        heap.pop();
        if (ele.second > coverage[ele.first])
        {
            ele.second = coverage[ele.first];
            heap.push(ele);
            continue;
        }

        maxInd = ele.first;
        vector<int>e = hg.hyperGB[maxInd];  //the edge influence
        influence += coverage[maxInd];
        RR_block_set.push_back(maxInd);
        nodeMark[maxInd] = false;
//        cout << 1.0 * coverage[maxInd] / usefulBR * totalINF << endl;     //saved nodes by each block node
        for (unsigned int j = 0; j < e.size(); ++j){
            if (edgeMark[e[j]])continue;

            vector<int>nList = hg.hyperGT[e[j]];
            for (unsigned int l = 0; l < nList.size(); ++l){
                if (nodeMark[nList[l]])coverage[nList[l]]--;
            }
            edgeMark[e[j]] = true;
        }
    }
    return 1.0 * influence / usefulBR * totalINF;

}

ii simpleBuild(Graph &g, HyperGraph &hg, long long num, double &inf, bool fflag, int lamda, float p)
{
    int numNodes = g.getNumNodes();
//	vvii hyperedges;
    vb visit(numNodes,false);
    vi visit_mark(numNodes,0);
    vi is_RB_set(num,0);
    vi set_dist(num,0);
    int ddl = g.getDeadline();
    int usefulTR = num;
    int usefulBR = 0;
    while((long long)hg.hyperGT.size() <= num)
        hg.hyperGT.push_back(vi());
    for(int i = 0; i< num; i++){
        auto ress = BuildHypergraphNode(g, hg, sfmt_genrand_uint32(&g.sfmt_seed) % numNodes, i, true, visit, visit_mark, fflag, lamda, p);
        is_RB_set[i] = ress.first;
        set_dist[i] = ress.second;
//        if(i % 5000 == 0)
//            cout << i << endl;
    }
    for(int i=0; i < num; i++){
        if (set_dist[i] <= ddl){
//            for(int t:hg.hyperGT[i])
//            {
//                hg.hyperG[t].push_back(i);
//            }
            if (is_RB_set[i]){
                for(int t:hg.hyperGT[i])
                    hg.hyperGB[t].push_back(i);
                usefulBR++;
            }
        }
//        else
//            usefulTR--;
    }
    inf = usefulBR * 1.0 / usefulTR * numNodes;
//    inf = CalSeedSetInf(g, hg, usefulTR);
//    cout << "The final inf = " << inf << endl;
    return make_pair(usefulTR,usefulBR);
}
bd max_cover(Graph &g, HyperGraph &hg,vi &RR_block_set, int budget, vi choose_neighbor, vi &final_coverage, double INF, int usefulBR, int run_time, vi &replace)
{
    int numNodes = g.getNumNodes();
    priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
    vector<int>coverage(numNodes, 0);
    for (int i = 0; i < choose_neighbor.size(); i++)
    {
        if(g.isFakeSeed(choose_neighbor[i]))
            continue;
        coverage[choose_neighbor[i]] = (int)hg.hyperGB[choose_neighbor[i]].size();
        if(run_time == 1){
            pair<int, int>tep(make_pair(choose_neighbor[i], (int)hg.hyperGB[choose_neighbor[i]].size()));
            heap.push(tep);
        }
    }
    int maxInd;
    int influence = 0;
    bool flag = false;
    long long numEdge = hg.hyperGT.size();

    // check if an edge is removed
    vector<bool> edgeMark(numEdge, false);
    // check if an node is remained in the heap
    vector<bool> nodeMark(numNodes + 1, true);
    if(run_time == 1){
        while ((int)RR_block_set.size() < budget)
        {
            pair<int, int>ele = heap.top();
            heap.pop();
            if (ele.second > coverage[ele.first])
            {
                ele.second = coverage[ele.first];
                heap.push(ele);
                continue;
            }

            maxInd = ele.first;
            vector<int>e = hg.hyperGB[maxInd];  //the edge influence
            influence += coverage[maxInd];
            RR_block_set.push_back(maxInd);
            final_coverage.emplace_back(coverage[maxInd]);
            nodeMark[maxInd] = false;
//        cout << 1.0 * coverage[maxInd] / usefulBR * totalINF << endl;     //saved nodes by each block node
            for (unsigned int j = 0; j < e.size(); ++j){
                if (edgeMark[e[j]])continue;

                vector<int>nList = hg.hyperGT[e[j]];
                for (unsigned int l = 0; l < nList.size(); ++l){
                    if (nodeMark[nList[l]])coverage[nList[l]]--;
                }
                edgeMark[e[j]] = true;
            }
        }
    }
    else{
        int cnt = 0;
        while(1){
            if(cnt == 0){
                for(int i = 0; i < budget -1; i++){
                    vector<int>ee = hg.hyperGB[RR_block_set[i]];  //the edge influence
                    nodeMark[RR_block_set[i]] = false;
                    for (unsigned int j = 0; j < ee.size(); ++j){
                        if (edgeMark[ee[j]])continue;
                        vector<int>nList = hg.hyperGT[ee[j]];
                        for (unsigned int l = 0; l < nList.size(); ++l){
                            if (nodeMark[nList[l]])coverage[nList[l]]--;
                        }
                        edgeMark[ee[j]] = true;
                    }
                }
                for (int i = 0; i < choose_neighbor.size(); i++) {
                    pair<int, int> tep(make_pair(choose_neighbor[i], coverage[choose_neighbor[i]]));
                    heap.push(tep);
                }
                cnt++;
            }
            pair<int, int>ele = heap.top();
            heap.pop();
            maxInd = ele.first;
            vector<int>e = hg.hyperGB[maxInd];  //the edge influence
            if(coverage[maxInd] > final_coverage[budget-1]){
                replace[run_time]++;
                flag = true;
                RR_block_set.pop_back();
                final_coverage.pop_back();
                RR_block_set.push_back(maxInd);
                final_coverage.emplace_back(coverage[maxInd]);
                nodeMark[maxInd] = false;
                for (unsigned int j = 0; j < e.size(); ++j){
                    if (edgeMark[e[j]])continue;
                    vector<int>nList = hg.hyperGT[e[j]];
                    for (unsigned int l = 0; l < nList.size(); ++l){
                        if (nodeMark[nList[l]])coverage[nList[l]]--;
                    }
                    edgeMark[e[j]] = true;
                }
                int nowID = budget - 1;
                while (final_coverage[nowID] > final_coverage[nowID-1]){
                    swap(final_coverage[nowID],final_coverage[nowID-1]);
                    swap(RR_block_set[nowID],RR_block_set[nowID-1]);
                    nowID--;
                }
            }
            else{
                for(auto cover : final_coverage)
                    influence += cover;
                break;
            }
        }

    }
    return bd(flag,1.0 * influence / usefulBR * INF);
}

double CalSaveNumberToBaseline(Graph &g, HyperGraph &hg,vi &RR_block_set, double INF, int usefulBR)
{
    int numNodes = g.getNumNodes();
    long long numEdge = hg.hyperGT.size();
    int k = RR_block_set.size();
    int influence = 0;
    vector<int>coverage(numNodes, 0);
    for (int i = 0; i < k; i++)
    {
        coverage[RR_block_set[i]] = (int)hg.hyperGB[RR_block_set[i]].size();
//        cout << coverage[RR_block_set[i]] << " ";
    }
//    cout << endl;
    // check if an edge is removed
    vector<bool> edgeMark(numEdge, false);
    // check if an node is remained in the heap
    vector<bool> nodeMark(numNodes + 1, true);
    for(int i = 0; i < k; i++){
        vector<int>e = hg.hyperGB[RR_block_set[i]];  //the edge influence
        influence += coverage[RR_block_set[i]];
//        cout << coverage[RR_block_set[i]] << " ";
        nodeMark[RR_block_set[i]] = false;
        for (unsigned int j = 0; j < e.size(); ++j){
            if (edgeMark[e[j]])continue;
            vector<int>nList = hg.hyperGT[e[j]];
            for (unsigned int l = 0; l < nList.size(); ++l){
                if (nodeMark[nList[l]])coverage[nList[l]]--;
            }
            edgeMark[e[j]] = true;
        }
    }
    cout << influence << "  " << usefulBR << "  " << INF << endl;
    return 1.0 * influence / usefulBR * INF;
}

// for baseline IMM
double getEPT(Graph &g, HyperGraph &coll_one,vi &block_set, int usefulBR, int k)
{
    int n = g.getNumNodes();
    priority_queue<pair<int, int>, vector<pair<int, int>>, CompareBySecond>heap;
    vector<int>coverage(n, 0);
    for (int i = 0; i < n; i++)
    {
        if(g.isFakeSeed(i))
            continue;
        pair<int, int>tep(make_pair(i, (int)coll_one.hyperGB[i].size()));
        heap.push(tep);
        coverage[i] = (int)coll_one.hyperGB[i].size();

    }
    int maxInd;
    long long influence = 0;
    long long numEdge = coll_one.hyperGT.size();
    vector<bool> edgeMark(numEdge, false);
    vector<bool> nodeMark(n + 1, true);
    while ((int)block_set.size() < k)
    {
        pair<int, int>ele = heap.top();
        heap.pop();
        if (ele.second > coverage[ele.first])
        {
            ele.second = coverage[ele.first];
            heap.push(ele);
            continue;
        }
        maxInd = ele.first;
        vector<int>e = coll_one.hyperGB[maxInd];  //the edge influence
        influence += coverage[maxInd];
        block_set.push_back(maxInd);
        nodeMark[maxInd] = false;
        for (unsigned int j = 0; j < e.size(); ++j){
            if (edgeMark[e[j]])continue;
            vector<int>nList = coll_one.hyperGT[e[j]];
            for (unsigned int l = 0; l < nList.size(); ++l){
                if (nodeMark[nList[l]])coverage[nList[l]]--;
            }
            edgeMark[e[j]] = true;
        }
    }
    return influence * 1.0 / usefulBR;
}