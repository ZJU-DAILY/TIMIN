#include "graph.h"
#include <algorithm>
#include <utility>
#include <vector>
#include <iostream>
#include <fstream>
#include <queue>
#include <cstdlib>
#include <sstream>
#include <queue>
#include <climits>
#include <cmath>
#include <unistd.h>
#include <ctime>
#include <unordered_map>
#include <map>

using namespace std;

typedef pair<int,int> ii;
typedef vector<int> vi;
typedef vector<bool> vb;
typedef vector<float> vf;
typedef vector<double> vd;
typedef vector<ii> vii;
typedef vector<vi> vvi;
typedef vector<vii> vvii;

Graph::Graph(int edge_time_max, int deadline)
{
    sfmt_init_gen_rand(&sfmt_seed, rand());
    _t=edge_time_max;
    _T=deadline;
}

void Graph::init_visit_thresh_hold()
{
    // assign all thresh hold to > 1 denote unvisited
    visit_thresh_hold = vector<double> (num_nodes, 2.0);
}

double decayFunc(int len, float alpha){
    if(len == 0)
        return 0;
    return  pow(len,alpha);
}

const vector<int> & Graph::getOutNeighbours (int u) const
{
    return adj_list[u];
}

const vector<int> & Graph::getInNeighbours (int u) const
{
    return rev_adj_list[u];
}

const vector<float> & Graph::getOutWeights (int u) const
{
    return weights[u];
}

const vector<float> & Graph::getInWeights (int u) const
{
    return rev_weights[u];
}

const vector<int> & Graph::getFakeSeeds () const
{
    return fake_seeds_index;
}

const vector<double> & Graph::getThreshold() const
{
    return visit_thresh_hold;
}

unsigned int Graph::getOutDegree(int u) const
{
    return adj_list[u].size();
}

unsigned int Graph::getInDegree(int u) const
{
    return rev_adj_list[u].size();
}

/*
* get the number of nodes
*/
unsigned int Graph::getNumNodes() const
{
    return num_nodes;
}

/*
* get the number of edges
*/
unsigned int Graph::getNumEdges() const
{
    return num_edges;
}

/*
 *get t
 */
int Graph::getMaxEdgeDelay()
{
    return _t;
}

/*
 *get T
 */
int Graph::getDeadline()
{
    return _T;
}
/*
 * set t
 */
void Graph::setMaxEdgeDelay(int t)
{
    _t = t;
}
/*
 * set T
 */
void Graph::setDeadline(int T)
{
    _T = T;
}
/*
* get the number of fake seeds
*/
unsigned int Graph::getNumFakeSeeds() const
{
    return num_fake_seeds;
}

/*
* determine if node v is a fake seed
*/
bool Graph::isFakeSeed(int v) const
{
    return fake_seed[v];
}

/*
* read input graph
*/
void Graph::readGraph(const char* filename, bool fixed, float w)
{
    FILE * pFile1;
    FILE * pFile2;
    string filename1 = filename;
    filename1.append(".bin");
    string filename2 = filename;
    filename2.append("_rev.bin");
    pFile1 = fopen(filename1.c_str(), "rb");
    pFile2 = fopen(filename2.c_str(), "rb");
    fread(&num_nodes, sizeof(int), 1, pFile1);
    fread(&num_edges, sizeof(long long), 1, pFile1);
    node_deg = vi(num_nodes);
    fread(&node_deg[0], sizeof(int), num_nodes, pFile1);
    rev_node_deg = vi(num_nodes);
    fread(&rev_node_deg[0], sizeof(int), num_nodes, pFile2);

    for (unsigned int i = 0; i < num_nodes; i++){
        vi tmp1(node_deg[i]);
        fread(&tmp1[0], sizeof(int), node_deg[i], pFile1);
        adj_list.push_back(tmp1);

        vi tmp2(rev_node_deg[i]);
        fread(&tmp2[0], sizeof(int), rev_node_deg[i], pFile2);
        rev_adj_list.push_back(tmp2);
    }

    for (unsigned int i = 0; i < num_nodes; i++){
        vf tmp1(node_deg[i], w);
        if (!fixed) fread(&tmp1[0], sizeof(float), node_deg[i], pFile1);
        weights.push_back(tmp1);

        vf tmp2(rev_node_deg[i], w);
        if (!fixed) fread(&tmp2[0], sizeof(float), rev_node_deg[i], pFile2);
        rev_weights.push_back(tmp2);
    }
}

/*
* read fake seeds
*/
void Graph::readFakeSeeds(const char* filename)
{
    int fs;
    ifstream in(filename);
    in >> num_fake_seeds;
    fake_seeds_index = vi(num_fake_seeds);
    fake_seed = vb(getNumNodes(), false);
    for (unsigned int i = 0; i < num_fake_seeds; i++){
        in >> fs;
        fake_seeds_index[i] = fs;
        fake_seed[fs] = true;
    }
    in.close();
}

// compute mitigation lower bound based on top k nodes from a depth 1 MIA
double Graph::computeMitigationLowerBound(unsigned int n, unsigned int k)
{
    int i, j;
    int num_seen = 0;
    vi seen_index(n,0);
    vf ap(n,0);

    int num_fs = getNumFakeSeeds();
    const vi &fs = getFakeSeeds();
    for (i = 0; i < num_fs; i++) {
        const vf &w = getOutWeights(fs[i]);
        const vi &neigh = getOutNeighbours(fs[i]);
        for (j = 0; j < node_deg[fs[i]]; j++) {
            if (ap[neigh[j]] < w[j]) {
                if (ap[neigh[j]] == 0) seen_index[num_seen++] = neigh[j];
                ap[neigh[j]] = w[j];
            }
        }
    }

    vf sorted_ap;
    for (i = 0; i < num_seen; i++) {
        sorted_ap.push_back(ap[seen_index[i]]);
    }
    sort(sorted_ap.begin(), sorted_ap.end(), greater<int>());

    double sum = 0.0;
    int len = (sorted_ap.size() < k) ? sorted_ap.size() : k;
    for (i = 0; i < len; i++) {
        sum += sorted_ap[i];
    }

    return sum * 2.0;
}

// generate a single forward monte carlo estimate of the influence of input seed
int Graph::generateInfluenceSample(vb &visit, vi &visit_index, vd &threshold, vi &dist, int root)
{
    int i, cur;
    float flip;

    int curPos = 0;
    int num_marked = 1;
    int edge_max_delay = getMaxEdgeDelay();
    int deadLine = getDeadline();
    visit[root] = true;
    visit_index[0] = root;
    dist[root] = 0;
    threshold[root] = -1.0;

    while(curPos < num_marked) {
        cur = visit_index[curPos];
        const vf &w = getOutWeights(cur);
        const vi &neigh = getOutNeighbours(cur);
        Time_list expand_time = poisson_time(node_deg[cur]);
        for (i = 0; i < node_deg[cur]; i++) {
            if(threshold[neigh[i]] < 0)
                continue;
            if(threshold[neigh[i]] > 1.0){
                flip = sfmt_genrand_uint32(&sfmt_seed) / (float)UI_MAX;
                threshold[neigh[i]] = flip;
            }
            else
                flip = threshold[neigh[i]];
            if(expand_time[i] <= edge_max_delay){
                const vi &inNeigh = getInNeighbours(neigh[i]);
                unsigned int inD = getInDegree(neigh[i]);
                const vf &inw = getInWeights(neigh[i]);
                for(int k = 0 ; k < inD ; k++){
                    int o = inNeigh[k];
                    // for activated node
                    if (threshold[o] < 0.0){
//                    cout<<"now visit "<<neigh[i]<<" his active neighbor is "<<o<<endl;
                        flip -= inw[k];
                    }
                }
                if (flip <= 0.0 && dist[cur] + expand_time[i] < dist[neigh[i]] && dist[cur] + expand_time[i] <= deadLine) {
                    dist[neigh[i]] = dist[cur] + expand_time[i];
                    threshold[neigh[i]] = -1.0;
                    if (!visit[neigh[i]]) {
                        visit_index[num_marked] = neigh[i];
                        num_marked++;
                        visit[neigh[i]] = true;
                    }
                }
            }
        }
        curPos++;
    }

    for(i = 0; i < num_marked; i++) {
        visit[visit_index[i]] = false;
    }
    for(i = 0; i < (int)threshold.size(); i++) {
        threshold[i] = 2.0;
    }
    for(i = 0; i < (int)dist.size(); i++) {
        dist[i] = 100000;
    }
    return num_marked;
}

// generate a single forward monte carlo estimate of the outward influence of F
int Graph::generateFakeInfluenceSample(vb &visit, vi &visit_index, vi &dist, priority_queue<ii, vii, greater<ii> > &pq)
{
    int edge_max_delay = getMaxEdgeDelay();
    int deadLine = getDeadline();
    auto n_nodes = getNumNodes();
    vector<vector<int> >liveChild(n_nodes,vector<int>());
    vector<vector<int> >delayChild(n_nodes,vector<int>());
    for(int i=0; i < n_nodes; i++){
        auto flip = sfmt_genrand_real1(&sfmt_seed);
        const vi &inNeigh = getInNeighbours(i);
        unsigned int inD = getInDegree(i);
        const vf &inw = getInWeights(i);
        Time_list delay_time = poisson_time(inD);
        for(int k = 0 ; k < inD ; k++){
            flip -= inw[k];
            if(flip > 0)
                continue;
            liveChild[inNeigh[k]].emplace_back(i);
            delayChild[inNeigh[k]].emplace_back(delay_time[k]);
            break;
        }
    }
    int i, cur;
    int num_marked = getNumFakeSeeds();
    const vi &fs = getFakeSeeds();
    auto fake_num = getNumFakeSeeds();
    for (i = 0; i < fake_num; i++) {
        pq.push(make_pair(0,fs[i]));
        visit_index[i]=fs[i];
        dist[fs[i]]=0;
    }

    while(!pq.empty()) {
        cur = pq.top().second; pq.pop();
        if (visit[cur]) continue; // duplicate entry in pq
        visit[cur] = true;
        if(!isFakeSeed(cur)){
            visit_index[num_marked] = cur;
            num_marked++;
        }
        auto degree = (int)liveChild[cur].size();
        for (i = 0; i < degree; i++) {
            int now_dist = dist[cur] + delayChild[cur][i];
            if (now_dist < dist[liveChild[cur][i]] && now_dist <= 32) {
                dist[liveChild[cur][i]] = now_dist;
                pq.push(make_pair(dist[liveChild[cur][i]],liveChild[cur][i]));
            }
        }
    }

    for(i = 0; i < num_marked; i++) {
        visit[visit_index[i]] = false;
    }
    for(i = 0; i < (int)dist.size(); i++) {
        dist[i] = 100000;
    }
//    cout<<"num_marked = " <<num_marked<<"\n";
    return num_marked - getNumFakeSeeds();
}
//one mc using live edge
int Graph::liveEdgeSampleGreedy(vb &visit, vi &visit_index, vi &dist, const vi &block, priority_queue<ii, vii, greater<ii> > &pq)
{
//    int edge_max_delay = getMaxEdgeDelay();
    int deadLine = getDeadline();
    auto n_nodes = getNumNodes();
    vector<vector<int> >liveChild(n_nodes,vector<int>());
    vector<vector<int> >delayChild(n_nodes,vector<int>());
    for(int i=0; i < n_nodes; i++){
        auto flip = sfmt_genrand_real1(&sfmt_seed);
        const vi &inNeigh = getInNeighbours(i);
        unsigned int inD = getInDegree(i);
        const vf &inw = getInWeights(i);
        Time_list delay_time = poisson_time(inD);
        for(int k = 0 ; k < inD ; k++){
            flip -= inw[k];
            if(flip > 0)
                continue;
            liveChild[inNeigh[k]].emplace_back(i);
            delayChild[inNeigh[k]].emplace_back(delay_time[k]);
            break;
        }
    }
    int i, cur;
    int num_marked = getNumFakeSeeds();
    unordered_map<int,bool> mm;
    for (i = 0; i< (int)block.size() ;i++)
        mm[block[i]] = true;
    const vi &fs = getFakeSeeds();
    auto fake_num = getNumFakeSeeds();
    for (i = 0; i < fake_num; i++) {
        pq.push(make_pair(0,fs[i]));
        visit_index[i]=fs[i];
        dist[fs[i]]=0;
    }

    while(!pq.empty()) {
        cur = pq.top().second; pq.pop();
        if (visit[cur]) continue; // duplicate entry in pq
        visit[cur] = true;
        if(!isFakeSeed(cur)){
            visit_index[num_marked] = cur;
            num_marked++;
        }
        auto degree = (int)liveChild[cur].size();
        for (i = 0; i < degree; i++) {
            if(mm.count(liveChild[cur][i]))
                continue;
            int now_dist = dist[cur] + delayChild[cur][i];
            if (now_dist < dist[liveChild[cur][i]] && now_dist <= deadLine) {
                dist[liveChild[cur][i]] = now_dist;
                pq.push(make_pair(dist[liveChild[cur][i]],liveChild[cur][i]));
            }
        }
    }

    for(i = 0; i < num_marked; i++) {
        visit[visit_index[i]] = false;
    }
    for(i = 0; i < (int)dist.size(); i++) {
        dist[i] = 100000;
    }
//    cout<<"num_marked = " <<num_marked<<"\n";
    return num_marked - getNumFakeSeeds();
}

//generate a forward of block graph
int Graph::generateBlockSampleGreedy(vb &visit, vi &visit_index, vd &threshold, vi &dist, const vi &block, priority_queue<ii, vii, greater<ii> > &pq)
{
    int edge_max_delay = getMaxEdgeDelay();
    int deadLine = getDeadline();
    int i, cur;
    float flip;
    int num_marked = getNumFakeSeeds();
    unordered_map<int,bool> mm;
    for (i = 0; i< (int)block.size() ;i++)
        mm[block[i]] = true;
    const vi &fs = getFakeSeeds();
    for (i = 0; i < num_marked; i++) {
        pq.push(make_pair(0,fs[i]));
        visit_index[i]=fs[i];
        threshold[fs[i]] = -1.0;
        dist[fs[i]]=0;
    }

    while(!pq.empty()) {
        cur = pq.top().second; pq.pop();
        if (visit[cur]) continue; // duplicate entry in pq
        visit[cur] = true;
        if(dist[cur]>0){
            visit_index[num_marked] = cur;
            num_marked++;
        }
        const vf &w = getOutWeights(cur);
        const vi &neigh = getOutNeighbours(cur);
        Time_list expand_time = poisson_time(node_deg[cur]);
        for (i = 0; i < node_deg[cur]; i++) {
//            if(threshold[neigh[i]] < 0 || mm[neigh[i]] == true)
            if(threshold[neigh[i]] < 0 || mm.count(neigh[i]))
                continue;
            if(threshold[neigh[i]] > 1.0){
                flip = sfmt_genrand_real1(&sfmt_seed);
//                flip = sfmt_genrand_uint32(&sfmt_seed) / (float)UI_MAX;
                threshold[neigh[i]] = flip;
            }
            else
                flip = threshold[neigh[i]];
            if(expand_time[i] < edge_max_delay){
                const vi &inNeigh = getInNeighbours(neigh[i]);
                unsigned int inD = getInDegree(neigh[i]);
                const vf &inw = getInWeights(neigh[i]);
                for(int k = 0 ; k < inD ; k++){
                    int o = inNeigh[k];
                    // for activated node
                    if (threshold[o] < 0.0){
//                    cout<<"now visit "<<neigh[i]<<" his active neighbor is "<<o<<endl;
                        flip -= inw[k];
                    }
                }
                if (flip <= 0.0 && dist[cur] + expand_time[i] < dist[neigh[i]] && dist[cur] + expand_time[i] < deadLine) {
                    dist[neigh[i]] = dist[cur] + expand_time[i];
                    threshold[neigh[i]] = -1.0;
                    pq.push(make_pair(dist[neigh[i]],neigh[i]));
                }
            }
        }
    }

    for(i = 0; i < num_marked; i++) {
        visit[visit_index[i]] = false;
    }
    for(i = 0; i < (int)threshold.size(); i++) {
        threshold[i] = 2.0;
    }
    for(i = 0; i < (int)dist.size(); i++) {
        dist[i] = 100000;
    }

    return num_marked - getNumFakeSeeds();
}

HyperGraph::HyperGraph(unsigned int n)
{
    sfmt_init_gen_rand(&sfmtSeed, rand());
    node_hyperedges = vvii(n);
    node_hyperedge_weight = vi(n);
    EstPro1=1.0, EstPro2=0.1, eps1=0., eps2=0.;
}

/*
* Add a hyperedge into the hypergraph
*/
void HyperGraph::addEdge(vii &edge)
{
    hyperedges.push_back(edge);
    unsigned int index = hyperedges.size() - 1;
    for (unsigned int i = 0; i < edge.size(); i++) {
        node_hyperedges[edge[i].first].push_back(make_pair(index, edge[i].second));
        node_hyperedge_weight[edge[i].first] += edge[i].second;
    }
}

/*
* get an edge from the hypergraph
*/
const vector<pair<int,int> > & HyperGraph::getEdge(int e) const
{
    return hyperedges[e];
}

/*
* get the list of hyperedges incident to node u
*/
const vector<pair<int,int> > & HyperGraph::getNode(int u) const
{
    return node_hyperedges[u];
}

/*
* get the list of hyperedges incident to node u
*/
int HyperGraph::getNodeWeight(int u) const
{
    return node_hyperedge_weight[u];
}

/*
* get the number of hyperedges
*/
int HyperGraph::getNumEdge() const
{
    return hyperedges.size();
}

/*
* remove all the hyperedges
*/
void HyperGraph::clearEdges()
{
    //hyperedges.clear();
    vvii().swap(hyperedges);
    node_hyperedges.clear();
    node_hyperedge_weight.clear();
    cout << "clear edges!" << endl;
}

// generating reachability set of fake campaign
//bool HyperGraph::phaseOne(Graph &g, ii &root_data, ii &traversal_data, priority_queue<ii, vii, greater<ii> > &pq, vb &visit, vi &visit_index, vb &dead_visit, vi &dead_visit_index, vi &dist,  vvi &dead_parents)
//{
//    int cur, root_index;
//    unsigned int i;
//    int num_marked = 0;
//    int dead_visit_num = 0;
//    int t=g.getMaxEdgeDelay();
//    int T=g.getDeadline();
//    unsigned int num_fs = g.getNumFakeSeeds();
//    const vi &fs = g.getFakeSeeds();
//    auto threshold = g.getThreshold();
//    for (i = 0; i < num_fs; i++) {
//        dist[fs[i]] = 0;
//        pq.push(make_pair(0,fs[i]));
//        threshold[fs[i]] = -1.0;
//    }
//
//    while ( !pq.empty() ) {
//        cur = pq.top().second; pq.pop();
//        if (visit[cur]) continue; // duplicate entry in pq
//        visit[cur] = true;
//        visit_index[num_marked] = cur;
//        num_marked++;
//        const vf &w = g.getOutWeights(cur);
//        const vi &neigh = g.getOutNeighbours(cur);
//        int len=g.node_deg[cur];
//        Time_list expand_time = poisson_time(len);
//        for (i = 0; i < len; i++) {
//            if (!g.fake_seed[neigh[i]]) { // do not consider in-neighbours of nodes in S_F
//                if(threshold[neigh[i]]<0)
//                    continue;
//                double randDouble;// generate threshold if not generated (node unvisited)
//                if(threshold[neigh[i]]>1.0){
//                    randDouble = sfmt_genrand_uint32(&sfmtSeed) / (double)(g.UI_MAX);
//                    cout<<"now visit "<<cur<<"'s neighbor "<<neigh[i]<<", his threshold is "<<randDouble<<endl;
//                    cout<<"now visit neighbor "<<neigh[i]<<", his delay is "<<expand_time[i]<<endl;
//                    threshold[neigh[i]] = randDouble;
//                } else
//                    randDouble = threshold[neigh[i]];
//                if (expand_time[i] < t) {
//                    const vi &inNeigh = g.getInNeighbours(neigh[i]);
//                    unsigned int inD = g.getInDegree(neigh[i]);
//                    const vf &inw = g.getInWeights(neigh[i]);
//                    for(int k = 0 ; k < inD ; k++){
//                        int o = inNeigh[k];
//                        // for activated node
//                        if (threshold[o] < 0.0){
//                            cout<<"now visit "<<neigh[i]<<" his active neighbor is "<<o<<endl;
//                            randDouble -= inw[k] * f_t[expand_time[i]];
//                        }
//
//                    }
//                    if (randDouble <= 0.0 && dist[cur] + expand_time[i] < dist[neigh[i]] && dist[cur] + expand_time[i] < T) {
//                        dist[neigh[i]] = dist[cur] + expand_time[i];
//                        threshold[neigh[i]] = -1.0;
//                        pq.push(make_pair(dist[neigh[i]],neigh[i]));
//                    }
//                }
//            }
//        }
//    }
//    cout << "\n*******************" << endl;
//    cout << "\tnum_marked="<<num_marked<< endl;
//    cout << "*******************" << endl;
//    // if nodes reachable from F
//    if (num_marked > num_fs) {
//        // select root uniformly at random from nodes reached by F
//        root_index = sfmt_genrand_uint32(&sfmtSeed)%(num_marked - num_fs);
//        root_data.first = visit_index[root_index + num_fs];
//        root_data.second = dist[root_data.first];
//    }
//
//    // reset local data structures
//    for(i = 0; i < num_marked; i++) {
//        dist[visit_index[i]] = INT_MAX;
//    }
//
//
//    traversal_data.first = num_marked;
//
//
//    if (num_marked <= num_fs) return true; // F failed to activate any nodes
//
//
//
//    return false;
//}

/*
* convert from an integer to a string
*/
string intToStr(int i) {
    stringstream ss;
    ss << i;
    return ss.str();
}

/*
* convert from a strong to an integer
*/
unsigned int strToInt(string s) {
    unsigned int i;
    istringstream myStream(s);

    if (myStream>>i) {
        return i;
    } else {
        cout << "String " << s << " is not a number." << endl;
        return atoi(s.c_str());
    }
    return i;
}

/*
* measure the consumed memory
*/
float getCurrentMemoryUsage() {
    string pid = intToStr(unsigned(getpid()));
    string outfile = "tmp_" + pid + ".txt";
    string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
    system(command.c_str());

    string mem_str;
    ifstream ifs(outfile.c_str());
    std::getline(ifs, mem_str);
    ifs.close();

    mem_str = mem_str.substr(0, mem_str.size()-1);
    float mem = (float)strToInt(mem_str);

    command = "rm " + outfile;
    system(command.c_str());

    return mem/1024;
}
//generate time delay
Time_list poisson_time(int len){
    static std::default_random_engine e(time(0));
    static std::poisson_distribution<UI> distribution(2);
    Time_list a;
    for(size_t i=0;i<len;i++){
        UI number = distribution(e)+1;
        a.emplace_back(number);
    }
    return a;
}