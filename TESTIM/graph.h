#ifndef _GRAPH_H
#define _GRAPH_H
#include <vector>
#include <unordered_set>
#include <random>
#include <functional>
#include <queue>
#include "sfmt/SFMT.h"

typedef uint32_t UI;
typedef std::vector<uint32_t> Time_list;
class Graph
{
	friend class HyperGraph;
	private:
		UI UI_MAX = 4294967295U;
		// number of nodes
		unsigned int num_nodes;
		// number of edges
		unsigned int num_edges;
		// fake seeds
		std::vector<int> fake_seeds_index;
		unsigned int num_fake_seeds;
		std::vector<bool> fake_seed;
		// dynamics parameters
		int _t; //the max time of each edge
        int _T; //the deadline
		// adjacency lists
		std::vector<std::vector<int> > adj_list;
		std::vector<std::vector<int> > rev_adj_list;
		std::vector<int> node_deg;
		std::vector<int> rev_node_deg;
		std::vector<std::vector<float> > weights;
		std::vector<std::vector<float> > rev_weights;
        std::vector<double> visit_thresh_hold;
	
	public:
        sfmt_t sfmt_seed;
		Graph(int t,int T);
		// get a vector of neighbours of node u in G
		const std::vector<int> & getOutNeighbours(int u) const;
		// get a vector of neighbours of node u in G^T
		const std::vector<int> & getInNeighbours(int u) const;
		// return weights of neighbours of node u in G
		const std::vector<float> & getOutWeights(int u) const;
		// return weights of neighbours of node u in G^T
		const std::vector<float> & getInWeights(int u) const;
		// return a vector of the fake seeds of G
		const std::vector<int> & getFakeSeeds() const;
		// get threshold of node u
        const std::vector<double> & getThreshold() const;
        //get out degree of node u
		unsigned int getOutDegree(int u) const;
        //init threshold for LT model
        void init_visit_thresh_hold();
		// get in degree of node u
		unsigned int getInDegree(int u) const;
		// get size of the graph
		unsigned int getNumNodes() const;
		// get number of edges
		unsigned int getNumEdges() const;
		// get number of fake seeds
		unsigned int getNumFakeSeeds() const;
		// determine if node v is a fake seed
		bool isFakeSeed(int v) const;
		// read graph from a file
		void readGraph(const char * filename, bool fixed, float w);
		// read fake seeds from a file
		void readFakeSeeds(const char* filename);
        //get t
        int getMaxEdgeDelay();
        //get T
        int getDeadline();
        // set t
        void setMaxEdgeDelay(int t);
        //set T
        void setDeadline(int T);
		// compute mitigation lower bound via depth 1 MIA
		double computeMitigationLowerBound(unsigned int n, unsigned int k);
		// generate a single forward monte carlo estimate of the influence of input seed
		int generateInfluenceSample(std::vector<bool> &visit, std::vector<int> &visit_index, std::vector<double> &threshold, std::vector<int> &dist, int root);
		// generate a single forward monte carlo estimate of the influence of F
		int generateFakeInfluenceSample(std::vector<bool> &visit, std::vector<int> &visit_index, std::vector<int> &dist, std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > &pq);
        // generate a single forward monte carlo estimate of the influence of F/B
        int generateBlockSampleGreedy(std::vector<bool> &visit, std::vector<int> &visit_index, std::vector<double> &threshold, std::vector<int> &dist, const std::vector<int> &block, std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > &pq);
        int liveEdgeSampleGreedy(std::vector<bool> &visit, std::vector<int> &visit_index, std::vector<int> &dist, const std::vector<int> &block, std::priority_queue<std::pair<int,int>, std::vector<std::pair<int,int>>, std::greater<std::pair<int,int>> > &pq);
};

class HyperGraph
{
	private:
        sfmt_t sfmtSeed;
		// store the hyperedges that a node is incident on together with the corresponding coverage
		std::vector<std::vector<std::pair<int,int> > > node_hyperedges;
		// store hyperedges
		std::vector<std::vector<std::pair<int,int> > > hyperedges;
		// store the weight of hyperedges that a node is covered by
		std::vector<int> node_hyperedge_weight;

		
	public:
        std::vector<std::vector<int > > hyperG;
        std::vector<std::vector<int > > hyperGB;
        std::vector<std::vector<int > > hyperGT;
        double EstPro1, EstPro2, eps1, eps2;
        //time possibility f(t) from f(1) to  f(16) under default setting
        // f(t) = t^(-0.02)
//        std::vector<double> f_t = {0, 1, 0.9862327, 0.9782674, 0.9726550, 0.9682380, 0.9647993, 0.9618294,
//                               0.9592641,0.9570071, 0.9549926, 0.9531739, 0.9515166, 0.9499946, 0.9485876, 0.9472796, 0.9460577};
        // f(t) = t^(-0.2)
        std::vector<double> f_t = {0, 1, 0.8705505, 0.80274156, 0.757858283, 0.724779663, 0.698827118772, 0.677610913,
                                   0.6597539553,0.6443940149};
        //unsigned int cur_hyperedge;
		HyperGraph(unsigned int n);
		void addEdge(std::vector<std::pair<int,int> > & edge);
		const std::vector<std::pair<int,int> > & getEdge(int e) const;
		const std::vector<std::pair<int,int> > & getNode(int u) const;
		int getNodeWeight(int u) const;
        int getNumEdge() const;
		void clearEdges();


};

float getCurrentMemoryUsage();
//generate time_delay
Time_list poisson_time(int len);
// decay function caused by time delay
double decayFunc(int len, float alpha = 0.0);

#endif
