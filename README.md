<h3>
	<center>Time-aware Influence Minimization via Blocking Social Networks</center>
</h3>

### Information

Version 1.0: Implementation of Algorithm for  Time-aware Influence Minimization in Social Networks. For more details about our code, please read our paper: "Xueqin C., Jiajie F., Qing L., Yunjun G., Baihua Z., Lu C., Time-aware Influence Minimization via Blocking Social Networks"

### Introduction

1. This repository contains the full version of our paper.
2. This repository contains the codes and datasets used in our paper.
3. **Time-aware Influence Minimization via Blocking Social Networks**.

Abstract: We study the problem of Time-aware Influence Minimization (Timin) in social networks, aiming to minimize the negative influence concerning a critical deadline by temporarily blocking some nodes of the given social network. To this end, first, we introduce a novel Time-delayed Linear Threshold (TLT) model by considering the time delay of influence, i.e., when a node is active, its out-neighbors receive the influence weight after a certain time delay. Building on the TLT model, we formally define the Timin problem, and prove that it is NP-hard, monotone, and supermodular. To tackle the Timin problem, we initially devise a basic greedy algorithm, Timin-Greedy, achieving ($1-1/e$) approximation. Since it is #P-hard to compute the exact negative influence spread for any node set in Timin-Greedy, we devise a Temporal Reverse Influence Sampling technique to estimate the expected negative influence spread, and propose a more efficient algorithm TESTIM, maintaining ($1-1/e-\epsilon$) approximation. To further improve the efficiency, we propose a heuristic algorithm NeighborReplace based on an important observation that potential blocking nodes are often located near the negative source. Furthermore, we investigate two variants of the Timin problem, which consider additional constraints. Finally, our extensive experiments demonstrate that (1) TESTIM is up to 10$\times$ faster than the baselines, yielding 30%-50% more negative influence spread reductions, and (2) compared with TESTIM, NeighborReplace exhibits 5$\times$ speedup while having comparable negative influence spread reductions.
### Datasets

We use six publicly available real-world road networks, including EmailCore, Epinions, Amazon, Youtube, FaceBook and LiveJournal datasets. 

All of them can be obtained from [1].

[1] Jure Leskovec and Andrej Krevl. 2014. SNAP Datasets: Stanford large network dataset collection. http://snap.stanford.edu/data.

### Algorithms

The following files are the codes for our proposed algorithms. We implemented all the codes using C++ with CLion 2022.3.2.

1. First we use el2bin.cpp$^{[2]}$ (can be found in genSeed folder) to convert from a graph file in weighted edge list to two binary files that encode the graph and its transpose;

```shell
 ./el2bin <input file> <output file> <transpose output file>
```

2. Then we use fake_seeds.cpp to generate fake_seeds (random or influential). Specifically, when generating influential seeds, we set -m top, and the -f parameter means that the fake seeds will be randomly drawn from the top f-th fraction for generating influential fake seeds, the orders of the nodes are determined by the  singleton influence file. When choosing random seeds, just use four parameters (-n -o -k and -m), while setting -m random. The usages are listed as follows (the first line is for influential seeds and the second is for random seeds):

```shell
./fake_seeds -n <number of nodes> -o <seed output file>  -k <number of seeds> -m <top> -f <fraction> -s <singleton influence file>
```

```shell
./fake_seeds -n <number of nodes> -o <seed output file>  -k <number of seeds> -m <random> 
```

3. Then use **algorithm** (can be found in TESTIM directory) to tackle our problem, **algorithm** includes:

- TLT:  A TM-loss Greedy function for finding the best seeds set for temporal influence minimization;
- Advanced_tlt: The heurist way for solving the problem;
- deadline_solution: For solving variant problem DSTIMIN;
- Minimal_block_solution: For solving variant problem BSTIMIN.

   The usages are listed as follows:

```shell
 ./TLT -i <input networkFile> -o <result output file> -fakeseeds <fakeSeed file> -k <budget of blockSet> -epsilon <xx> -t <max edge delay> -T <deadline> -lamda <edge delay parameters> -delta <xx>
```

```shell
 ./Advanced_tlt -i <input networkFile> -o <result output file> -fakeseeds <fakeSeed file> -k <budget of blockSet> -t <max edge delay> -T <deadline> -lamda <edge delay parameters> 
```

```shell
./deadline_solution -i <input networkFile> -o <result output file> -fakeseeds <fakeSeed file> -k <budget of blockSet> -t <max edge delay> -T <deadline> -lamda <edge delay parameters> -preFile <file for determine the max_T and max_alpha> -alpha <the percentage of users influenced>
```

```shell
./Minimal_block_solution -i <input networkFile> -o <result output file> -fakeseeds <fakeSeed file> -t <max edge delay> -T <deadline> -lamda <edge delay parameters> -alpha <the percentage of users influenced>
```

[2] Michael Simpson, Farnoosh Hashemi, and Laks VS Lakshmanan. 2022. Misinformation mitigation under differential propagation rates and temporal penalties. Proceedings of the VLDB Endowment 15, 10 (2022), 2216–2229.

### Running Environment

A 64-bit Linux-based OS. 

GCC 4.7.2 and later.