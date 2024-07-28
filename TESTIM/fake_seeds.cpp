#include "option.h"
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

using namespace std;

void generateSeedsTop(int n, int k, float f, const char* infFile, vector<int> &seeds) {
	ifstream in(infFile);
	int i, v1, rand_pos;
	float v2;
	int num = n * f;
	vector<int> candidates(num,0);
	srand (time(NULL));

	for (i = 0; i < num; i++) {
		in >> v1 >> v2;
		rand_pos = rand() % (i+1);
		if (rand_pos != i) candidates[i] = candidates[rand_pos];
		candidates[rand_pos] = v1;
	}
	in.close();

	for (i = 0; i < k; i++) {
		seeds[i] = candidates[i];
	}
}

void generateSeedsTopAgain(int n, int k, float f, const char* infFile, char* infile2, vector<int> &seeds) {
    ifstream in(infFile);
    ifstream in2(infile2);
    int generated_num, id;
    map<int,bool>already_choose; //mark the nodes that are already generate
    in2 >> generated_num;
    int i, v1, rand_pos, j;
    for (j = 0; j < generated_num; j++){
        in2 >> id;
        already_choose[id] = true;
        seeds[j] = id;
    }

    float v2;
    int num = n * f;
    vector<int> candidates(num,0);
    srand (time(NULL));

    for (i = 0; i < num; i++) {
        in >> v1 >> v2;
        rand_pos = rand() % (i+1);
        if (rand_pos != i) candidates[i] = candidates[rand_pos];
        candidates[rand_pos] = v1;
    }
    in.close();

    for (i = generated_num, j = 0; i < k; j++) {
        if (already_choose[candidates[j]] == false){
            seeds[i] = candidates[j];
            i++;
        }
    }
}

void generateSeedsRandom(int n, int k, vector<int> &seeds) {
	int i, v, rand_pos;
	vector<int> candidates(n,0);
	srand (time(NULL));

	for (i = 0; i < n; i++) {
		candidates[i] = i;
	}

	for (i = n-1; i > 0; i--) {
		rand_pos = rand() % (i+1);
		v = candidates[rand_pos];
		candidates[rand_pos] = candidates[i];
		candidates[i] = v;
	}

	for (i = 0; i < k; i++) {
		seeds[i] = candidates[i];
	}
}

void generateSeedsRandomAgain(int n, int k, char* infile2, vector<int> &seeds) {
    ifstream in2(infile2);
    int generated_num, id;
    map<int,bool>already_choose; //mark the nodes that are already generate
    in2 >> generated_num;
    int i, v, rand_pos, j;
    for (j = 0; j < generated_num; j++){
        in2 >> id;
        already_choose[id] = true;
        seeds[j] = id;
    }

    float v2;
    vector<int> candidates(n,0);
    srand (time(NULL));
    for (i = 0; i < n; i++) {
        candidates[i] = i;
    }

    for (i = n-1; i > 0; i--) {
        rand_pos = rand() % (i+1);
        v = candidates[rand_pos];
        candidates[rand_pos] = candidates[i];
        candidates[i] = v;
    }

    for (i = generated_num, j = 0; i < k; j++) {
        if (already_choose[candidates[j]] == false){
            seeds[i] = candidates[j];
            i++;
        }
    }
}

void generateSeedsComm(Graph &g, int n, int k, const char* commFile, vector<int> &seeds) {
	return;
}

int main(int argc, char ** argv)
{
	srand(time(NULL));
	
	OptionParser op(argc, argv);
	if (!op.validCheck()){
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}

	char * outFile = op.getPara("-o");
	if (outFile == NULL){
		outFile = (char*)"fake.seeds";
	}

	int n = -1;
	char * tmp = op.getPara("-n");
	if (tmp != NULL){
		n = atoi(tmp);
	}
	if (n < 0) {
		printf("Must provide graph size!");
		return -1;
	}

	int k = 1;
	tmp = op.getPara("-k");
	if (tmp != NULL){
		k = atoi(tmp);
	}

	float f = 0.1;
	tmp = op.getPara("-f");
	if (tmp != NULL){
		f = atof(tmp);
	}

	char * method = op.getPara("-m");
	if (method == NULL){
		method = (char*)"top";
	}

	char * infFile = op.getPara("-s");
	if (infFile == NULL){
		infFile = (char*)"singleton.inf";
	}

    char * infFile2 = op.getPara("-a");
    if (infFile == NULL){
        infFile = (char*)"already.seeds";
    }

	vector<int> seeds(k,0);

	if (strcmp(method, "top") == 0) {
		generateSeedsTop(n, k, f, infFile, seeds);
	} else if (strcmp(method, "random") == 0) {
		generateSeedsRandom(n, k, seeds);
	} else if (strcmp(method, "ta") == 0) {
        generateSeedsTopAgain(n, k, f, infFile, infFile2, seeds);
    } else if (strcmp(method, "ra") == 0) {
        generateSeedsRandomAgain(n, k, infFile2, seeds);
    } else {
		printf("Incorrect method option!");
		return -1;
	}

    ofstream of(outFile);
    of << k << endl;
	for (int i = 0; i < k; i++){
		of << seeds[i] << endl;
	}
	of.close();

	//
	// ALL DONE
	//
}