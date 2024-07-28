#include <cstdio>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, char ** argv)
{
	ifstream in(argv[1]);
	unsigned int flag = atoi(argv[3]);

	unsigned long long n,m,i;
	in >> n >> m;
	printf("%lld, %lld\n", n, m);
	vector<unsigned long long> node(n,0);
	vector<unsigned long long> v1, v2;
	v1.reserve(m);
	v2.reserve(m);
	
	unsigned long long t1,t2;
	in >> t1 >> t2;
	v1.push_back(t1);
	v2.push_back(t2);
	node[t2]++;
	if (flag == 0) node[t1]++;
	printf("Reading the graph!\n");
	for (i = 1; i < m; i++) {
		in >> t1 >> t2;
		node[t2]++;
		if (flag == 0) node[t1]++;

		v1.push_back(t1);
		v2.push_back(t2);
		if (i %100000 == 0) printf("%lld\n", i);
	}
	in.close();
	ofstream out(argv[2]);
	printf("Writing down to file!\n");
	out << n << " ";

	if (flag == 0) {
		out << v1.size()*2 << endl;
	} else{
		out << v1.size() << endl;
	}
    for (i = 0; i < v1.size(); i++) {
        out << v1[i] << " " << v2[i] << " " << (float) 1.0/(float)node[v2[i]] << endl;
        if (flag == 0) out << v2[i] << " " << v1[i] << " " << (float) 1.0/(float)node[v1[i]] << endl;
    }

	out.close();
}
