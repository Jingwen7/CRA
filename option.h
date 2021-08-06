#ifndef OPTIONS_H_
#define OPTIONS_H_

class idxopt_t {
public:
	uint32_t k; 
	uint32_t w;

	idxopt_t (uint32_t K, uint32_t W) : k(K), w(W) {};

};

class fragopt_t {
public:
	int CleanMaxDiag;
	int DiagMinCluster;
	int clusterMaxDiag;
	int clusterMaxDist;
	int clusterMinLength;
	int clusterTrimedge;
	uint32_t code;
	bool debug;
	int freq;


	fragopt_t () 
	{
		CleanMaxDiag = 100;
		DiagMinCluster = 100;
		clusterMaxDiag = 500;
		clusterMaxDist = 1000;
		clusterMinLength = 1000;
		clusterTrimedge = 200;
		code = 0;
		debug = true;
		freq = 50;


	};

	
};

#endif
