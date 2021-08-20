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
	double ovpfrac;


	fragopt_t () 
	{
		CleanMaxDiag = 100;
		DiagMinCluster = 500;
		clusterMaxDiag = 20;
		clusterMaxDist = 200; // 800
		clusterMinLength = 1000;
		// CleanMaxDiag = 100;
		// DiagMinCluster = 50;
		// clusterMaxDiag = 1;
		// clusterMaxDist = 1000;
		// clusterMinLength = 500;
		clusterTrimedge = 400;
		code = 0;
		debug = true;
		freq = 50;
		ovpfrac = 0.6;
	};

	
};

#endif
