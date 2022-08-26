#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <tuple>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mt19937ar.h"
#include "MemoryOperation.h"
#include "include/stats.hpp"

using namespace std;

//#define MeasureTime	1
#define MeasureTime	0

string EdgeFile;
int NodeNum;
double EpsT;	// epsilon
double EpsL;	// epsilon in the local randomizer
double EpsD;	// epsilon for degree (Alg = 4)
double Delta;	// delta
string EpsT_s, Delta_s;
int PairNum;
int ItrNum;
int Alg;
double AlgPrm;
int NumericalBound;
string Alg_s;
int Bip;

// Parameters in the 2-rounds local algorithm [Imola+, USENIX22]
double TClip;
double EClip;
double EpsNsDeg;
double Eps1st;
double Eps2ndTrSt;
double Mu;

// Initialization of statslib
stats::rand_engine_t engine(1776);

FILE *FileOpen(string filename, const char *mode) {
	FILE *fp;

	if ((fp = fopen(filename.c_str(), mode)) == NULL) {
		cout << "cannot open " << filename << endl;
		exit(-1);
	}
	return fp;
}

bool checkFileExistence(const std::string& str) {
    std::ifstream ifs(str);
    return ifs.is_open();
}

// Randomly generate 0, 1, 2, ..., size-1, and store the first num values into rndperm
void MakeRndPerm(int *rndperm, int size, int num) {
	int rnd;
	int *ordperm;
	int i, j;

	// 0, 1, 2, ..., size-1 --> ordperm
	ordperm = (int *)malloc(size * sizeof(int));
	for (i = 0; i < size; i++) {
		ordperm[i] = i;
	}

	for (i = 0; i < num; i++) {
		rnd = genrand_int32() % (size - i);
		rndperm[i] = ordperm[rnd];
		for (j = rnd + 1; j < size - i; j++) {
			ordperm[j - 1] = ordperm[j];
		}
	}

	free(ordperm);
}

// Read edges from the edge file
void ReadEdges(map<int, int> *a_mat, int *node_order){
	int node1, node2;
	int i;
	char s[1025];
	char *tok;
	FILE *fp;
	int type1, type2;

	fp = FileOpen(EdgeFile, "r");
	for(i=0;i<3;i++) fgets(s, 1024, fp);
	while(fgets(s, 1024, fp) != NULL){
		// 1st node --> node1
		tok = strtok(s, ",");
		node1 = atoi(tok);
		// 2nd node --> node2
		tok = strtok(NULL, ",");
		node2 = atoi(tok);
		if(node1 == node2) continue;

		// Make a bipatite graph (type = 1 (odd) or 0 (even))
		if(Bip == 1){
			type1 = node_order[node1] % 2;
			type2 = node_order[node2] % 2;
			if(type1 == type2) continue;
		}

		// If both nodes exist, add the edge
		if(node_order[node1] < NodeNum && node_order[node2] < NodeNum){
			a_mat[node_order[node1]][node_order[node2]] = 1;
			a_mat[node_order[node2]][node_order[node1]] = 1;
		}
	}
	fclose(fp);
}

// Read the numerical upper-bound on epsilon in the local randomizer [Feldman+, FOCS21]
double ReadNumericalBound(int n, double eps, string delta_s){
	double epsl;
	double eps_tmp;
	string indir;
	string infile;
	int i;
	char s[1025];
	char *tok;
	FILE *fp;

	i = EdgeFile.find_last_of("/");
	indir = EdgeFile.substr(0, i+1);
	infile = indir + "numerical-bound_n" + to_string(n) + "_d10-" + delta_s + ".csv";
	if((fp = FileOpen(infile, "r")) == NULL){
		cout << "Cannot open " + infile << endl;
		exit(-1);
	}
	fgets(s, 1024, fp);
	while(fgets(s, 1024, fp) != NULL){
		// epsL --> epsl
		tok = strtok(s, ",");
		epsl = atof(tok);
		// eps_upper --> eps_tmp
		tok = strtok(NULL, ",");
		tok = strtok(NULL, ",");
		eps_tmp = atof(tok);
		if(eps_tmp <= eps) break;
	}
	fclose(fp);

	return epsl;
}

// Calculate the closed-form upper bound on epsilon in the local randomizer [Feldman+, FOCS21]
double CalcEpsL(int n, double eps, double delta){
	double epsl_min = 0.001;
	double epsl_max = 100;
	double epsl_sht = 0.001;
	double epsl;
	double alpha;
	double x1, x2, x3;
	double eps_tmp;

	epsl_max = log((double)n / (16.0 * log(2.0 / pow(10, -delta))));

	for(epsl=epsl_max; epsl>=epsl_min; epsl-=epsl_sht){
		alpha = exp(epsl);
		x1 = (alpha - 1.0) / (alpha + 1.0);
		x2 = 8.0 * sqrt(alpha * log(4.0/pow(10, -delta))) / sqrt((double)n);
		x3 = 8.0 * alpha / (double)n;
		eps_tmp = log(1.0 + x1 * (x2 + x3));
		if(eps_tmp <= eps) break;
	}

	if(epsl > epsl_max) epsl = epsl_max;
//	cout << epsl << " " << epsl_max << endl;

	return epsl;
}

// Calculate #triangles in the centralized model (sensitivity = maximum degree)
void CalcCentTri(long long tri_num, int *deg, double &tri_num_ns, double &sen_tri){
	int max_deg;
	int i;
	FILE *fp;

	// Initialization
	tri_num_ns = tri_num;

	// Sensitivity using max degree --> sen_tri
	// max(deg) --> max_deg
	max_deg = 0;
	for(i=0;i<NodeNum;i++){
		if(max_deg < deg[i]) max_deg = deg[i];
	}
	sen_tri = (double)max_deg;

	// Add Lap(sen_tri/Eps) --> tri_num_ns
	tri_num_ns += stats::rlaplace(0.0, sen_tri/EpsT, engine);
}

// Calculate #2-stars in the centralized model (sensitivity = maximum degree)
void CalcCentSt(long long st2_num, int *deg, double &st2_num_ns, double &sen_st2){
	int max_deg;
	int i;

    // Initialization
    st2_num_ns = st2_num;

	// Sensitivity using max degree --> sen_st2
	sen_st2 = 0.0;
	// max(deg) --> max_deg
	max_deg = 0;
	for(i=0;i<NodeNum;i++){
		if(max_deg < deg[i]) max_deg = deg[i];
	}
	sen_st2 = 2.0 * (double)max_deg;

	// Add Lap(sen_tri/Eps) --> st2_num_ns
	st2_num_ns += stats::rlaplace(0.0, sen_st2/EpsT, engine);
}

// Calculate #triangles in the shuffle model
void CalcShuffleTri(map<int, int> *a_mat, double &tri_num_ns, int *deg, int alg){
	double *deg_ns;
	double avg_deg, avg_deg_ns;
	double deg_thr;
	int pair;
	int *node_order;
	int node1, node2;
	double p2, q2, p1w, q1w;
	int edge12_ns, edge21_ns;
	int wedge, wedge_ns;
	long long tri_num, ed2_num, ed1_num, non_num;
	double c1, e1, e3;
	int deg_sum;
	double edge_prob;
	double rnd;
	int i, j;
	double *a_mat_node1, *a_mat_node2;
	map<int, int>::iterator aitr;
	int pair_num;

	int wedge_num;
	double numerator1, numerator2, denominator;

	// Initialization
	tri_num_ns = 0.0;

	// malloc
	malloc1D(&deg_ns, NodeNum);
	malloc1D(&node_order, NodeNum);
	malloc1D(&a_mat_node1, NodeNum);
	malloc1D(&a_mat_node2, NodeNum);

	if (alg == 3){
		// average degree --> avg_deg
		avg_deg = 0.0;
		for(i=0;i<NodeNum;i++) avg_deg += (double)deg[i];
		avg_deg /= (double)NodeNum;
		// Degree threshold --> deg_thr
		deg_thr = avg_deg * AlgPrm;
	}
	else if (alg == 4){
		// Noisy degree --> deg_ns
		for(i=0;i<NodeNum;i++){
			deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, 1.0/EpsD, engine);
		}
		// Noisy average degree --> avg_deg_ns
		avg_deg_ns = 0.0;
		for(i=0;i<NodeNum;i++) avg_deg_ns += (double)deg_ns[i];
		avg_deg_ns /= (double)NodeNum;
		// Noisy degree threshold --> deg_thr
		deg_thr = avg_deg_ns * AlgPrm;
	}

	// shuffle model
	if (alg == 2 || alg == 3 || alg == 4){
		// Flip probability (wedge) --> q1w
		q1w = 1.0 / (exp(EpsL) + 1.0);
		p1w = 1 - q1w;
	}
	// local model
	else if (alg == 5){
		// Flip probability (wedge) --> q1w
		q1w = 1.0 / (exp(EpsT) + 1.0);
		p1w = 1 - q1w;
	}

	// Flip probability (edge) --> q2
    q2 = 1.0 / (exp(EpsT) + 1.0);
	p2 = 1 - q2;

	// Randomly generate 0, 1, 2, ..., NodeNum-1
	MakeRndPerm(node_order, NodeNum, NodeNum);

	// For each pair of two nodes
	pair = 0;
	for(i=0;i<NodeNum;i+=2){
		// Initialization
		tri_num = ed2_num = ed1_num = non_num = 0;

		// Two nodes --> node1, node2
		node1 = node_order[i];
		node2 = node_order[i+1];

		// Noisy edge between node1 & node2 (RR) --> edge12_ns, edge21_ns
		rnd = genrand_real2();
		// 0 --> 1 (flip)
		if(rnd < q2 && a_mat[node1].count(node2) == 0) edge12_ns = 1;
		// 1 --> 1 (not flip)
		else if(rnd >= q2 && a_mat[node1].count(node2) == 1) edge12_ns = 1;
		else edge12_ns = 0;

		rnd = genrand_real2();
		// 0 --> 1 (flip)
		if(rnd < q2 && a_mat[node2].count(node1) == 0) edge21_ns = 1;
		// 1 --> 1 (not flip)
		else if(rnd >= q2 && a_mat[node2].count(node1) == 1) edge21_ns = 1;
		else edge21_ns = 0;

		// a_mat[node1] --> a_mat_node1
		for(j=0;j<NodeNum;j++) a_mat_node1[j] = 0;
		for (aitr = a_mat[node1].begin(); aitr != a_mat[node1].end(); aitr++) {
			a_mat_node1[aitr->first] = 1;
		}
		// a_mat[node2] --> a_mat_node2
		for(j=0;j<NodeNum;j++) a_mat_node2[j] = 0;
		for (aitr = a_mat[node2].begin(); aitr != a_mat[node2].end(); aitr++) {
			a_mat_node2[aitr->first] = 1;
		}

		// For each node
		wedge_num = 0;
		for(j=0;j<NodeNum;j++){
			if(j == node1 || j == node2) continue;

			// Original wedge --> wedge
			if (a_mat_node1[j] == 1 && a_mat_node2[j] == 1) wedge = 1;
			else wedge = 0;

			// Noisy wedge --> wedge_ns
			rnd = genrand_real2();
			// 0 --> 1 (flip)
			if(rnd < q1w && wedge == 0) wedge_ns = 1;
			// 1 --> 1 (not flip)
			else if(rnd >= q1w && wedge == 1) wedge_ns = 1;
			else wedge_ns = 0;

			// Update #wedges --> wedge_num
			if(wedge_ns == 1) wedge_num++;

			/*
			// *** another way of estimating #triangles (produce the same results) *** //
			// Update #noisy triangles --> tri_num
			if(edge12_ns == 1 && wedge_ns == 1) tri_num++;
			// Update #2-edges --> ed2_num
			if(edge12_ns == 0 && wedge_ns == 1) ed2_num++;
			// Update #1-edges --> ed1_num
			if(edge12_ns == 1 && wedge_ns == 0) ed1_num++;
			// Update #no-edges --> ed1_num
			if(edge12_ns == 0 && wedge_ns == 0) non_num++;

			// Update #noisy triangles --> tri_num
			if(edge21_ns == 1 && wedge_ns == 1) tri_num++;
			// Update #2-edges --> ed2_num
			if(edge21_ns == 0 && wedge_ns == 1) ed2_num++;
			// Update #1-edges --> ed1_num
			if(edge21_ns == 1 && wedge_ns == 0) ed1_num++;
			// Update #no-edges --> ed1_num
			if(edge21_ns == 0 && wedge_ns == 0) non_num++;
			*/
		}
		// The former part of the numerator --> numerator1
		numerator1 = (double)edge12_ns + (double)edge21_ns - 2.0 * q2;
		// The latter part of the numerator --> numerator2
		numerator2 = (double)wedge_num - (double)(NodeNum - 2) * q1w;
		// The denominator --> denominator
		denominator = 2.0 * (2.0 * p2 - 1.0) * (2.0 * p1w - 1.0);
		// Calculate the unbiased estimate of #triangles --> c1
		c1 = numerator1 * numerator2 / denominator;

		/*
		// *** another way of estimating #triangles (produce the same results) *** //
		// Calculate the unbiased estimate of #triangles --> c1
		e1 = (p2 * (double)tri_num + (p2 - 1.0) * (double)ed2_num) / (2.0 * p2 - 1.0);
		e3 = (p2 * (double)ed1_num + (p2 - 1.0) * (double)non_num) / (2.0 * p2 - 1.0);
		c1 = (p1w * (double)e1 + (p1w - 1.0) * (double)e3) / (2.0 * p1w - 1.0);
		*/

		// Shuffle (wedge)
		if(alg == 2 || alg == 5){
			// Increase the unbiased estimate of #triangles --> tri_num_ns
			tri_num_ns += c1;
		}
		// Shuffle (wedge + ignore, degree threshold)
		else if(alg == 3){
			// Increase the unbiased estimate of #triangles only when node1 & node2 are dense --> tri_num_ns
			if(deg[node1] > deg_thr && deg[node2] > deg_thr) tri_num_ns += c1;
		}
		// Shuffle (wedge + ignore, noisy degree threshold)
		else if(alg == 4){
			// Increase the unbiased estimate of #triangles only when node1 & node2 are dense --> tri_num_ns
			if(deg_ns[node1] > deg_thr && deg_ns[node2] > deg_thr) tri_num_ns += c1;
		}

		pair++;
		if(PairNum != -1 && PairNum == pair) break;
	}

	// Number of pairs --> pairnum
	if(PairNum != -1) pair_num = PairNum;
	else pair_num = (int)(NodeNum / 2);

	// Calculate the unbiased estimate of #triangles --> tri_num_ns
	tri_num_ns = tri_num_ns * (double)NodeNum * (double)(NodeNum - 1) / (6.0 * pair_num);

	/*
	// *** another way of estimating #triangles (produce the same results) *** //
	tri_num_ns = tri_num_ns * (double)NodeNum * (double)(NodeNum - 1) / (12.0 * pair_num);
	*/

	free1D(deg_ns);
	free1D(node_order);
	free1D(a_mat_node1);
	free1D(a_mat_node2);
}

// Calculate #triangles in the local model (ARR)
void CalcLocTriARR(map<int, int> *a_mat, double &tri_num_ns){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	int *deg_ns;									// noisy degree
	long long tot_edge_num_ns;
	long long tri_num, st2_num, ed2_num, ed1_num, non_num;
	double tri_num_bs, ed2_num_bs, ed1_num_bs, non_num_bs;
	double q;
	double alp, alp_1_3, q_inv_11, q_inv_21, q_inv_31, q_inv_41;
	double rnd;
	int i, j, k;
	double mu, murho, p1, p2, q2;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	malloc1D(&deg_ns, NodeNum);

	// Parameter in RR --> p1
	p1 = exp(EpsT) / (exp(EpsT) + 1.0);
	// Sampling rate --> p2 (p2^3 = (1/n)^1/3 * Sampling Weight)
	p2 = pow((1.0/(double)NodeNum), 1.0/3.0) * AlgPrm;

	// Parameters in ARR (Asymmetric RR) --> mu (1 --> 1), murho (0 --> 1)
	mu = p1 * p2;
	murho = mu / exp(EpsT);

	// Flip 0/1 in a_mat with probability q --> a_mat_ns
	for(i=0;i<NodeNum;i++){
		for(j=i+1;j<NodeNum;j++){
			rnd = genrand_real2();
			// 0 --> 1 (flip)
			if(rnd < murho && a_mat[i].count(j) == 0){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
			// 1 --> 1 (not flip)
			else if(rnd < mu && a_mat[i].count(j) == 1){
				a_mat_ns[i][j] = 1;
				a_mat_ns[j][i] = 1;
			}
		}
	}

	// Degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) deg_ns[i] += 1;
	}

	// Total number of edges --> tot_edge_num_ns
	tot_edge_num_ns = 0;
	for(i=0;i<NodeNum;i++) tot_edge_num_ns += (long long)deg_ns[i];
	tot_edge_num_ns /= 2;

	// #triangles --> tri_num
	tri_num = 0;
	for(i=0;i<NodeNum;i++){
		for (aitr = a_mat_ns[i].begin(); aitr != a_mat_ns[i].end(); aitr++) {
			j = aitr->first;
			if (i >= j) continue;
			for (aitr2 = a_mat_ns[i].begin(); aitr2 != a_mat_ns[i].end(); aitr2++) {
				k = aitr2->first;
				if (j >= k) continue;
				if(a_mat_ns[j].count(k) > 0) tri_num++;
			}
		}
	}

	// empirical estimation
	// #2-stars --> st2_num
	st2_num = 0;
	for(i=0;i<NodeNum;i++){
		st2_num += ((long long)deg_ns[i] * ((long long)deg_ns[i]-1)) / 2;
	}

	// #2-edges --> ed2_num
	ed2_num = st2_num - 3*tri_num;
	// #1-edge --> ed1_num
	ed1_num = (long long)tot_edge_num_ns*(NodeNum-2) - 2*ed2_num - 3*tri_num;

	// Calculate #triangles, #2-edges, #1-edge before sampling --> tri_num_bs, ed2_num_bs, ed1_num_bs
	q2 = 1.0 - p2;
	tri_num_bs = (double)tri_num / (p2 * p2 * p2);
	ed2_num_bs = (double)ed2_num / (p2 * p2) - 3.0 * q2 * tri_num_bs;
	ed1_num_bs = (double)ed1_num / p2 - 2.0 * q2 * ed2_num_bs - 3.0 * q2 * q2 * tri_num_bs;

	// #none --> non_num_bs
	non_num_bs = (double)NodeNum*(NodeNum-1)*(NodeNum-2)/6 - tri_num_bs - ed2_num_bs - ed1_num_bs;

	alp = exp(EpsT);
	alp_1_3 = (alp-1.0)*(alp-1.0)*(alp-1.0);
	q_inv_11 = (alp*alp*alp) / alp_1_3;
	q_inv_21 = - alp*alp / alp_1_3;
	q_inv_31 = alp / alp_1_3;
	q_inv_41 = - 1.0 / alp_1_3;

	tri_num_ns = tri_num_bs * q_inv_11 + ed2_num_bs * q_inv_21 + ed1_num_bs * q_inv_31 + non_num_bs * q_inv_41;

	delete[] a_mat_ns;
	free1D(deg_ns);
}

double CalcKL(double p1, double p2){
	return p1 * log(p1 / p2) + (1.0 - p1) * log((1.0 - p1) / (1.0 - p2));
}

double CalcRedSen(double deg_ns, int alg){
	double kappa, kappa_min, kappa_max;
	double mu_pow;
	double tclip_prob, tclip_prob_thr;
	double kappa_deg_ns, kl;

	if (alg <= 0 || alg >= 4){
		printf("Error: incorrect alg @ CalcRedSen\n");
		exit(-1);
	}

	if(deg_ns > 0.0){
		if (alg == 1) mu_pow = Mu;
		else mu_pow = Mu * Mu;

		if (alg == 1 || alg == 2) kappa_min = mu_pow * deg_ns;
		else kappa_min = Mu * Mu * Mu * deg_ns;
		kappa_max = deg_ns;
		tclip_prob_thr = pow(10, -TClip);

		// Find kappa s.t. the triangle clipping probability is small enough
		for(kappa = kappa_min; kappa <= kappa_max; kappa += kappa_min){
			// triangle clipping probability --> tclip_prob
			if (alg == 3 && kappa < Mu * Mu * deg_ns) tclip_prob = Mu;
			else{
				kappa_deg_ns = kappa / deg_ns;
				kl = CalcKL(kappa_deg_ns, mu_pow);
				tclip_prob = exp(- deg_ns * kl);
				if(alg == 3) tclip_prob = Mu * tclip_prob;
			}

			if(tclip_prob <= tclip_prob_thr) break;
		}
	}
	else{
		kappa = 0.0;
	}

	return kappa;
}

// Calculate #triangles and #2-stars in the interactive local model (efficient algorithm II)
void CalcLocTri2R(map<int, int> *a_mat, int *deg, double &tri_num_ns, double &sen_tri, double &eclip_sum, double &tclip_sum, int &eclip_num, int &tclip_num){
	map<int, int> *a_mat_ns;			// noisy adjacency matrix
	map<int, int> *a_mat_del;			// deleted adjacency matrix after projection
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<tuple<int, int>, int> a2_mat;
	map<tuple<int, int>, int>::iterator a2_itr;
	double murho;
	double *tri_num_u, *st2_num_u, *trist2_num_u;
	int max_deg;
	int sen_st2_ij;
	int del_num;
	int *rndperm;
	double rnd;
	int i, j, k, x;
	FILE *fp;

	double *deg_ns;
	double *red_sen;
	double tri_num_u_ij;
	double sen;
	int deg_ns_floor;

	// Initialization
	a_mat_ns = new map<int, int>[NodeNum];
	a_mat_del = new map<int, int>[NodeNum];
	malloc1D(&tri_num_u, NodeNum);
	malloc1D(&st2_num_u, NodeNum);
	malloc1D(&trist2_num_u, NodeNum);
	malloc1D(&deg_ns, NodeNum);
	malloc1D(&red_sen, NodeNum);

	// Parameters in asymmetric RR --> Mu (1 --> 1), murho (0 --> 1)
	murho = Mu / exp(Eps1st);

	// Noisy degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, 1.0/EpsNsDeg, engine);
	}
	// Add positive bias (EClip) --> deg_ns
	for(i=0;i<NodeNum;i++) deg_ns[i] += EClip;
	for(i=0;i<NodeNum;i++) deg_ns[i] = max(deg_ns[i], 0.0);

	// Reduced sensitivity --> red_sen
	for(i=0;i<NodeNum;i++) red_sen[i] = CalcRedSen(deg_ns[i], 2);

	// Graph projection (edge clipping) for each user
	for(i=0;i<NodeNum;i++){
		// If deg[i] exceeds deg_ns[i], then perform graph projection (edge clipping)
		if((double)deg[i] > deg_ns[i]){
			eclip_sum += ((double)deg[i] - floor(deg_ns[i]));
			eclip_num += 1;

			// Randomly generate 0, 1, ..., deg[i]-1 --> rndperm
			malloc1D(&rndperm, deg[i]);
			MakeRndPerm(rndperm, deg[i], deg[i]);

			// Randomly delete (deg[i] - floor(deg_ns[i])) edges from a_mat[i]
			deg_ns_floor = (int)floor(deg_ns[i]);
			x = 0;
			for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
				if(rndperm[x] >= deg_ns_floor){
					j = aitr->first;
					// Deleted edge --> a_mat_del[i][j]
					a_mat_del[i][j] = 1;
				}
				x++;
			}
			free1D(rndperm);
		}

		// Count #noisy triangles and #noisy 2-stars --> tri_num_u, st2_num_u
		tri_num_u[i] = 0.0;
		st2_num_u[i] = 0.0;
		for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
			j = aitr->first;
			if (j >= i) continue;
			// Continue if the edge is deleted
			if(a_mat_del[i].count(j) == 1) continue;
			tri_num_u_ij = 0;
			for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
				k = aitr2->first;
				// j < k < i
				if (j >= k || k >= i) continue;
				// Continue if the edge is deleted
				if(a_mat_del[i].count(k) == 1) continue;
				st2_num_u[i] += 1.0;

				// If a_mat_ns[k][i] does not exist
				if(a_mat_ns[k].count(i) == 0){
					// Flip 0/1 in a_mat[k][i] --> a_mat_ns[k][i]
					rnd = genrand_real2();
					// 0 --> 1 (flip)
					if(rnd < murho && a_mat[k].count(i) == 0){
						a_mat_ns[k][i] = 1;
					}
					// 1 --> 1 (not flip)
					else if(rnd < Mu && a_mat[k].count(i) == 1){
						a_mat_ns[k][i] = 1;
					}
					// 1 --> 0 (flip) or 0 --> 0 (not flip)
					else{
						a_mat_ns[k][i] = 0;
					}
				}

				// If a_mat_ns[j][k] does not exist
				if(a_mat_ns[j].count(k) == 0){
					// Flip 0/1 in a_mat[j][k] --> a_mat_ns[j][k]
					rnd = genrand_real2();
					// 0 --> 1 (flip)
					if(rnd < murho && a_mat[j].count(k) == 0){
						a_mat_ns[j][k] = 1;
					}
					// 1 --> 1 (not flip)
					else if(rnd < Mu && a_mat[j].count(k) == 1){
						a_mat_ns[j][k] = 1;
					}
					// 1 --> 0 (flip) or 0 --> 0 (not flip)
					else{
						a_mat_ns[j][k] = 0;
					}
				}

				if(a_mat_ns[k][i] == 1 && a_mat_ns[j][k] == 1) tri_num_u_ij += 1.0;
			}
			// Triangle clipping
			if(tri_num_u_ij > red_sen[i]){
				tclip_sum += (tri_num_u_ij - red_sen[i]);
				tclip_num += 1;
				tri_num_u_ij = red_sen[i];
			}
			tri_num_u[i] += tri_num_u_ij;
		}
		// #triangles - Mu * murho * #2-stars --> trist2_num_u
		trist2_num_u[i] = tri_num_u[i] - Mu * murho * st2_num_u[i];
	}

	// Add Lap for each user
	sen_tri = 0.0;
	for(i=0;i<NodeNum;i++){
		sen_tri += red_sen[i];
		trist2_num_u[i] += stats::rlaplace(0.0, red_sen[i]/Eps2ndTrSt, engine);
	}
	sen_tri /= (double)NodeNum;


    // Empirical estimate --> tri_num_ns
    tri_num_ns = 0;
	for(i=0;i<NodeNum;i++) tri_num_ns += trist2_num_u[i];
	// Divide #triangles by Mu * (Mu - murho)
	tri_num_ns /= (Mu * (Mu - murho));

	delete[] a_mat_ns;
	delete[] a_mat_del;
	free1D(tri_num_u);
	free1D(st2_num_u);
	free1D(trist2_num_u);
	free1D(deg_ns);
	free1D(red_sen);
}

// Calculate #4-cycles in the shuffle model
void CalcShuffleCy4(map<int, int> *a_mat, double &cy4_num_ns, int alg){
	int pair;
	int *node_order;
	int node1, node2;
	double p1w, q1w;
	int wedge, wedge_ns;
	double rnd;
	int i, j;
	double *a_mat_node1, *a_mat_node2;
	map<int, int>::iterator aitr;
	int pair_num;
	int wedge_num;
	double c1;
	double term1, term2;

	// Initialization
	cy4_num_ns = 0.0;

	// malloc
	malloc1D(&node_order, NodeNum);
	malloc1D(&a_mat_node1, NodeNum);
	malloc1D(&a_mat_node2, NodeNum);

	// shuffle model
	if (alg == 8){
		// Flip probability (wedge) --> q1w
		q1w = 1.0 / (exp(EpsL) + 1.0);
		p1w = 1 - q1w;
	}
	// local model
	else if (alg == 9){
		// Flip probability (wedge) --> q1w
		q1w = 1.0 / (exp(EpsT) + 1.0);
		p1w = 1 - q1w;
	}

	// Randomly generate 0, 1, 2, ..., NodeNum-1
	MakeRndPerm(node_order, NodeNum, NodeNum);

	// For each pair of two nodes
	pair = 0;
	for(i=0;i<NodeNum;i+=2){
		// Initialization
		wedge_num = 0;

		// Two nodes --> node1, node2
		node1 = node_order[i];
		node2 = node_order[i+1];

		// a_mat[node1] --> a_mat_node1
		for(j=0;j<NodeNum;j++) a_mat_node1[j] = 0;
		for (aitr = a_mat[node1].begin(); aitr != a_mat[node1].end(); aitr++) {
			a_mat_node1[aitr->first] = 1;
		}
		// a_mat[node2] --> a_mat_node2
		for(j=0;j<NodeNum;j++) a_mat_node2[j] = 0;
		for (aitr = a_mat[node2].begin(); aitr != a_mat[node2].end(); aitr++) {
			a_mat_node2[aitr->first] = 1;
		}

		// For each node
		c1 = 0.0;
		for(j=0;j<NodeNum;j++){
			if(j == node1 || j == node2) continue;

			// Original wedge --> wedge
			if (a_mat_node1[j] == 1 && a_mat_node2[j] == 1) wedge = 1;
			else wedge = 0;

			// Noisy wedge --> wedge_ns
			rnd = genrand_real2();
			// 0 --> 1 (flip)
			if(rnd < q1w && wedge == 0) wedge_ns = 1;
			// 1 --> 1 (not flip)
			else if(rnd >= q1w && wedge == 1) wedge_ns = 1;
			else wedge_ns = 0;

			// Update the unbiased estimate of #wedges --> c1
			c1 += ((double)wedge_ns - (1.0 - p1w)) / (2.0 * p1w - 1.0);
		}

		// Update the unbiased estimate of #4-cycles --> cy4_num_ns
		term1 = (double)c1 * ((double)c1 - 1.0) / 2.0;
		term2 = (((double)NodeNum - 2.0) * p1w * (1.0 - p1w)) / (2.0 * (2.0 * p1w - 1.0) * (2.0 * p1w - 1.0));
		cy4_num_ns += (term1 - term2);

		pair++;
		if(PairNum != -1 && PairNum == pair) break;
	}

	// Number of pairs --> pairnum
	if(PairNum != -1) pair_num = PairNum;
	else pair_num = (int)(NodeNum / 2);

	// Update the unbiased estimate of #4-cycles --> cy4_num_ns
	cy4_num_ns = ((double)NodeNum * ((double)NodeNum - 1.0) / (4.0 * (double)pair_num)) * cy4_num_ns;

	free1D(node_order);
	free1D(a_mat_node1);
	free1D(a_mat_node2);
}

// Calculate #2-stars in the local model [Imola+, USENIX22]
void CalcLocSt(long long st2_num, int *deg, double &st2_num_ns, double &sen_st2){
	int max_deg;
	int i;
	FILE *fp;

	double sen;
	int deg_ns_floor;
	double *deg_ns;

	double eps_nsdeg = EpsT / 10;		// epsilon for calculating the noisy degree
	double eps_allbut_nsdeg = EpsT - eps_nsdeg;
	double eclip = 150.0;

    // Initialization
    st2_num_ns = st2_num;
	malloc1D(&deg_ns, NodeNum);

	// Noisy degree --> deg_ns
	for(i=0;i<NodeNum;i++){
		deg_ns[i] = (double)deg[i] + stats::rlaplace(0.0, 1.0/eps_nsdeg, engine);
	}
	// Add positive bias (eclip) --> deg_ns
	for(i=0;i<NodeNum;i++) deg_ns[i] += eclip;
	for(i=0;i<NodeNum;i++) deg_ns[i] = max(deg_ns[i], 0.0);

	// Average sensitivity using each user's degree --> sen_st2
	sen_st2 = 0.0;
	// Graph projection for each user
	st2_num_ns = 0;
	for(i=0;i<NodeNum;i++){
		deg_ns_floor = (int)floor(deg_ns[i]);
		// If deg[i] exceeds floor(deg_ns[i]), then perform graph projection
		if(deg[i] > deg_ns_floor){
			st2_num_ns += ((long long)deg_ns_floor * ((long long)deg_ns_floor-1)) / 2;
		}
		else{
			st2_num_ns += ((long long)deg[i] * ((long long)deg[i]-1)) / 2;
		}
		// Sensitivity for #2-stars --> sen
		sen = (double)deg_ns_floor;
		// Add Lap(sen/Eps) --> st2_num_ns
		st2_num_ns += stats::rlaplace(0.0, sen/eps_allbut_nsdeg, engine);
		sen_st2 += sen;
	}
	sen_st2 /= (double)NodeNum;

	free1D(deg_ns);
}

// Calculate the clustering-coefficient
double CalcClstCoef(double tri_num_ns, double st2_num_ns){
	double clst_ns;

    if(tri_num_ns < 0) tri_num_ns = 0;
    if(st2_num_ns < 0) st2_num_ns = 0;
    if(st2_num_ns == 0) clst_ns = 1.0;
    else clst_ns = 3.0 * tri_num_ns / st2_num_ns;
    if(clst_ns > 1.0) clst_ns = 1.0;
    else if(clst_ns < 0.0) clst_ns = 0.0;

	return clst_ns;
}

int main(int argc, char *argv[])
{
	int all_node_num;
	int **node_order;
	map<int, int> *a_mat;			// adjacency matrix
	map<int, int>::iterator aitr;
	map<int, int>::iterator aitr2;
	map<int, int>::iterator aitr3;
	int *deg;									// degree
	int *deg_lower;								// degree in the lower-triangular part of a_max
	int max_deg;
	long long tri_num, st2_num;
	long long cy4_num;
	long long pa2_pow2;
	map<int, int> pa2_cnt;
	double clst;
	double tri_num_ns, sen_tri;
	double st2_num_ns, sen_st2;
	double clst_ns;
	double tri_re_ns, tri_l2_ns;
	double tri_re_ns_avg, tri_l2_ns_avg;
	double st2_re_ns, st2_l2_ns;
	double st2_re_ns_avg, st2_l2_ns_avg;
	double clst_re_ns, clst_l2_ns;
	double clst_re_ns_avg, clst_l2_ns_avg;
	double cy4_num_ns;
	double cy4_re_ns, cy4_l2_ns;
	double cy4_re_ns_avg, cy4_l2_ns_avg;
	int fix_perm;
	int itr;
	int i, j, k, l;
	string outdir;
	string outfile;
	char s[1025];
	char *tok;
	FILE *fp;
	clock_t start_clock, end_clock;
	double calc_time, calc_time_avg;
	double eclip_sum, tclip_sum;
	int eclip_num, tclip_num;

	// Initialization of Mersennne Twister
	unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
	init_by_array(init, length);

	if (argc < 2) {
	    printf("Usage: %s [EdgeFile (in)] ([#nodes (default: -1)] [epsilon-delta (default: 1-8)] [#pairs (default:-1)] [#itr(-1) (default: 1)] [alg (default: 0)] ([bi (default: 0)]))]))\n\n", argv[0]);
		printf("[EdgeFile]: Edge file\n");
		printf("[#nodes]: Number of nodes (-1: all)\n");
		printf("[epsilon-delta]: Privacy budgets (delta: -exponent)\n");
		printf("[#pairs]: Number of pairs (-1: all)\n");
		printf("[#itr(-1)]: Number of iterations (set #itr-1 to fix the permutation of nodes)\n");
		printf("[alg]: Algorithm (1: central triangle, 2(n): shuffle triangle (wedge) (n: numerical), 3(n)(-thr (default: 1)): shuffle triangle (wedge + ignore (thr * d_avg)) (n: numerical), 4(n)(-thr (default: 1)): shuffle triangle (wedge + ignore (thr * noisy d_avg) (n: numerical), 5: local triangle (wedge), 6(-smpl weight (default: 1)): local triangle (ARR), 7(-mu* (default: 1)): local triangle (2-rounds), 8(n): shuffle 4-cycle (n: numerical), 9: local 4-cycle)\n");
		printf("[bi]: Make a bipartite graph (1: yes, 0: no)\n");
		return -1;
	}

	EdgeFile = argv[1];

	NodeNum = -1;
	if (argc >= 3) NodeNum = atoi(argv[2]);

	EpsT = 1.0;
	Delta = 8.0;
	EpsT_s = "1";
	Delta_s = "8";
	if (argc >= 4){
		if((tok = strtok(argv[3], "-")) == NULL){
			printf("Error: incorrect [epsilon1-epsilon2]\n");
			exit(-1);
		}
		EpsT = atof(tok);
		EpsT_s = tok;

		if((tok = strtok(NULL, "-")) == NULL){
			printf("Error: incorrect [epsilon1-epsilon2]\n");
			exit(-1);
		}
		Delta = atof(tok);
		Delta_s = tok;
	}

	PairNum = -1;
	if (argc >= 5) PairNum = atoi(argv[4]);

	ItrNum = 1;
	fix_perm = 0;
	if (argc >= 6){
		tok  = strtok(argv[5], "-");
		ItrNum = atoi(tok);
		if((tok  = strtok(NULL, "-")) != NULL){
			if (strcmp(tok, "1") != 0){
				printf("Error: incorrect [#itr(-1)]\n");
				exit(-1);
			}
			else fix_perm = 1;
		}
	}

	Alg = 0;
	AlgPrm = 1;
	NumericalBound = 0;
	Alg_s = "0";
	if (argc >= 7){
		Alg_s = string(argv[6]);
		tok  = strtok(argv[6], "-");
		if(strlen(tok) == 1){
			Alg = atoi(tok);
		}
		else if(strlen(tok) == 2 && tok[1] == 'n'){
			Alg = tok[0] - '0';
			NumericalBound = 1;
		}
		else{
			printf("Error: incorrect [alg]\n");
			exit(-1);
		}

		if((tok  = strtok(NULL, "-")) != NULL){
			AlgPrm = atof(tok);
		}
	}

	Bip = 0;
	if (argc >= 8) Bip = atoi(argv[7]);

	// Privacy budget allocation in shuffle (wedge + ignore (thr * noisy d_avg))
	if (Alg == 4){
		EpsD = 0.1 * EpsT;
		EpsT = EpsT - EpsD;
	}

	// Total number of nodes --> all_node_num
	fp = FileOpen(EdgeFile, "r");
	for(i=0;i<2;i++) fgets(s, 1024, fp);
	all_node_num = atoi(s);
	fclose(fp);

	// malloc
	malloc2D(&node_order, ItrNum, all_node_num);

	// Use all nodes
	if (NodeNum == -1){
		NodeNum = all_node_num;
		for(j=0;j<NodeNum;j++) node_order[0][j] = j;
	}
	// Randomly generate the order of nodes --> node_order
	else{
		i = EdgeFile.find_last_of("/");
		outdir = EdgeFile.substr(0, i+1);
		outfile = outdir + "node-order_itr" + to_string(ItrNum) + ".csv";
		if(checkFileExistence(outfile)){
			fp = FileOpen(outfile, "r");
			for(j=0;j<all_node_num;j++){
				fgets(s, 1024, fp);
				strtok(s, ",");
				for(i=0;i<ItrNum;i++){
					node_order[i][j] = atoi(strtok(NULL, ","));
				}
			}
			fclose(fp);
		}
		else{
			for(i=0;i<ItrNum;i++){
				MakeRndPerm(node_order[i], all_node_num, all_node_num);
			}
			fp = FileOpen(outfile, "w");
			for(j=0;j<all_node_num;j++){
				fprintf(fp, "%d,", j);
				for(i=0;i<ItrNum;i++) fprintf(fp, "%d,", node_order[i][j]);
				fprintf(fp, "\n");
			}
			fclose(fp);
		}

		// Use only the first permutation
		if (fix_perm){
			for(j=0;j<all_node_num;j++){
				for(i=1;i<ItrNum;i++) node_order[i][j] = node_order[0][j];
			}
		}
	}

	// Calculate epsilon in the local randomizer from EpsT [Feldman+, arXiv21] --> EpsL
	if(NumericalBound == 1){
		EpsL = ReadNumericalBound(NodeNum, EpsT, Delta_s);
	}
	else{
		EpsL = CalcEpsL((double)NodeNum-2, EpsT, Delta);
	}
	if (EpsL < EpsT) EpsL = EpsT;

	// Initialization
	malloc1D(&deg, NodeNum);
	malloc1D(&deg_lower, NodeNum);
	tri_re_ns_avg = tri_l2_ns_avg = 0.0;
	st2_re_ns_avg = st2_l2_ns_avg = 0.0;
	clst_re_ns_avg = clst_l2_ns_avg = 0.0;
	cy4_re_ns_avg = cy4_l2_ns_avg = 0.0;
	calc_time_avg = 0.0;

	// Output the header
	i = EdgeFile.find_last_of("/");
	outdir = EdgeFile.substr(0, i+1);
	if(Bip == 1){
		if(fix_perm) outfile = outdir + "res_n" + to_string(NodeNum) + "b_alg" + Alg_s + "_eps" + EpsT_s + "-" + Delta_s + "_pair" + to_string(PairNum) + "_itr" + to_string(ItrNum) + "-1.csv";
		else outfile = outdir + "res_n" + to_string(NodeNum) + "b_alg" + Alg_s + "_eps" + EpsT_s + "-" + Delta_s + "_pair" + to_string(PairNum) + "_itr" + to_string(ItrNum) + ".csv";
	}
	else{
		if(fix_perm) outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + Alg_s + "_eps" + EpsT_s + "-" + Delta_s + "_pair" + to_string(PairNum) + "_itr" + to_string(ItrNum) + "-1.csv";
		else outfile = outdir + "res_n" + to_string(NodeNum) + "_alg" + Alg_s + "_eps" + EpsT_s + "-" + Delta_s + "_pair" + to_string(PairNum) + "_itr" + to_string(ItrNum) + ".csv";
	}
	fp = FileOpen(outfile, "w");
	// triangles
	if (Alg <= 7){
		if(MeasureTime) fprintf(fp, "#tri(true),#tri(est),#tri(rel-err),#tri(l2-loss),#2st(true),#2st(est),#2st(rel-err),#2st(l2-loss),clst(true),clst(est),clst(rel-err),clst(l2-loss),EpsL,max_deg,calctime\n");
		else fprintf(fp, "#tri(true),#tri(est),#tri(rel-err),#tri(l2-loss),#2st(true),#2st(est),#2st(rel-err),#2st(l2-loss),clst(true),clst(est),clst(rel-err),clst(l2-loss),EpsL,max_deg\n");
	}
	// 4-cycles
	else{
		fprintf(fp, "#4cyc(true),#4cyc(est),#4cyc(rel-err),#4cyc(l2-loss),EpsL,max_deg\n");
	}
	fclose(fp);

	// For each iteration
	for(itr=0;itr<ItrNum;itr++){
		// Read edges for each iteration when NodeNum < all_node_num
		if(NodeNum < all_node_num || itr == 0){
			// Initialization
			a_mat = new map<int, int>[NodeNum];

			// Read edges from the edge file --> a_mat
			ReadEdges(a_mat, node_order[itr]);

			// Degree --> deg
			for(i=0;i<NodeNum;i++) deg[i] = 0;
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) deg[i] += 1;
			}

			// Degree --> deg_lower
			for(i=0;i<NodeNum;i++) deg_lower[i] = 0;
			for(i=0;i<NodeNum;i++){
				for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++){
					if(aitr->first < i) deg_lower[i] += 1;
				}
			}

			// max(deg) --> max_deg
			max_deg = 0;
			for(i=0;i<NodeNum;i++){
				if(max_deg < deg[i]) max_deg = deg[i];
			}

			// #2-stars --> st2_num
			st2_num = 0;
			for(i=0;i<NodeNum;i++){
				st2_num += ((long long)deg[i] * ((long long)deg[i]-1)) / 2;
			}

			// Count triangles/4-cycles in the original graph
			if (fix_perm == 0 || itr == 0){
				// triangles
				if (Alg <= 7){
					// #triangles --> tri_num
					tri_num = 0;
					for(i=0;i<NodeNum;i++){
						for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
							j = aitr->first;
							if (i >= j) continue;
							for (aitr2 = a_mat[i].begin(); aitr2 != a_mat[i].end(); aitr2++) {
								k = aitr2->first;
								if (j >= k) continue;
								if(a_mat[j].count(k) > 0) tri_num++;
							}
						}
					}

					// clustering coefficient --> clst
					if(st2_num != 0) clst = 3.0 * (double)tri_num / (double)st2_num;
					else clst = 1.0;
				}
				// 4-cycles
				else{
					// #4-cycles --> cy4_num
					/*
					// Simple method (inefficient)
					// Count i--j--k--l--i (i,j,k,l are distinct)
					cy4_num = 0;
					for(i=0;i<NodeNum;i++){
						for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
							j = aitr->first;
							if (i==j) continue;
							for (aitr2 = a_mat[j].begin(); aitr2 != a_mat[j].end(); aitr2++) {
								k = aitr2->first;
								if (i==k || j==k) continue;
								for (aitr3 = a_mat[k].begin(); aitr3 != a_mat[k].end(); aitr3++) {
									l = aitr3->first;
									if (i==l || j==l || k==l) continue;
									if (a_mat[l].count(i) == 1){
	//									cout << cy4_num<< ":" << i << "-" << j << "-" << k << "-" << l << endl;
										cy4_num++;
									}
								}
							}
						}
					}
					// We count each 4-cycle 8 times (we have 4 start nodes and 2 directions), so divide cy4_num by 8
					cy4_num /= 8;
					cout << "#4-cycles (simple counting):" << cy4_num << endl;
					*/

					// 2-path method (efficient and output the same value as the simple method)
					// #2-path from i to k --> pa2_cnt[k]
					cy4_num = 0;
					pa2_pow2 = 0;
					// i--j--k (i,j,k are distinct)
					for(i=0;i<NodeNum;i++){
						for (aitr = a_mat[i].begin(); aitr != a_mat[i].end(); aitr++) {
							j = aitr->first;
							if (i==j) continue;
							for (aitr2 = a_mat[j].begin(); aitr2 != a_mat[j].end(); aitr2++) {
								k = aitr2->first;
								if (i==k || j==k) continue;
								// Update #2-path from i to k --> pa2_cnt[k]
								if(pa2_cnt.count(k) == 0) pa2_cnt[k] = 1;
								else pa2_cnt[k] += 1;
							}
						}
						// Update the square of #2-paths from i to k (i < k) --> pa2_pow2
						for (aitr = pa2_cnt.begin(); aitr != pa2_cnt.end(); aitr++) {
							k = aitr->first;
							if (i>=k) continue;
							pa2_pow2 += pa2_cnt[k] * pa2_cnt[k];
						}
						pa2_cnt.clear();
					}
					// We use the fact that (#2-paths)^2 = #2-stars + 4 * #4-cycles
					cy4_num = (pa2_pow2 - st2_num) / 4;
				}
			}
		}

		/************************ Calculate sub-graph counts ************************/
		// Central
		if (Alg == 1){
        	// Calculate #triangles
			CalcCentTri(tri_num, deg_lower, tri_num_ns, sen_tri);
			// Calculate #2-stars
        	CalcCentSt(st2_num, deg, st2_num_ns, sen_st2);
		}
		// Shuffle (wedge)
		else if (Alg == 2){
        	// Calculate #triangles
			start_clock = clock();
			CalcShuffleTri(a_mat, tri_num_ns, deg, Alg);
			end_clock = clock();
			calc_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
			// Calculate #2-stars
        	CalcLocSt(st2_num, deg, st2_num_ns, sen_st2);
		}
		// Shuffle (wedge + ignore, thr: degree threshold)
		else if (Alg == 3){
        	// Calculate #triangles
			start_clock = clock();
			CalcShuffleTri(a_mat, tri_num_ns, deg, Alg);
			end_clock = clock();
			calc_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
			// Calculate #2-stars
        	CalcLocSt(st2_num, deg, st2_num_ns, sen_st2);
		}
		// Shuffle (wedge + ignore, thr: noisy degree threshold)
		else if (Alg == 4){
        	// Calculate #triangles
			start_clock = clock();
			CalcShuffleTri(a_mat, tri_num_ns, deg, Alg);
			end_clock = clock();
			calc_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
			// Calculate #2-stars
        	CalcLocSt(st2_num, deg, st2_num_ns, sen_st2);
		}
		// Local (wedge)
		else if (Alg == 5){
        	// Calculate #triangles
			start_clock = clock();
			CalcShuffleTri(a_mat, tri_num_ns, deg, Alg);
			end_clock = clock();
			calc_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
			// Calculate #2-stars
        	CalcLocSt(st2_num, deg, st2_num_ns, sen_st2);
		}
		// Local (ARR)
		else if (Alg == 6){
        	// Calculate #triangles
			start_clock = clock();
			CalcLocTriARR(a_mat, tri_num_ns);
			end_clock = clock();
			calc_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
			// Calculate #2-stars
        	CalcLocSt(st2_num, deg, st2_num_ns, sen_st2);
		}
		// Local (2-rounds)
		else if (Alg == 7){
			// Epsilon
			EpsNsDeg = 0.1 * EpsT;
			Eps1st = 0.45 * EpsT;
			Eps2ndTrSt = 0.45 * EpsT;

			// Clipping parameters
			TClip = 6.0;
			EClip = 150.0;

			// assign mu*(=mu^2) --> Mu
			Mu = AlgPrm;
			// mu^2 = exp(Eps1st)/(exp(Eps1st)+1) when it is set to 1 --> Mu
			if(Mu == 1.0) Mu = exp(Eps1st) / (exp(Eps1st) + 1);
			// Calculate mu from mu^2 --> Mu
			Mu = pow(Mu, 1.0/2.0);

			// Calculate #2-stars
        	CalcLocSt(st2_num, deg, st2_num_ns, sen_st2);
        	// Calculate #triangles (efficient algorithm II)
			CalcLocTri2R(a_mat, deg_lower, tri_num_ns, sen_tri, eclip_sum, tclip_sum, eclip_num, tclip_num);
		}
		// Shuffle 4-cycle
		else if (Alg == 8){
        	// Calculate #4-cycles
			start_clock = clock();
			CalcShuffleCy4(a_mat, cy4_num_ns, Alg);
			end_clock = clock();
			calc_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
		}
		// Local 4-cycle
		else if (Alg == 9){
        	// Calculate #4-cycles
			start_clock = clock();
			CalcShuffleCy4(a_mat, cy4_num_ns, Alg);
			end_clock = clock();
			calc_time = (double)(end_clock - start_clock) / CLOCKS_PER_SEC;
		}

		// triangles
		if (Alg <= 7){
			/******************** Calculate the cluster coefficient *********************/
			clst_ns = CalcClstCoef(tri_num_ns, st2_num_ns);

			/**************************** Evaluate the loss *****************************/
			// relative error --> tri_re_ns
			tri_re_ns = fabs(tri_num_ns - (double)tri_num) / max((double)tri_num, 0.001 * NodeNum);
			tri_re_ns_avg += tri_re_ns;
			// l2_loss --> tri_l2_ns
			tri_l2_ns = (tri_num_ns - (double)tri_num)*(tri_num_ns - (double)tri_num);
			tri_l2_ns_avg += tri_l2_ns;

			// relative error --> st2_re_ns
			st2_re_ns = fabs(st2_num_ns - (double)st2_num) / max((double)st2_num, 0.001 * NodeNum);
			st2_re_ns_avg += st2_re_ns;
			// l2_loss --> st2_l2_ns
			st2_l2_ns = (st2_num_ns - (double)st2_num)*(st2_num_ns - (double)st2_num);
			st2_l2_ns_avg += st2_l2_ns;

			// relative error --> clst_re_ns
			clst_re_ns = fabs(clst_ns - clst) / (double)clst;
			clst_re_ns_avg += clst_re_ns;
			// l2_loss --> clst_l2_ns
			clst_l2_ns = (clst_ns - clst)*(clst_ns - clst);
			clst_l2_ns_avg += clst_l2_ns;

			// Calculation time
			calc_time_avg += calc_time;

			/**************************** Output the results ****************************/
			fp = FileOpen(outfile, "a");
			if(MeasureTime){
				fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%f,%d,%e\n", 
				(double)tri_num, tri_num_ns, tri_re_ns, tri_l2_ns, 
				(double)st2_num, st2_num_ns, st2_re_ns, st2_l2_ns, 
				clst, clst_ns, clst_re_ns, clst_l2_ns, 
				EpsL, max_deg, calc_time);
			}
			else{
				fprintf(fp, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%f,%d\n", 
				(double)tri_num, tri_num_ns, tri_re_ns, tri_l2_ns, 
				(double)st2_num, st2_num_ns, st2_re_ns, st2_l2_ns, 
				clst, clst_ns, clst_re_ns, clst_l2_ns, 
				EpsL, max_deg);
			}
			fclose(fp);
		}
		// 4-cycles
		else{
			/**************************** Evaluate the loss *****************************/
			// relative error --> cy4_re_ns
			cy4_re_ns = fabs(cy4_num_ns - (double)cy4_num) / max((double)cy4_num, 0.001 * NodeNum);
			cy4_re_ns_avg += cy4_re_ns;
			// l2_loss --> cy4_l2_ns
			cy4_l2_ns = (cy4_num_ns - (double)cy4_num)*(cy4_num_ns - (double)cy4_num);
			cy4_l2_ns_avg += cy4_l2_ns;

			// Calculation time
			calc_time_avg += calc_time;

			/**************************** Output the results ****************************/
			fp = FileOpen(outfile, "a");
			if(MeasureTime){
				fprintf(fp, "%e,%e,%e,%e,%f,%d,%e\n", 
				(double)cy4_num, cy4_num_ns, cy4_re_ns, cy4_l2_ns, 
				EpsL, max_deg, calc_time);
			}
			else{
				fprintf(fp, "%e,%e,%e,%e,%f,%d\n", 
				(double)cy4_num, cy4_num_ns, cy4_re_ns, cy4_l2_ns, 
				EpsL, max_deg);
			}
			fclose(fp);
		}

		if(NodeNum < all_node_num || itr == ItrNum - 1){
			delete[] a_mat;
		}
	}

	/************************* Output the results (AVG) *************************/
	// triangles
	if (Alg <= 7){
		tri_re_ns_avg /= (double)ItrNum;
		tri_l2_ns_avg /= (double)ItrNum;
		st2_re_ns_avg /= (double)ItrNum;
		st2_l2_ns_avg /= (double)ItrNum;
		clst_re_ns_avg /= (double)ItrNum;
		clst_l2_ns_avg /= (double)ItrNum;
		calc_time_avg /= (double)ItrNum;

		fp = FileOpen(outfile, "a");
		fprintf(fp, "function,AVG(rel-err),AVG(l2-loss)\n");
		fprintf(fp, "Triangles,%e,%e\n", tri_re_ns_avg, tri_l2_ns_avg);
		fprintf(fp, "2-stars,%e,%e\n", st2_re_ns_avg, st2_l2_ns_avg);
		fprintf(fp, "Clst,%e,%e\n", clst_re_ns_avg, clst_l2_ns_avg);
		if(MeasureTime) fprintf(fp, "CalcTime,%e\n", calc_time_avg);
		fclose(fp);
	}
	// 4-cycles
	else{
		cy4_re_ns_avg /= (double)ItrNum;
		cy4_l2_ns_avg /= (double)ItrNum;
		calc_time_avg /= (double)ItrNum;

		fp = FileOpen(outfile, "a");
		fprintf(fp, "function,AVG(rel-err),AVG(l2-loss)\n");
		fprintf(fp, "4-cycles,%e,%e\n", cy4_re_ns_avg, cy4_l2_ns_avg);
		if(MeasureTime) fprintf(fp, "CalcTime,%e\n", calc_time_avg);
		fclose(fp);
	}

	// free
	free2D(node_order, ItrNum);
	free1D(deg);
	free1D(deg_lower);

	return 0;
}
