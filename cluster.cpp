#include "rfpriv.h"
#include "assert.h"
#include "hit.h"
#include "cluster.h"
#include <set>
#include <tuple>
#include <utility>
#include <map>

bool clustDiagonalSort (const cluster &a, const cluster &b) {
	int64_t aDiag = (int64_t) a.xStart - (int64_t) a.yStart;
	int64_t bDiag = (int64_t) b.xStart - (int64_t) b.yStart;
	if (aDiag != bDiag)
		return aDiag < bDiag;
	else
		return a.xStart < b.xStart;
}

bool clustAntiDiagonalSortOp(const cluster &a, const cluster &b) {
	int64_t aDiag = (int64_t) a.xStart + (int64_t) a.yStart;
	int64_t bDiag = (int64_t) b.xStart + (int64_t) b.yStart;
	if (aDiag != bDiag)
		return aDiag < bDiag;
	else
		return a.xStart < b.xStart;
}

void storeDiagCluster (vector<hit> &hits, clusters &clust, bool st, idxopt_t &iopt, fragopt_t &fopts) 
{
	if (st == 0)
		sort(hits.begin(), hits.end(), hitDiagonalSort);
	else
		sort(hits.begin(), hits.end(), hitAntiDiagonalSort);

	uint32_t n = hits.size();
	uint32_t cs = s, ce = s;
	while (cs < n) {
		ce = cs + 1;
		uint32_t  xStart = hits[cs].x, 
				  xEnd = hits[cs].x + iopt.k, 
			      yStart = hits[cs].y, 
			      yEnd = hits[cs].y + iopt.k;

		while (ce < n and abs(DiagonalDifference(hits[ce], hits[ce - 1], st)) < fopts.clusterMaxDiag) {
			xStart = min(xStart, hits[ce].x);
			xEnd   = max(xEnd, hits[ce].x + iopt.k);
			yStart = min(yStart, hits[ce].y);
			yEnd   = max(yEnd, hits[ce].y + iopt.k);
			ce++;
		}	

		if (ce - cs >= fopts.DiagMinCluster and xEnd - xStart >= fopts.clusterMinLength and yEnd - yStart >= fopts.clusterMinLength) 
			clust.push_back(cluster(cs, ce, xStart, xEnd, yStart, yEnd, st));
		cs=ce;
	}	

	if (fopts.debug) {
		ofstream rclust("diagcluster.dots");
		for (int m = 0; m < clust.size(); m++) {
			for (int c = clust[m].s; c < clust[m].e; ++c) {
				rclust << hits[c].x << "\t" << hits[c].y << "\t" << hits[c].x + iopt.k << "\t"
					   << hits[c].y + iopt.k << "\t" << m << "\t" << clust[m].strand << "\t" << iopt.k << endl;				
			}
		}
		rclust.close();		
	}
}

// void mergeDiagCluster(clusters &clusts, clusters &reclusts, int64_t mergeDiag, bool st)
// {
// 	uint32_t n = clusts.size();
// 	uint32_t cs = s, ce = s;
// 	while (cs < n) {
// 		ce = cs + 1;
// 		uint32_t  xStart = clusts[cs].xStart, 
// 				  xEnd = clusts[cs].xEnd, 
// 			      yStart = clusts[cs].yStart, 
// 			      yEnd = clusts[cs].yEnd;

// 		while (ce < n and abs(DiagonalDifference(clusts[ce], clusts[ce - 1], st)) < mergeDiag) {
// 			xStart = min(xStart, clusts[cs].xStart);
// 			xEnd   = max(xEnd, clusts[cs].xEnd);
// 			yStart = min(yStart, clusts[cs].yStart);
// 			yEnd   = max(yEnd, clusts[cs].yEnd);
// 			ce++;
// 		}	

// 		reclusts.push_back(cluster(cs, ce, xStart, xEnd, yStart, yEnd, st));
// 		cs=ce;
// 	}	
// 	// clusts.clear();
// 	// clusts = reclusts;
// 	// clusts.resize(reclusts.size());
// }

// void clustersContained(cluster &a, cluster &b, const fragopt_t &fopt, state &s) 
// {
// 	if ((a.xstart >= b.xstart and a.xEnd <= b.xEnd) or (b.xstart >= a.xstart and b.xEnd <= a.xEnd)) // contained on x
// 		return cg_x;
// 	else if ((a.ystart >= b.ystart and a.yEnd <= b.yEnd) or (b.ystart >= a.ystart and b.yEnd <= a.yEnd)) // contained on y
// 		return cg_y;
// 	else if ((int64_t) a.xstart >= (int64_t) b.xstart - (int64_t) fopt.clusterTrimedge and (int64_t) a.xEnd <= (int64_t) b.xEnd + (int64_t) fopt.clusterTrimedge) 
// 	     or ((int64_t) b.xstart >= (int64_t) a.xstart - (int64_t) fopt.clusterTrimedge and (int64_t) b.xEnd <= (int64_t) a.xEnd + (int64_t) fopt.clusterTrimedge)  // shifted contained on x
// 	    return sft_x;
// 	else if (((int64_t) a.ystart >= (int64_t) b.ystart - (int64_t) fopt.clusterTrimedge and (int64_t) a.yEnd <= (int64_t) b.yEnd + (int64_t) fopt.clusterTrimedge) 
// 	      or ((int64_t) b.ystart >= (int64_t) a.ystart - (int64_t) fopt.clusterTrimedge and (int64_t) b.yEnd <= (int64_t) a.yEnd + (int64_t) fopt.clusterTrimedge)) // shifted contained on y
// 	    return sft_y;
// 	return None;

// }

bool sortPoint (const pointInfo &a, const pointInfo &b)
{
	if ((get<0>(a) == s and get<0>(b) == s) or (get<0>(a) == e and get<0>(b) == e)) {
		return get<1>(a) < get<1>(b);
	}
	else if (get<0>(a) == s and get<0>(b) == e) {
		return 1;
	}
	else if (get<0>(a) == e and get<0>(b) == s) {
		return 0;
	}
	return 0;
}

bool sortFragmatch (const fragmatch &a, const fragmatch &b) 
{
	if (a.s != b.s) {
		return a.s < b.s;
	}
	else {
		cerr << " Two fragmatch has the same starts!";
		return a.e < b.e;		
	}
}	 

void trim(clusters &clusts, const vector<pointInfo> &Spoints, const vector<pointInfo> &Epoints, int clusterTrimedge) {
	uint32_t i, j, k;
	uint32_t diff;
	i = 1;
	while (i < Spoints.size()) {
		uint32_t range = get<1>(Spoints[i - 1]) + clusterTrimedge;
		j = i;
		while (j < Spoints.size()) {
			if (get<1>(Spoints[j]) <= range) {
				j += 1;
			}
		}
		if (j > i) { // [i - 1, j) are close
			for (k = i - 1; k < j - 1; ++k) {
				diff = clusts[j - 1].yStart > clusts[k].yStart ? clusts[j - 1].yStart - clusts[k].yStart : 0;
				clusts[k].yStart = clusts[j - 1].yStart;
				clusts[k].xStart += diff;
			}
			i = j + 1;
		}
		else {i += 1;}
	}	

	i = 1;
	while (i < Epoints.size()) {
		uint32_t range = get<1>(Epoints[i - 1]) + clusterTrimedge;
		j = i;
		while (j < Epoints.size()) {
			if (get<1>(Epoints[j]) <= range) {
				j += 1;
			}
		}
		if (j > i) { // [i - 1, j) are close
			for (k = i; k < j; ++k) {
				diff = clusts[i - 1].yEnd < clusts[k].yEnd ? clusts[k].yEnd - clusts[i - 1].yEnd : 0;
				clusts[k].yEnd = clusts[i - 1].yEnd;
				clusts[k].xEnd -= diff;
			}
			i = j + 1;
		}
		else {i += 1;}
	}	
}

void trim_ovpClusters(clusters &clusts, int clusterTrimedge)
{
	uint32_t i, j, k = 0;
	int diff = 0;

	vector<pointInfo> ySPoints, yEPoints;

	// trim based on ystart/yEnd points
	for (i = 0; i < clusts.size(); ++i) {
		ySPoints.push_back(make_tuple(s, clusts[i].yStart, i));
		yEPoints.push_back(make_tuple(e, clusts[i].yEnd, i));
	}
	sort(ySPoints.begin(), ySPoints.end(), sortPoint);
	sort(yEPoints.begin(), yEPoints.end(), sortPoint);
	trim(clusts, ySPoints, yEPoints, clusterTrimedge);

	vector<pointInfo> xSPoints, xEPoints;

	// trim based on xstart/xEnd points
	for (i = 0; i < clusts.size(); ++i) {
		xSPoints.push_back(make_tuple(s, clusts[i].xStart, i));
		xEPoints.push_back(make_tuple(e, clusts[i].xEnd, i));
	}
	sort(xSPoints.begin(), xSPoints.end(), sortPoint);
	sort(xEPoints.begin(), xEPoints.end(), sortPoint);
	trim(clusts, xSPoints, xEPoints, clusterTrimedge);

	// trim based on xstart/xEnd points again
	xSPoints.clear(); xEPoints.clear();
	for (i = 0; i < clusts.size(); ++i) {
		xSPoints.push_back(make_tuple(s, clusts[i].xStart, i));
		xEPoints.push_back(make_tuple(e, clusts[i].xEnd, i));
	}
	sort(xSPoints.begin(), xSPoints.end(), sortPoint);
	sort(xEPoints.begin(), xEPoints.end(), sortPoint);
	trim(clusts, xSPoints, xEPoints, clusterTrimedge);	
}

void splitClusters(clusters & clusts, clusters & splitclusts) 
{
	set<uint32_t> qSet, tSet;
	//
	// insert x/y coordinates of each cluster into qSet/tSet;
	//
	uint32_t m = 0;
	for (m = 0; m < clusts.size(); ++m) {
		qSet.insert(clusts[m].xStart);
		qSet.insert(clusts[m].xEnd);
		tSet.insert(clusts[m].yStart);
		tSet.insert(clusts[m].yEnd);				
	} 

	//
	// Find what coordinates appear in the interval of each cluster
	//
	for (m = 0; m < clusts.size(); m++) {
	
		IntervalSet itlSet(clusts[m]);
		set<uint32_t>::iterator its, ite;
		its = qSet.upper_bound(clusts[m].xStart);
		ite = qSet.lower_bound(clusts[m].xEnd); // this points to clusts[m].yEnd

		for (set<uint32_t>::iterator it = its; it != ite; it++) {
			itlSet.Set.push_back(make_pair(*it, 0)); 
		}

		its = tSet.upper_bound(clusts[m].yStart);
		ite = tSet.lower_bound(clusts[m].yEnd);

		for (set<uint32_t>::iterator it = its; it != ite; it++) {
			itlSet.Set.push_back(make_pair(*it, 1)); 
		}		

		itlSet.Sort();

		//
		// Split clusts[m] 
		//
		pair<uint32_t, uint32_t> prev;
		if (clusts[m].strand == 0) prev = make_pair(clusts[m].xStart, clusts[m].yStart);
		else prev = make_pair(clusts[m].xStart, clusts[m].yEnd); 

		vector<pair<uint32_t, bool>>::iterator it = itlSet.Set.begin();
		for (; it < itlSet.Set.end(); it++) {

			if (it->second == 0) { // split on x coord.
				uint32_t t = (uint32_t) ceil(itlSet.slope * it->first + itlSet.intercept);

				if (prev.first < it->first) { 
					if (clusts[m].strand == 0 and it->first >= prev.first + 3 and t >= prev.second + 3)
						splitclusts.push_back(cluster(prev.first, it->first, prev.second, t, clusts[m].strand, m)); // initialize coarse to specify the index of original index
					else if (clusts[m].strand ==  1 and it->first >= prev.first + 3 and prev.second >= t + 3)
						splitclusts.push_back(cluster(prev.first, it->first, t, prev.second, clusts[m].strand, m));	
					// cerr << "splitclusts: " << splitclusts.size() << " " << prev.first << " " << it->first << endl;				
				}
				else continue;
				prev = make_pair(it->first, t);					
			}
			else { // split on y coord.
				uint32_t q = (uint32_t) ceil((it->first - itlSet.intercept) / itlSet.slope);
				
				if (prev.first < q) {
					if (clusts[m].strand == 0 and q >= prev.first + 3 and it->first >= prev.second + 3) 
						splitclusts.push_back(cluster(prev.first, q, prev.second, it->first, clusts[m].strand, m));
					else if (clusts[m].strand ==  1 and q >= prev.first + 3 and prev.second >= it->first + 3)
						splitclusts.push_back(cluster(prev.first, q, it->first, prev.second, clusts[m].strand, m));
					// cerr << "splitclusts: " << splitclusts.size() << " " << prev.first << " " << q << endl;
				}
				else continue;
				prev = make_pair(q, it->first);					
			}

		} 

		if (prev.first < clusts[m].xEnd) {
			if (clusts[m].strand == 0 and clusts[m].xEnd >= prev.first + 3 and clusts[m].yEnd >= prev.second + 3) {
				splitclusts.push_back(cluster(prev.first, clusts[m].xEnd, prev.second, clusts[m].yEnd, clusts[m].strand, m));
			}
			else if (clusts[m].strand ==  1 and clusts[m].xEnd >= prev.first + 3 and prev.second >= clusts[m].yStart + 3){
				splitclusts.push_back(cluster(prev.first, clusts[m].xEnd, clusts[m].yStart, prev.second, clusts[m].strand, m));	
			}		
			// cerr << "splitclusts: " << splitclusts.size() << " " << prev.first << " " << clusters[m].xEnd << endl;
		}
	}
}

void fragLabel(const clusters &clusts, fragopt_t &fopts)
{
	uint32_t i;
	vector<fragmatch> fms;
	map<pair<uint32_t, uint32_t>, uint32_t> coords; // (s, e) -> code
	fragmatch f1, f2;

	for (i = 0; i < clusts.size(); ++i) {
		f1 = {clusts[i].xStart, clusts[i].xEnd, fopts.code, i};
		f2 = {clusts[i].yStart, clusts[i].yEnd, fopts.code, i};
		pair<uint32_t, uint32_t> p1 = make_pair(f1.s, f1.e);
		pair<uint32_t, uint32_t> p2 = make_pair(f2.s, f2.e);

		if (!coords.empty()) {
			auto it1 = coords.find(p1);
			auto it2 = coords.find(p2);
			if (it1 == coords.end() and it2 == coords.end()) {
				coords[p1] = fopts.code;
				coords[p2] = fopts.code;
				fms.push_back(f1); fms.push_back(f2);
				fopts.code += 1;
			}
			else if (it1 != coords.end()) {
				f2.code = it1->second;
				fms.push_back(f2);
			}
			else if (it2 != coords.end()) {
				f1.code = it2->second;
				fms.push_back(f1);
			}
		}
		else {
			coords[p1] = fopts.code;
			coords[p2] = fopts.code;
			fms.push_back(f1); fms.push_back(f2);
			fopts.code += 1;
		}
	}

	sort(fms.begin(), fms.end(), sortFragmatch);

	if (fopts.debug) {
		ofstream fclust("Fragmatch.bed");
		for (uint32_t m = 0; m < fms.size(); ++m) {
			fclust << fms[m].s << "\t" << fms[m].e << "\t" << fms[m].code << "\t"
					<< fms[m].cidx << "\t" << m << endl;
		}
		fclust.close();			
	}
}


