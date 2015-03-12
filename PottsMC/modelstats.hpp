#ifndef MODELSTATS_H
#define MODELSTATS_H

#include <vector>
#include <valarray>
#include <cstdio>

#include "lattice.hpp"

using namespace std;

//structure to compute and store statistics and compute convergence
class modelValues{
public:
  modelValues(int nmodels, int ncorr);    //general constructor
  void populate(vector<lattice> &models); //calculate values from models
  void setLimit(double limit);            //set convergence limit
  bool converged(FILE *outfile=stdout);   //test for convergence and output results
  void printResults(double temp, FILE *outfile=stdout); //print current stored values
  void print(FILE *outfile=stdout);
  void printHeader(FILE *outfile=stdout); //print output file header
  //statistic storage variables
  valarray<double> magnetization;
  vector<valarray<double> > correlation;
private:
  double _limit;
  //zerolimit denotes minimum value considered non-zero and prevents returning negative numbers
  double _zerolimit;
};

#endif
