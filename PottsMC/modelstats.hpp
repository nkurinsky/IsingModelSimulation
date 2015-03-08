#ifndef MODELSTATS_H
#define MODELSTATS_H

#include <vector>
#include <valarray>
#include <cstdio>

#include "lattice.hpp"

using namespace std;

class modelValues{
public:
  modelValues(int nmodels, int ncorr);
  void populate(vector<lattice> &models);
  void setLimit(double limit);
  bool converged(FILE *outfile=stdout);
  void printResults(double temp, FILE *outfile=stdout);
  void print(FILE *outfile=stdout);
  void printHeader(FILE *outfile=stdout);
  valarray<double> magnetization;
  vector<valarray<double> > correlation;
private:
  double _limit;
  double _zerolimit;
};

#endif
