#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <valarray>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <unordered_map>
#include <random>
#include <memory>
#include <functional>
#include <stack>
#include <cmath>

#define SITE_LIMIT 1000
#define CALCMEAN true

using namespace std;

typedef short spin;

struct probe{
  probe();
  vector<int> zeros;
  vector<vector<int> > points;
  double totalPoints;
  double correlationMean;
  vector<double> correlation;
};

class location{
public:
  location(int ndim, int size, int index=-1);
  location(const location &obj);
  void randomize();
  void set_index(int index);
  void move(int dim=0, int distance=0);
  int index(int dim=0, int distance=0);
  void indices(vector<int> array);
  void print();
  void neighbor(location &nb, int dim, int distance);
  const int ndim() const { return _ndim; };
 private:
  int _ndim;
  int _size;
  int _index_max;
  vector<int> _indices;
  vector<int> _dim_step;
}; 

class lattice{
public:
  lattice(int ndim, int size, short q=2);
  void randomize();
  void setTemp(double T);
  spin randomSpin();
  bool addBond();
  spin get(location site) const;
  spin get(int index) const;
  void flip(location site, spin newValue);
  void flip(int index, spin newValue);
  double magnetization(bool mean=false);
  double correlation(int length, bool mean=false);
  bool correlated(location site1, location site2);
  int flipCluster(int Nflips=1);
  unsigned long flips() const{ return _total_flips;};
  void display(FILE *outfile=stdout) const;
private:
  double chainMean(vector<double> chain);
  double cos_LatticeAngle(spin value);
  double sin_LatticeAngle(spin value);
  valarray<spin> _lattice;
  double _T;
  float _pBond;
  short _ndim;
  short _q;
  short _size;
  unsigned long _total_flips;
  mutable unsigned long _display_calls;
  //magnetization variables
  double _spinAngle;
  vector<double> magValues;
  unsigned long _magLastCheckSize;
  double _lastMagMean;
  //correlation variables
  unordered_map<int,probe> probes;
};

#endif
