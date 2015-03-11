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

typedef char spin;
typedef vector<int> coords;

struct probe{
  probe();
  vector<coords> zeros;
  vector<vector<coords> > points;
  double totalPoints;
  double correlationMean;
  vector<double> correlation;
};

class location{
public:
  location(int ndim, int size);
  location(const location &obj);
  void randomize();
  void move(coords array);
  void move(location &obj);
  void move(int dim=0, int distance=0);
  coords get(int dim=0, int distance=0);
  void indices(vector<int> array);
  void print();
  int ndim() const { return _ndim; };
 private:
  int _ndim;
  int _size;
  coords _indices;
}; 

class lattice{
public:
  lattice(int ndim, int size, char q=2);
  void randomize();
  void setTemp(double T);
  spin randomSpin();
  bool addBond();
  spin get(location site);
  //spin get(int index) const;
  void flip(location site, spin newValue);
  //void flip(int index, spin newValue);
  double magnetization(bool mean=false);
  double correlation(int length, bool mean=false);
  bool correlated(location site1, location site2);
  int flipCluster(int Nflips=1);
  unsigned long flips() const{ return _total_flips;};
  void display(FILE *outfile=stdout);
private:
  double chainMean(vector<double> chain);
  double cos_LatticeAngle(spin value);
  double sin_LatticeAngle(spin value);
  map<coords,spin> _lattice;
  double _T;
  float _pBond;
  short _ndim;
  char _q;
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
