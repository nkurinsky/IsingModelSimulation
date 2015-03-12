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

//limit on correlation probe number
#define SITE_LIMIT 1000
//define for clarity whether to calculate means
#define CALCMEAN true

using namespace std;

typedef char spin;

//data structure to store location of correlation probes
struct probe{
  probe(); //basic constructor to initialize totalPoints and correlationMean
  vector<int> zeros;           //location of r=0
  vector<vector<int> > points; //location of r=r point for each distance
  double totalPoints;          //probes per dimension times dimensions
  double correlationMean;      //last calculated mean of correlation vector
  vector<double> correlation;  //stores previously calculated values
};

//data structure to handle indexing throuh multidimensional array efficiently
class location{
public:
  location(int ndim, int size, int index=-1); //general constructor
  location(const location &obj);              //copy constructor
  void randomize();                           //set all indices randomly between 0 and size-1
  void set_index(int index);                  //set location based on 1d index
  void move(int dim=0, int distance=0);       //move relative to current position
  int index(int dim=0, int distance=0);       //get index relative to current position
  void indices(vector<int> array);            //obtain index array
  void print();                               //output information about position
  void neighbor(location &nb, int dim, int distance); //get neighboring location structure
  int ndim() const { return _ndim; };
 private:
  int _ndim;
  int _size;
  int _index_max;
  vector<int> _indices;    //ndimensional length vector of indices
  vector<int> _dim_step;   //1d base for each dimension, conversion to ndim from 1d
}; 

//data structure for storing and simulating spin lattice for arbitrary dimension, size and spin
class lattice{
public:
  lattice(int ndim, int size, char q=2); //general constructor
  void randomize();                      //randomly set each lattice site to randomSpin()
  void setTemp(double T);                //set temperature for clustering algorithm
  spin randomSpin();                     //generate random integer in (1,q)
  bool addBond();                        //wolf algorithm bond generation probability
  spin get(location site) const;         //get spin from location structure
  spin get(int index) const;             //get spin from 1d index in _lattice
  void flip(location site, spin newValue); //flip location to new spin value
  void flip(int index, spin newValue);     //flip index to new spin value
  double magnetization(bool mean=false);   //calculate lattice order parameter, optionally calc mean
  //calculate lattice correlation, optionally calculating correlation mean
  double correlation(int length, bool mean=false);
  bool correlated(location site1, location site2);    //return whether spins are equal
  int flipCluster(int Nflips=1);                      //run wolf algorith for N minimum spin flips
  unsigned long flips() const{ return _total_flips;}; //total number of spins flipped since initialization
  void display(FILE *outfile=stdout) const;           //output current lattice to file or terminal
private:
  double chainMean(vector<double> chain); //calculate mean of latter half of values (throw away earlier length)
  double cos_LatticeAngle(spin value);    //calculate cos(2pi*s/q)
  double sin_LatticeAngle(spin value);    //calculate sin(2pi*s/q)
  valarray<spin> _lattice;                //lattice storage (lightweight array class)
  double _T;                              //simulation temperature
  float _pBond;                           //bond creation probability
  short _ndim;                            
  char _q;
  short _size;
  unsigned long _total_flips;
  mutable unsigned long _display_calls;   //number of times lattice has been displayed (to differentiate)
  //magnetization variables
  double _spinAngle;                      //2*pi/q
  vector<double> magValues;               //vector of calculated magnetizations
  unsigned long _magLastCheckSize;        //number of spins flipped at last calculation
  double _lastMagMean;                    //last calculated magnetization mean
  //correlation variables
  unordered_map<int,probe> probes;        //storage for correlation probes
};

#endif
