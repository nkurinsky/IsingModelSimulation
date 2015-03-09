#include <cstdio>
#include <cstdlib>
#include <thread>
#include <iostream>
#include <unistd.h>
#include <valarray>

#include "lattice.hpp"
#include "modelstats.hpp"

#define NMODELS 5
#define CORR_MAX 10

#define N4DMODELS 3
#define N5DMODELS 3

#define MAX_THREADS 5
#define MAX_ITERATIONS 10000

using namespace std;

void outputLattices(vector<lattice> &models,FILE * outfile){
  for(unsigned int i=0;i<models.size();i++)
    models[i].display(outfile);
}

int main(int argc, char *argv[]){

  if(argc < 5){
    printf("Invalid Argument Number\nCalling Sequence: %s Dim Size Q T [outfile] \n",argv[0]);
    exit(1);
  }

  int dimension=atoi(argv[1]);
  int size=atoi(argv[2]);
  short q=static_cast<short>(atoi(argv[3]));
  double T=atof(argv[4]);
  bool filewrite=false;
  FILE * outfile=stdout;
  if(argc > 5){
    filewrite=true;
    outfile=fopen(argv[5],"w+");
    if(outfile == NULL){
      printf("Could not open output file, exiting\n");
      exit(2);
    }
  }
  
  int modelnum=NMODELS;
  if(dimension > 3){
    if(dimension == 4){
      modelnum=N4DMODELS
	}
    else if(dimension == 5){
      modelnum=N5DMODELS
	}
    else{
      fprintf(stderr,"Dimension number %i too large, will max out memory\n");
      exit(3);
    }
  }

  int N=pow(size,dimension);
  int ncorr = (size/2 > CORR_MAX) ? CORR_MAX : size/2 ;
  modelValues vals(modelnum,ncorr);
  vals.setLimit(0.01);

  fprintf(outfile,"Simulation Parameters:\n");
  fprintf(outfile,"\tDimensions:   %i\n",dimension);
  fprintf(outfile,"\tLattice Size: %i\n",size);
  fprintf(outfile,"\tPotts States: %i\n",q);
  fprintf(outfile,"\tTemperature:  %f\n",T);
  fprintf(outfile,"\tLattices:     %i\n\n",modelnum);

  vector<lattice> models;
  models.clear();
  for(int i=0;i<modelnum;i++){
    models.push_back(lattice(dimension,size,q));
    models[i].setTemp(T);
    models[i].flipCluster(10*N);
  }    
  
  vals.printHeader(outfile);
  vals.populate(models);
  vals.print(outfile);

  int iter=0;
  while((not vals.converged(outfile)) and (iter < MAX_ITERATIONS)){
    iter++;
    for(unsigned int i=0;i<models.size();i++){
      models[i].flipCluster(N);
    }
    vals.populate(models);
    vals.print(outfile);
  }
  if(iter < MAX_ITERATIONS)
    fprintf(outfile,"\nConverged!\n");

  vals.print(outfile);
  vals.printResults(T,outfile);

  if(filewrite){
    fclose(outfile);
  }
  
  return 0;
}
