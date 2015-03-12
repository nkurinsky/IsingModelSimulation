#include <cstdio>
#include <cstdlib>
#include <thread>
#include <iostream>
#include <unistd.h>
#include <valarray>

#include "lattice.hpp"
#include "modelstats.hpp"

//set number of models to use for convergence testing
#define NMODELS 5
#define N4DMODELS 3
#define N5DMODELS 3
#define CORR_MAX 10

//upper run limit for convergence monitoring loop
#define MAX_ITERATIONS 100000

//variable used to enable/disable lattice displat
#define DISPLAYGRID

using namespace std;

//function to call display function for each lattice
#ifdef DISPLAYGRID
void outputLattices(vector<lattice> &models,FILE * outfile){
  for(unsigned int i=0;i<models.size();i++)
    models[i].display(outfile);
}
#endif /* DISPLAYGRID */

int main(int argc, char *argv[]){

  //I/O code for reading in size, dimension, q and temperature, as well as optional output file
  if(argc < 5){
    printf("Invalid Argument Number\nCalling Sequence: %s Dim Size Q T [outfile] \n",argv[0]);
    exit(1);
  }

  int dimension=atoi(argv[1]);
  int size=atoi(argv[2]);
  char q=static_cast<char>(atoi(argv[3]));
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
  if(dimension < 1){
    printf("Invalid number of dimensions (Entered %i, range 1-5)\n",dimension);
    exit(3);
  }

  //set number of models based on dimension (save memory and time in higher dimensions)
  if(dimension > 3){
    if(dimension == 4){
      modelnum=N4DMODELS;
	}
    else if(dimension == 5){
      modelnum=N5DMODELS;
	}
    else{
      fprintf(stderr,"Dimension number %i too large, will max out memory\n",dimension);
      exit(3);
    }
  }

  int N=pow(size,dimension); //number of sites in lattice
  //ncorr set by overall lattice size, do not compute more than size/2
  int ncorr = (size/2 > CORR_MAX) ? CORR_MAX : size/2 ; //number of correlation values to compute
  modelValues vals(modelnum,ncorr); //data structure to store and print statistics
  vals.setLimit(0.01); //set convergence criterion

  fprintf(outfile,"Simulation Parameters:\n");
  fprintf(outfile,"\tDimensions:   %i\n",dimension);
  fprintf(outfile,"\tLattice Size: %i\n",size);
  fprintf(outfile,"\tPotts States: %i\n",q);
  fprintf(outfile,"\tTemperature:  %f\n",T);
  fprintf(outfile,"\tLattices:     %i\n\n",modelnum);

  //initial lattices by randomizing and running 10 full monte carlo steps
  vector<lattice> models;
  models.clear();
  for(int i=0;i<modelnum;i++){
    models.push_back(lattice(dimension,size,q));
    models[i].setTemp(T);
    models[i].flipCluster(10*N);
  }    
  
  //output initial lattice statistics
  vals.printHeader(outfile);
  vals.populate(models);
  vals.print(outfile);

  //begin simulation loop
  int iter=0;
  while((not vals.converged(outfile)) and (iter < MAX_ITERATIONS)){
    iter++;
    for(unsigned int i=0;i<models.size();i++){
      //flip clusters until at least N spins have flipped
      //here N is the number of lattice sites in the model
      models[i].flipCluster(N);
    }
    vals.populate(models); //get correlation and magnetization
    vals.print(outfile); //output current values to stdout or file
  }
  //end while loop

  if(iter < MAX_ITERATIONS) //end condition was convergence
    fprintf(outfile,"\nConverged!\n");
  else //end condition was timeout
    fprintf(outfile,"\nDid not converge in %i iterations\n",MAX_ITERATIONS);
  
  //output final values
  vals.print(outfile);
  vals.printResults(T,outfile);

  //displace final lattice
#ifdef DISPLAYGRID
  outputLattices(models,outfile);
#endif /* DISPLAYGRID */
  
  if(filewrite){
    fclose(outfile);
  }
  
  return 0;
}
