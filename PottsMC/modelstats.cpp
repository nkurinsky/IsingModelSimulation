#include "modelstats.hpp"

modelValues::modelValues(int nmodels, int ncorr){
  magnetization.resize(nmodels);
  correlation.resize(ncorr,valarray<double>(nmodels));
  _limit=0.01;
}

void modelValues::populate(vector<lattice> &models){
  for(unsigned long i=0;i<models.size();i++){
    magnetization[i]=models[i].magnetization(true);
    for(unsigned long j=0;j<correlation.size();j++){
      correlation[j][i]=models[i].correlation(j,true);
    }
  }
}

bool modelValues::converged(FILE *outfile){

  bool retval=true;
  double mean = magnetization.sum()/static_cast<double>(magnetization.size());
  fprintf(outfile,"Mean %f\t",mean);
  for(unsigned long i=0;i<correlation.size();i++){
    mean = correlation[i].sum()/static_cast<double>(correlation[i].size());
    fprintf(outfile,"%f\t",mean);
  }
  fprintf(outfile,"\n");


  double CI = (magnetization.max()-magnetization.min());
  if(CI > _limit)
    retval=false;
  fprintf(outfile,"CI %f\t",CI);
  for(unsigned long i=0;i<correlation.size();i++){
    CI = (correlation[i].max()-correlation[i].min());
    fprintf(outfile,"%f\t",CI);
    if(CI > _limit)
      retval=false;
  }
  fprintf(outfile,"\n");

  return retval;
}

void modelValues::printResults(double temp, FILE *outfile){

  double mean = magnetization.sum()/static_cast<double>(magnetization.size());
  double CI = (magnetization.max()-magnetization.min());

  fprintf(outfile,"T \tM eM \tCorrelation eC\n");
  fprintf(outfile,"%f\t%f %f\t",temp,mean,CI);
  for(unsigned long i=0;i<correlation.size();i++){
    mean = correlation[i].sum()/static_cast<double>(correlation[i].size());
    CI = (correlation[i].max()-correlation[i].min());
    fprintf(outfile,"%f %f\t",mean,CI);
  }
  fprintf(outfile,"\n");
}

void modelValues::print(FILE *outfile){
  for(unsigned long i=0;i<magnetization.size();i++){
    fprintf(outfile,"%lu %f\t",i,magnetization[i]);
    for(unsigned long j=0;j<correlation.size();j++)
      fprintf(outfile,"%f\t",correlation[j][i]);
    fprintf(outfile,"\n");
  }
}

void modelValues::printHeader(FILE *outfile){
  fprintf(outfile,"Magnetization\tCorrelation (1-%lu)\n",correlation.size()); 
}

