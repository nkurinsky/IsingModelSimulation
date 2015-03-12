#include "modelstats.hpp"

//constructor implementation
// -- ncorr = distance out to which to calculate correlation
// -- nmodels = number of models beinbg simulated
modelValues::modelValues(int nmodels, int ncorr){
  magnetization.resize(nmodels);
  correlation.resize(ncorr,valarray<double>(nmodels));

  //default values for convergence/zeroing
  _limit=0.01;
  _zerolimit=0.0001;
}

void modelValues::populate(vector<lattice> &models){
  for(unsigned long i=0;i<models.size();i++){ //for each model
    //calculate mean magnetization, 
    //set to nonzero value if returned value greater than limit
    magnetization[i]=(models[i].magnetization(true) > _zerolimit) ? models[i].magnetization(true) : 0;
    for(unsigned long j=0;j<correlation.size();j++){
      //calculate mean correlation,
      //set to nonzero value if returned value greater than limit
      correlation[j][i]=(models[i].correlation(j,true) > _zerolimit) ? models[i].correlation(j,true) : 0;
    }
  }
}

void modelValues::setLimit(double limit){
  _limit=limit;
  _zerolimit=0.01*limit;
}

bool modelValues::converged(FILE *outfile){

  bool retval=true;
  //compute average magnetization across all models
  double mean = magnetization.sum()/static_cast<double>(magnetization.size());
  fprintf(outfile,"Mean %f\t",mean);
  for(unsigned long i=0;i<correlation.size();i++){
    //compute average correlation across all models
    mean = correlation[i].sum()/static_cast<double>(correlation[i].size());
    fprintf(outfile,"%f\t",mean);
  }
  fprintf(outfile,"\n");

  //compute total range of magnetization values across models
  double CI = (magnetization.max()-magnetization.min());
  if(CI > _limit) //not converged if above limit
    retval=false;
  fprintf(outfile,"CI %f\t",CI);
  for(unsigned long i=0;i<correlation.size();i++){
    //compute total range of correlation values across models
    CI = (correlation[i].max()-correlation[i].min());
    fprintf(outfile,"%f\t",CI);
    if(CI > _limit) //not converged if above limit
      retval=false;
  }
  fprintf(outfile,"\n");

  //will be true is all parameters vary by less than convergence limit across models
  return retval;
}

//calculate cross-chain values without computing convergence, printing temperature
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

