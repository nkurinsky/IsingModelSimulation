#include "lattice.hpp"

probe::probe(){
  totalPoints=0;
  correlationMean=0;
}

location::location(int ndim, int size){
  static bool called=false;
  if(not called){
    srand(time(NULL));
    called=true;
  }
  
  _ndim=ndim;
  _size=size;
  
  _indices.resize(ndim,0);
  
  for(int i=0;i<ndim;i++)
    _indices[i]=0;
  
  randomize();

}

location::location(const location &obj){
  _ndim=obj._ndim;
  _size=obj._size;
  _indices=obj._indices;
}

void location::randomize(){
  for(int i=0;i<_ndim;i++){
    _indices[i]=rand() % _size;
  }
}

void location::move(location &obj){
  _ndim=obj._ndim;
  _size=obj._size;
  _indices=obj._indices;
}

void location::move(coords array){
  if(array.size() == _indices.size()){
    _indices=array;
  }
  else{
    fprintf(stderr,"Invalid array passed to location::move\n");
    exit(10);
  }
}

void location::move(int dim, int distance){
  if(dim < 0 or dim >= _ndim){
    printf("Invalid dimension passed to location::move");
    exit(10);
  }

  if (distance != 0){
    _indices[dim]+=distance;
    while(_indices[dim]>_size-1)
      _indices[dim]-=_size;
    while(_indices[dim]<0)
      _indices[dim]+=_size;
  }    
}

coords location::get(int dim, int distance){
  static coords tc(_ndim,0);

  if(dim < 0 or dim >= _ndim){
    printf("Invalid dimension passed to location::get");
    exit(10);
  }

  if(distance == 0){
    return _indices;
  }    

  tc=_indices;
  tc[dim]+=distance;
  while(tc[dim]>_size-1)
    tc[dim]-=_size;
  while(tc[dim]<0)
    tc[dim]+=_size;

  return tc;
}

void location::indices(vector<int> array){
  array.resize(_ndim);
  for(int i=0;i<_ndim;i++)
    array[i]=_indices[i];
}

void location::print(){
  printf("Location Details:\n");
  printf("\tDimensions: %i\n",_ndim);
  printf("\tLattice Size: %i\n",_size);
  printf("\tND Indices:");
  for(int i=0;i<_ndim;i++)
    printf(" %i",_indices[i]);
}

lattice::lattice(int ndim, int size, char q){
  _ndim=ndim;
  _size=size;
  _q=q;

  _T=0;
  _total_flips=0;
  _display_calls=0;

  _magLastCheckSize=0;
  _lastMagMean=0;
  _spinAngle=(2*M_PI)/static_cast<double>(_q);
  magValues.clear();

  randomize();
}

void lattice::randomize(){
  //needs to size lattice as well
  coords ct(_ndim,0);
  
  while(ct[0] != _size){
    _lattice[ct]=randomSpin();
    ct.back()++;
    for(int i=_ndim-1;i>0;i--){
      if(ct[i] == _size){
	ct[i]=0;
	ct[i-1]++;
      }
    }
  }
}

void lattice::setTemp(double T){
  _T=T;
  _pBond=1-exp(-1.0/_T);
}

spin lattice::randomSpin(){
  return static_cast<spin>((rand() % _q) + 1);
}

bool lattice::addBond(){
  static double norm=1.0/static_cast<double>(RAND_MAX);
  return (norm*static_cast<double>(rand()) < _pBond);
}

spin lattice::get(location site){
  return _lattice[site.get()];
}

void lattice::flip(location site, spin newValue){
  if(site.ndim() == _ndim){
    if((newValue <= _q) and (newValue > 0)){
      _lattice[site.get()]=newValue;
      _total_flips++;
    }
    else{
      fprintf(stderr,"ERROR: Invalid q value %i passed to lattice::flip\n",newValue);
      exit(10);
    }
  }
  else
    fprintf(stderr,"ERROR: Invalid location passed to lattice::flip\n");
}

double lattice::magnetization(bool mean){
  if(_total_flips > _magLastCheckSize){
    double re=0;
    double im=0;
    for(auto it = _lattice.begin();it != _lattice.end();it++){
      re+=cos_LatticeAngle(it->second);
      im+=sin_LatticeAngle(it->second);
    }

    magValues.push_back(sqrt(pow(re,2)+pow(im,2))/static_cast<double>(_lattice.size()));
    _magLastCheckSize=_total_flips;
    _lastMagMean=chainMean(magValues);
  }

  if(mean)
    return _lastMagMean;
  else
    return magValues.back(); 
}

double lattice::correlation(int length, bool mean){
  static double latticeSize=static_cast<double>(pow(_size,_ndim));
  static double totalPoints=((latticeSize/4) > SITE_LIMIT) ? SITE_LIMIT : latticeSize/4;
  
  if(probes.count(length) == 0){
    probes[length].zeros.resize(totalPoints);
    probes[length].points.resize(totalPoints);
    probe*tprobe=&probes[length];
    location site0(_ndim,_size);
    for(int i=0;i<static_cast<int>(totalPoints);i++){
      site0.randomize();
      tprobe->zeros[i]=site0.get();
      tprobe->points[i].resize(_ndim);
      for(int j=0;j<_ndim;j++){
	tprobe->points[i][j]=site0.get(j,length);
      }
    }
    tprobe->totalPoints=totalPoints*static_cast<double>(_ndim);
    tprobe->correlationMean=0;
  }

  probe * lprobe=&probes[length];
  static double q=static_cast<double>(_q);
  static double norm=q/(q-1)/(lprobe->totalPoints);
  static double uncorrvalue=-1/q;
  static double corrvalue=1+uncorrvalue;

  double val=0;
  for(unsigned long i=0;i<lprobe->zeros.size();i++){
    for(unsigned long j=0;j<lprobe->points[0].size();j++){
      val += (_lattice[lprobe->zeros[i]]==_lattice[lprobe->points[i][j]]) ? corrvalue : uncorrvalue;
    }
  }
  lprobe->correlation.push_back(norm*val);

  if(mean){
    lprobe->correlationMean=chainMean(lprobe->correlation);
    return lprobe->correlationMean-pow(magnetization(CALCMEAN),2);
  }
  else
    return lprobe->correlation.back()-pow(magnetization(),2);
}

bool lattice::correlated(location site1, location site2){
  return (_lattice[site1.get()] == _lattice[site2.get()]);
}

//Wolf Clustering Algorithm
int lattice::flipCluster(int Nflips){
  static location site(_ndim,_size);
  static location tc(_ndim,_size);
  stack<coords> clusterSites;  
  int flips=0;

  while(flips < Nflips){
    site.randomize();
    int newspin=randomSpin();
    int oldspin=get(site);
    while(newspin == oldspin){
      newspin=randomSpin();
    }
    flip(site,newspin);
    flips++;
    clusterSites.push(site.get());
    while(not clusterSites.empty()){
      site.move(clusterSites.top());
      clusterSites.pop();
      for(int dim=0;dim<_ndim;dim++){
	for(int dir=-1;dir<2;dir+=2){
	  tc.move(site.get(dim,dir));
	  if((get(tc) == oldspin) and (addBond())){
	    flip(tc,newspin);
	    flips++;
	    clusterSites.push(tc.get());
	  }
	}
      }
    } 
  }
  
  _total_flips+=flips;
  return flips;
}

void lattice::display(FILE *outfile){
  int twoDsize=_size*_size;
  string separator(_size*2,'-');

  if(outfile != NULL){
    _display_calls++;
    fprintf(outfile,"- Display Call %lu Start -\n",_display_calls);

    coords ct(_ndim,0);
    int i=0;

    while(ct[0] != _size){
      
      if(i % twoDsize == 0){
	if(i > 0)
	  fprintf(outfile,"\n");
	fprintf(outfile,"%s\n ",separator.c_str());
      }
      else if(i % _size == 0)
	fprintf(outfile,"\n ");
      fprintf(outfile,"%i ",_lattice[ct]);
      
      ct.back()++;
      for(int j=_ndim-1;j>0;j--){
	if(ct[j] == _size){
	  ct[j]=0;
	  ct[j-1]++;
	}
      }
      i++;
    }

    fprintf(outfile,"\n%s\n",separator.c_str());
    fprintf(outfile,"-  Display Call %lu End  -\n",_display_calls);
  }
  else{
    printf("ERROR: invalid file pointer passed to lattice::displayLattice\n");
    throw(1);
  }
}

double lattice::chainMean(vector<double> chain){
  double retval=0;
  if(chain.size() == 1){
    return chain[0];
  }
  for(unsigned long i=chain.size()/2;i<chain.size();i++){
    retval+=chain[i];
  }
  retval/=static_cast<double>(chain.size()-chain.size()/2);
  return retval;
}

double lattice::cos_LatticeAngle(spin value){
  static vector<double> cosVals;
  static double angle = (2*M_PI)/_q;

  while(value > _q)
    value-=_q;
  while(value < 1)
    value+=_q;

  //initialization loop
  if(cosVals.size() == 0){
    cosVals.resize(_q+1,1);
    for(int i=1;i<=_q;i++){
      cosVals[i] = cos(angle*static_cast<double>(i));
    }
  }

  return cosVals[value];
}

double lattice::sin_LatticeAngle(spin value){
  static vector<double> sinVals;
  static double angle = (2*M_PI)/_q;

  while(value > _q)
    value-=_q;
  while(value < 1)
    value+=_q;

  //initialization loop
  if(sinVals.size() == 0){
    sinVals.resize(_q+1,0);
    for(int i=1;i<=_q;i++){
      sinVals[i] = sin(angle*static_cast<double>(i));
    }
  }

  return sinVals[value];
}
