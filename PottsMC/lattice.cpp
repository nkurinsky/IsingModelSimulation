#include "lattice.hpp"

//general constructor implementation
probe::probe(){
  totalPoints=0;
  correlationMean=0;
}

//general constructor implementattion
location::location(int ndim, int size, int index){
  static bool called=false;

  //initialiaze random number generator
  if(not called){
    srand(time(NULL));
    called=true;
  }
  
  _ndim=ndim;
  _size=size;
  _index_max=pow(size,ndim);
  
  _indices.resize(ndim,0);
  _dim_step.resize(ndim,0);
  
  //calculate base for each dimension
  for(int i=0;i<ndim;i++){
    _dim_step[i]=pow(size,i);
  }
  
  for(int i=0;i<ndim;i++)
    _indices[i]=0;
  
  //if specific index not passed, randomize
  if(index < 0)
    randomize();
  else 
    set_index(index);
}

//copy constructor implementation
location::location(const location &obj){
  _ndim=obj._ndim;
  _size=obj._size;
  _index_max=obj._index_max;
  _indices=obj._indices; //deep copy
  _dim_step=obj._dim_step;
}

void location::set_index(int index){
  if (index < _index_max){ //error checking
    int index_temp=index;
    //loop to convert index to coordinate vector
    for(int i=0;i<_ndim;i++){
      _indices[i]=index_temp % _size;
      index_temp=index_temp/_size;
    }
  }
  else{
    printf("Invalid index passed to location::set_index\n");
    throw(1);
  }
}

void location::move(int dim, int distance){
  set_index(index(dim,distance));
}

void location::randomize(){
  for(int i=0;i<_ndim;i++){
    _indices[i]=rand() % _size; //random number in range 0,size-1
  }
}

//get index within two dimensional lattice array
int location::index(int dim, int distance){
  //making this a static data structure ensures that only one exists across all classes
  static map< vector<int>, int> _1d_indices;
  static bool initialized=false;
  
  if(not initialized){
    vector<int> _tindices;
    _tindices.resize(_ndim,0);
    //initialize lookup table
    int j;
    bool done=false;
    //loop through all possible index vectors
    while(not done){
      _1d_indices[_tindices]=0;
      //calculate 1d index for index vector
      for(j=0;j<_ndim;j++){
	_1d_indices[_tindices]+=_dim_step[j]*_tindices[j];
      }

      //find next vector to instantiate
      _tindices[0]++;
      for(j=0;j<_ndim;j++){
	if(_tindices[j] == _size){
	  _tindices[j]=0;
	  if(j+1<_ndim)
	    _tindices[j+1]++;
	  else
	    done=true;
	}
      }
    }
    initialized=true;
  }

  //if no displacment just return index for current location
  if(distance == 0)
    return _1d_indices[_indices];
  else{
    //create temporary vector at current location
    vector<int> irel(_indices);
    //move temporary vector
    irel[dim]+=distance;
    //apply periodic boundary conditions
    while(irel[dim]>_size-1)
      irel[dim]-=_size;
    while(irel[dim]<0)
      irel[dim]+=_size;

    return _1d_indices[irel];
  }
}

void location::indices(vector<int> array){
  //resize and fill array
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
  printf("\n\t1D Index: %i\n",index());
}

void location::neighbor(location &nb, int dim, int distance){
  nb.set_index(index(dim,distance));
}

//general constructor implementation
lattice::lattice(int ndim, int size, char q){
  _ndim=ndim;
  _size=size;
  _q=q;
  //resize to total number of lattice sites
  _lattice.resize(pow(size,ndim));

  _T=0;
  _total_flips=0;
  _display_calls=0;

  _magLastCheckSize=0;
  _lastMagMean=0;
  _spinAngle=(2*M_PI)/static_cast<double>(_q);
  magValues.clear();

  //set lattice values to random spins
  randomize();
}

void lattice::randomize(){
  for(unsigned long i=0;i<_lattice.size();i++){
    _lattice[i]=randomSpin();
  }
}

//set temperature and bond generation probability (simple boltzmann factor with k=1)
void lattice::setTemp(double T){
  _T=T;
  _pBond=1-exp(-1.0/_T);
}

//generate integer in (1,q)
spin lattice::randomSpin(){
  return static_cast<spin>((rand() % _q) + 1);
}

//decide whether to add bond based on precomputed probability compared to random variate
bool lattice::addBond(){
  static double norm=1.0/static_cast<double>(RAND_MAX);
  return (norm*static_cast<double>(rand()) < _pBond);
}

spin lattice::get(location site) const{
  return get(site.index());
}

spin lattice::get(int index) const{
  //error checking
  if((index < static_cast<int>(_lattice.size())) and (index >= 0)){
    return _lattice[index];
  }
  else{
    printf("ERROR: Invalid index %i passed to lattice::flip\n",index);
    exit(1);
  }
}

void lattice::flip(location site, spin newValue){
  flip(site.index(),newValue);
}

void lattice::flip(int index, spin newValue){
  //lattice vector range error checking
  if((index < static_cast<int>(_lattice.size())) and (index >= 0)){
    if((newValue <= _q) and (newValue > 0)){
      //perform flip and increment counter
      _lattice[index]=newValue;
      _total_flips++;
    }
    else
      printf("ERROR: Invalid q value %i passed to lattice::flip\n",newValue);
  }
  else
    printf("ERROR: Invalid index %i passed to lattice::flip\n",index);
}

double lattice::magnetization(bool mean){
  //recalculate if update has been performed since last calculation
  if(_total_flips > _magLastCheckSize){
    //calculate complex magnitude of sum of spins
    double re=0;
    double im=0;
    for(unsigned long i=0;i<_lattice.size();i++){
      re+=cos_LatticeAngle(_lattice[i]);
      im+=sin_LatticeAngle(_lattice[i]);
    }
    //add complex magnitude per spin to magnetization array
    magValues.push_back(sqrt(pow(re,2)+pow(im,2))/static_cast<double>(_lattice.size()));
    _magLastCheckSize=_total_flips;
    //calculate new vector mean for magnetization
    _lastMagMean=chainMean(magValues);
  }

  //return mean or last value depending on input
  if(mean)
    return _lastMagMean;
  else
    return magValues.back(); 
}

double lattice::correlation(int length, bool mean){
  static double latticeSize=static_cast<double>(pow(_size,_ndim));
  //number of correlation probes to generate, standard is 1/4 of lattice size
  static double totalPoints=((latticeSize/4) > SITE_LIMIT) ? SITE_LIMIT : latticeSize/4;
  
  //check if probe for this length has been initialized
  if(probes.count(length) == 0){
    probes[length].zeros.resize(totalPoints);
    probes[length].points.resize(totalPoints);
    probe*tprobe=&probes[length];
    location site0(_ndim,_size);
    //for totalPoints probes
    for(int i=0;i<static_cast<int>(totalPoints);i++){
      //generate random zero point
      site0.randomize();
      tprobe->zeros[i]=site0.index();
      tprobe->points[i].resize(_ndim);
      //generate one point per dimension for each 0 at requested length
      for(int j=0;j<_ndim;j++){
	tprobe->points[i][j]=site0.index(j,length);
      }
    }
    tprobe->totalPoints=totalPoints*static_cast<double>(_ndim);
    tprobe->correlationMean=0;
  }

  //get probe structure for requested length
  probe * lprobe=&probes[length];
  static double q=static_cast<double>(_q);
  static double norm=q/(q-1)/(lprobe->totalPoints);
  //values to add for correlated or uncorrelated spins
  static double uncorrvalue=-1/q;
  static double corrvalue=1+uncorrvalue;

  //loop over lattice sites, checking correlation at each point
  double val=0;
  for(unsigned long i=0;i<lprobe->zeros.size();i++){
    for(unsigned long j=0;j<lprobe->points[0].size();j++){
      //lambda expression to determine which value to add, based on whether spins are equal
      val += (_lattice[lprobe->zeros[i]]==_lattice[lprobe->points[i][j]]) ? corrvalue : uncorrvalue;
    }
  }
  //normalize and store correlation value
  lprobe->correlation.push_back(norm*val);

  //return mean or last value depending on input
  if(mean){
    lprobe->correlationMean=chainMean(lprobe->correlation);
    return lprobe->correlationMean-pow(magnetization(CALCMEAN),2);
  }
  else
    return lprobe->correlation.back()-pow(magnetization(),2);
}

bool lattice::correlated(location site1, location site2){
  return (_lattice[site1.index()] == _lattice[site2.index()]);
}

//Wolf Clustering Algorithm
int lattice::flipCluster(int Nflips){
  static location site(_ndim,_size);
  stack<int> clusterSites;
  int tind;
  int flips=0;

  //continute to flip clusters until desired number of spins flipped
  while(flips < Nflips){
    //start at random site
    site.randomize();
    //find spin at base site
    int newspin=randomSpin();
    //generate new spin, not the same as old spin
    int oldspin=get(site);
    while(newspin == oldspin){
      newspin=randomSpin();
    }
    //flip site, increment counter
    flip(site,newspin);
    flips++;
    //push initial site to stack
    clusterSites.push(site.index());
    //loop while stack is not empty (more spins to visit)
    while(not clusterSites.empty()){
      //set site to top of stack and delete from stack
      site.set_index(clusterSites.top());
      clusterSites.pop();
      //loop over all neighboring sites
      for(int dim=0;dim<_ndim;dim++){
	for(int dir=-1;dir<2;dir+=2){
	  //set location to neighbor
	  tind=site.index(dim,dir);
	  //if neighbor has same spin and new bond is requested
	  if((get(tind) == oldspin) and (addBond())){
	    //flip spin and add to stack
	    flip(tind,newspin);
	    flips++;
	    clusterSites.push(tind);
	  }
	}
      }
    } 
  }
  
  _total_flips+=flips;
  return flips;
}

//output lattice, in blocks of 2d lattices
void lattice::display(FILE *outfile) const{
  //size of base 2d lattice (2d slice of ndim lattice)
  int twoDsize=_size*_size;
  string separator(_size*2,'-');

  if(outfile != NULL){
    _display_calls++;
    fprintf(outfile,"- Display Call %lu Start -\n",_display_calls);
    //loop over all spins
    for(unsigned long i=0;i<_lattice.size();i++){
      //print separators between slices
      if(i % twoDsize == 0){
	if(i > 0)
	  fprintf(outfile,"\n");
	fprintf(outfile,"%s\n ",separator.c_str());
      }
      //print newlines at the end of each row
      else if(i % _size == 0)
	fprintf(outfile,"\n ");
      //print lattice value
      fprintf(outfile,"%i ",_lattice[i]);
    }
    fprintf(outfile,"\n%s\n",separator.c_str());
    fprintf(outfile,"-  Display Call %lu End  -\n",_display_calls);
  }
  else{
    printf("ERROR: invalid file pointer passed to lattice::displayLattice\n");
    throw(1);
  }
}

//calculate average of latter half of values in an array
//the latter half is used to throw out "thermalization" period
double lattice::chainMean(vector<double> chain){
  double retval=0;
  //if chain only one link, return link
  if(chain.size() == 1){
    return chain[0];
  }
  //add latter half values
  for(unsigned long i=chain.size()/2;i<chain.size();i++){
    retval+=chain[i];
  }
  //divide by number of values
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

  //initialization loop for value lookup table
  if(cosVals.size() == 0){
    cosVals.resize(_q+1,1);
    for(int i=1;i<=_q;i++){
      cosVals[i] = cos(angle*static_cast<double>(i));
    }
  }

  //get lookup table value
  return cosVals[value];
}

double lattice::sin_LatticeAngle(spin value){
  static vector<double> sinVals;
  static double angle = (2*M_PI)/_q;

  while(value > _q)
    value-=_q;
  while(value < 1)
    value+=_q;

  //initialization loop for value lookup table
  if(sinVals.size() == 0){
    sinVals.resize(_q+1,0);
    for(int i=1;i<=_q;i++){
      sinVals[i] = sin(angle*static_cast<double>(i));
    }
  }

  //get lookup table value
  return sinVals[value];
}
