#include <vector>
#include <random>
#include <math.h>
#include <iostream>

class RandSeed{
public:
  RandSeed(int seed);
  void reinit(int _N,
	      int _N_fields,
	      double _val,
	      double _r,
	      double _R);
  double eval(const double x,
	      const double y,
	      const int f) const;

  int N, N_fields;
  double val, r, R;

  //std::random_device rd;
  std::mt19937 gen;

  std::vector<std::vector<double> > centers; // Assume 2D for now
  std::vector<int> eta_i;
  std::vector<double> value;
  
};

RandSeed::RandSeed(int seed)
  :
  gen(seed){}

void RandSeed::reinit(int _N,
		      int _N_fields,
		      double _val,
		      double _r,
		      double _R){
  
  N = _N;//Number of seeds
  N_fields = _N_fields; //Number of order parameters
  val = _val;
  r = _r;
  R = _R;
  
  centers.resize(N,std::vector<double>(2));
  eta_i.resize(N);
  value.resize(N);

  std::uniform_real_distribution<double> dist(-R,R);
  std::uniform_int_distribution<int> dist_int(0,2*N_fields-1);

  double x,y;
  for (int i=0; i<N; ++i){
    while(true){
      x = dist(gen);
      y = dist(gen);
      if (sqrt(x*x + y*y) <= R){
	break;
      }
    }
    centers[i] = {x,y};
    int variant = dist_int(gen);
    eta_i[i] = variant/2;
    double sign = 2.*((variant%2) - 0.5);
    value[i] = sign*val;
  }

}

double RandSeed::eval(const double x,
		      const double y,
		      const int f) const{

  double val = 0.;
  for (int i=0; i<N; ++i){
    if ((pow(x-centers[i][0],2) + pow(y-centers[i][1],2) < r*r)){
      if ((eta_i[i] == f) && (pow(x-centers[i][0],2) + pow(y-centers[i][1],2) < 0.9*0.9*r*r)){
	val = value[i];
      }
      break;
    }
  }

  return val;
}

