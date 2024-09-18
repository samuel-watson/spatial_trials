#include <RcppEigen.h>
#include <random>
#include <vector>
#include <stdexcept>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(glmmrBase)]]
#define R_BUILD
#include "optim.h"

using namespace Eigen;
using namespace Rcpp;

template<typename T1, typename T2> SEXP wrap( const std::pair<T1,T2>& _v ) {
  return Rcpp::List::create(
    Rcpp::Named("first")  = Rcpp::wrap<T1>( _v.first ),
    Rcpp::Named("second") = Rcpp::wrap<T2>( _v.second )
  );
};

// ys test statistic
inline VectorXd ys(const VectorXd& y,
                   const MatrixXd& B){
  VectorXd By = B*y;
  By *= 1.0/By.norm();
  return By;
}

// LogSumExp function
inline ArrayXd fn(const ArrayXd& dist,
                  const double l,
                  const double kappa,
                  const double nu,
                  const std::vector<double>& del,
                  bool misspec = false){
  ArrayXd result(dist);
  if(!misspec){
    if(del.size()==1){
      result *= -l/del[0];
      result = result.exp() + exp(-1.0*l);
      result = result.log();
      result *= (-1.0/l);
      result = 1.0 - result.pow(kappa);
    } else {
      ArrayXd signd(dist.size());
      for(int i = 0; i < signd.size(); i++)signd(i) = dist(i) >= 0 ? 1.0 : -1.0;
      result -= del[1];
      result *= -l/(del[0] - del[1]);
      result *= signd;
      result = result.exp() + (-0.5*l*(signd+1.0)).exp();
      result = result.log();
      result *= (-1.0/l);
      result *= signd;
      result = 1.0 - result.pow(kappa);
    }
    return result.pow(nu);
  } else {
    if(del.size()==1){
      result *= M_PI/(2.0*del[0]);
      result = result.cos();
      for(int i = 0; i < result.size(); i++)if(dist(i) > del[0])result(i) = 0.0;
    } else {
      ArrayXd signd(dist.size());
      for(int i = 0; i < signd.size(); i++)signd(i) = dist(i) >= 0 ? 1.0 : -1.0;
      result -= del[1];
      result *= -l/(del[0] - del[1]);
      result *= signd;
      result = result.exp() + (-0.5*l*(signd+1.0)).exp();
      result = result.log();
      result *= (-1.0/l);
      result *= signd;
      result = 1.0 - result.pow(kappa);
    }
    return result;
  }
  
  
}

inline double xsr(const ArrayXd& dists,
                  const VectorXd& ystat,
                  const MatrixXd& B,
                  const double l,
                  const double kappa,
                  const double nu,
                  const std::vector<double>& del){
  ArrayXd fnv = fn(dists,l,kappa,nu,del);
  VectorXd xs = ys(fnv.matrix(),B);
  double R = (xs.transpose() * ystat)(0);
  return R;
}

inline double xsr(const ArrayXd& dists,
                  const VectorXd& ystat,
                  const MatrixXd& B,
                  const double l,
                  const double kappa,
                  const double nu,
                  const std::vector<double>& del,
                  const ArrayXd& pred_fn){
  ArrayXd fnv = fn(dists,l,kappa,nu,del);
  fnv -= pred_fn;
  VectorXd xs = ys(fnv.matrix(),B);
  double R = (xs.transpose() * ystat)(0);
  return R;
}

inline int random_int(const int lower, const int upper){
  std::random_device rd;     // Only used once to initialise (seed) engine
  std::mt19937 rng(rd());    // Random-number engine used (Mersenne-Twister in this case)
  std::uniform_int_distribution<> uni(lower,upper); // Guaranteed unbiased
  auto random_integer = uni(rng);
  return random_integer;
}

inline ArrayXd new_dist_permute(const MatrixXd& D,
                                const int n,
                                const double min_dist){
  int dim = (int)D.cols();
  std::vector<int> idx;
  int idx_new = random_int(0,dim-1);
  idx.push_back(idx_new);
  int iter = 0;
  double min_dist_curr;
  
  if(D.rows()==D.cols()){
    while(idx.size() < n){
      idx_new = random_int(0,dim-1);
      min_dist_curr = 10.0;
      for(int i = 0; i < idx.size(); i++){
        if(D(idx[i],idx_new) < min_dist_curr) min_dist_curr = D(idx[i],idx_new);
      }
      if(min_dist_curr > min_dist) idx.push_back(idx_new);
      iter++;
      if(iter > 100)break;
    }
  } else {
    for(int i = 0; i < n; i++){
      while(idx.size() < n){
        idx_new = random_int(0,dim-1);
        if(std::find(idx.begin(), idx.end(), idx_new) == idx.end()) {
          idx.push_back(idx_new);
        } 
      }
    }
  }
  
  ArrayXd dists(D.rows());
  if(iter > 100){
    dists.setZero();
    throw std::runtime_error("Iterations exceed limit");
  } else {
    for(int i = 0; i < dists.size(); i++){
      min_dist_curr = 10.0;
      for(int j = 0; j < idx.size(); j++){
        if(D(i,idx[j]) < min_dist_curr) min_dist_curr = D(i,idx[j]);
      }
      dists(i) = min_dist_curr;
    }
  }
  return dists;
}

class Rstat {
public:
  VectorXd y;
  VectorXd ysv;
  ArrayXd dists;
  ArrayXd dists_alt;
  MatrixXd B; // orthogonal matrix
  MatrixXd D; // distance matrix
  double l = 50.0;
  double kappa = 4.0;
  double nu = 8.0;
  std::vector<double> lower_bounds; 
  std::vector<double> upper_bounds;
  
  Rstat(const VectorXd& y_, const MatrixXd& B_, const ArrayXd& dists_, const MatrixXd& D_) : y(y_), ysv(y.size()), 
    dists(dists_), dists_alt(dists_), B(B_), D(D_), stored_fn(dists_.size()) {
    ysv = ys(y,B_);
    lower_bounds.push_back(0.01);
    upper_bounds.push_back(0.44);
    if((dists < 0).any()){
      lower_bounds.push_back(-0.2);
      upper_bounds.push_back(0.0);
    }
  };
  
  void update_y(const VectorXd& y_){
    y = y_;
    ysv = ys(y,B);
  }
  
  void update_ysv(){
    ysv = ys(y,B);
  }
  
  ArrayXd calc_fn(const std::vector<double>& del){
    ArrayXd newfn = fn(dists,l,kappa,nu,del,use_misspec);
    return newfn;
  }
  
  void update_ysv(const double beta){
    ysv = ys(y - beta*stored_fn.matrix(),B);
  }
  
  void update_stored_fn(const std::vector<double>& d0){
    stored_fn = calc_fn(d0);
  }
  
  double rstat(const std::vector<double>& x){
    if(x.size()==2 & (dists < 0).any() != 1 ) throw std::runtime_error("No negative distances");
    last_del = x[0];
    if(x.size()>1)last_del_i = x[1];
    double R;
    if(use_alt_dists){
      if(use_sub_pred){
        R = xsr(dists_alt,ysv,B,l,kappa,nu,x,stored_fn);
      } else {
        R = xsr(dists_alt,ysv,B,l,kappa,nu,x);
      }
    } else {
      if(use_sub_pred){
        R = xsr(dists,ysv,B,l,kappa,nu,x,stored_fn);
      } else {
        R = xsr(dists,ysv,B,l,kappa,nu,x);
      }
    }
    last_r = abs(R);
    return -1.0*abs(R);
  }
  
  void update_bounds(const std::vector<double>& lower, 
                     const std::vector<double>& upper){
    lower_bounds = lower;
    upper_bounds = upper;
  }
  
  void update_shapes(const double l_, const double kappa_, const double nu_){
    l = l_;
    kappa = kappa_;
    nu = nu_;
  }
  
  void update_B(const MatrixXd& B_){
    B = B_;
  }
  
  void update_dists(const ArrayXd& newdists, bool alt = true){
    if(alt){
      dists_alt = newdists; 
    } else {
      dists = newdists;
    }
  }
  
  void permute_dists(const double n, const double min_dist){
    dists_alt = new_dist_permute(D,n,min_dist);
  }
  
  std::vector<double> values(){
    std::vector<double> result;
    result.push_back(last_del);
    result.push_back(last_r);
    if((dists < 0).any())result.push_back(last_del_i);
    return result;
  }
  
  void use_alt(bool alt, bool sub_pred = false){
    use_alt_dists = alt;
    use_sub_pred = sub_pred;
  }
  
  void set_misspec(bool mis){
    use_misspec = mis;
  }
  
private:
  double last_del = 0.1;
  double last_del_i = 0.1;
  double last_r = 0.0;
  bool use_alt_dists = false;
  bool use_sub_pred = false;
  ArrayXd stored_fn;
  bool use_misspec = false;
};

