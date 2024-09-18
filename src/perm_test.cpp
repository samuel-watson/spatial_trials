#include "perm.h"

// test an optimiser
// [[Rcpp::export]]
SEXP new_r_stat(const VectorXd& y, 
                const MatrixXd& B_, 
                const ArrayXd& dists_,
                const MatrixXd& all_dists_,
                const std::vector<double>& lower,
                const std::vector<double>& upper,
                const double l,
                const double kappa,
                const double nu){
  XPtr<Rstat> ptr(new Rstat(y,B_,dists_,all_dists_));
  ptr->update_bounds(lower,upper);
  ptr->update_shapes(l,kappa,nu);
  return wrap(ptr);
}

// [[Rcpp::export]]
SEXP optim_r(SEXP xp, const std::vector<double>& del){
  XPtr<Rstat> ptr(xp);
  std::vector<double> start(del);
  optim<double(const std::vector<double>&),BOBYQA> op(start);
  op.fn<&Rstat::rstat, Rstat >(ptr);
  op.set_bounds(ptr->lower_bounds,ptr->upper_bounds);
  op.minimise();
  return wrap(ptr->values());
}

// [[Rcpp::export]]
SEXP fun(const ArrayXd& dist,
   const double l,
   const double kappa,
   const double nu,
   const std::vector<double>& del,
   const bool mis = false){
  
  ArrayXd result = fn(dist,l,kappa,nu,del,mis);
  return wrap(result);
}

// [[Rcpp::export]]
SEXP get_r_stat(SEXP xp,
                const std::vector<double>& del){
  XPtr<Rstat> ptr(xp);
  return wrap(ptr->rstat(del));
}

// [[Rcpp::export]]
void update_y(SEXP xp,
              const VectorXd& y,
              const ArrayXd& dists_){
  XPtr<Rstat> ptr(xp);
  ptr->update_y(y);
  ptr->update_dists(dists_,false);
}

// [[Rcpp::export]]
void update_b(SEXP xp,
              const MatrixXd& B){
  XPtr<Rstat> ptr(xp);
  ptr->update_B(B);
  ptr->update_ysv();
  
}

// [[Rcpp::export]]
void use_misspec(SEXP xp){
  XPtr<Rstat> ptr(xp);
  ptr->set_misspec(true);
}

// [[Rcpp::export]]
SEXP permute_distances(const MatrixXd& D,
                       const int n,
                       const double min_dist){
  ArrayXd dists = new_dist_permute(D,n,min_dist);
  return wrap(dists);
}

// [[Rcpp::export]]
double permute_p_value_b(SEXP xp,
                         const int n,
                         const double min_dist,
                         const int n_iter,
                         const std::vector<double>& del){
  XPtr<Rstat> ptr(xp);
  ArrayXd newdists(ptr->dists);
  std::vector<double> result;
  double test_stat;
  
  std::vector<double> start(del);
  optim<double(const std::vector<double>&),BOBYQA> op(start);
  op.fn<&Rstat::rstat, Rstat >(ptr);
  op.set_bounds(ptr->lower_bounds,ptr->upper_bounds);
  op.minimise();
  result = ptr->values();
  test_stat = abs(result[1]);
  // Rcpp::Rcout << "\nTest stat: " << test_stat;
  
  int exceeds = 0;
  ptr->use_alt(true);
  for(int i = 0; i < n_iter; i++){
    ptr->permute_dists(n,min_dist);
    op.minimise();
    result = ptr->values();
    // Rcpp::Rcout << "\nIter: " << i << " new stat " << abs(result.second);
    if(abs(result[1]) >= test_stat) exceeds++;
  }
  ptr->use_alt(false);
  return (1.0/(n_iter+1.0))*(exceeds + 1);
}

// [[Rcpp::export]]
SEXP calc_fn(SEXP xp,const std::vector<double>& d0){
  XPtr<Rstat> ptr(xp);
  ArrayXd newdists(ptr->dists.size());
  newdists = ptr->calc_fn(d0);
  return wrap(newdists);
}

// [[Rcpp::export]]
double confint_b(SEXP xp,
                 const int n,
                 const double min_dist,
                 const int n_iter,
                 const double bound_start,
                 const double b0,
                 const std::vector<double>& d0){
  
  XPtr<Rstat> ptr(xp);
  ArrayXd dists(ptr->dists);
  ArrayXd newdists(dists);
  std::vector<double> result;
  
  std::vector<double> start(d0);
  optim<double(const std::vector<double>&),BOBYQA> op(start);
  op.fn<&Rstat::rstat, Rstat >(ptr);
  op.set_bounds(ptr->lower_bounds,ptr->upper_bounds);
  
  double bound = bound_start;
  double step, alt_stat, test_stat;
  const double k_tmp = 0.8352199;
  double k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));  
  
  ptr->update_stored_fn(d0);
    
  for(int i = 0; i < n_iter; i++){
    ptr->update_ysv(bound);
    ptr->use_alt(false);
    op.minimise();
    result = ptr->values();
    test_stat = abs(result[1]);
    
    ptr->use_alt(true);
    ptr->permute_dists(n,min_dist);
    op.minimise();
    result = ptr->values();
    alt_stat = abs(result[1]);
    
    step = k*(b0 - bound);
    //Rcpp::Rcout << "\nIter: " << i << " | test: " << test_stat << " alt: " << alt_stat;
    if(test_stat > alt_stat){
      bound += step*0.05/(i+1);
    } else {
      bound -= step*0.95/(i+1);
    }
    Rcpp::Rcout << "\rIter: " << i << " | bound: " << bound;
  }
  ptr->update_ysv(0.0);
  ptr->use_alt(false);
  return bound;
}

// [[Rcpp::export]]
SEXP confint_del(SEXP xp,
                 const int n,
                 const double min_dist,
                 const int n_iter,
                 const std::vector<double>& bound_start,
                 const double b0,
                 const std::vector<double>& d0){
  
  XPtr<Rstat> ptr(xp);
  ArrayXd newdists(ptr->dists);
  std::vector<double> result;
  
  std::vector<double> start(d0);
  optim<double(const std::vector<double>&),BOBYQA> op(start);
  op.fn<&Rstat::rstat, Rstat >(ptr);
  op.set_bounds(ptr->lower_bounds,ptr->upper_bounds);
  
  std::vector<double> bound = bound_start;
  std::vector<double> step(bound);
  double alt_stat, test_stat;
  const double k_tmp = 0.8352199;
  double k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));  
  
  for(int i = 0; i < n_iter; i++){
    ptr->update_stored_fn(bound);
    ptr->update_ysv(b0);
    ptr->use_alt(false,true);
    op.minimise();
    result = ptr->values();
    test_stat = abs(result[1]);
    
    ptr->use_alt(true,true);
    ptr->permute_dists(n,min_dist);
    op.minimise();
    result = ptr->values();
    alt_stat = abs(result[1]);
    
    for(int j = 0; j < step.size(); j++)step[j] = k*(d0[j] - bound[j]);
    
    // Rcpp::Rcout << "\nIter: " << i << " | test: " << test_stat << " alt: " << alt_stat;
    Rcpp::Rcout << "\rIter: " << i << " | bound:";
    for(int j = 0; j < step.size(); j++){
      if(test_stat > alt_stat){
        bound[j] += step[j]*0.05/(i+1);
      } else {
        bound[j] -= step[j]*0.95/(i+1);
      }
      if(j==0){
        bound[j] = bound[j] < ptr->lower_bounds[j] ? ptr->lower_bounds[j] : bound[j];
      } else {
        bound[j] = bound[j] > ptr->upper_bounds[j] ? ptr->upper_bounds[j] : bound[j];
      }
      Rcpp::Rcout << " " << bound[j];
    }
    
  }
  ptr->update_ysv(0.0);
  ptr->use_alt(false);
  return wrap(bound);
}