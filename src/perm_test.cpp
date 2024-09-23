#include "perm.h"

// test an optimiser
// [[Rcpp::export]]
SEXP new_r_stat(const VectorXd& y, 
                const MatrixXd& B_, 
                const ArrayXd& dists_,
                const MatrixXd& all_dists_,
                const std::vector<double>& lower_e,
                const std::vector<double>& upper_e,
                const std::vector<double>& lower_i,
                const std::vector<double>& upper_i,
                const double l,
                const double kappa,
                const double nu,
                const ArrayXi& time){
  XPtr<Rstat> ptr(new Rstat(y,B_,dists_,all_dists_,time));
  ptr->update_bounds(lower_e,upper_e,lower_i,upper_i);
  ptr->update_shapes(l,kappa,nu);
  return wrap(ptr);
}

// [[Rcpp::export]]
SEXP optim_r(SEXP xp, const std::vector<double>& del){
  XPtr<Rstat> ptr(xp);
  std::vector<double> start(del);
  std::vector<double> lower(ptr->lower_bounds_e);
  std::vector<double> upper(ptr->upper_bounds_e);
  if((ptr->dists < 0).any()){
    for(int i = 0; i < ptr->lower_bounds_i.size(); i++){
      lower.push_back(ptr->lower_bounds_i[i]);
      upper.push_back(ptr->upper_bounds_i[i]);
    }
  }
  optim<double(const std::vector<double>&),BOBYQA> op(start);
  // op.control.trace = 1;
  op.fn<&Rstat::rstat, Rstat >(ptr);
  op.set_bounds(lower,upper);
  op.minimise();
  return wrap(ptr->values());
}

// [[Rcpp::export]]
SEXP fun(const ArrayXd& dist,
         const double l,
         const double kappa,
         const double nu,
         const std::vector<double>& del,
         int nT,
         const ArrayXi& time,
         const bool mis = false){
  ArrayXd result(dist.size());
  if(nT == 1){
    result = fn(dist,l,kappa,nu,del,mis);
  } else {
    result = fn(dist,l,kappa,nu,del,nT,time,mis);
  }
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
  std::vector<double> lower(ptr->lower_bounds_e);
  std::vector<double> upper(ptr->upper_bounds_e);

  if((ptr->dists < 0).any()){
    for(int i = 0; i < ptr->lower_bounds_i.size(); i++){
      lower.push_back(ptr->lower_bounds_i[i]);
      upper.push_back(ptr->upper_bounds_i[i]);
    }
  }
  if(del.size() != lower.size())throw std::runtime_error("Del size != lower size");
  op.set_bounds(lower,upper);
  op.minimise();
  result = ptr->values();
  test_stat = abs(result[ptr->nT]);
  
  int exceeds = 0;
  ptr->use_alt(true);
  Rcpp::Rcout << "\nPermutations: \n";
  for(int i = 0; i < n_iter; i++){
    ptr->permute_dists(n,min_dist);
    op.minimise();
    result = ptr->values();
    if(abs(result[ptr->nT]) >= test_stat) exceeds++;
    Rcpp::Rcout << "\rIter: " << i+1 << " of " << n_iter;
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
SEXP confint_b(SEXP xp,
                 const int n,
                 const double min_dist,
                 const int n_iter,
                 const std::vector<double>& bound_start,
                 const std::vector<double>& b0,
                 const std::vector<double>& d0,
                 bool joint = false){
  
  XPtr<Rstat> ptr(xp);
  ArrayXd dists(ptr->dists);
  ArrayXd newdists(dists);
  std::vector<double> start(d0);
  std::vector<double> bound = bound_start;
  std::vector<double> step(bound);
  std::vector<double> b0bound(b0);
  double alt_stat, test_stat;
  const double k_tmp = 0.8352199;
  double k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));  
  
  ptr->update_stored_fn(d0);
  
  for(int i = 0; i < n_iter; i++){
    if(!joint){
      Rcpp::Rcout << "\rIter: " << i << " | bound:";
      for(int t = 0; t < bound.size(); t++){
        b0bound = b0;
        b0bound[t] = bound[t];
        ptr->update_ysv(b0bound);
        ptr->use_alt(false);
        double result0 = ptr->rstat(d0);
        test_stat = abs(result0);
        
        ptr->use_alt(true);
        ptr->permute_dists(n,min_dist);
        result0 = ptr->rstat(d0);
        alt_stat = abs(result0);
        step[t] = k*(b0[t] - bound[t]);
        if(test_stat > alt_stat){
          bound[t] += step[t]*0.05/(i+1);
        } else {
          bound[t] -= step[t]*0.95/(i+1);
        }
        Rcpp::Rcout << " " << bound[t];
      }
    } else {
      ptr->update_ysv(bound);
      ptr->use_alt(false);
      double result0 = ptr->rstat(d0);
      test_stat = abs(result0);
      
      ptr->use_alt(true);
      ptr->permute_dists(n,min_dist);
      result0 = ptr->rstat(d0);
      alt_stat = abs(result0);
      
      for(int j = 0; j < step.size(); j++)step[j] = k*(b0[j] - bound[j]);
      Rcpp::Rcout << "\rIter: " << i << " | bound:";
      for(int j = 0; j < step.size(); j++){
        if(test_stat > alt_stat){
          bound[j] += step[j]*0.05/(i+1);
        } else {
          bound[j] -= step[j]*0.95/(i+1);
        }
        Rcpp::Rcout << " " << bound[j];
      }
    }
  }
  
  ptr->update_ysv(0.0);
  ptr->use_alt(false);
  return wrap(bound);
}

// [[Rcpp::export]]
SEXP confint_del(SEXP xp,
                 const int n,
                 const double min_dist,
                 const int n_iter,
                 const std::vector<double>& bound_start,
                 const std::vector<double>& b0,
                 const std::vector<double>& d0,
                 bool joint = false){
  
  XPtr<Rstat> ptr(xp);
  ArrayXd newdists(ptr->dists);
  std::vector<double> start(d0);
  std::vector<double> bound = bound_start;
  std::vector<double> step(bound);
  std::vector<double> d0bound(bound);
  double alt_stat, test_stat;
  const double k_tmp = 0.8352199;
  const double k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));  
  ptr->update_b0(b0);
  
  for(int i = 0; i < n_iter; i++){
    if(!joint){
      Rcpp::Rcout << "\rIter: " << i << " | bound:";
      for(int t = 0; t < bound.size(); t++){
        d0bound = d0;
        d0bound[t] = bound[t];
        ptr->update_stored_fn(d0bound);
        ptr->update_ysv(b0);
        ptr->use_alt(false);
        double result0 = ptr->rstat(d0);
        test_stat = abs(result0);

        ptr->use_alt(true);
        ptr->permute_dists(n,min_dist);
        result0 = ptr->rstat(d0);
        alt_stat = abs(result0);

        step[t] = k*(d0[t] - bound[t]);
        if(test_stat > alt_stat){
          bound[t] += step[t]*0.05/(i+1);
        } else {
          bound[t] -= step[t]*0.95/(i+1);
        }

        if(t < ptr->nT){
          bound[t] = bound[t] < ptr->lower_bounds_e[t] ? ptr->lower_bounds_e[t] : bound[t];
        } else {
          bound[t] = bound[t] > ptr->upper_bounds_i[t] ? ptr->upper_bounds_i[t] : bound[t];
        }

        Rcpp::Rcout << " " << bound[t];
      }
    } else {
      ptr->update_stored_fn(bound);
      ptr->update_ysv(b0);
      ptr->use_alt(false);
      double result0 = ptr->rstat(d0);
      test_stat = abs(result0);
      
      ptr->use_alt(true);
      ptr->permute_dists(n,min_dist);
      result0 = ptr->rstat(d0);
      alt_stat = abs(result0);
      
      for(int j = 0; j < step.size(); j++)step[j] = k*(d0[j] - bound[j]);
      Rcpp::Rcout << "\rIter: " << i << " | bound:";
      for(int j = 0; j < step.size(); j++){
        if(test_stat > alt_stat){
          bound[j] += step[j]*0.05/(i+1);
        } else {
          bound[j] -= step[j]*0.95/(i+1);
        }
        if(j < ptr->nT){
          bound[j] = bound[j] < ptr->lower_bounds_e[j] ? ptr->lower_bounds_e[j] : bound[j];
        } else {
          bound[j] = bound[j] > ptr->upper_bounds_i[j] ? ptr->upper_bounds_i[j] : bound[j];
        }
        Rcpp::Rcout << " " << bound[j];
      }
    }
  }
  
  ptr->update_ysv(0.0);
  ptr->use_alt(false);
  return wrap(bound);
}

// [[Rcpp::export]]
SEXP confint_both(SEXP xp,
                 const int n,
                 const double min_dist,
                 const int n_iter,
                 const std::vector<double>& bound_start_b,
                 const std::vector<double>& bound_start_del,
                 const std::vector<double>& d0,
                 const std::vector<double>& b0){
  
  XPtr<Rstat> ptr(xp);
  ArrayXd newdists(ptr->dists);
  
  std::vector<double> start(d0);
  
  std::vector<double> bound_b = bound_start_b;
  std::vector<double> bound_del = bound_start_del;
  std::vector<double> step_b(bound_b);
  std::vector<double> step_del(bound_del);
  double alt_stat, test_stat;
  const double k_tmp = 0.8352199;
  const double k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));  
  ptr->update_b0(b0);
  
  for(int i = 0; i < n_iter; i++){
    ptr->update_stored_fn(bound_del);
    ptr->update_ysv(bound_b);
    // ptr->y_sub_pred();
    ptr->use_alt(false);
    double result0 = ptr->rstat(d0);
    test_stat = abs(result0);
    
    ptr->use_alt(true);
    ptr->permute_dists(n,min_dist);
    result0 = ptr->rstat(d0);
    alt_stat = abs(result0);
    
    for(int j = 0; j < step_b.size(); j++)step_b[j] = k*(b0[j] - bound_b[j]);
    Rcpp::Rcout << "\rIter: " << i << " | bound_b:";
    for(int j = 0; j < step_b.size(); j++){
      if(test_stat > alt_stat){
        bound_b[j] += step_b[j]*0.05/(i+1);
      } else {
        bound_b[j] -= step_b[j]*0.95/(i+1);
      }
      Rcpp::Rcout << " " << bound_b[j];
    }
    
    for(int j = 0; j < step_del.size(); j++)step_del[j] = k*(d0[j] - bound_del[j]);
    for(int j = 0; j < step_del.size(); j++){
      if(test_stat > alt_stat){
        bound_del[j] += step_del[j]*0.05/(i+1);
      } else {
        bound_del[j] -= step_del[j]*0.95/(i+1);
      }
      if(j < ptr->nT){
        bound_del[j] = bound_del[j] < ptr->lower_bounds_e[j] ? ptr->lower_bounds_e[j] : bound_del[j];
      } else {
        bound_del[j] = bound_del[j] > ptr->upper_bounds_i[j] ? ptr->upper_bounds_i[j] : bound_del[j];
      }
      Rcpp::Rcout << " " << bound_del[j];
    }
    
    // if(i % 2 == 0){
    //   for(int j = 0; j < step_b.size(); j++)step_b[j] = k*(b0[j] - bound_b[j]);
    //   Rcpp::Rcout << "\rIter: " << i << " | bound_b:";
    //   for(int j = 0; j < step_b.size(); j++){
    //     if(test_stat > alt_stat){
    //       bound_b[j] += step_b[j]*0.05/(i+1);
    //     } else {
    //       bound_b[j] -= step_b[j]*0.95/(i+1);
    //     }
    //     Rcpp::Rcout << " " << bound_b[j];
    //   }
    // } else {
    //   for(int j = 0; j < step_del.size(); j++)step_del[j] = k*(d0[j] - bound_del[j]);
    //   for(int j = 0; j < step_del.size(); j++){
    //     if(test_stat > alt_stat){
    //       bound_del[j] += step_del[j]*0.05/(i);
    //     } else {
    //       bound_del[j] -= step_del[j]*0.95/(i);
    //     }
    //     if(j < ptr->nT){
    //       bound_del[j] = bound_del[j] < ptr->lower_bounds_e[j] ? ptr->lower_bounds_e[j] : bound_del[j];
    //     } else {
    //       bound_del[j] = bound_del[j] > ptr->upper_bounds_i[j] ? ptr->upper_bounds_i[j] : bound_del[j];
    //     }
    //     Rcpp::Rcout << " " << bound_del[j];
    //   }
    // }
    
    
    
  }
  
  ptr->update_ysv(0.0);
  ptr->use_alt(false);
  
  bound_b.insert(bound_b.end(),bound_del.begin(),bound_del.end());
  
  return wrap(bound_b);
}