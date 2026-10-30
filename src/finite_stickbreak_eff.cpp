
// This script is mainly to be used for aiding in calculating Cross validated LPPL for
// a validation set. Also we plan to replace sampler for C by some efficient
// way. Apart from that the change is in saving sum of loglikelihoods of holdout entries
//Also we need to provide input of the holdout entries. But if the data already has missuing
// entries, what we do is that we feed -Inf for those entries and our created l_pdf_x will return 0;
// thus nullifying its effect in the sum,

// We have to create an R wrapper which if user specifies val=T, will choose
// a validation set on its own and get LPPL for that set, without user 
// worrying too much about creating those things. The wrapper will also put -Inf
// for those entries which are already missing say

// Note that due to change in one argument, we have to change most of 
// our earlier code if we use this script to get LPPL as a measure of fit 
// in earlier simulations. What it looks like, later our oruginal cpp file
// have to be replaced by this one since we will probably want to keep 
// the option of getting LPPL in artifcally induced held out dataset

#include <RcppArmadillo.h>
#include <RcppTN.h>
#include <Rmath.h>  
#include <cmath>

// [[Rcpp::depends(RcppTN)]]

// [[Rcpp::depends(RcppArmadillo)]]


// The following C++ implementation of rcat has been borrowed from hesim package's github:- https://github.com/hesim-dev/hesim/blob/main/inst/include/hesim/stats/distributions.h#L48
// ***************************************************************************//**
// * Random number generation for the categorical distribution.
// * @param probs A vector of length @c K specifying the probability of each of the 
// * @c K categories. Internally normalized to sum to 1.
// * @return A random sample from the categorical distribution.
// ******************************************************************************/ 


inline int rcat(arma::rowvec probs) {
  int k = probs.n_elem;
  double probs_sum = accu(probs);
  probs = probs/probs_sum;
  Rcpp::IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  int max = which_max(ans);
  return(max);
}


// Calculating log of sum of exponential of a vector using max trick
// ***************************************************************************//**
// log(sum(exp(x)))=m+log(sum(exp(x-m))) where m=max(x)
// ******************************************************************************/ 


// [[Rcpp::export]]

double log_sum_exp(arma::vec x){
  double m= arma::max(x);
  arma::vec x_new=x-m*arma::ones<arma::vec>(x.n_elem);
  
  return m+std::log(arma::accu(arma::exp(x_new)));
}

// ***************************************************************************//**
// Calculating log likelihood of an observed entry (i.e. whose true value is known)
// based on the mixture model we have.
// If we have truely unobserved data, i.e. whose value we don't know at all,
// then we ask the function to return 0. We will make those entries Inf/-inf
// when fitting the model. Hence the if condition below.
// ******************************************************************************/ 


double l_pdf_x_DP(double x, arma::vec ws, arma::vec mus, arma::vec sds){
  
  if(std::isinf(x))
    return 0.0;
  
  arma::vec u= arma::log(ws) + arma::log_normpdf(x,mus, sds);
  
  return log_sum_exp(u);
  
}



// [[Rcpp::export]]

Rcpp::List mcmcDraws_fsbp(arma::mat X, int K, int R, arma::vec mus_init, arma::vec sigma2s_init, arma::umat C_ij_init, arma::vec all_mis_elem, int M=5000, int burn_in=2000, int thin=5, double tau2= 10000 ){
  
  // X is the matrix with missing entries, after introducing user induced artifical missing elements
  //note. X can already have its own truely missing entries. 
  
  //all_mis_elem= the true values of the missing entries (ideal simulation setup), includes -Inf/Inf if something is already missing
  //-------it is in the same order as how missing entries's posotion was determined using sorted column major ordering----
  // Basically X.elem(mis)=all_mis_elem
  
  
  //tau2 is the variance of prior for mu_k's i.e. muk~N(0,tau2)
  // last three arguments are initial values supp;ied for mu_k's sigma2_k's and C_ij's
  // We can get these initial estimates based on k-means say. We assume our X is scaled,
  // so that we initialize missing entry as zeroes (Refer to BAMITA paper as well) to get full X
  // Or we can initialize missing entries to mean of all other observarions while doing knn
  
  // Then we get C_ij's, mu_s etc for all (i,j)'s. 
  //Be CAREFUL that C++ indexing starts from Zero. So, We have to provide initial C_ij as
  // per C++ ordering, that is, if C_ij=1 as obtained from k-means then it should be C_ij=0 when we supply here
  
  // This ensures that we have to have X as scaled when provided as input
  // and we have to initialize missing observations at zero; or initialize at mean of all other non missing observation
  
  //Later all these things can be wrapped up in an R function where first few initial
  // steps will be the things we described above and then it calls the actual cpp code, 
  // thus giving user much less thing to worry about
  
  arma::mat X_dummy=X; //X_dummy is just the completed version of X, with truely missing entrues as Inf/-Inf and artifically hold out entries set to their true values
  // We need X_dummy to calculate the cross validated log likelihood
  
  arma::uword N=X.n_rows;
  arma::uword Q=X.n_cols;
  //save the indices of missing entries of X in column major ordering:-
  arma::uvec mis= arma::find_nan(X); // Assuming that the missing entries are coded as NA in our input in R
  
  X_dummy.elem(mis)=all_mis_elem;
  
  arma::uvec non_mis= arma::find_nonnan(X); //non missing elements
  
  arma::cube C3(N,Q*K,1); 
  
  //initialize [[Z1,Z2,Z3]] first, then create Z based on tensor decomposition
  
  
  arma::field<arma::mat> Z_comp(3); //save Z1 (N*R),Z2(Q*R),Z3(K*R) in a field
  
  // Below the method of construction of Z as a tensor CP factorization of Z1,Z2,Z3 is shown
  // However, we are initializing Z1,Z2,Z3 by zeros. So, all the things done below is not needed.
  // since we can just initialize Z as a cube with fill::zeros (Z=O if Z1,Z2,Z3 are O )
  
  //You can also use fill::randn where the construction will come handy.
  
  Z_comp(0)=arma::mat(N,R, arma::fill::zeros);                //arma::fill::randn
  Z_comp(1)=arma::mat(Q,R, arma::fill::zeros);
  Z_comp(2)=arma::mat(K,R, arma::fill::zeros);
  arma::vec a3(N*Q*K,arma::fill::zeros);
  for(int r=0;r<R;r++){
    a3=a3+arma::kron(arma::kron((Z_comp(2)).col(r),(Z_comp(1)).col(r)),(Z_comp(0)).col(r));
  }
  C3.slice(0)=arma::reshape(a3,N,Q*K);
  arma::cube Z=arma::reshape(C3,N,Q,K);
  
  // Then given Z and initialized C_ij, draw Z_star 
  
  arma::cube Z_star(N,Q,K); 
  
  arma::umat C_ij = C_ij_init ;
  
  
  arma::mat Sigma_e= arma::eye(K,K) ;
  
  arma::uword k0;
  for(arma::uword i=0;i<N;++i){
    for(arma::uword j=0;j<Q;++j){
      k0=C_ij(i,j);
      if( k0==0 ){
        Z_star(i,j,0)=RcppTN::rtn1(Z(i,j,0),1.0,0.0,std::numeric_limits<double>::infinity());
        for(arma::uword k=1;k<K;++k){
          Z_star(i,j,k)=arma::randn(arma::distr_param(Z(i,j,k),1.0));
        }
        
      }
      
      if(  (k0 >0) &&  (k0 < K-1)        ){
        Z_star(i,j,k0)=RcppTN::rtn1(Z(i,j,k0),1.0,0.0,std::numeric_limits<double>::infinity());
        for(arma::uword k=0; k<k0; ++k){
          Z_star(i,j,k)=RcppTN::rtn1(Z(i,j,k),1.0,-std::numeric_limits<double>::infinity(),0.0);
        }
        for(arma::uword k=k0+1; k<K; ++k){
          Z_star(i,j,k)=arma::randn(arma::distr_param(Z(i,j,k),1.0));
        }
        
      }
      
      if(k0==K-1){
        for(arma::uword k=0; k<k0; ++k){
          Z_star(i,j,k)=RcppTN::rtn1(Z(i,j,k),1.0,-std::numeric_limits<double>::infinity(),0.0);
          
        }
        
        Z_star(i,j,k0)=arma::randn(arma::distr_param(Z(i,j,k0),1.0));
      }
      
    }
  }
  
  arma::vec mus=  mus_init ;
  arma::vec sigma2s=  sigma2s_init ;
  arma::vec sds=arma::sqrt(sigma2s);
  
  
  
  
  
  arma::uword ind;
  
  
  // initialize missing entries by mean of all non-missing entries
  
  double mean_nonmis=arma::mean(arma::omit_nan(arma::vectorise(X)));
  
  for (arma::uword i : mis){
    ind=C_ij(i);
    X(i)= mean_nonmis;
  }
  
  
  //start from here
  
  
  
  int numsave =(M-burn_in)/thin-1;
  int save_counter =0;
  double c;
  double p;
  double u;
  int nk;
  double a_post;
  double b_post;
  double post_mu_mean;
  double post_mu_var;
  
  
  
  arma::field<arma::mat>Z_star_mat(3);
  arma::mat Z23(Q*K,R);
  arma::mat Z13(N*K,R);
  arma::mat Z12(N*Q,R);
  
  arma::vec a1;
  arma::vec a2;
  arma::mat cov_inv;
  arma::cube C1(Q,N,K); 
  arma::cube C2(K,N,Q); 
  
  
  
  arma::cube V(N,Q,K);
  arma::cube P(N,Q,K);
  
  arma::vec sub_V;
  arma::vec Pi(K);
  arma::vec logPi(K);
  
  //save the samples
  arma::mat S_mu(K,numsave);
  //arma::mat S_sd(K,numsave);
  arma::mat S_mis(mis.n_elem,numsave);
  
  arma::vec S_Pt(numsave); //samples of summed CV LPPL for heldout entries
  
  double Pt=0.0;
  
  //arma::mat S_P(K, numsave); // save P[i,j,] for some specific (i,j)
  
  //arma::cube S_Zstar(N,K,numsave);
  //arma::cube S_Cij(N,Q,numsave);
  //arma::vec S_C(numsave);
  
  //S_Zstar=array(dim=c(N,K,numsave)) #save only Z_star[,j,] for fixed j
  
  
  arma::vec V_tube(K);
  arma::vec P_tube(K);
  arma::uword i0;
  arma::uword j0;
  
  for (int t=0;t<M;t++){
    
    if(t % 200 ==0){
      Rcpp::Rcout << "Sample No: " << t << std::endl; // print every 50th sample
    }
    
    // We update the parameters in terms of the order we initialized. So we initialized C_ij,mu,sigma and misisng entries
    // Now, we set Z=0 and then generated Z_star giveb Z, C_ij. So, now we start with updating Z first
    // and then going to calculating pi_ijk based on Z and the order follows
    
    //(1) generate from Z=[Z1,Z2,Z3]|Z_star   #entries of Z1,Z2,Z3 are given N(0,1) prior.
    
    //Rcpp::Rcout << "Update Z1,Z2,Z3 " << std::endl;
    
    //i=0 i.e update Z1 given Z2,Z3 & others
    
    Z_star_mat(0)=arma::reshape(Z_star,N,Q*K,1);
    for(int r=0;r<R;r++){
      a1=(Z_comp(1)).col(r);
      a2=(Z_comp(2)).col(r);
      Z23.col(r)=arma::kron(a2,a1);
    }
    
    // Rcpp::Rcout << "Check1:"  << std::endl;
    
    
    cov_inv=arma::inv_sympd( Z23.t() * Z23 + arma::eye(R,R)  );
    
    
    
    Z_comp(0)=Z_star_mat(0) * Z23 * cov_inv;
    
    
    
    
    
    Z_comp(0)=Z_comp(0)+(arma::mvnrnd(arma::zeros(R),cov_inv,N)).t();
    
    
    
    
    //i=1
    for(arma::uword s=0;s<K;++s){
      C1.slice(s)=(Z_star.slice(s)).t();
    }
    Z_star_mat(1) =arma::reshape(C1,Q,(N*K),1);
    for(int r=0;r<R;r++){
      a1=(Z_comp(0)).col(r);
      a2=(Z_comp(2)).col(r);
      Z13.col(r)=arma::kron(a2,a1);
    }
    
    //Rcpp::Rcout << "Check2:"  << std::endl;
    
    
    
    cov_inv=arma::inv_sympd(Z13.t() * Z13 + arma::eye(R,R));
    
    
    Z_comp(1)=Z_star_mat(1) * Z13 * cov_inv;
    
    Z_comp(1)=Z_comp(1)+(arma::mvnrnd(arma::zeros(R),cov_inv,Q)).t();
    
    
    
    
    //i=2
    for(arma::uword j=0;j<Q;++j){
      C2.slice(j)=(Z_star.col_as_mat(j)).t();
    }
    Z_star_mat(2)=arma::reshape(C2,K,(N*Q),1);
    for(int r=0;r<R;r++){
      a1=(Z_comp(0)).col(r);
      a2=(Z_comp(1)).col(r);
      Z12.col(r)=arma::kron(a2,a1);
    }
    
    //Rcpp::Rcout << "Check3:"  << std::endl;
    
    
    cov_inv= arma::inv_sympd( Z12.t()*Z12 + arma::eye(R,R));
    
    
    
    
    
    Z_comp(2)=Z_star_mat(2) * Z12 * cov_inv;
    Z_comp(2)=Z_comp(2)+(arma::mvnrnd(arma::zeros(R),cov_inv,K)).t();
    
    
    
    
    //finally calculate Z
    arma::vec a3(N*Q*K,arma::fill::zeros);
    for(int r=0;r<R;r++){
      a3=a3+arma::kron(arma::kron((Z_comp(2)).col(r),(Z_comp(1)).col(r)),(Z_comp(0)).col(r));
    }
    C3.slice(0)=arma::reshape(a3,N,Q*K);
    Z=arma::reshape(C3,N,Q,K);
    
    
    
    
    //(2) Update from full conditional of C_ij,Z_star,P|Z,mu,sigma, complete data
    
    // (2a) Update pi_ijk=P[i,j,k] given Z using probit stick breaking
    // Note, Z determines V and V in-turn determine P
    
    for(arma::uword k=0; k<K-1; ++k){
      V.slice(k)=arma::normcdf(Z.slice(k));
    }
    
    V.slice(K-1)=arma::ones(N,Q);
    
    for(arma::uword i=0; i<N; ++i){
      for(arma::uword j=0; j<Q; ++j){
        
        V_tube=V.tube(i,j);
        sub_V=arma::ones<arma::vec>(K)-V_tube;
        
        P(i,j,0)=V(i,j,0);
        for(arma::uword k=1; k<K; ++k){
          P(i,j,k)=V(i,j,k)* arma::prod(sub_V.subvec(0,k-1));
        }
        
      }
    }
    
    
    // [(**)Now update C_ij,Z_star|P,Z,mu,sigma, complete data]
    
    
    // [2b] Update C_ij|P,mu,sigma,complete data 
    
    
    //Rcpp::Rcout << "Update C: " << std::endl;
    
    for(arma::uword i=0; i<N; ++i){
      for(arma::uword j=0; j<Q; ++j){
        
        //logPi=P.tube(i,j);
        
        Pi=P.tube(i,j);
        logPi=arma::log(Pi)+ arma::log_normpdf(X(i,j), mus, sds);
        
        //for(arma::uword k=0; k<K; ++k){
        
        //logPi(k)=std::log(P(i,j,k))+arma::log_normpdf(X(i,j), mus(k), sds(k));
        //prolly here i can use vectorized version?
        
        //Pi(k)=P(i,j,k)*arma::normpdf(X(i,j), mus(k), sds(k) );
        
        //}
        
        Pi=arma::exp(logPi-arma::max(logPi)*arma::ones(K));
        
        C_ij(i,j)=rcat(Pi.t());
        
      }
    }
    
    
    // [2c] Update Z_star|C_ij, Z
    
    //Rcpp::Rcout << "Update Z_star: " << std::endl;
    
    for(arma::uword i=0;i<N;++i){
      for(arma::uword j=0;j<Q;++j){
        k0=C_ij(i,j);
        if( k0==0 ){
          Z_star(i,j,0)=RcppTN::rtn1(Z(i,j,0),1.0,0.0,std::numeric_limits<double>::infinity());
          for(arma::uword k=1;k<K;++k){
            Z_star(i,j,k)=arma::randn(arma::distr_param(Z(i,j,k),1.0));
          }
          
        }
        
        if(  (k0 >0) &&  (k0 < K-1)        ){
          Z_star(i,j,k0)=RcppTN::rtn1(Z(i,j,k0),1.0,0.0,std::numeric_limits<double>::infinity());
          for(arma::uword k=0; k<k0; ++k){
            Z_star(i,j,k)=RcppTN::rtn1(Z(i,j,k),1.0,-std::numeric_limits<double>::infinity(),0.0);
          }
          for(arma::uword k=k0+1; k<K; ++k){
            Z_star(i,j,k)=arma::randn(arma::distr_param(Z(i,j,k),1.0));
          }
          
        }
        
        if(k0==K-1){
          for(arma::uword k=0; k<k0; ++k){
            Z_star(i,j,k)=RcppTN::rtn1(Z(i,j,k),1.0,-std::numeric_limits<double>::infinity(),0.0);
            
          }
          
          Z_star(i,j,k0)=arma::randn(arma::distr_param(Z(i,j,k0),1.0));
        }
        
      }
    }
    
    
    
    //(3) ###generate mu_k, sigma2_k from normal-inverse-gamma posterior --muk~N(0,tau2),sigma_2k~InvGamma(2,1)
    // given Data and C_ij
    
    //Rcpp::Rcout << "Update mu_sds: " << std::endl;
    
    
    for(arma::uword k=0;k<K;++k){
      arma::uvec q1=arma::find(C_ij==k);
      nk=q1.n_elem;
      
      if(nk >0 ){
        
        arma::vec x_k = X.elem(q1);
        
        a_post = 2.0 + 0.5 * nk;
        
        b_post =
          1.0 +
          0.5 * arma::accu(arma::square(x_k - mus(k)));
        
        sigma2s(k) =
          1.0 /
            arma::randg<double>(
              arma::distr_param(a_post, 1.0 / b_post)
            );
        post_mu_mean=(tau2*arma::accu(X.elem(q1)))/(nk*tau2+sigma2s(k));
        post_mu_var=sigma2s(k)*tau2/(nk*tau2+sigma2s(k));
        
        
        mus(k)=arma::randn<double>(arma::distr_param(post_mu_mean,sqrt(post_mu_var)));

        
      } else{
        
        
        // if no observations for this atom, draw from the prior
        
        mus(k)=arma::randn<double>(arma::distr_param(0.0,sqrt(tau2)));
        sigma2s(k)=1.0/arma::randg<double>(arma::distr_param(2.0,1.0));
        
        
      }
      
      
      
    } 
    
    sds=arma::sqrt(sigma2s);
    
    
    
    
    
    
    
    //(4) #update missing entries given parameters
    // and also save summed crossvalidated loglikelihoods of artifically heldout
    // missing entries. This will help us to calculate LPPL later
    
    //initialize the sum LPPL to be zero first
    Pt=0.0;
    
    
    //Rcpp::Rcout << "Update missing entry: " << std::endl;
    
    for (arma::uword i : mis){
      //ind=C_ij(i);
      
      i0 = i % N;
      j0 = i / N;
      P_tube=P.tube(i0,j0);
      
      ind=rcat(P_tube.t());
      
      X(i)=arma::randn<double>( arma::distr_param(mus(ind),sds(ind)) );
      
      
      Pt=Pt+l_pdf_x_DP(X_dummy(i0,j0), P_tube, mus, sds    );
      
    }
    
    
    
    
    if(t >burn_in && t % thin ==0){
      S_mu.col(save_counter) =mus;
      //S_sd.col(save_counter) =sds;
      S_mis.col(save_counter)= X.elem(mis);
      
      S_Pt(save_counter)=Pt;
      
      //P_tube=P.tube(5,0);
      //S_P.col(save_counter)=P_tube;
      //S_C(save_counter)=C_ij(0);
      //S_Zstar.slice(save_counter)=Z_star.col(2);
      //S_Cij.slice(save_counter)=C_ij;
      //S_Zstar[,,save_counter]=Z_star[,3,] #saving Z_star[i,3,k] for all i,ks
      save_counter=save_counter+1;
      
    }
    
    
  }
  
  //std::cout << "Time_Cij:--" << total_time_Cij << std::endl;
  //std::cout << "Time_WStar:--" << total_time_WStar << std::endl;
  //std::cout << "Time_W:--" << total_time_W << std::endl;
  
  Rcpp::List L(3);
  
  L[0]=S_mis;
  L[1]=S_mu;
  L[2]=S_Pt;
  //L[2]=S_P;
  //L[3]=S_C;
  //L[1]=S_sd;
  //L[2]=S_mis;
  //L[3]=mis;
  //L[3]=S_Zstar;
  //L[4]=S_Cij;
  
  return L;
  
  
  
}