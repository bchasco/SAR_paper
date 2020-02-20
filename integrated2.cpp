// Separable covariance on lattice with AR1 structure in each direction.
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(yShift); //the offset for the environmental and pitTotal start years
  DATA_IVECTOR(yr);
  DATA_IVECTOR(j);
  DATA_VECTOR(s_n); //total number of years across env. and pit data
  DATA_VECTOR(s_k); //total number of years across env. and pit data
  DATA_IVECTOR(k);
  DATA_INTEGER(nvar); //total number of env. vars.
  DATA_ARRAY(env);
  DATA_ARRAY(timing);
  DATA_MATRIX(xm);
  DATA_IVECTOR(marVars);
  DATA_SCALAR(sd);
  DATA_INTEGER(re_j); //day effect flag
  DATA_INTEGER(re_t); //year effect flag
  DATA_INTEGER(re_jt); //day X year effect flag
  DATA_INTEGER(cov_pars); //env. covariance flag
  DATA_INTEGER(retro);
  DATA_VECTOR(env_mu);
  DATA_VECTOR(env_sc);
  DATA_VECTOR(env_mu_2000_2015);
  DATA_VECTOR(env_sc_2000_2015);
  DATA_INTEGER(calibration_flag);
  DATA_MATRIX(calibration_timing);
  DATA_IVECTOR(calibration_years);
  
  PARAMETER_VECTOR(mu_s);
  PARAMETER_VECTOR(frho_j); //temporal correlation for each env. variable
  PARAMETER_VECTOR(frho_t); //temporal correlation for each env. variable
  PARAMETER_VECTOR(frho1_jt); //temporal correlation for each env. variable
  PARAMETER_VECTOR(frho2_jt); //temporal correlation for each env. variable
  PARAMETER_VECTOR(fpsi_j); //temporal correlation for each env. variable
  PARAMETER_VECTOR(fpsi_t); //temporal correlation for each env. variable
  PARAMETER_VECTOR(fpsi_jt); //temporal correlation for each env. variable
  
  PARAMETER_MATRIX(eps_j);    //random effects for each variable, nt * ni
  PARAMETER_MATRIX(eps_t);    //random effects for each variable, nt * ni
  PARAMETER_ARRAY(eps_jt);    //random effects for each variable, nt * ni
  
  PARAMETER_MATRIX(beta_mar);


  PARAMETER(frho_x); //temporal correlation for env. variables
  PARAMETER_VECTOR(frho_Rx); //unstructured correlations
  PARAMETER_VECTOR(fpsi_x); //process error
  PARAMETER_ARRAY(eps_x);    //random effects for each variable, nt * ni
  
  // 
  using namespace density;
  
  Type ff = 0.;  //joint likelihood

  int nk = mu_s.size();
  int ni = s_k.size();
  int nt = eps_t.rows();
  int nte = env.dim[0];
  int nj = eps_j.rows();
  int rb = beta_mar.rows();
  
  //Correlation coefficient and process variances
  vector<Type> psi_x = exp(fpsi_x);

  Type rho_x = 1/(1+exp(-frho_x));
  UNSTRUCTURED_CORR_t<Type> Sigma_x(frho_Rx);
  if(cov_pars==1){
	  ff += AR1(rho_x,Sigma_x)(eps_x);
//     for(int t=1;t<nte;t++){
// 		  ff += Sigma_x(eps_x.col(t)-rho_x*eps_x.col(t));
// 	  }	
  }

  matrix<Type> env_hat(nte,nvar);
  if(cov_pars==1){
    for(int j=0;j<nvar;j++){
      for(int t=0;t<nte;t++){
        env_hat(t,j) = eps_x(j,t)*psi_x(j);
        if(env(t,j)>-100 & cov_pars==1){
          ff -= dnorm(env_hat(t,j), env(t,j), sd, true);
        }
      }
    }
  }
  vector<Type> rho_j = atan(frho_j)*2./3.154;//1/(1+exp(-frho_j));//
  vector<Type> rho_t = atan(frho_t)*2./3.154;//1/(1+exp(-frho_t));//
  vector<Type> rho1_jt = atan(frho1_jt)*2./3.154;//1/(1+exp(-frho1_jt));//
  vector<Type> rho2_jt = atan(frho2_jt)*2./3.154;//1/(1+exp(-frho2_jt));//
  vector<Type> psi_j = exp(fpsi_j);
  vector<Type> psi_t = exp(fpsi_t);
  vector<Type> psi_jt = exp(fpsi_jt);

  matrix<Type> eMar(ni,nk);
  
  // Survival random effects
  for(int kk=0;kk<nk;kk++){
    array<Type> tmp_eps(nj,nt);
    for(int tt=0;tt<nt;tt++){
      for(int jj=0;jj<nj;jj++){
        tmp_eps(jj,tt) = eps_jt(kk,jj,tt);
      }
    }
    if(re_jt==1){
      ff += AR1(rho1_jt(kk),AR1(rho2_jt(kk)))(tmp_eps);		
    }
  }

  vector<Type> s_hat(ni);
  
  //Get the subset of the marine variable that you want to model
  matrix<Type> xm_tmp(ni,rb);
  for(int rr=0;rr<rb;rr++){
    xm_tmp.col(rr) = xm.col(marVars(rr));
  }
  //Model the temporal processes and get the marine effects
  for(int kk=0;kk<nk;kk++){
    if(re_j==1){
      ff += AR1(rho_j(kk))(eps_j.col(kk));		
    }
    if(re_t==1){
      ff += AR1(rho_t(kk))(eps_t.col(kk));
    }
    
    //Predicted survival
    eMar.col(kk) = xm_tmp*beta_mar.col(kk);
  }
  
  for(int i=0;i<s_k.size();i++){
    Type nu = mu_s(k(i)) +        //Mean
      eMar(i) +         //Annual environmental effects
      eps_t(yr(i) - yShift,k(i)) * psi_t(k(i)) +        //1D AR1 year
      eps_j(j(i),k(i)) * psi_j(k(i)) +       //1D AR1 day
      eps_jt(k(i),j(i),yr(i)-yShift) * psi_jt(k(i)); //2D AR1XAR1 yearXday

    s_hat(i) = exp(nu)/(1+exp(nu)); //Logit
    ff -= dbinom(s_k(i),s_n(i),s_hat(i),true);
  }


  vector<Type> s_t(nt); //Retrospective survival
  matrix<Type> s_tj(nt,nj); //Retrospective survival
  matrix<Type> xm_tmp2(nt,rb); //Retrospective environmental variables
  for(int rr=0;rr<rb;rr++){
    for(int t=0;t<nt;t++){
      //Un-zscore
      xm_tmp2(t,rr) = env(t+yShift,marVars(rr)) * env_sc(marVars(rr)) + env_mu(marVars(rr));
      //Zscore for the 2000 to 2015 mean and sd
      xm_tmp2(t,rr) = (xm_tmp2(t,rr) - env_mu_2000_2015(marVars(rr)))/env_sc_2000_2015(marVars(rr));
    }
  }
  
  vector<Type> eMar2(nt); //Marine effects for retrospective
  eMar2 = xm_tmp2*beta_mar.col(0);
  s_t = 0.;
  for(int tt=0;tt<nt;tt++){
    for(int jj=0;jj<nj;jj++){
      Type nu = mu_s(0) +        //Mean
        eMar2(tt) +         //Annual environmental effects
        eps_t(tt,0) * psi_t(0) +        //1D AR1 year
        eps_j(jj,0) * psi_j(0) +       //1D AR1 day
        eps_jt(0,jj,tt) * psi_jt(0); //2D AR1XAR1 yearXday
      
      s_t(tt) += nu * timing(jj,tt);
      s_tj(tt,jj) = nu;
    }
  }

  vector<Type> s_j = eps_j.col(0) * psi_j(0);
  s_j += mu_s(0);

  REPORT(mu_s);
  REPORT(s_t);
  REPORT(frho_j); //temporal correlation for each env. variable
  REPORT(frho_t); //temporal correlation for each env. variable
  REPORT(frho1_jt); //temporal correlation for each env. variable
  REPORT(frho2_jt); //temporal correlation for each env. variable
  REPORT(fpsi_j); //temporal correlation for each env. variable
  REPORT(fpsi_t); //temporal correlation for each env. variable
  REPORT(fpsi_jt); //temporal correlation for each env. variable

  REPORT(eps_j);    //random effects for each variable, nt * ni
  REPORT(eps_t);    //random effects for each variable, nt * ni
  REPORT(eps_jt);    //random effects for each variable, nt * ni
  REPORT(beta_mar);
  REPORT(s_hat);
  REPORT(nj);
  REPORT(nt);
  REPORT(nte);
  REPORT(nk);
  REPORT(xm);
  REPORT(eMar2);
  REPORT(xm_tmp2);
  REPORT(calibration_timing);
  REPORT(psi_x);
  REPORT(rho_x);
  REPORT(Sigma_x.cov());
  REPORT(ff);
  REPORT(eps_x);
  REPORT(env);
  REPORT(env_hat);
  ADREPORT(s_t);
  ADREPORT(s_tj);
  ADREPORT(s_j);
  
  return ff;
}
