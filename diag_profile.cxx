#include<iostream>
#include<cmath>
#include<cstring>
#include<algorithm>
#include "stdlib.h"
using namespace std;
#include<vector>
#include<map>
#include "WCPLEEANA/TOsc.h"
#include "WCPLEEANA/Configure_Osc.h"
#include "TROOT.h"
#include "TDecompChol.h"
// diag_profile: at fixed (dm2, sin2_2theta_mue) dump the chi2 profile over the
// LOWER-branch t24 family on a fine grid, and report what a 50-point log scan
// would have found. Answers whether profiling density can explain an 8% limit
// shift without running any toys.
static inline double qf(const double*v,const double*M,int n){
  double s=0; for(int i=0;i<n;i++){ double vi=v[i]; if(vi==0) continue;
    const double*r=M+(Long64_t)i*n; double a=0; for(int j=0;j<n;j++) a+=r[j]*v[j]; s+=vi*a; } return s; }
static inline double qfd(const double*a,const double*b,const double*M,int n,double*buf){
  for(int j=0;j<n;j++) buf[j]=a[j]-b[j]; return qf(buf,M,n); }
static TMatrixD pinv(const TMatrixD&P,const TMatrixD&S){
  int n=P.GetNcols(); TMatrixDSym C(n);
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) C(i,j)=0.5*(S(i,j)+S(j,i));
  for(int i=0;i<n;i++){ double st=P(0,i); if(st==0) st=1e-6; C(i,i)+=st; }
  TDecompChol ch(C); ch.Decompose(); ch.Invert(C); return TMatrixD(C); }
int main(int argc,char**argv){
  gROOT->SetBatch(1);
  double DM2=(argc>1)?atof(argv[1]):2.371;
  double X  =(argc>2)?atof(argv[2]):1.0e-3;
  int NFINE =(argc>3)?atoi(argv[3]):400;
  TOsc *osc_test = new TOsc();

  ///////

  osc_test->tosc_scaleF_POT_BNB  = 1;
  osc_test->tosc_scaleF_POT_NuMI = 1;
     
  ///////

  osc_test->flag_apply_oscillation_BNB  = Configure_Osc::flag_apply_oscillation_BNB;
  osc_test->flag_apply_oscillation_NuMI = Configure_Osc::flag_apply_oscillation_NuMI;

  osc_test->flag_goodness_of_fit_CNP    = Configure_Osc::flag_goodness_of_fit_CNP;

  ///////

  osc_test->flag_syst_dirt   = Configure_Osc::flag_syst_dirt;
  osc_test->flag_syst_mcstat = Configure_Osc::flag_syst_mcstat;
  osc_test->flag_syst_flux   = Configure_Osc::flag_syst_flux;
  osc_test->flag_syst_geant  = Configure_Osc::flag_syst_geant;
  osc_test->flag_syst_Xs     = Configure_Osc::flag_syst_Xs;
  osc_test->flag_syst_det    = Configure_Osc::flag_syst_det;
  
  ///////

  osc_test->flag_NuMI_nueCC_from_intnue         = Configure_Osc::flag_NuMI_nueCC_from_intnue;
  osc_test->flag_NuMI_nueCC_from_overlaynumu    = Configure_Osc::flag_NuMI_nueCC_from_overlaynumu;
  osc_test->flag_NuMI_nueCC_from_appnue         = Configure_Osc::flag_NuMI_nueCC_from_appnue;
  osc_test->flag_NuMI_nueCC_from_appnumu        = Configure_Osc::flag_NuMI_nueCC_from_appnumu;
  osc_test->flag_NuMI_nueCC_from_overlaynueNC   = Configure_Osc::flag_NuMI_nueCC_from_overlaynueNC;
  osc_test->flag_NuMI_nueCC_from_overlaynumuNC  = Configure_Osc::flag_NuMI_nueCC_from_overlaynumuNC;
  
  osc_test->flag_NuMI_numuCC_from_overlaynumu   = Configure_Osc::flag_NuMI_numuCC_from_overlaynumu;
  osc_test->flag_NuMI_numuCC_from_overlaynue    = Configure_Osc::flag_NuMI_numuCC_from_overlaynue;
  osc_test->flag_NuMI_numuCC_from_appnue        = Configure_Osc::flag_NuMI_numuCC_from_appnue;
  osc_test->flag_NuMI_numuCC_from_appnumu       = Configure_Osc::flag_NuMI_numuCC_from_appnumu;
  osc_test->flag_NuMI_numuCC_from_overlaynumuNC = Configure_Osc::flag_NuMI_numuCC_from_overlaynumuNC;
  osc_test->flag_NuMI_numuCC_from_overlaynueNC  = Configure_Osc::flag_NuMI_numuCC_from_overlaynueNC;  
  
  osc_test->flag_NuMI_CCpi0_from_overlaynumu    = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumu;
  osc_test->flag_NuMI_CCpi0_from_appnue         = Configure_Osc::flag_NuMI_CCpi0_from_appnue;
  osc_test->flag_NuMI_CCpi0_from_overlaynumuNC  = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_CCpi0_from_overlaynueNC   = Configure_Osc::flag_NuMI_CCpi0_from_overlaynueNC;
  
  osc_test->flag_NuMI_NCpi0_from_overlaynumu    = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumu;
  osc_test->flag_NuMI_NCpi0_from_appnue         = Configure_Osc::flag_NuMI_NCpi0_from_appnue;
  osc_test->flag_NuMI_NCpi0_from_overlaynumuNC  = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_NCpi0_from_overlaynueNC   = Configure_Osc::flag_NuMI_NCpi0_from_overlaynueNC;


  ///////
  
  osc_test->flag_BNB_nueCC_from_intnue         = Configure_Osc::flag_BNB_nueCC_from_intnue;
  osc_test->flag_BNB_nueCC_from_overlaynumu    = Configure_Osc::flag_BNB_nueCC_from_overlaynumu;
  osc_test->flag_BNB_nueCC_from_appnue         = Configure_Osc::flag_BNB_nueCC_from_appnue;
  osc_test->flag_BNB_nueCC_from_appnumu        = Configure_Osc::flag_BNB_nueCC_from_appnumu;
  osc_test->flag_BNB_nueCC_from_overlaynueNC   = Configure_Osc::flag_BNB_nueCC_from_overlaynueNC;
  osc_test->flag_BNB_nueCC_from_overlaynumuNC  = Configure_Osc::flag_BNB_nueCC_from_overlaynumuNC;
  
  osc_test->flag_BNB_numuCC_from_overlaynumu   = Configure_Osc::flag_BNB_numuCC_from_overlaynumu;
  osc_test->flag_BNB_numuCC_from_overlaynue    = Configure_Osc::flag_BNB_numuCC_from_overlaynue;
  osc_test->flag_BNB_numuCC_from_appnue        = Configure_Osc::flag_BNB_numuCC_from_appnue;
  osc_test->flag_BNB_numuCC_from_appnumu       = Configure_Osc::flag_BNB_numuCC_from_appnumu;
  osc_test->flag_BNB_numuCC_from_overlaynumuNC = Configure_Osc::flag_BNB_numuCC_from_overlaynumuNC;
  osc_test->flag_BNB_numuCC_from_overlaynueNC  = Configure_Osc::flag_BNB_numuCC_from_overlaynueNC;  
  
  osc_test->flag_BNB_CCpi0_from_overlaynumu    = Configure_Osc::flag_BNB_CCpi0_from_overlaynumu;
  osc_test->flag_BNB_CCpi0_from_appnue         = Configure_Osc::flag_BNB_CCpi0_from_appnue;
  osc_test->flag_BNB_CCpi0_from_overlaynumuNC  = Configure_Osc::flag_BNB_CCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_CCpi0_from_overlaynueNC   = Configure_Osc::flag_BNB_CCpi0_from_overlaynueNC;
  
  osc_test->flag_BNB_NCpi0_from_overlaynumu    = Configure_Osc::flag_BNB_NCpi0_from_overlaynumu;
  osc_test->flag_BNB_NCpi0_from_appnue         = Configure_Osc::flag_BNB_NCpi0_from_appnue;
  osc_test->flag_BNB_NCpi0_from_overlaynumuNC  = Configure_Osc::flag_BNB_NCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_NCpi0_from_overlaynueNC   = Configure_Osc::flag_BNB_NCpi0_from_overlaynueNC;
  
  /////// set only one time

  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,Configure_Osc::default_dirtadd_file,
                               Configure_Osc::default_mcstat_file,Configure_Osc::default_fluxXs_dir,
                               Configure_Osc::default_detector_dir);
  osc_test->Set_oscillation_base(Configure_Osc::default_eventlist_dir);
  osc_test->flag_osc_appearance_only=false; osc_test->flag_toy_clip0=true;
  osc_test->Set_oscillation_pars(0,0,0,0); osc_test->Apply_oscillation(); osc_test->Set_apply_POT();
  osc_test->Set_meas2fitdata();
  TMatrixD DATA=osc_test->matrix_tosc_fitdata_newworld;
  TMatrixD P3=osc_test->matrix_tosc_eff_newworld_pred;
  int NB=P3.GetNcols(); vector<double> buf(NB);
  const double*p3=P3.GetMatrixArray(), *dat=DATA.GetMatrixArray();
  double lg_lo=log10(X);
  printf("# dm2=%.4g  s22tue=%.4g   LOWER branch only\n",DM2,X);
  printf("#   t24        t14        chi2_4nu(data)     D(Asimov sep)\n");
  double bestD_fine=1e99, bestC_fine=1e99, t24_fine=0;
  double bestD_50=1e99,  bestC_50=1e99,  t24_50=0;
  for(int k=0;k<NFINE;k++){
    double t24=pow(10, lg_lo+(0-lg_lo)*double(k)/(NFINE-1)); if(t24<X) t24=X;
    double r=X/t24; if(r>1) r=1; double root=sqrt(1-r);
    double t14=(1-root)/2;
    osc_test->Set_oscillation_pars(DM2,t14,t24,0); osc_test->Apply_oscillation(); osc_test->Set_apply_POT();
    TMatrixD P4=osc_test->matrix_tosc_eff_newworld_pred,S4=osc_test->matrix_tosc_eff_newworld_abs_syst_total;
    TMatrixD C4=pinv(P4,S4);
    double c2=qfd(dat,P4.GetMatrixArray(),C4.GetMatrixArray(),NB,&buf[0]);
    double D =qfd(p3, P4.GetMatrixArray(),C4.GetMatrixArray(),NB,&buf[0]);
    if(c2<bestC_fine){ bestC_fine=c2; bestD_fine=D; t24_fine=t24; }
    if(k%10==0) printf(" %10.4e %10.5f %14.4f %14.4f\n", t24,t14,c2,D);
  }
  for(int k=0;k<50;k++){
    double t24=pow(10, lg_lo+(0-lg_lo)*double(k)/49.0); if(t24<X) t24=X;
    double r=X/t24; if(r>1) r=1; double root=sqrt(1-r);
    double t14=(1-root)/2;
    osc_test->Set_oscillation_pars(DM2,t14,t24,0); osc_test->Apply_oscillation(); osc_test->Set_apply_POT();
    TMatrixD P4=osc_test->matrix_tosc_eff_newworld_pred,S4=osc_test->matrix_tosc_eff_newworld_abs_syst_total;
    TMatrixD C4=pinv(P4,S4);
    double c2=qfd(dat,P4.GetMatrixArray(),C4.GetMatrixArray(),NB,&buf[0]);
    double D =qfd(p3, P4.GetMatrixArray(),C4.GetMatrixArray(),NB,&buf[0]);
    if(c2<bestC_50){ bestC_50=c2; bestD_50=D; t24_50=t24; }
  }
  printf("\n fine scan (%d pts): chi2_4nu_min=%.4f at t24=%.4e,  D there = %.4f\n",NFINE,bestC_fine,t24_fine,bestD_fine);
  printf(" 50-pt log scan   : chi2_4nu_min=%.4f at t24=%.4e,  D there = %.4f\n",bestC_50,t24_50,bestD_50);
  printf(" => d(chi2_4nu) = %.4f ,  d(D) = %.4f  (D drives the sensitivity; x ~ D^(1/2))\n",
     bestC_50-bestC_fine, bestD_50-bestD_fine);
  return 0;
}
