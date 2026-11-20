#include<iostream>
#include<sstream>
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
// diag_dscale: measure how the Asimov Dchi2 separation D(x) scales with the
// appearance amplitude x at fixed dm2, using the SAME profiled family as the
// scanner. Gaussian-CLs theory: 1-sigma band width = (1/n)*log10(7.438/1.994)
// where D ~ x^n. Published band 0.296 dex implies n = 2.
static inline double qf(const double *v,const double *M,int n){
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
  double DM2 = (argc>1)? atof(argv[1]) : 1.0;
  int n_t24 = 50;
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
  TMatrixD P3=osc_test->matrix_tosc_eff_newworld_pred;
  int NB=P3.GetNcols(); vector<double> buf(NB);
  const double*p3=P3.GetMatrixArray();
  printf("# dm2 = %.4g   [D = Asimov Dchi2 at the profiled 4nu hypothesis]\n", DM2);
  printf("#    s22tue        D        local n=dlnD/dlnx   implied band width (dex)\n");
  vector<double> xs, Ds;
  for(int k=0;k<22;k++){
    double x = pow(10, -3.6 + 0.1*k);
    double best=1e99;
    double lg_lo=log10(x), lg_hi=0;
    for(int t=0;t<n_t24;t++){
      double t24=pow(10, lg_lo+(lg_hi-lg_lo)*double(t)/(n_t24-1)); if(t24<x) t24=x;
      double r=x/t24; if(r>1) r=1; double root=sqrt(1-r);
      for(int br=0;br<2;br++){
        if(br==1&&root==0) continue;
        double t14=(br==0)?(1-root)/2:(1+root)/2;
        osc_test->Set_oscillation_pars(DM2,t14,t24,0); osc_test->Apply_oscillation(); osc_test->Set_apply_POT();
        TMatrixD P4=osc_test->matrix_tosc_eff_newworld_pred;
        TMatrixD S4=osc_test->matrix_tosc_eff_newworld_abs_syst_total;
        TMatrixD C4=pinv(P4,S4);
        double D=qfd(p3,P4.GetMatrixArray(),C4.GetMatrixArray(),NB,&buf[0]);
        if(D<best) best=D;
      }
    }
    xs.push_back(x); Ds.push_back(best);
  }
  for(size_t k=1;k+1<xs.size();k++){
    double n=(log(Ds[k+1])-log(Ds[k-1]))/(log(xs[k+1])-log(xs[k-1]));
    printf(" %12.4e %10.3f      %6.2f            %6.3f%s\n", xs[k], Ds[k], n, (1.0/n)*log10(7.438/1.994),
      (Ds[k]>2.5&&Ds[k]<6.0)?"   <-- near the 95% crossing":"");
  }
  return 0;
}
