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
// diag_band2: at one (dm2, x) point, compare the Dchi2 distribution of the SM
// toys (band ensemble) with that of the 3nu pdf pseudos. They are thrown by the
// same generator at the same hypothesis, so they MUST agree; CL_b over SM toys
// must then be uniform(0,1). Any deviation localises the band compression.
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
static double quant(vector<double> v,double q){ sort(v.begin(),v.end()); return v[(int)(q*v.size())]; }
int main(int argc,char**argv){
  gROOT->SetBatch(1);
  double DM2=(argc>1)?atof(argv[1]):1.0;
  double XAMP=(argc>2)?atof(argv[2]):2.0e-3;
  int NTOY=10000, NSM=1000, n_t24=50;
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
  TMatrixD P3=osc_test->matrix_tosc_eff_newworld_pred, S3=osc_test->matrix_tosc_eff_newworld_abs_syst_total;
  TMatrixD C3=pinv(P3,S3);
  int NB=P3.GetNcols(); vector<double> buf(NB);
  const double *p3=P3.GetMatrixArray(), *c3=C3.GetMatrixArray();
  // pdf 3nu toys
  osc_test->Set_toy_variations(NTOY);
  vector<double> x3((size_t)NTOY*NB), chi3_x3(NTOY);
  for(int i=1;i<=NTOY;i++){ const double*s=osc_test->map_matrix_tosc_toy_pred[i].GetMatrixArray();
    copy(s,s+NB,x3.begin()+(size_t)(i-1)*NB); }
  for(int i=0;i<NTOY;i++) chi3_x3[i]=qfd(&x3[(size_t)i*NB],p3,c3,NB,&buf[0]);
  // SM toys
  osc_test->Set_toy_variations(NSM);
  vector<double> xs((size_t)NSM*NB), chi3_xs(NSM);
  for(int i=1;i<=NSM;i++){ const double*s=osc_test->map_matrix_tosc_toy_pred[i].GetMatrixArray();
    copy(s,s+NB,xs.begin()+(size_t)(i-1)*NB); }
  for(int i=0;i<NSM;i++) chi3_xs[i]=qfd(&xs[(size_t)i*NB],p3,c3,NB,&buf[0]);
  // profiled 4nu hypothesis at (DM2, XAMP): minimise Asimov separation like the scanner does on data
  double best=1e99,bt14=0,bt24=0;
  double lg_lo=log10(XAMP);
  for(int t=0;t<n_t24;t++){
    double t24=pow(10,lg_lo+(0-lg_lo)*double(t)/(n_t24-1)); if(t24<XAMP) t24=XAMP;
    double r=XAMP/t24; if(r>1)r=1; double root=sqrt(1-r);
    for(int br=0;br<2;br++){ if(br==1&&root==0) continue;
      double t14=(br==0)?(1-root)/2:(1+root)/2;
      osc_test->Set_oscillation_pars(DM2,t14,t24,0); osc_test->Apply_oscillation(); osc_test->Set_apply_POT();
      TMatrixD P4=osc_test->matrix_tosc_eff_newworld_pred,S4=osc_test->matrix_tosc_eff_newworld_abs_syst_total;
      TMatrixD C4=pinv(P4,S4);
      double D=qfd(p3,P4.GetMatrixArray(),C4.GetMatrixArray(),NB,&buf[0]);
      if(D<best){best=D;bt14=t14;bt24=t24;} } }
  osc_test->Set_oscillation_pars(DM2,bt14,bt24,0); osc_test->Apply_oscillation(); osc_test->Set_apply_POT();
  TMatrixD P4=osc_test->matrix_tosc_eff_newworld_pred,S4=osc_test->matrix_tosc_eff_newworld_abs_syst_total;
  TMatrixD C4=pinv(P4,S4);
  const double *p4=P4.GetMatrixArray(), *c4=C4.GetMatrixArray();
  printf("\n=== dm2=%.3g  x=%.3e   profiled t14=%.4f t24=%.4g   Asimov D=%.3f ===\n",DM2,XAMP,bt14,bt24,best);
  // Dchi2 of both ensembles
  vector<double> Dpdf(NTOY), Dsm(NSM);
  for(int i=0;i<NTOY;i++) Dpdf[i]=qfd(&x3[(size_t)i*NB],p4,c4,NB,&buf[0])-chi3_x3[i];
  for(int i=0;i<NSM;i++)  Dsm[i] =qfd(&xs[(size_t)i*NB],p4,c4,NB,&buf[0])-chi3_xs[i];
  printf(" Dchi2 quantiles      16%%        50%%        84%%     (84-16 spread)\n");
  printf("   pdf 3nu toys  %10.3f %10.3f %10.3f   %8.3f\n",
     quant(Dpdf,.1587),quant(Dpdf,.5),quant(Dpdf,.8413),quant(Dpdf,.8413)-quant(Dpdf,.1587));
  printf("   SM band toys  %10.3f %10.3f %10.3f   %8.3f\n",
     quant(Dsm,.1587),quant(Dsm,.5),quant(Dsm,.8413),quant(Dsm,.8413)-quant(Dsm,.1587));
  printf("   theory N(+D,4D): %7.3f %10.3f %10.3f   %8.3f\n",
     best-2*sqrt(best), best, best+2*sqrt(best), 4*sqrt(best));
  // CL_b over SM toys (should be uniform)
  sort(Dpdf.begin(),Dpdf.end());
  vector<double> clb(NSM);
  for(int m=0;m<NSM;m++){ int nb=Dpdf.end()-upper_bound(Dpdf.begin(),Dpdf.end(),Dsm[m]); clb[m]=nb*1./NTOY; }
  printf(" CL_b over SM toys: 16%%=%.3f 50%%=%.3f 84%%=%.3f  (uniform expects 0.159/0.500/0.841)\n",
     quant(clb,.1587),quant(clb,.5),quant(clb,.8413));
  return 0;
}
