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
// diag_band3: at fixed dm2, scan amplitude and track for each SM toy its
// standardised position z within the 3nu pdf. Gaussian-CLs theory assumes z is
// constant in x; if z regresses towards 0 as x grows, per-toy crossings bunch
// together and the +-1 sigma band narrows. Also reports the band width from the
// empirical CLs crossings.
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
  double DM2=(argc>1)?atof(argv[1]):1.0;
  int NTOY=10000, NSM=1000, n_t24=50, NX=26;
  int PROFMODE=(argc>2)?atoi(argv[2]):0;// 0 = profile on 3nu Asimov, 1 = profile on DATA
  double LX0=-3.4, LX1=-2.2;
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
  TMatrixD P3=osc_test->matrix_tosc_eff_newworld_pred, S3=osc_test->matrix_tosc_eff_newworld_abs_syst_total;
  TMatrixD C3=pinv(P3,S3);
  int NB=P3.GetNcols(); vector<double> buf(NB);
  const double *p3=P3.GetMatrixArray(), *c3=C3.GetMatrixArray();
  osc_test->Set_toy_variations(NTOY);
  vector<double> x3((size_t)NTOY*NB), chi3_x3(NTOY);
  for(int i=1;i<=NTOY;i++){ const double*s=osc_test->map_matrix_tosc_toy_pred[i].GetMatrixArray();
    copy(s,s+NB,x3.begin()+(size_t)(i-1)*NB); }
  for(int i=0;i<NTOY;i++) chi3_x3[i]=qfd(&x3[(size_t)i*NB],p3,c3,NB,&buf[0]);
  osc_test->Set_toy_variations(NSM);
  vector<double> xs((size_t)NSM*NB), chi3_xs(NSM);
  for(int i=1;i<=NSM;i++){ const double*s=osc_test->map_matrix_tosc_toy_pred[i].GetMatrixArray();
    copy(s,s+NB,xs.begin()+(size_t)(i-1)*NB); }
  for(int i=0;i<NSM;i++) chi3_xs[i]=qfd(&xs[(size_t)i*NB],p3,c3,NB,&buf[0]);
  vector<vector<double>> Z(NX, vector<double>(NSM)), CLS(NX, vector<double>(NSM));
  vector<double> XS(NX), DS(NX);
  printf("# dm2=%.3g  PROFMODE=%d (0=Asimov-profiled, 1=data-profiled)\n",DM2,PROFMODE);
  printf("#   x            D      corr(z,z_ref)   sd(z)\n");
  for(int k=0;k<NX;k++){
    double lx = BINCENTRE ? (LX0+(LX1-LX0)*(k+0.5)/NX) : (LX0+(LX1-LX0)*k/(NX-1));
    double x=pow(10,lx); XS[k]=lx;
    double best=1e99,bt14=0,bt24=0;
    for(int t=0;t<n_t24;t++){
      double t24=pow(10,lx+(0-lx)*double(t)/(n_t24-1)); if(t24<x) t24=x;
      double r=x/t24; if(r>1)r=1; double root=sqrt(1-r);
      for(int br=0;br<2;br++){ if(br==1&&root==0) continue;
        double t14=(br==0)?(1-root)/2:(1+root)/2;
        osc_test->Set_oscillation_pars(DM2,t14,t24,0); osc_test->Apply_oscillation(); osc_test->Set_apply_POT();
        TMatrixD Pp=osc_test->matrix_tosc_eff_newworld_pred,Ss=osc_test->matrix_tosc_eff_newworld_abs_syst_total;
        TMatrixD Cc=pinv(Pp,Ss);
        double D   = qfd(p3,Pp.GetMatrixArray(),Cc.GetMatrixArray(),NB,&buf[0]);
        double Dsel= (PROFMODE==0)? D
                   : qfd(DATA.GetMatrixArray(),Pp.GetMatrixArray(),Cc.GetMatrixArray(),NB,&buf[0]);
        if(Dsel<best){best=Dsel;bt14=t14;bt24=t24;} } }
    osc_test->Set_oscillation_pars(DM2,bt14,bt24,0); osc_test->Apply_oscillation(); osc_test->Set_apply_POT();
    TMatrixD P4=osc_test->matrix_tosc_eff_newworld_pred,S4=osc_test->matrix_tosc_eff_newworld_abs_syst_total;
    TMatrixD C4=pinv(P4,S4);
    const double *p4=P4.GetMatrixArray(), *c4=C4.GetMatrixArray();
    DS[k]=best;
    vector<double> T3(NTOY);
    for(int i=0;i<NTOY;i++) T3[i]=qfd(&x3[(size_t)i*NB],p4,c4,NB,&buf[0])-chi3_x3[i];
    double mu=0; for(double v:T3) mu+=v; mu/=NTOY;
    double sd=0; for(double v:T3) sd+=(v-mu)*(v-mu); sd=sqrt(sd/(NTOY-1));
    osc_test->Set_toy_variations(NTOY);
    vector<double> T4(NTOY);
    for(int i=0;i<NTOY;i++){ const double*x4=osc_test->map_matrix_tosc_toy_pred[i+1].GetMatrixArray();
      T4[i]=qfd(x4,p4,c4,NB,&buf[0])-qfd(x4,p3,c3,NB,&buf[0]); }
    sort(T3.begin(),T3.end()); sort(T4.begin(),T4.end());
    for(int m=0;m<NSM;m++){
      double T=qfd(&xs[(size_t)m*NB],p4,c4,NB,&buf[0])-chi3_xs[m];
      Z[k][m]=(T-mu)/sd;
      int nb=T3.end()-upper_bound(T3.begin(),T3.end(),T);
      int ns=T4.end()-upper_bound(T4.begin(),T4.end(),T);
      double CLb=nb*1./NTOY, CLsb=ns*1./NTOY;
      CLS[k][m]= (CLb==0)? 1.0 : CLsb/CLb;
    }
    double sx=0,sy=0,sxx=0,syy=0,sxy=0;
    for(int m=0;m<NSM;m++){ double a=Z[0][m],b=Z[k][m]; sx+=a;sy+=b;sxx+=a*a;syy+=b*b;sxy+=a*b; }
    double corr=(NSM*sxy-sx*sy)/sqrt((NSM*sxx-sx*sx)*(NSM*syy-sy*sy));
    double zsd=sqrt((NSM*syy-sy*sy))/NSM;
    { vector<double> cc(NSM); for(int m=0;m<NSM;m++) cc[m]=CLS[k][m];
      sort(cc.begin(),cc.end());
      printf(" %10.3e %9.3f    %7.4f      %6.3f   CLs q16/50/84 = %8.4f %8.4f %8.4f\n",
        x, best, corr, zsd, cc[(int)(.1587*NSM)], cc[NSM/2], cc[(int)(.8413*NSM)]); }
  }
  vector<double> L;
  for(int m=0;m<NSM;m++) for(int k=1;k<NX;k++){
    double c0=CLS[k-1][m],c1=CLS[k][m];
    if(c0>0.05&&c1<=0.05){ L.push_back(XS[k-1]+(XS[k]-XS[k-1])*(0.05-c0)/(c1-c0)); break; } }
  sort(L.begin(),L.end()); int n=L.size();
  if(n>20) printf("\n band from empirical CLs crossings: width=%.3f dex  median=%.3e  (N=%d)\n",
     L[(int)(.8413*n)]-L[(int)(.1587*n)], pow(10,L[n/2]), n);
  return 0;
}
