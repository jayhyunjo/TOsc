#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstring>
#include "stdlib.h"
using namespace std;
#include<vector>
#include<map>
#include "WCPLEEANA/TOsc.h"
#include "WCPLEEANA/Configure_Osc.h"
#include "TROOT.h"
#include "TDecompChol.h"

// check_dchi2_map: point-level validation of the TOsc chi2 against the OFFICIAL
// published Delta-chi2 map (Zenodo 10.5281/zenodo.17161263, NuMI-only, full 3+1,
// Pearson). Reads lines "dm2 s2t14 s2t24 dchi2_official", evaluates our
// chi2(data | 4nu point) on the NuMI subset, and reports chi2_ours - dchi2_official
// (should be CONSTANT = chi2_globalmin if the machinery matches).
// Usage: ./check_dchi2_map <points.txt>

static double quadform(const TMatrixD &v, const TMatrixD &Minv)
{
  int n = v.GetNcols(); double sum=0;
  for(int i=0;i<n;i++){ double vi=v(0,i); if(vi==0) continue;
    const double *row = Minv.GetMatrixArray() + (Long64_t)i*n;
    double acc=0; for(int j=0;j<n;j++) acc += row[j]*v(0,j);
    sum += vi*acc; }
  return sum;
}
static TMatrixD pearson_inv(const TMatrixD &P, const TMatrixD &S)
{
  int n=P.GetNcols(); TMatrixDSym C(n);
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) C(i,j)=0.5*(S(i,j)+S(j,i));
  for(int i=0;i<n;i++){ double st=P(0,i); if(st==0) st=1e-6; C(i,i)+=st; }
  TDecompChol chol(C);
  if(!chol.Decompose()){ cerr<<" ERROR: Cholesky failed"<<endl; exit(1); }
  chol.Invert(C);
  return TMatrixD(C);
}

int main(int argc, char** argv)
{
  if(argc<2){ cerr<<"usage: check_dchi2_map <points.txt>"<<endl; return 1; }
  gROOT->SetBatch(1);

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


  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
                               Configure_Osc::default_dirtadd_file,
                               Configure_Osc::default_mcstat_file,
                               Configure_Osc::default_fluxXs_dir,
                               Configure_Osc::default_detector_dir);
  osc_test->Set_oscillation_base(Configure_Osc::default_eventlist_dir);

  osc_test->flag_osc_appearance_only = false;// FULL 3+1

  int lo = 26*7, hi = 26*14 - 1;// NuMI subset

  // data (independent of pars; need one Apply first to size matrices)
  osc_test->Set_oscillation_pars(0,0,0,0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();
  osc_test->Set_meas2fitdata();
  TMatrixD D = osc_test->matrix_tosc_fitdata_newworld.GetSub(0,0, lo,hi);

  ifstream fin(argv[1]);
  double dm2, t14, t24, dco;
  vector<double> diffs; vector<double> v_dm2, v_dco, v_ours;
  while( fin>>dm2>>t14>>t24>>dco ){
    osc_test->Set_oscillation_pars(dm2, t14, t24, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();
    TMatrixD P = osc_test->matrix_tosc_eff_newworld_pred.GetSub(0,0, lo,hi);
    TMatrixD S = osc_test->matrix_tosc_eff_newworld_abs_syst_total.GetSub(lo,hi, lo,hi);
    TMatrixD Cinv = pearson_inv(P, S);
    double chi2 = quadform(D - P, Cinv);
    diffs.push_back(chi2 - dco); v_dm2.push_back(dm2); v_dco.push_back(dco); v_ours.push_back(chi2);
  }

  int n = diffs.size();
  double mean=0; for(double d : diffs) mean += d; mean /= n;
  double rms=0, mx=0; int imx=0;
  for(int i=0;i<n;i++){ double r=diffs[i]-mean; rms+=r*r; if(fabs(r)>mx){mx=fabs(r); imx=i;} }
  rms = sqrt(rms/n);

  printf("\n ---> %d points: (chi2_ours - dchi2_official) mean = %.4f (= chi2_globalmin if consistent)\n", n, mean);
  printf(" ---> residual about mean: rms = %.4f, max|res| = %.4f (at dm2 %.4g, official dchi2 %.4f, ours %.4f)\n",
         rms, mx, v_dm2[imx], v_dco[imx], v_ours[imx]);
  printf("\n     first 10 points:  dm2        dchi2_official   chi2_ours   diff\n");
  for(int i=0;i<10 && i<n;i++)
    printf("     %18.6g  %12.4f  %12.4f  %10.4f\n", v_dm2[i], v_dco[i], v_ours[i], v_ours[i]-v_dco[i]);
  return 0;
}
