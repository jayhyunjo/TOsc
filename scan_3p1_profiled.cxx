#include<iostream>
#include<fstream>
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
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TDecompChol.h"
#include "Math/ProbFuncMathCore.h"

// scan_3p1_profiled: FULL 3+1 scan for the Nature Fig. 3a sin2(2theta_mue)
// exclusion, following the published procedure (arXiv:2512.07159 Methods):
//   - chi2 = (N-M)^T Sigma^-1 (N-M), Pearson stat + full syst (as TOsc::FCN)
//   - 2D grid in (dm2_41, sin2_2theta_mue); at each 2D point the third freedom
//     is PROFILED: hypotheses (t14, t24) with 4*t14*(1-t14)*t24 = s22tue,
//     t24 in [s22tue, 1] (log grid), BOTH t14 branches (<1/2 and >1/2),
//     pick min chi2 -- on data for the data curve, on the 3nu Asimov for the
//     sensitivity curve.
//   - toy CLs at the profiled hypotheses: T = chi2_3nu - chi2_4nu;
//     CLs = P(T<=Tobs|4nu)/P(T<=Tobs|3nu); sens at median 3nu-toy T.
//
// Usage: ./scan_3p1_profiled -j <ijob> -n <njobs> [-o outdir] [-t ntoys] [-k n_t24]

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
  int ijob = 0, njobs = 1, ntoys = 2000, n_t24 = 20, gate = 0;
  int num_dm2 = 60, num_amp = 60;
  TString outdir = "output_scan_3p1";

  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-j")==0 ) { stringstream c(argv[i+1]); c>>ijob; }
    if( strcmp(argv[i],"-n")==0 ) { stringstream c(argv[i+1]); c>>njobs; }
    if( strcmp(argv[i],"-t")==0 ) { stringstream c(argv[i+1]); c>>ntoys; }
    if( strcmp(argv[i],"-k")==0 ) { stringstream c(argv[i+1]); c>>n_t24; }
    if( strcmp(argv[i],"-o")==0 ) { outdir = argv[i+1]; }
    if( strcmp(argv[i],"-gdm2")==0 ) { stringstream c(argv[i+1]); c>>num_dm2; }
    if( strcmp(argv[i],"-gamp")==0 ) { stringstream c(argv[i+1]); c>>num_amp; }
    if( strcmp(argv[i],"-G")==0 )    { stringstream c(argv[i+1]); c>>gate; }
  }

  gROOT->SetBatch(1);

  const int    NUM_dm2 = num_dm2;  const double DM2_LO = -2, DM2_HI = 2;
  const int    NUM_amp = num_amp;  const double AMP_LO = -4, AMP_HI = 0;
  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "", NUM_dm2, DM2_LO, DM2_HI);
  TH1D *h1d_amp = new TH1D("h1d_amp", "", NUM_amp, AMP_LO, AMP_HI);

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

  //////////////////////////////////////////// H0 = 3nu

  osc_test->Set_oscillation_pars(0, 0, 0, 0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();
  TMatrixD P3 = osc_test->matrix_tosc_eff_newworld_pred;
  TMatrixD S3 = osc_test->matrix_tosc_eff_newworld_abs_syst_total;
  TMatrixD C3inv = pearson_inv(P3, S3);

  osc_test->Set_toy_variations(ntoys);
  vector<TMatrixD> toys3; toys3.reserve(ntoys);
  for(int i=1;i<=ntoys;i++) toys3.push_back( osc_test->map_matrix_tosc_toy_pred[i] );

  vector<double> chi2_3_toys3(ntoys);
  for(int i=0;i<ntoys;i++) chi2_3_toys3[i] = quadform(toys3[i] - P3, C3inv);

  osc_test->Set_meas2fitdata();
  TMatrixD D = osc_test->matrix_tosc_fitdata_newworld;
  double chi2_3_data = quadform(D - P3, C3inv);
  cout<<endl<<" ---> chi2_3nu(data) = "<<chi2_3_data<<endl;

  //////////////////////////////////////////// output

  system( TString::Format("mkdir -p %s", outdir.Data()) );
  TFile *outfile = new TFile(TString::Format("%s/scan_3p1_%03d.root", outdir.Data(), ijob), "recreate");
  TTree *outtree = new TTree("tree", "full 3+1 profiled toy-CLs grid");
  int    t_idm2, t_iamp;
  double t_dm2, t_amp;
  double t_chi2d_min, t_t14_d, t_t24_d;     // data-profiled hypothesis
  double t_chi2a_min, t_t14_a, t_t24_a;     // asimov-profiled hypothesis
  double t_Tobs, t_p3obs, t_p4obs, t_cls_data;
  double t_Tmed3, t_p4med, t_cls_sens;
  outtree->Branch("idm2",&t_idm2,"idm2/I");             outtree->Branch("iamp",&t_iamp,"iamp/I");
  outtree->Branch("dm2",&t_dm2,"dm2/D");                outtree->Branch("s22tue",&t_amp,"s22tue/D");
  outtree->Branch("chi2d_min",&t_chi2d_min,"chi2d_min/D");
  outtree->Branch("t14_d",&t_t14_d,"t14_d/D");          outtree->Branch("t24_d",&t_t24_d,"t24_d/D");
  outtree->Branch("chi2a_min",&t_chi2a_min,"chi2a_min/D");
  outtree->Branch("t14_a",&t_t14_a,"t14_a/D");          outtree->Branch("t24_a",&t_t24_a,"t24_a/D");
  outtree->Branch("T_obs",&t_Tobs,"T_obs/D");
  outtree->Branch("p3_obs",&t_p3obs,"p3_obs/D");        outtree->Branch("p4_obs",&t_p4obs,"p4_obs/D");
  outtree->Branch("cls_data",&t_cls_data,"cls_data/D");
  outtree->Branch("T_med3",&t_Tmed3,"T_med3/D");        outtree->Branch("p4_med",&t_p4med,"p4_med/D");
  outtree->Branch("cls_sens",&t_cls_sens,"cls_sens/D");
  int t_gated_d, t_gated_s; double t_clsg_data, t_clsg_sens, t_dchi2a_d, t_dchi2a_s;
  outtree->Branch("gated_data",&t_gated_d,"gated_data/I");
  outtree->Branch("gated_sens",&t_gated_s,"gated_sens/I");
  outtree->Branch("clsg_data",&t_clsg_data,"clsg_data/D");
  outtree->Branch("clsg_sens",&t_clsg_sens,"clsg_sens/D");
  outtree->Branch("dchi2a_d",&t_dchi2a_d,"dchi2a_d/D");
  outtree->Branch("dchi2a_s",&t_dchi2a_s,"dchi2a_s/D");

  //////////////////////////////////////////// grid loop

  auto apply_point=[&](double dm2, double t14, double t24, TMatrixD &P4, TMatrixD &S4){
    osc_test->Set_oscillation_pars(dm2, t14, t24, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();
    P4.ResizeTo(osc_test->matrix_tosc_eff_newworld_pred);     P4 = osc_test->matrix_tosc_eff_newworld_pred;
    S4.ResizeTo(osc_test->matrix_tosc_eff_newworld_abs_syst_total); S4 = osc_test->matrix_tosc_eff_newworld_abs_syst_total;
  };

  int ndone = 0;
  for(int ipt=0; ipt<NUM_dm2*NUM_amp; ipt++) {
    if( ipt % njobs != ijob ) continue;
    t_idm2 = ipt / NUM_amp;  t_iamp = ipt % NUM_amp;
    t_dm2 = pow(10, h1d_dm2->GetBinCenter(t_idm2+1));
    t_amp = pow(10, h1d_amp->GetBinCenter(t_iamp+1));

    /////////////// profiling: hypotheses family 4*t14*(1-t14)*t24 = s22tue
    t_chi2d_min = 1e99; t_chi2a_min = 1e99;
    double lg_lo = log10(t_amp), lg_hi = 0;// t24 from s22tue to 1
    for(int k=0;k<n_t24;k++){
      double t24 = pow(10, lg_lo + (lg_hi-lg_lo)*( n_t24==1 ? 1.0 : double(k)/(n_t24-1) ));
      if( t24 < t_amp ) t24 = t_amp;
      double r = t_amp/t24; if( r>1 ) r = 1;
      double root = sqrt(1-r);
      for(int br=0; br<2; br++){
        double t14 = (br==0) ? (1-root)/2 : (1+root)/2;
        if( br==1 && root==0 ) continue;// branches coincide
        TMatrixD P4, S4;
        apply_point(t_dm2, t14, t24, P4, S4);
        TMatrixD C4inv = pearson_inv(P4, S4);
        double c2d = quadform(D  - P4, C4inv);
        double c2a = quadform(P3 - P4, C4inv);
        if( c2d < t_chi2d_min ){ t_chi2d_min = c2d; t_t14_d = t14; t_t24_d = t24; }
        if( c2a < t_chi2a_min ){ t_chi2a_min = c2a; t_t14_a = t14; t_t24_a = t24; }
      }
    }

    /////////////// toys at the DATA-profiled hypothesis -> cls_data
    {
      TMatrixD P4, S4;
      apply_point(t_dm2, t_t14_d, t_t24_d, P4, S4);
      TMatrixD C4inv = pearson_inv(P4, S4);
      t_Tobs = chi2_3_data - quadform(D - P4, C4inv);
      t_dchi2a_d = quadform(P3 - P4, C4inv);
      // Gaussian-CLs estimate (arXiv:1407.5052): T~N(+-D, 4D)
      { double Dv = t_dchi2a_d>1e-9 ? t_dchi2a_d : 1e-9; double sg = 2*sqrt(Dv);
        t_clsg_data = ROOT::Math::normal_cdf((t_Tobs - Dv)/sg) / ROOT::Math::normal_cdf((t_Tobs + Dv)/sg); }
      t_gated_d = ( gate && (t_clsg_data < 0.001 || t_clsg_data > 0.85) ) ? 1 : 0;
      if( t_gated_d ){
        t_p3obs = -1; t_p4obs = -1; t_cls_data = t_clsg_data;
      } else {
      osc_test->Set_toy_variations(ntoys);
      vector<double> T3(ntoys), T4(ntoys);
      for(int i=0;i<ntoys;i++){
        const TMatrixD &x4 = osc_test->map_matrix_tosc_toy_pred[i+1];
        T4[i] = quadform(x4 - P3, C3inv) - quadform(x4 - P4, C4inv);
        const TMatrixD &x3 = toys3[i];
        T3[i] = chi2_3_toys3[i] - quadform(x3 - P4, C4inv);
      }
      sort(T3.begin(), T3.end()); sort(T4.begin(), T4.end());
      auto tailfrac=[&](vector<double>&v, double t)->double{
        double c = upper_bound(v.begin(), v.end(), t) - v.begin();
        if( c < 0.5 ) c = 0.5; if( c > ntoys-0.5 ) c = ntoys;
        return c/ntoys; };
      t_p3obs = tailfrac(T3, t_Tobs);
      t_p4obs = tailfrac(T4, t_Tobs);
      t_cls_data = t_p4obs / t_p3obs;
      }
    }

    /////////////// toys at the ASIMOV-profiled hypothesis -> cls_sens
    {
      TMatrixD P4, S4;
      apply_point(t_dm2, t_t14_a, t_t24_a, P4, S4);
      TMatrixD C4inv = pearson_inv(P4, S4);
      t_dchi2a_s = quadform(P3 - P4, C4inv);
      { double Dv = t_dchi2a_s>1e-9 ? t_dchi2a_s : 1e-9;
        t_clsg_sens = 2*ROOT::Math::normal_cdf(-sqrt(Dv)); }
      t_gated_s = ( gate && (t_clsg_sens < 0.001 || t_clsg_sens > 0.85) ) ? 1 : 0;
      if( t_gated_s ){
        t_Tmed3 = 0; t_p4med = -1; t_cls_sens = t_clsg_sens;
      } else {
      osc_test->Set_toy_variations(ntoys);
      vector<double> T3(ntoys), T4(ntoys);
      for(int i=0;i<ntoys;i++){
        const TMatrixD &x4 = osc_test->map_matrix_tosc_toy_pred[i+1];
        T4[i] = quadform(x4 - P3, C3inv) - quadform(x4 - P4, C4inv);
        const TMatrixD &x3 = toys3[i];
        T3[i] = chi2_3_toys3[i] - quadform(x3 - P4, C4inv);
      }
      sort(T3.begin(), T3.end()); sort(T4.begin(), T4.end());
      auto tailfrac=[&](vector<double>&v, double t)->double{
        double c = upper_bound(v.begin(), v.end(), t) - v.begin();
        if( c < 0.5 ) c = 0.5; if( c > ntoys-0.5 ) c = ntoys;
        return c/ntoys; };
      t_Tmed3 = 0.5*( T3[ntoys/2] + T3[(ntoys-1)/2] );
      t_p4med = tailfrac(T4, t_Tmed3);
      t_cls_sens = t_p4med / tailfrac(T3, t_Tmed3);
      }
    }

    outtree->Fill();
    ndone++;
    if( ndone % 10 == 0 ) cout<<TString::Format("     done %5d pts (ipt %d: cls_data %.4f, cls_sens %.4f, t24hat_d %.3g)", ndone, ipt, t_cls_data, t_cls_sens, t_t24_d)<<endl;
  }

  outtree->Write();
  outfile->Close();
  cout<<endl<<" ---> job "<<ijob<<"/"<<njobs<<" finished, "<<ndone<<" points"<<endl;
  return 0;
}
