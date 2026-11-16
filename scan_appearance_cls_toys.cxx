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

// scan_appearance_cls_toys: nue appearance-only chi2 grid with TOY-BASED CLs,
// following the published procedure (Nature 648, 64 (2025) / arXiv:2512.07159):
//   pseudo-experiments generated with the analysis' own recipe
//   (TOsc::Set_toy_variations: Gaussian fluctuation along eigenvectors of the
//   systematic covariance + bin-wise Poisson), under both H0=3nu and H1=4nu.
//   Test statistic T(x) = chi2_3nu(x) - chi2_4nu(x) (Pearson, as in TOsc::FCN).
//   CLs = P(T <= T_obs | 4nu) / P(T <= T_obs | 3nu); exclusion at CLs < 0.05.
//   Sensitivity: CLs evaluated at the median T of the 3nu toys.
//
// Usage: ./scan_appearance_cls_toys -j <ijob> -n <njobs> [-o outdir] [-b beam] [-t ntoys]
//        beam: 0 = BNB+NuMI (default), 1 = BNB only, 2 = NuMI only
//        ntoys: pseudo-experiments per hypothesis per grid point (default 1000)

static double quadform(const TMatrixD &v, const TMatrixD &Minv)// v: 1xN
{
  int n = v.GetNcols();
  double sum = 0;
  for(int i=0;i<n;i++){
    double vi = v(0,i);
    if( vi==0 ) continue;
    const double *row = Minv.GetMatrixArray() + (Long64_t)i*n;
    double acc = 0;
    for(int j=0;j<n;j++) acc += row[j]*v(0,j);
    sum += vi*acc;
  }
  return sum;
}

static TMatrixD pearson_inv(const TMatrixD &P, const TMatrixD &S)// P: 1xN pred, S: NxN syst
{
  int n = P.GetNcols();
  TMatrixDSym C(n);
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) C(i,j) = 0.5*( S(i,j) + S(j,i) );
  for(int i=0;i<n;i++){
    double stat = P(0,i);
    if( stat==0 ) stat = 1e-6;// as in TOsc::FCN
    C(i,i) += stat;
  }
  TDecompChol chol(C);
  if( !chol.Decompose() ) { cerr<<" ERROR: Cholesky failed"<<endl; exit(1); }
  chol.Invert(C);
  return TMatrixD(C);
}

int main(int argc, char** argv)
{
  int ijob = 0, njobs = 1, beam = 0, ntoys = 1000;
  TString outdir = "output_scan_app_toys";

  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-j")==0 ) { stringstream c(argv[i+1]); c>>ijob; }
    if( strcmp(argv[i],"-n")==0 ) { stringstream c(argv[i+1]); c>>njobs; }
    if( strcmp(argv[i],"-b")==0 ) { stringstream c(argv[i+1]); c>>beam; }
    if( strcmp(argv[i],"-t")==0 ) { stringstream c(argv[i+1]); c>>ntoys; }
    if( strcmp(argv[i],"-o")==0 ) { outdir = argv[i+1]; }
  }

  gROOT->SetBatch(1);

  const int    NUM_dm2 = 60;  const double DM2_LO = -2, DM2_HI = 2;
  const int    NUM_amp = 60;  const double AMP_LO = -4, AMP_HI = 0;
  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "", NUM_dm2, DM2_LO, DM2_HI);
  TH1D *h1d_amp = new TH1D("h1d_amp", "", NUM_amp, AMP_LO, AMP_HI);

  //////////////////////////////////////////// TOsc initialization

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

  osc_test->flag_osc_appearance_only = true;

  //////////////////////////////////////////// beam sub-range (newworld: 26*7 BNB + 26*7 NuMI)

  int lo = 0, hi = 26*14 - 1;
  if( beam==1 ) { lo = 0;     hi = 26*7 - 1;  }
  if( beam==2 ) { lo = 26*7;  hi = 26*14 - 1; }
  int nsub = hi - lo + 1;
  cout<<endl<<TString::Format(" ---> beam mode %d: bins [%d, %d], ntoys %d/hypothesis", beam, lo, hi, ntoys)<<endl;

  //////////////////////////////////////////// H0 = 3nu: pred, syst, toys, data

  osc_test->Set_oscillation_pars(0, 0, 0, 0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();

  TMatrixD P3 = osc_test->matrix_tosc_eff_newworld_pred.GetSub(0,0, lo,hi);
  TMatrixD S3 = osc_test->matrix_tosc_eff_newworld_abs_syst_total.GetSub(lo,hi, lo,hi);
  TMatrixD C3inv = pearson_inv(P3, S3);

  osc_test->Set_toy_variations(ntoys);
  vector<TMatrixD> toys3; toys3.reserve(ntoys);
  for(int i=1;i<=ntoys;i++) toys3.push_back( osc_test->map_matrix_tosc_toy_pred[i].GetSub(0,0, lo,hi) );

  osc_test->Set_meas2fitdata();
  TMatrixD D = osc_test->matrix_tosc_fitdata_newworld.GetSub(0,0, lo,hi);

  double chi2_3_data = quadform(D - P3, C3inv);
  cout<<" ---> chi2_3nu(data) = "<<chi2_3_data<<"  ("<<nsub<<" bins)"<<endl;

  //////////////////////////////////////////// output

  system( TString::Format("mkdir -p %s", outdir.Data()) );
  TFile *outfile = new TFile(TString::Format("%s/scan_app_toys_%03d.root", outdir.Data(), ijob), "recreate");
  TTree *outtree = new TTree("tree", "appearance-only toy-CLs grid");
  int    t_idm2, t_iamp;
  double t_dm2, t_amp, t_dchi2_asimov, t_Tobs, t_p3obs, t_p4obs, t_cls_data, t_Tmed3, t_p4med, t_cls_sens;
  outtree->Branch("idm2",&t_idm2,"idm2/I");           outtree->Branch("iamp",&t_iamp,"iamp/I");
  outtree->Branch("dm2",&t_dm2,"dm2/D");              outtree->Branch("s22tue",&t_amp,"s22tue/D");
  outtree->Branch("dchi2_asimov",&t_dchi2_asimov,"dchi2_asimov/D");
  outtree->Branch("T_obs",&t_Tobs,"T_obs/D");
  outtree->Branch("p3_obs",&t_p3obs,"p3_obs/D");      outtree->Branch("p4_obs",&t_p4obs,"p4_obs/D");
  outtree->Branch("cls_data",&t_cls_data,"cls_data/D");
  outtree->Branch("T_med3",&t_Tmed3,"T_med3/D");      outtree->Branch("p4_med",&t_p4med,"p4_med/D");
  outtree->Branch("cls_sens",&t_cls_sens,"cls_sens/D");

  //////////////////////////////////////////// grid loop

  int ndone = 0;
  for(int ipt=0; ipt<NUM_dm2*NUM_amp; ipt++) {
    if( ipt % njobs != ijob ) continue;
    t_idm2 = ipt / NUM_amp;  t_iamp = ipt % NUM_amp;
    t_dm2 = pow(10, h1d_dm2->GetBinCenter(t_idm2+1));
    t_amp = pow(10, h1d_amp->GetBinCenter(t_iamp+1));

    // H1 = 4nu at the grid point
    osc_test->Set_oscillation_pars(t_dm2, t_amp, 0, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();
    TMatrixD P4 = osc_test->matrix_tosc_eff_newworld_pred.GetSub(0,0, lo,hi);
    TMatrixD S4 = osc_test->matrix_tosc_eff_newworld_abs_syst_total.GetSub(lo,hi, lo,hi);
    TMatrixD C4inv = pearson_inv(P4, S4);

    osc_test->Set_toy_variations(ntoys);

    t_dchi2_asimov = quadform(P3 - P4, C4inv);
    t_Tobs = chi2_3_data - quadform(D - P4, C4inv);

    vector<double> T3(ntoys), T4(ntoys);
    for(int i=0;i<ntoys;i++){
      TMatrixD x4 = osc_test->map_matrix_tosc_toy_pred[i+1].GetSub(0,0, lo,hi);
      T4[i] = quadform(x4 - P3, C3inv) - quadform(x4 - P4, C4inv);
      const TMatrixD &x3 = toys3[i];
      T3[i] = quadform(x3 - P3, C3inv) - quadform(x3 - P4, C4inv);
    }
    sort(T3.begin(), T3.end());
    sort(T4.begin(), T4.end());

    auto tailfrac=[&](vector<double>&v, double t)->double{
      double c = upper_bound(v.begin(), v.end(), t) - v.begin();
      if( c < 0.5 ) c = 0.5;                 // floor: avoid 0/0, half-a-toy convention
      if( c > ntoys-0.5 ) c = ntoys;         // cap
      return c/ntoys;
    };

    t_p3obs = tailfrac(T3, t_Tobs);
    t_p4obs = tailfrac(T4, t_Tobs);
    t_cls_data = t_p4obs / t_p3obs;

    t_Tmed3 = 0.5*( T3[ntoys/2] + T3[(ntoys-1)/2] );
    t_p4med = tailfrac(T4, t_Tmed3);
    t_cls_sens = t_p4med / tailfrac(T3, t_Tmed3);

    outtree->Fill();
    ndone++;
    if( ndone % 30 == 0 ) cout<<TString::Format("     done %5d pts (last: ipt %d, cls_data %.4f, cls_sens %.4f)", ndone, ipt, t_cls_data, t_cls_sens)<<endl;
  }

  outtree->Write();
  outfile->Close();
  cout<<endl<<" ---> job "<<ijob<<"/"<<njobs<<" (beam "<<beam<<") finished, "<<ndone<<" points"<<endl;
  return 0;
}
