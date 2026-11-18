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

// scan_3p1_exact: full 3+1 scan reproducing the PUBLISHED construction exactly,
// as implemented in Xiangpan's production CLs stage
// (winxji/Frequentist_CL_calculation, ea_frequentis_CLs.cc):
//
//   * per grid point, two pseudo-experiment ensembles at FIXED hypotheses give
//     the Dchi2 = chi2_4nu - chi2_3nu distributions (func_cal_CLs):
//         CL_sb = #{4nu pseudo : Dchi2_obs < Dchi2_pseudo}/N
//         CL_b  = #{3nu pseudo : Dchi2_obs < Dchi2_pseudo}/N
//         CLs   = CL_sb/CL_b          (exclusion at CLs < 0.05)
//     (fixed-hypothesis Dchi2 for the ensembles is confirmed by the derived
//      production code gen_chisquare_grid.cxx in jpmendez79/TOsc)
//   * DATA is profiled over the third parameter before entering the test
//     statistic (production: total_data_grids_dm2_t14_fit_t24.root carries
//     min_chi2 / min_sin2_theta_24_val per grid point)
//   * SENSITIVITY: a SEPARATE ensemble of 3nu "SM toys" is treated EXACTLY like
//     data -- each one profiled, each one given its own CLs value at every grid
//     point. Offline, each SM toy yields its own 95% contour and the median /
//     1sigma / 2sigma bands are QUANTILES ACROSS THAT ENSEMBLE OF CONTOURS
//     (func_CLs_toys + func_sensitivity + func_quantiles). This file therefore
//     stores the per-SM-toy CLs values; plot_3p1_exact.C does the contour and
//     quantile stage.
//
// NO approximations: no Gaussian CLs, no "CLs at the median test statistic".
//
// Usage: ./scan_3p1_exact -j <ijob> -n <njobs> [-o outdir] [-t ntoys] [-m nsm]
//                         [-k n_t24] [-gdm2 N] [-gamp N]

static inline double quadform_raw(const double *vv, const double *MM, int n)
{
  double sum=0;
  for(int i=0;i<n;i++){
    double vi = vv[i]; if(vi==0) continue;
    const double *row = MM + (Long64_t)i*n;
    double acc=0; for(int j=0;j<n;j++) acc += row[j]*vv[j];
    sum += vi*acc;
  }
  return sum;
}
// identical arithmetic to quadform_raw on (a-b); only avoids a TMatrixD temporary
static inline double quadform_diff(const double *a, const double *b, const double *MM, int n, double *buf)
{
  for(int j=0;j<n;j++) buf[j] = a[j]-b[j];
  return quadform_raw(buf, MM, n);
}
static double quadform(const TMatrixD &v, const TMatrixD &Minv)
{
  return quadform_raw(v.GetMatrixArray(), Minv.GetMatrixArray(), v.GetNcols());
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
  int ijob=0, njobs=1, ntoys=7000, nsm=500, n_t24=50, num_dm2=60, num_amp=60, onerow=-1, crn=0, fcmode=0, toyflip=0, toyclip=0;
  TString outdir = "output_scan_3p1_exact";

  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-j")==0 )    { stringstream c(argv[i+1]); c>>ijob; }
    if( strcmp(argv[i],"-n")==0 )    { stringstream c(argv[i+1]); c>>njobs; }
    if( strcmp(argv[i],"-t")==0 )    { stringstream c(argv[i+1]); c>>ntoys; }
    if( strcmp(argv[i],"-m")==0 )    { stringstream c(argv[i+1]); c>>nsm; }
    if( strcmp(argv[i],"-k")==0 )    { stringstream c(argv[i+1]); c>>n_t24; }
    if( strcmp(argv[i],"-gdm2")==0 ) { stringstream c(argv[i+1]); c>>num_dm2; }
    if( strcmp(argv[i],"-gamp")==0 ) { stringstream c(argv[i+1]); c>>num_amp; }
    if( strcmp(argv[i],"-o")==0 )    { outdir = argv[i+1]; }
    if( strcmp(argv[i],"-row")==0 )  { stringstream c(argv[i+1]); c>>onerow; }
    if( strcmp(argv[i],"-crn")==0 )  { stringstream c(argv[i+1]); c>>crn; }
    if( strcmp(argv[i],"-fc")==0 )   { stringstream c(argv[i+1]); c>>fcmode; }
    if( strcmp(argv[i],"-toyflip")==0 ){ stringstream c(argv[i+1]); c>>toyflip; }
    if( strcmp(argv[i],"-toyclip")==0 ){ stringstream c(argv[i+1]); c>>toyclip; }
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
  osc_test->flag_toy_flip_all = toyflip;
  osc_test->flag_toy_clip0 = toyclip;

  //////////////////////////////////////////// H0 = 3nu reference

  osc_test->Set_oscillation_pars(0, 0, 0, 0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();
  TMatrixD P3 = osc_test->matrix_tosc_eff_newworld_pred;
  TMatrixD S3 = osc_test->matrix_tosc_eff_newworld_abs_syst_total;
  TMatrixD C3inv = pearson_inv(P3, S3);

  //////////////////////////////////////////// ensemble A: 3nu pseudo for the p-value distribution

  osc_test->Set_toy_variations(ntoys);
  const int NB = P3.GetNcols();
  vector<double> buf(NB);
  const double *p3raw = P3.GetMatrixArray();
  const double *c3raw = C3inv.GetMatrixArray();

  vector<double> x3(( size_t)ntoys*NB);
  vector<double> chi2_3_on_x3(ntoys);
  for(int i=1;i<=ntoys;i++){
    const double *src = osc_test->map_matrix_tosc_toy_pred[i].GetMatrixArray();
    copy(src, src+NB, x3.begin()+(size_t)(i-1)*NB);
  }
  for(int i=0;i<ntoys;i++) chi2_3_on_x3[i] = quadform_diff(&x3[(size_t)i*NB], p3raw, c3raw, NB, &buf[0]);

  //////////////////////////////////////////// ensemble B: 3nu "SM toys" = pseudo-DATA for sensitivity
  //////////////////////////////////////////// (independent set; each treated exactly like data)

  osc_test->Set_toy_variations(nsm);
  vector<double> xsm((size_t)nsm*NB);
  vector<double> chi2_3_on_xsm(nsm);
  for(int i=1;i<=nsm;i++){
    const double *src = osc_test->map_matrix_tosc_toy_pred[i].GetMatrixArray();
    copy(src, src+NB, xsm.begin()+(size_t)(i-1)*NB);
  }
  for(int i=0;i<nsm;i++) chi2_3_on_xsm[i] = quadform_diff(&xsm[(size_t)i*NB], p3raw, c3raw, NB, &buf[0]);

  //////////////////////////////////////////// data

  osc_test->Set_meas2fitdata();
  TMatrixD D = osc_test->matrix_tosc_fitdata_newworld;
  double chi2_3_data = quadform_diff(D.GetMatrixArray(), p3raw, c3raw, NB, &buf[0]);
  cout<<endl<<" ---> chi2_3nu(data) = "<<chi2_3_data
      <<"   [ntoys "<<ntoys<<", nSMtoys "<<nsm<<", grid "<<NUM_dm2<<"x"<<NUM_amp<<", k "<<n_t24<<"]"<<endl;

  //////////////////////////////////////////// output

  system( TString::Format("mkdir -p %s", outdir.Data()) );
  TFile *outfile = new TFile(TString::Format("%s/scan_exact_%03d.root", outdir.Data(), ijob), "recreate");
  TTree *outtree = new TTree("tree", "exact published-construction 3+1 CLs grid");
  int    t_idm2, t_iamp;
  double t_dm2, t_amp, t_t14_d, t_t24_d, t_chi2d_min, t_Dchi2_data, t_cls_data;
  vector<double> *t_cls_sm = new vector<double>;      // consistent: SM toy at the grid-point hypothesis
  vector<double> *t_cls_sm_prof = new vector<double>; // diagnostic: SM toy profiled to its own best fit
  outtree->Branch("idm2",&t_idm2,"idm2/I");        outtree->Branch("iamp",&t_iamp,"iamp/I");
  outtree->Branch("dm2",&t_dm2,"dm2/D");           outtree->Branch("s22tue",&t_amp,"s22tue/D");
  outtree->Branch("t14_d",&t_t14_d,"t14_d/D");     outtree->Branch("t24_d",&t_t24_d,"t24_d/D");
  outtree->Branch("chi2d_min",&t_chi2d_min,"chi2d_min/D");
  outtree->Branch("Dchi2_data",&t_Dchi2_data,"Dchi2_data/D");
  outtree->Branch("cls_data",&t_cls_data,"cls_data/D");
  outtree->Branch("cls_sm", &t_cls_sm);
  outtree->Branch("cls_sm_prof", &t_cls_sm_prof);

  //////////////////////////////////////////// grid loop

  int ndone = 0;
  vector<double> chi2sm_min(nsm);
  for(int ipt=0; ipt<NUM_dm2*NUM_amp; ipt++) {
    if( ipt % njobs != ijob ) continue;
    t_idm2 = ipt / NUM_amp;  t_iamp = ipt % NUM_amp;
    if( onerow>=0 && t_idm2!=onerow ) continue;
    t_dm2 = pow(10, h1d_dm2->GetBinCenter(t_idm2+1));
    t_amp = pow(10, h1d_amp->GetBinCenter(t_iamp+1));

    /////////////// profile data AND every SM toy over the (t14,t24) family
    t_chi2d_min = 1e99;
    for(int m=0;m<nsm;m++) chi2sm_min[m] = 1e99;

    double lg_lo = log10(t_amp), lg_hi = 0;
    for(int k=0;k<n_t24;k++){
      double t24 = pow(10, lg_lo + (lg_hi-lg_lo)*( n_t24==1 ? 1.0 : double(k)/(n_t24-1) ));
      if( t24 < t_amp ) t24 = t_amp;
      double r = t_amp/t24; if( r>1 ) r = 1;
      double root = sqrt(1-r);
      for(int br=0; br<2; br++){
        if( br==1 && root==0 ) continue;
        double t14 = (br==0) ? (1-root)/2 : (1+root)/2;
        osc_test->Set_oscillation_pars(t_dm2, t14, t24, 0);
        osc_test->Apply_oscillation();
        osc_test->Set_apply_POT();
        TMatrixD P4 = osc_test->matrix_tosc_eff_newworld_pred;
        TMatrixD S4 = osc_test->matrix_tosc_eff_newworld_abs_syst_total;
        TMatrixD C4inv = pearson_inv(P4, S4);
        const double *p4raw = P4.GetMatrixArray();
        const double *c4raw = C4inv.GetMatrixArray();

        double c2d = quadform_diff(D.GetMatrixArray(), p4raw, c4raw, NB, &buf[0]);
        if( c2d < t_chi2d_min ){ t_chi2d_min = c2d; t_t14_d = t14; t_t24_d = t24; }

        for(int m=0;m<nsm;m++){
          double c = quadform_diff(&xsm[(size_t)m*NB], p4raw, c4raw, NB, &buf[0]);
          if( c < chi2sm_min[m] ) chi2sm_min[m] = c;
        }
      }
    }

    /////////////// p-value distributions at the grid point's 4nu hypothesis (FIXED, as in production)
    osc_test->Set_oscillation_pars(t_dm2, t_t14_d, t_t24_d, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();
    TMatrixD P4 = osc_test->matrix_tosc_eff_newworld_pred;
    TMatrixD S4 = osc_test->matrix_tosc_eff_newworld_abs_syst_total;
    TMatrixD C4inv = pearson_inv(P4, S4);

    // Common random numbers across grid points (variance reduction): the 4nu
    // ensemble at every grid point is drawn from the same random stream, so the
    // CL_sb surface varies smoothly in the oscillation parameters instead of
    // carrying independent Monte-Carlo noise at each point. Each point still gets
    // a valid 7000-sample ensemble from its own 4nu distribution, so the CLs
    // estimate at each point remains unbiased; only the point-to-point jitter
    // (which otherwise fragments the per-toy 95% contours and biases a
    // first-crossing limit low) is removed.
    if( crn ) osc_test->tosc_rand->SetSeed(4357);
    osc_test->Set_toy_variations(ntoys);// 4nu pseudo-experiments at this hypothesis

    const double *p4raw2 = P4.GetMatrixArray();
    const double *c4raw2 = C4inv.GetMatrixArray();
    vector<double> Dchi2_4vPseudo(ntoys), Dchi2_3vPseudo(ntoys);
    for(int i=0;i<ntoys;i++){
      const double *x4 = osc_test->map_matrix_tosc_toy_pred[i+1].GetMatrixArray();
      Dchi2_4vPseudo[i] = quadform_diff(x4, p4raw2, c4raw2, NB, &buf[0]) - quadform_diff(x4, p3raw, c3raw, NB, &buf[0]);
      Dchi2_3vPseudo[i] = quadform_diff(&x3[(size_t)i*NB], p4raw2, c4raw2, NB, &buf[0]) - chi2_3_on_x3[i];
    }
    sort(Dchi2_4vPseudo.begin(), Dchi2_4vPseudo.end());
    sort(Dchi2_3vPseudo.begin(), Dchi2_3vPseudo.end());

    // func_cal_CLs: CL_sb = #{Dchi2_obs < Dchi2_4vPseudo}/N, CL_b likewise on 3nu
    auto cls_of=[&](double dobs)->double{
      int n_sb = Dchi2_4vPseudo.end() - upper_bound(Dchi2_4vPseudo.begin(), Dchi2_4vPseudo.end(), dobs);
      int n_b  = Dchi2_3vPseudo.end() - upper_bound(Dchi2_3vPseudo.begin(), Dchi2_3vPseudo.end(), dobs);
      double CL_sb = n_sb*1./ntoys;
      double CL_b  = n_b *1./ntoys;
      if( CL_b==0 ) return 1.0;// production: CL_s stays 1
      return CL_sb/CL_b;
    };

    /////////////// data
    t_Dchi2_data = t_chi2d_min - chi2_3_data;
    t_cls_data   = cls_of( t_Dchi2_data );

    if( fcmode ){
      // PROFILED-EVERYWHERE variant (Feldman-Cousins-like): chi2_4nu of EVERY
      // tested dataset -- each 4nu pseudo, each 3nu pseudo, each SM toy, and the
      // data -- is min over the (t14,t24) family, exactly as the data is treated.
      // Tests whether production produced its toy chi2 vectors through the same
      // fit_t24 code path as the data grids.
      vector<double> min4(ntoys,1e99), min3(ntoys,1e99), minsm(nsm,1e99);
      vector<const double*> x4raw(ntoys);
      for(int i=0;i<ntoys;i++) x4raw[i] = osc_test->map_matrix_tosc_toy_pred[i+1].GetMatrixArray();
      for(int k=0;k<n_t24;k++){
        double t24 = pow(10, lg_lo + (lg_hi-lg_lo)*( n_t24==1 ? 1.0 : double(k)/(n_t24-1) ));
        if( t24 < t_amp ) t24 = t_amp;
        double r = t_amp/t24; if( r>1 ) r = 1;
        double root = sqrt(1-r);
        for(int br=0; br<2; br++){
          if( br==1 && root==0 ) continue;
          double t14 = (br==0) ? (1-root)/2 : (1+root)/2;
          osc_test->Set_oscillation_pars(t_dm2, t14, t24, 0);
          osc_test->Apply_oscillation();
          osc_test->Set_apply_POT();
          TMatrixD Ph = osc_test->matrix_tosc_eff_newworld_pred;
          TMatrixD Sh = osc_test->matrix_tosc_eff_newworld_abs_syst_total;
          TMatrixD Chinv = pearson_inv(Ph, Sh);
          const double *ph = Ph.GetMatrixArray();
          const double *chv = Chinv.GetMatrixArray();
          for(int i=0;i<ntoys;i++){
            double c4 = quadform_diff(x4raw[i], ph, chv, NB, &buf[0]); if(c4<min4[i]) min4[i]=c4;
            double c3 = quadform_diff(&x3[(size_t)i*NB], ph, chv, NB, &buf[0]); if(c3<min3[i]) min3[i]=c3;
          }
          for(int m=0;m<nsm;m++){
            double c = quadform_diff(&xsm[(size_t)m*NB], ph, chv, NB, &buf[0]); if(c<minsm[m]) minsm[m]=c;
          }
        }
      }
      for(int i=0;i<ntoys;i++){
        Dchi2_4vPseudo[i] = min4[i] - quadform_diff(x4raw[i], p3raw, c3raw, NB, &buf[0]);
        Dchi2_3vPseudo[i] = min3[i] - chi2_3_on_x3[i];
      }
      sort(Dchi2_4vPseudo.begin(), Dchi2_4vPseudo.end());
      sort(Dchi2_3vPseudo.begin(), Dchi2_3vPseudo.end());
      t_Dchi2_data = t_chi2d_min - chi2_3_data;   // data chi2_4 is already the profiled min
      t_cls_data   = cls_of( t_Dchi2_data );
      t_cls_sm->clear(); t_cls_sm_prof->clear();
      for(int m=0;m<nsm;m++){
        t_cls_sm->push_back( cls_of( minsm[m] - chi2_3_on_xsm[m] ) );      // profiled, consistent with ensembles
        t_cls_sm_prof->push_back( cls_of( chi2sm_min[m] - chi2_3_on_xsm[m] ) );// (identical here; kept for format)
      }
      outtree->Fill();
      ndone++;
      if( ndone % 5 == 0 ) cout<<TString::Format("     done %5d pts (FC mode, ipt %d: cls_data %.4f)", ndone, ipt, t_cls_data)<<endl;
      continue;
    }

    /////////////// every SM toy, evaluated at the SAME hypothesis as the ensembles
    // (the grid point's hypothesis, fixed from the data profile). Using each toy's
    // own profiled minimum instead would compare the observed statistic against
    // distributions built for a different hypothesis -- inconsistent, and it makes
    // the per-toy CLs surface noise-dominated.
    t_cls_sm->clear(); t_cls_sm->reserve(nsm);
    t_cls_sm_prof->clear(); t_cls_sm_prof->reserve(nsm);
    for(int m=0;m<nsm;m++){
      double chi2_at_hyp = quadform_diff(&xsm[(size_t)m*NB], p4raw2, c4raw2, NB, &buf[0]);
      t_cls_sm->push_back( cls_of( chi2_at_hyp - chi2_3_on_xsm[m] ) );
      t_cls_sm_prof->push_back( cls_of( chi2sm_min[m] - chi2_3_on_xsm[m] ) );
    }

    outtree->Fill();
    ndone++;
    if( ndone % 5 == 0 ) cout<<TString::Format("     done %5d pts (ipt %d: cls_data %.4f, t24hat %.3g)", ndone, ipt, t_cls_data, t_t24_d)<<endl;
  }

  outtree->Write();
  outfile->Close();
  cout<<endl<<" ---> job "<<ijob<<"/"<<njobs<<" finished, "<<ndone<<" points"<<endl;
  return 0;
}
