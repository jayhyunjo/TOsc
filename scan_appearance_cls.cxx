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
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

// scan_appearance_cls: nue appearance-only ("simple case") chi2 grid for
// Gaussian-CLs contours, following the published Wire-Cell 3+1 convention:
//   grid in (dm2_41, sin2_2theta_mue), log10 bin centers (Jesse-style TH1D grid)
//   H0 = 3nu (no oscillation), H1 = 4nu at the grid point
//   per point: chi2_4v on the 3nu-Asimov sample (-> Delta chi2_Asimov, since
//   chi2_3v(Asimov)=0 by closure) and chi2_4v on real data; chi2_3v(data) once.
// Gaussian CLs (Qian et al., arXiv:1407.5052):
//   sigma = 2*sqrt(dchi2_A); T_obs = chi2_3v(data) - chi2_4v(data)
//   CLs = (1 - Phi((T_obs - dchi2_A)/sigma)) / (1 - Phi((T_obs + dchi2_A)/sigma))
//   sensitivity contour: CLs(T_obs = -dchi2_A) = 0.05  <=>  dchi2_A = 3.84
//
// Usage: ./scan_appearance_cls -j <ijob> -n <njobs> [-o outdir]
//        grid is strided over jobs; each job writes outdir/scan_app_<ijob>.root

int main(int argc, char** argv)
{
  int ijob = 0, njobs = 1;
  TString outdir = "output_scan_app";

  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-j")==0 ) { stringstream c(argv[i+1]); c>>ijob; }
    if( strcmp(argv[i],"-n")==0 ) { stringstream c(argv[i+1]); c>>njobs; }
    if( strcmp(argv[i],"-o")==0 ) { outdir = argv[i+1]; }
  }

  gROOT->SetBatch(1);

  //////////////////////////////////////////// grid definition (log10 bin centers)

  const int    NUM_dm2 = 60;  const double DM2_LO = -2, DM2_HI = 2;   // dm2: 0.01 - 100 eV2
  const int    NUM_amp = 60;  const double AMP_LO = -4, AMP_HI = 0;   // sin2_2theta_mue: 1e-4 - 1

  TH1D *h1d_dm2 = new TH1D("h1d_dm2", "", NUM_dm2, DM2_LO, DM2_HI);
  TH1D *h1d_amp = new TH1D("h1d_amp", "", NUM_amp, AMP_LO, AMP_HI);

  //////////////////////////////////////////// TOsc initialization (as in read_oscillation_v01)

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

  //////////////////////////////////////////// appearance-only mode

  osc_test->flag_osc_appearance_only = true;

  //////////////////////////////////////////// H0 = 3nu: build Asimov, closure check

  double pars_3v[4] = {0, 0, 0, 0};

  osc_test->Set_oscillation_pars(0, 0, 0, 0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();
  osc_test->Set_asimov2fitdata();

  double chi2_3v_asimov = osc_test->FCN( pars_3v );
  cout<<endl<<" ---> closure chi2_3v(3nu Asimov) = "<<chi2_3v_asimov<<"  (expect 0)"<<endl;

  //////////////////////////////////////////// output

  system( TString::Format("mkdir -p %s", outdir.Data()) );
  TFile *outfile = new TFile(TString::Format("%s/scan_app_%03d.root", outdir.Data(), ijob), "recreate");
  TTree *outtree = new TTree("tree", "appearance-only chi2 grid");
  int    t_idm2, t_iamp;
  double t_dm2, t_amp, t_chi2_sens, t_chi2_4v_data, t_chi2_3v_data;
  outtree->Branch("idm2",         &t_idm2,         "idm2/I");
  outtree->Branch("iamp",         &t_iamp,         "iamp/I");
  outtree->Branch("dm2",          &t_dm2,          "dm2/D");
  outtree->Branch("s22tue",       &t_amp,          "s22tue/D");
  outtree->Branch("chi2_sens",    &t_chi2_sens,    "chi2_sens/D");     // chi2_4v on 3nu Asimov == dchi2_Asimov
  outtree->Branch("chi2_4v_data", &t_chi2_4v_data, "chi2_4v_data/D");
  outtree->Branch("chi2_3v_data", &t_chi2_3v_data, "chi2_3v_data/D");  // same value in every row

  //////////////////////////////////////////// pass 1: sensitivity (fitdata = 3nu Asimov)

  vector<double> v_sens(NUM_dm2*NUM_amp, -1);

  for(int ipt=0; ipt<NUM_dm2*NUM_amp; ipt++) {
    if( ipt % njobs != ijob ) continue;
    int idm2 = ipt / NUM_amp;
    int iamp = ipt % NUM_amp;
    double dm2 = pow(10, h1d_dm2->GetBinCenter(idm2+1));
    double amp = pow(10, h1d_amp->GetBinCenter(iamp+1));
    double p4[4] = {dm2, amp, 0, 0};
    v_sens[ipt] = osc_test->FCN( p4 );
    if( ipt % (10*njobs) == ijob ) cout<<TString::Format("     sens %4d / %d", ipt, NUM_dm2*NUM_amp)<<endl;
  }

  //////////////////////////////////////////// pass 2: data (fitdata = measurement)

  osc_test->Set_meas2fitdata();
  t_chi2_3v_data = osc_test->FCN( pars_3v );
  cout<<endl<<" ---> chi2_3v(data) = "<<t_chi2_3v_data<<endl;

  for(int ipt=0; ipt<NUM_dm2*NUM_amp; ipt++) {
    if( ipt % njobs != ijob ) continue;
    t_idm2 = ipt / NUM_amp;
    t_iamp = ipt % NUM_amp;
    t_dm2 = pow(10, h1d_dm2->GetBinCenter(t_idm2+1));
    t_amp = pow(10, h1d_amp->GetBinCenter(t_iamp+1));
    double p4[4] = {t_dm2, t_amp, 0, 0};
    t_chi2_4v_data = osc_test->FCN( p4 );
    t_chi2_sens    = v_sens[ipt];
    outtree->Fill();
    if( ipt % (10*njobs) == ijob ) cout<<TString::Format("     data %4d / %d", ipt, NUM_dm2*NUM_amp)<<endl;
  }

  outtree->Write();
  outfile->Close();

  cout<<endl<<" ---> job "<<ijob<<"/"<<njobs<<" finished"<<endl;
  return 0;
}
