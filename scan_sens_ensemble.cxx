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

// scan_sens_ensemble: TEST of the sensitivity construction hypothesis.
// Median expected sensitivity built as ENSEMBLE of 3nu pseudo-experiments,
// EACH profiled over the (t14,t24) family like real data (vs the shortcut of
// fixed Asimov-profiled H1 + median toy T). Gaussian CLs per pseudo-experiment
// (arXiv:1407.5052) -- fast; if this moves the high-dm2 sens onto the
// published gh_median_val, the construction is confirmed and a toy-grade
// version can follow.
//
// Usage: ./scan_sens_ensemble -j <ijob> -n <njobs> -irows "49,61,74,92"
//        [-e n_pseudo] [-k n_t24] [-gdm2 100] [-gamp 100] [-o outdir]

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
  int ijob=0, njobs=1, n_pseudo=500, n_t24=50, num_dm2=100, num_amp=100;
  TString outdir = "output_sens_ensemble";
  TString rowlist = "";

  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-j")==0 ) { stringstream c(argv[i+1]); c>>ijob; }
    if( strcmp(argv[i],"-n")==0 ) { stringstream c(argv[i+1]); c>>njobs; }
    if( strcmp(argv[i],"-e")==0 ) { stringstream c(argv[i+1]); c>>n_pseudo; }
    if( strcmp(argv[i],"-k")==0 ) { stringstream c(argv[i+1]); c>>n_t24; }
    if( strcmp(argv[i],"-gdm2")==0 ) { stringstream c(argv[i+1]); c>>num_dm2; }
    if( strcmp(argv[i],"-gamp")==0 ) { stringstream c(argv[i+1]); c>>num_amp; }
    if( strcmp(argv[i],"-o")==0 ) { outdir = argv[i+1]; }
    if( strcmp(argv[i],"-irows")==0 ) { rowlist = argv[i+1]; }
  }

  vector<int> rows;
  { stringstream ss(rowlist.Data()); string tok;
    while( getline(ss, tok, ',') ) if(!tok.empty()) rows.push_back(atoi(tok.c_str()));
  }
  if( rows.empty() ){ cerr<<" ERROR: give -irows"<<endl; return 1; }

  gROOT->SetBatch(1);
  TH1D *h1d_dm2 = new TH1D("h1d_dm2","",num_dm2,-2,2);
  TH1D *h1d_amp = new TH1D("h1d_amp","",num_amp,-4,0);

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
  osc_test->flag_osc_appearance_only = false;

  //////// 3nu: pred, cov, pseudo-experiment ensemble
  osc_test->Set_oscillation_pars(0,0,0,0);
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();
  TMatrixD P3 = osc_test->matrix_tosc_eff_newworld_pred;
  TMatrixD S3 = osc_test->matrix_tosc_eff_newworld_abs_syst_total;
  TMatrixD C3inv = pearson_inv(P3, S3);

  osc_test->Set_toy_variations(n_pseudo);
  vector<TMatrixD> xe; xe.reserve(n_pseudo);
  for(int i=1;i<=n_pseudo;i++) xe.push_back( osc_test->map_matrix_tosc_toy_pred[i] );
  vector<double> chi2_3_e(n_pseudo);
  for(int e=0;e<n_pseudo;e++) chi2_3_e[e] = quadform(xe[e]-P3, C3inv);

  //////// output
  system( TString::Format("mkdir -p %s", outdir.Data()) );
  TFile *outfile = new TFile(TString::Format("%s/sens_ens_%03d.root", outdir.Data(), ijob), "recreate");
  TTree *outtree = new TTree("tree","ensemble sensitivity test");
  int t_idm2, t_iamp; double t_dm2, t_amp, t_cls_ens_med, t_cls_old_gauss;
  outtree->Branch("idm2",&t_idm2,"idm2/I");  outtree->Branch("iamp",&t_iamp,"iamp/I");
  outtree->Branch("dm2",&t_dm2,"dm2/D");     outtree->Branch("s22tue",&t_amp,"s22tue/D");
  outtree->Branch("cls_ens_med",&t_cls_ens_med,"cls_ens_med/D");   // NEW construction
  outtree->Branch("cls_old_gauss",&t_cls_old_gauss,"cls_old_gauss/D"); // old shortcut (gaussian)

  //////// loop over selected rows
  int npts_total = (int)rows.size()*num_amp, idx_lin = -1, ndone=0;
  for(int ir=0; ir<(int)rows.size(); ir++){
    t_idm2 = rows[ir];
    t_dm2  = pow(10, h1d_dm2->GetBinCenter(t_idm2+1));
    for(int ja=0; ja<num_amp; ja++){
      idx_lin++;
      if( idx_lin % njobs != ijob ) continue;
      t_iamp = ja;
      t_amp  = pow(10, h1d_amp->GetBinCenter(ja+1));

      vector<double> best(n_pseudo, 1e99), bestD(n_pseudo, 0);
      double bestA = 1e99, bestA_D = 0;

      double lg_lo = log10(t_amp), lg_hi = 0;
      for(int k=0;k<n_t24;k++){
        double t24 = pow(10, lg_lo + (lg_hi-lg_lo)*( n_t24==1 ? 1.0 : double(k)/(n_t24-1) ));
        if( t24 < t_amp ) t24 = t_amp;
        double r = t_amp/t24; if(r>1) r=1;
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
          double DA = quadform(P3 - P4, C4inv);
          if( DA < bestA ){ bestA = DA; bestA_D = DA; }
          for(int e=0;e<n_pseudo;e++){
            double c = quadform(xe[e] - P4, C4inv);
            if( c < best[e] ){ best[e] = c; bestD[e] = DA; }
          }
        }
      }

      // per pseudo-experiment: Gaussian CLs with per-e profiled hypothesis
      vector<double> cls_e(n_pseudo);
      for(int e=0;e<n_pseudo;e++){
        double T  = chi2_3_e[e] - best[e];
        double Dv = bestD[e] > 1e-9 ? bestD[e] : 1e-9;
        double sg = 2*sqrt(Dv);
        cls_e[e]  = ROOT::Math::normal_cdf((T - Dv)/sg) / ROOT::Math::normal_cdf((T + Dv)/sg);
      }
      sort(cls_e.begin(), cls_e.end());
      t_cls_ens_med = 0.5*( cls_e[n_pseudo/2] + cls_e[(n_pseudo-1)/2] );

      // old shortcut, gaussian: CLs at median under 3nu with FIXED asimov-profiled D
      { double Dv = bestA_D > 1e-9 ? bestA_D : 1e-9;
        t_cls_old_gauss = 2*ROOT::Math::normal_cdf(-sqrt(Dv)); }

      outtree->Fill();
      ndone++;
      if( ndone % 10 == 0 ) cout<<TString::Format("     done %4d (row %d amp %d: ens %.4f old %.4f)", ndone, t_idm2, ja, t_cls_ens_med, t_cls_old_gauss)<<endl;
    }
  }

  outtree->Write(); outfile->Close();
  cout<<" ---> job "<<ijob<<" finished, "<<ndone<<" points of "<<npts_total<<endl;
  return 0;
}
