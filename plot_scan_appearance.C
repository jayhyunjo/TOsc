// plot_scan_appearance.C
// Build 95% Gaussian-CLs appearance-only contours (sensitivity + data) from the
// scan_appearance_cls output grid and overlay the published Nature curves.
//
//   T = chi2_3v(x) - chi2_4v(x)   (more 4nu-like => larger T)
//   under H_3nu: T ~ N(-D, 4D);  under H_4nu: T ~ N(+D, 4D);  D = chi2_sens (Asimov dchi2)
//   CLs = Phi((T_obs - D)/(2 sqrt(D))) / Phi((T_obs + D)/(2 sqrt(D)))   [arXiv:1407.5052]
//   sensitivity contour: CLs(T_obs=-D)=0.05  <=>  D = 3.84
//
// Usage: root -l -b -q plot_scan_appearance.C

void plot_scan_appearance()
{
  const int NUM_dm2 = 60, NUM_amp = 60;
  const double DM2_LO = -2, DM2_HI = 2, AMP_LO = -4, AMP_HI = 0;

  TH1D *h1d_dm2 = new TH1D("h1d_dm2","",NUM_dm2,DM2_LO,DM2_HI);
  TH1D *h1d_amp = new TH1D("h1d_amp","",NUM_amp,AMP_LO,AMP_HI);

  ////////////////////////////////////// read grid

  TChain ch("tree");
  ch.Add("output_scan_app/scan_app_*.root");
  int idm2, iamp; double dm2, s22tue, chi2_sens, chi2_4v_data, chi2_3v_data;
  ch.SetBranchAddress("idm2",&idm2); ch.SetBranchAddress("iamp",&iamp);
  ch.SetBranchAddress("dm2",&dm2);   ch.SetBranchAddress("s22tue",&s22tue);
  ch.SetBranchAddress("chi2_sens",&chi2_sens);
  ch.SetBranchAddress("chi2_4v_data",&chi2_4v_data);
  ch.SetBranchAddress("chi2_3v_data",&chi2_3v_data);

  std::vector<std::vector<double>> D(NUM_dm2, std::vector<double>(NUM_amp, -1));
  std::vector<std::vector<double>> CLS(NUM_dm2, std::vector<double>(NUM_amp, 1));
  int nent = ch.GetEntries();
  printf(" grid entries: %d (expect %d)\n", nent, NUM_dm2*NUM_amp);

  for(int i=0;i<nent;i++){
    ch.GetEntry(i);
    D[idm2][iamp] = chi2_sens;
    double d = chi2_sens;
    if( d < 1e-9 ) { CLS[idm2][iamp] = 1; continue; }
    double T_obs = chi2_3v_data - chi2_4v_data;
    double s = 2*sqrt(d);
    double cls_num = ROOT::Math::normal_cdf( (T_obs - d)/s );
    double cls_den = ROOT::Math::normal_cdf( (T_obs + d)/s );
    CLS[idm2][iamp] = (cls_den>0)? cls_num/cls_den : 1;
  }

  ////////////////////////////////////// extract 95% contours: first crossing in amp per dm2 row

  TGraph *g_sens = new TGraph(); TGraph *g_data = new TGraph();
  const double D_SENS = 3.84;   // CLs_sens = 2*Phi(-sqrt(D)) = 0.05
  const double CLS_CUT = 0.05;

  for(int i=0;i<NUM_dm2;i++){
    double y = pow(10, h1d_dm2->GetBinCenter(i+1));
    // sensitivity: D crosses 3.84
    for(int j=1;j<NUM_amp;j++){
      double d0=D[i][j-1], d1=D[i][j];
      if( d0>0 && d1>0 && d0<D_SENS && d1>=D_SENS ){
        double x0=h1d_amp->GetBinCenter(j), x1=h1d_amp->GetBinCenter(j+1);
        double x = x0 + (x1-x0)*(D_SENS-d0)/(d1-d0);
        g_sens->SetPoint(g_sens->GetN(), pow(10,x), y);
        break;
      }
    }
    // data: CLs crosses 0.05 downward
    for(int j=1;j<NUM_amp;j++){
      double c0=CLS[i][j-1], c1=CLS[i][j];
      if( c0>CLS_CUT && c1<=CLS_CUT ){
        double x0=h1d_amp->GetBinCenter(j), x1=h1d_amp->GetBinCenter(j+1);
        double x = x0 + (x1-x0)*(CLS_CUT-c0)/(c1-c0);
        g_data->SetPoint(g_data->GetN(), pow(10,x), y);
        break;
      }
    }
  }
  printf(" sens contour points: %d, data contour points: %d\n", g_sens->GetN(), g_data->GetN());

  ////////////////////////////////////// published curves

  TFile *fref = new TFile("/home/jjo/Work/MicroBooNE/natureReproducer/data/BNBNuMI_3p1_data_median_brazil_asimov_Apeparance.root","read");
  TGraph *gh_meas   = (TGraph*)fref->Get("gh_meas");
  TGraph *gh_asimov = (TGraph*)fref->Get("gh_asimov");

  TFile *flsnd = new TFile("/home/jjo/Work/MicroBooNE/natureReproducer/data/BNBNuMI_3p1_App_data_and_lsnd_and_mb.root","read");
  TGraph *lsnd1 = (TGraph*)flsnd->Get("lsnd1");
  TGraph *lsnd2 = (TGraph*)flsnd->Get("lsnd2");

  ////////////////////////////////////// draw

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1",800,700);
  c1->SetLogx(); c1->SetLogy();
  c1->SetLeftMargin(0.13); c1->SetBottomMargin(0.12);

  TH2D *frame = new TH2D("frame",";sin^{2}2#theta_{#mue};#Deltam^{2}_{41} (eV^{2})",
                         10, 1e-4, 1, 10, 1e-2, 1e2);
  frame->Draw();

  if(lsnd1){ lsnd1->SetFillColorAlpha(kBlue-7,0.5); lsnd1->Draw("f same"); }
  if(lsnd2){ lsnd2->SetFillColorAlpha(kBlue-9,0.5); lsnd2->Draw("f same"); }

  if(gh_asimov){ gh_asimov->SetLineColor(kCyan+1);  gh_asimov->SetLineWidth(3); gh_asimov->SetLineStyle(2); gh_asimov->Draw("l same"); }
  if(gh_meas)  { gh_meas->SetLineColor(kGray+2);    gh_meas->SetLineWidth(3);   gh_meas->SetLineStyle(2);   gh_meas->Draw("l same"); }

  g_sens->SetLineColor(kAzure+2); g_sens->SetLineWidth(2); g_sens->Draw("l same");
  g_data->SetLineColor(kRed);     g_data->SetLineWidth(2); g_data->Draw("l same");

  TLegend *leg = new TLegend(0.15,0.15,0.62,0.38);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(g_sens,     "This work: sens (Gaussian CLs 95%)","l");
  leg->AddEntry(g_data,     "This work: data (Gaussian CLs 95%)","l");
  if(gh_asimov) leg->AddEntry(gh_asimov,"Published: Asimov sens","l");
  if(gh_meas)   leg->AddEntry(gh_meas,  "Published: data","l");
  if(lsnd1)     leg->AddEntry(lsnd1,    "LSND allowed","f");
  leg->Draw();

  c1->SaveAs("contour_appearance_cls.png");
  c1->SaveAs("contour_appearance_cls.pdf");

  TFile *fout = new TFile("contour_appearance_cls.root","recreate");
  g_sens->Write("g_sens_thiswork");
  g_data->Write("g_data_thiswork");
  if(gh_asimov) gh_asimov->Write("gh_asimov_published");
  if(gh_meas)   gh_meas->Write("gh_meas_published");
  fout->Close();

  // quick numeric comparison at a few dm2 values
  printf("\n   dm2 (eV2)  |  this-sens  published-sens |  this-data  published-data\n");
  for(double ydm : {0.1, 0.3, 1.0, 3.0, 10.0, 50.0}) {
    auto evalx=[&](TGraph*g, double y)->double{
      if(!g) return -1;
      for(int k=0;k<g->GetN()-1;k++){
        double x0,y0,x1,y1; g->GetPoint(k,x0,y0); g->GetPoint(k+1,x1,y1);
        if((y0-y)*(y1-y)<=0 && y0!=y1) return exp(log(x0)+(log(x1)-log(x0))*(log(y/y0))/(log(y1/y0)));
      }
      return -1;
    };
    printf("   %8.2f   |  %9.2e  %9.2e      |  %9.2e  %9.2e\n", ydm,
           evalx(g_sens,ydm), evalx(gh_asimov,ydm), evalx(g_data,ydm), evalx(gh_meas,ydm));
  }
}
