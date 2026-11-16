// plot_scan_appearance_toys.C
// Extract 95% toy-CLs contours (sensitivity + data) from scan_appearance_cls_toys
// output and overlay the corresponding published curves.
//
// Usage:
//   root -l -b -q 'plot_scan_appearance_toys.C("output_scan_app_toys_bnb", 1)'
//   beam: 0 = BNB+NuMI, 1 = BNB only, 2 = NuMI only

void plot_scan_appearance_toys(TString indir = "output_scan_app_toys", int beam = 0, int NUM_dm2 = 60, int NUM_amp = 60, TString tagsuffix = "")
{
  const double DM2_LO = -2, DM2_HI = 2, AMP_LO = -4, AMP_HI = 0;
  TH1D *h1d_dm2 = new TH1D("h1d_dm2","",NUM_dm2,DM2_LO,DM2_HI);
  TH1D *h1d_amp = new TH1D("h1d_amp","",NUM_amp,AMP_LO,AMP_HI);

  ////////////////////////////////////// read grid

  TChain ch("tree");
  ch.Add(indir + "/*.root");
  int idm2, iamp; double cls_data, cls_sens;
  ch.SetBranchAddress("idm2",&idm2); ch.SetBranchAddress("iamp",&iamp);
  ch.SetBranchAddress("cls_data",&cls_data); ch.SetBranchAddress("cls_sens",&cls_sens);

  std::vector<std::vector<double>> CD(NUM_dm2, std::vector<double>(NUM_amp,-1));
  std::vector<std::vector<double>> CS(NUM_dm2, std::vector<double>(NUM_amp,-1));
  int nent = ch.GetEntries();
  printf(" grid entries: %d (expect %d)\n", nent, NUM_dm2*NUM_amp);
  for(int i=0;i<nent;i++){ ch.GetEntry(i); CD[idm2][iamp]=cls_data; CS[idm2][iamp]=cls_sens; }

  ////////////////////////////////////// extract CLs = 0.05 crossings per dm2 row

  const double CUT = 0.05;
  auto extract=[&](std::vector<std::vector<double>> &M)->TGraph*{
    TGraph *g = new TGraph();
    for(int i=0;i<NUM_dm2;i++){
      double y = pow(10, h1d_dm2->GetBinCenter(i+1));
      for(int j=1;j<NUM_amp;j++){
        double c0=M[i][j-1], c1=M[i][j];
        if( c0<0 || c1<0 ) continue;
        if( c0>CUT && c1<=CUT ){
          double x0=h1d_amp->GetBinCenter(j), x1=h1d_amp->GetBinCenter(j+1);
          double x = x0 + (x1-x0)*(CUT-c0)/(c1-c0);
          g->SetPoint(g->GetN(), pow(10,x), y);
          break;
        }
      }
    }
    return g;
  };
  TGraph *g_data = extract(CD);
  TGraph *g_sens = extract(CS);
  printf(" contour points: data %d, sens %d\n", g_data->GetN(), g_sens->GetN());

  ////////////////////////////////////// published references per beam

  TString refdir = "/home/jjo/Work/MicroBooNE/natureReproducer/data/";
  TGraph *ref_data = 0, *ref_sens = 0;
  TString ref_label, tag;
  if( beam==1 ){
    TFile *f = new TFile(refdir+"BNB_Only_App_Only_data_and_sens.root");
    ref_data = (TGraph*)f->Get("gh_appONLY_data");
    ref_sens = (TGraph*)f->Get("gh_appONLY_sensitivity");
    ref_label = "Published app-only (BNB)"; tag = "bnbonly"+tagsuffix;
  } else if( beam==2 ){
    TFile *f = new TFile(refdir+"NuMI_only_app_data.root");
    ref_data = (TGraph*)f->Get("NuMIonly_data");
    ref_label = "Published (NuMI)"; tag = "numionly"+tagsuffix;
  } else {
    TFile *f = new TFile(refdir+"BNBNuMI_3p1_data_median_brazil_asimov_Apeparance.root");
    ref_data = (TGraph*)f->Get("gh_meas");
    ref_sens = (TGraph*)f->Get("gh_asimov");
    ref_label = "Published 3+1 profiled (BNB+NuMI)"; tag = "bnbnumi"+tagsuffix;
  }

  ////////////////////////////////////// draw

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1",800,700);
  c1->SetLogx(); c1->SetLogy();
  c1->SetLeftMargin(0.13); c1->SetBottomMargin(0.12);
  TH2D *frame = new TH2D("frame",";sin^{2}2#theta_{#mue};#Deltam^{2}_{41} (eV^{2})",
                         10, 1e-4, 1, 10, 1e-2, 1e2);
  frame->Draw();

  if(ref_sens){ ref_sens->SetLineColor(kCyan+1); ref_sens->SetLineWidth(3); ref_sens->SetLineStyle(2); ref_sens->Draw("l same"); }
  if(ref_data){ ref_data->SetLineColor(kGray+2); ref_data->SetLineWidth(3); ref_data->SetLineStyle(2); ref_data->Draw("l same"); }
  g_sens->SetLineColor(kAzure+2); g_sens->SetLineWidth(2); g_sens->Draw("l same");
  g_data->SetLineColor(kRed);     g_data->SetLineWidth(2); g_data->Draw("l same");

  TLegend *leg = new TLegend(0.15,0.15,0.68,0.36);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(g_sens, "This work: sens (toy CLs 95%)","l");
  leg->AddEntry(g_data, "This work: data (toy CLs 95%)","l");
  if(ref_sens) leg->AddEntry(ref_sens, ref_label+" sens","l");
  if(ref_data) leg->AddEntry(ref_data, ref_label+" data","l");
  leg->Draw();

  c1->SaveAs("contour_app_toys_"+tag+".png");
  c1->SaveAs("contour_app_toys_"+tag+".pdf");

  TFile *fout = new TFile("contour_app_toys_"+tag+".root","recreate");
  g_sens->Write("g_sens_thiswork");
  g_data->Write("g_data_thiswork");
  if(ref_sens) ref_sens->Write("ref_sens");
  if(ref_data) ref_data->Write("ref_data");
  fout->Close();

  ////////////////////////////////////// numeric comparison

  auto evalx=[&](TGraph*g, double y)->double{
    if(!g) return -1;
    for(int k=0;k<g->GetN()-1;k++){
      double x0,y0,x1,y1; g->GetPoint(k,x0,y0); g->GetPoint(k+1,x1,y1);
      if((y0-y)*(y1-y)<=0 && y0!=y1) return exp(log(x0)+(log(x1)-log(x0))*(log(y/y0))/(log(y1/y0)));
    }
    return -1;
  };
  printf("\n   dm2 (eV2)  |  this-sens   ref-sens  ratio |  this-data   ref-data  ratio\n");
  for(double ydm : {0.1, 0.3, 1.0, 3.0, 10.0, 50.0}) {
    double ts=evalx(g_sens,ydm), rs=evalx(ref_sens,ydm), td=evalx(g_data,ydm), rd=evalx(ref_data,ydm);
    printf("   %8.2f   |  %9.2e  %9.2e  %5.2f |  %9.2e  %9.2e  %5.2f\n",
           ydm, ts, rs, (rs>0&&ts>0)?ts/rs:-1, td, rd, (rd>0&&td>0)?td/rd:-1);
  }
}
