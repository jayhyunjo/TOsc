// plot_3p1_exact.C
// Contour + sensitivity stage, replicating Xiangpan's production CLs code
// (winxji/Frequentist_CL_calculation, ea_frequentis_CLs.cc) step for step:
//
//   func_data_CLs    : fill the (theta, dm2) map with 1 - CLs for DATA, take the
//                      0.95 contour with ROOT "cont z list"      -> data exclusion
//   func_CLs_toys    : the same, once per 3nu SM toy             -> one contour per toy
//   func_sensitivity : per dm2, quantiles ACROSS the toy contours
//                      q = 0.5           -> median sensitivity
//                      q = 0.5 +- 0.6827/2 -> 1 sigma band
//                      q = 0.5 +- 0.9545/2 -> 2 sigma band
//   func_quantiles   : "exact value" branch (sorted list, first idx with idx+1 >= N*q)
//
// Usage: root -l -b -q 'plot_3p1_exact.C("output_scan_3p1_exact", 60, 60)'

#include <vector>
#include <algorithm>

static TGraph* get_contour_95(TH2D *h2_CL, int index)
{
  // production func_get_contours: contours at 0.9 / 0.95 / 0.9973, take the 0.95 one
  gROOT->SetBatch(1);
  const int Ncontour = 3;
  double contours[Ncontour] = {0.9, 0.95, 0.9973};

  TString roostr = TString::Format("canv_h2_CL_%d", index);
  TCanvas *canv = new TCanvas(roostr, roostr, 800, 600);
  h2_CL->SetStats(0);
  h2_CL->SetContour(Ncontour, contours);
  h2_CL->Draw("cont z list");
  canv->Update();

  TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  TGraph *out = new TGraph();
  out->SetName(TString::Format("gh_95_%d", index));
  if( conts && conts->GetSize()>1 ){
    TList *lvl = (TList*)conts->At(1);// index 1 -> 0.95
    if( lvl && lvl->GetSize()>0 ){
      TGraph *g = (TGraph*)lvl->First();
      for(int i=0;i<g->GetN();i++){
        double x,y; g->GetPoint(i,x,y);
        out->SetPoint(out->GetN(), pow(10.,x), pow(10.,y));
      }
    }
  }
  delete canv;
  return out;
}

// production func_quantiles, "exact value" branch
static double quantile_of(std::vector<TGraph*> &gT, double qq, double dm2_val)
{
  std::vector<double> vc;
  for(size_t i=0;i<gT.size();i++){
    int nps = gT[i]->GetN();
    if( nps<2 ) continue;
    double d1,t1,d2,t2;
    gT[i]->GetPoint(0,d1,t1);
    gT[i]->GetPoint(nps-1,d2,t2);
    bool in = (dm2_val>=d1 && dm2_val<=d2) || (dm2_val>=d2 && dm2_val<=d1);
    if( !in ) continue;
    vc.push_back( gT[i]->Eval(dm2_val) );
  }
  if( vc.empty() ) return -1;
  std::sort(vc.begin(), vc.end());
  int n = vc.size();
  for(int i=0;i<n;i++) if( i+1 >= n*qq ) return vc[i];
  return vc[n-1];
}

void plot_3p1_exact(TString indir = "output_scan_3p1_exact", int NUM_dm2 = 60, int NUM_amp = 60)
{
  const double DM2_LO = -2, DM2_HI = 2, AMP_LO = -4, AMP_HI = 0;

  ////////////////////////////////////// read the grid

  TChain ch("tree");
  ch.Add(indir + "/*.root");
  int idm2, iamp; double cls_data;
  std::vector<double> *cls_sm = 0;
  ch.SetBranchAddress("idm2",&idm2); ch.SetBranchAddress("iamp",&iamp);
  ch.SetBranchAddress("cls_data",&cls_data);
  ch.SetBranchAddress("cls_sm",&cls_sm);

  int nent = ch.GetEntries();
  ch.GetEntry(0);
  int NSM = cls_sm->size();
  printf(" grid entries: %d (expect %d),  SM toys: %d\n", nent, NUM_dm2*NUM_amp, NSM);

  TH2D *h2_data = new TH2D("h2_data","",NUM_amp,AMP_LO,AMP_HI,NUM_dm2,DM2_LO,DM2_HI);
  std::vector<TH2D*> h2_sm(NSM);
  for(int m=0;m<NSM;m++) h2_sm[m] = new TH2D(TString::Format("h2_sm_%d",m),"",NUM_amp,AMP_LO,AMP_HI,NUM_dm2,DM2_LO,DM2_HI);

  for(int i=0;i<nent;i++){
    ch.GetEntry(i);
    h2_data->SetBinContent(iamp+1, idm2+1, 1 - cls_data);// production stores 1 - CLs
    for(int m=0;m<NSM;m++) h2_sm[m]->SetBinContent(iamp+1, idm2+1, 1 - cls_sm->at(m));
  }

  ////////////////////////////////////// data contour

  TGraph *g_data = get_contour_95(h2_data, 999999);
  g_data->SetName("g_data_thiswork");
  printf(" data contour points: %d\n", g_data->GetN());

  ////////////////////////////////////// one contour per SM toy, then quantiles

  std::vector<TGraph*> gT;// transposed (dm2, theta), as in production
  for(int m=0;m<NSM;m++){
    TGraph *g = get_contour_95(h2_sm[m], m);
    TGraph *t = new TGraph();
    for(int i=0;i<g->GetN();i++){ double x,y; g->GetPoint(i,x,y); t->SetPoint(i,y,x); }
    if( t->GetN()>1 ) gT.push_back(t);
    if( (m+1)%100==0 ) printf("   ... %d/%d toy contours\n", m+1, NSM);
  }
  printf(" usable toy contours: %d\n", (int)gT.size());

  const double val_1sigma = 0.6827, val_2sigma = 0.9545;
  TGraph *g_med = new TGraph(); g_med->SetName("g_sens_median");
  TGraph *g_1L  = new TGraph(); TGraph *g_1R = new TGraph();
  TGraph *g_2L  = new TGraph(); TGraph *g_2R = new TGraph();

  for(int i=1;i<=10000;i++){
    double dm2_val = i*0.05;// production dm2_step
    if( dm2_val < 0.01 ) continue;
    if( dm2_val > 94 ) break;
    double v;
    v = quantile_of(gT, 0.5, dm2_val);                 if(v>0) g_med->SetPoint(g_med->GetN(), v, dm2_val);
    v = quantile_of(gT, 0.5-val_1sigma/2, dm2_val);    if(v>0) g_1L->SetPoint(g_1L->GetN(), v, dm2_val);
    v = quantile_of(gT, 0.5+val_1sigma/2, dm2_val);    if(v>0) g_1R->SetPoint(g_1R->GetN(), v, dm2_val);
    v = quantile_of(gT, 0.5-val_2sigma/2, dm2_val);    if(v>0) g_2L->SetPoint(g_2L->GetN(), v, dm2_val);
    v = quantile_of(gT, 0.5+val_2sigma/2, dm2_val);    if(v>0) g_2R->SetPoint(g_2R->GetN(), v, dm2_val);
  }
  printf(" median sens points: %d\n", g_med->GetN());

  // bands as closed polygons
  TGraph *band1 = new TGraph(); TGraph *band2 = new TGraph();
  for(int i=0;i<g_1L->GetN();i++){ double x,y; g_1L->GetPoint(i,x,y); band1->SetPoint(band1->GetN(),x,y); }
  for(int i=g_1R->GetN()-1;i>=0;i--){ double x,y; g_1R->GetPoint(i,x,y); band1->SetPoint(band1->GetN(),x,y); }
  for(int i=0;i<g_2L->GetN();i++){ double x,y; g_2L->GetPoint(i,x,y); band2->SetPoint(band2->GetN(),x,y); }
  for(int i=g_2R->GetN()-1;i>=0;i--){ double x,y; g_2R->GetPoint(i,x,y); band2->SetPoint(band2->GetN(),x,y); }

  ////////////////////////////////////// published references

  TString refdir = "/home/jjo/Work/MicroBooNE/natureReproducer/data/";
  TFile *fref = new TFile(refdir+"BNBNuMI_3p1_data_median_brazil_asimov_Apeparance.root");
  TGraph *ref_data = (TGraph*)fref->Get("gh_meas");
  TGraph *ref_med  = (TGraph*)fref->Get("gh_median_val");

  ////////////////////////////////////// draw

  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1",850,720);
  c1->SetLogx(); c1->SetLogy();
  c1->SetLeftMargin(0.14); c1->SetBottomMargin(0.13);
  TH2D *frame = new TH2D("frame",";sin^{2}2#theta_{#mue};#Deltam^{2}_{41} (eV^{2})",10,1e-4,1,10,1e-2,1e2);
  frame->Draw();

  band2->SetFillColor(kYellow);      band2->SetLineWidth(0); band2->Draw("f same");
  band1->SetFillColor(kGreen);       band1->SetLineWidth(0); band1->Draw("f same");
  if(ref_med) { ref_med->SetLineColor(kCyan+1); ref_med->SetLineWidth(3); ref_med->SetLineStyle(2); ref_med->Draw("l same"); }
  if(ref_data){ ref_data->SetLineColor(kGray+2); ref_data->SetLineWidth(3); ref_data->SetLineStyle(2); ref_data->Draw("l same"); }
  g_med->SetLineColor(kAzure+2); g_med->SetLineWidth(2); g_med->Draw("l same");
  g_data->SetLineColor(kRed);    g_data->SetLineWidth(2); g_data->Draw("l same");

  TLegend *leg = new TLegend(0.16,0.15,0.70,0.37);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(g_data,"This work: data (toy CLs 95%)","l");
  leg->AddEntry(g_med, "This work: median sens (toy ensemble)","l");
  leg->AddEntry(band1, "This work: #pm1#sigma / #pm2#sigma","f");
  if(ref_data) leg->AddEntry(ref_data,"Published: data","l");
  if(ref_med)  leg->AddEntry(ref_med, "Published: median sens","l");
  leg->Draw();

  c1->SaveAs("contour_3p1_exact.png");
  c1->SaveAs("contour_3p1_exact.pdf");

  TFile *fout = new TFile("contour_3p1_exact.root","recreate");
  g_data->Write("g_data_thiswork");
  g_med->Write("g_sens_median_thiswork");
  band1->Write("g_band_1sigma"); band2->Write("g_band_2sigma");
  g_1L->Write("g_1s_LL"); g_1R->Write("g_1s_RR");
  g_2L->Write("g_2s_LL"); g_2R->Write("g_2s_RR");
  if(ref_data) ref_data->Write("ref_data");
  if(ref_med)  ref_med->Write("ref_median");
  fout->Close();

  ////////////////////////////////////// numeric comparison

  auto evalx=[&](TGraph*g,double y)->double{
    if(!g) return -1;
    for(int k=0;k<g->GetN()-1;k++){
      double x0,y0,x1,y1; g->GetPoint(k,x0,y0); g->GetPoint(k+1,x1,y1);
      if((y0-y)*(y1-y)<=0 && y0!=y1) return exp(log(x0)+(log(x1)-log(x0))*(log(y/y0))/(log(y1/y0)));
    }
    return -1;
  };
  printf("\n   dm2 (eV2)  |  this-data    pub-data  ratio |  this-med     pub-med  ratio\n");
  for(double y : {0.1,0.3,1.0,3.0,10.0,50.0}){
    double td=evalx(g_data,y), rd=evalx(ref_data,y), tm=evalx(g_med,y), rm=evalx(ref_med,y);
    printf("   %8.2f   |  %9.3e  %9.3e  %5.2f |  %9.3e  %9.3e  %5.2f\n",
           y, td, rd, (rd>0&&td>0)?td/rd:-1, tm, rm, (rm>0&&tm>0)?tm/rm:-1);
  }
}
