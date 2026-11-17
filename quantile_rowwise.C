// quantile_rowwise.C
// Sensitivity quantiles computed row-wise: for each dm2 row and each SM toy,
// locate the CLs = 0.05 crossing directly in amplitude (log-linear interpolation
// between adjacent grid points), then take quantiles across toys at that row.
// Same mathematical object as production's per-toy 95% contour + func_quantiles,
// but immune to the contour fragmentation (18-39 disconnected pieces per toy)
// and to Eval() on a non-monotonic transposed graph.
void quantile_rowwise(TString indir="output_scan_3p1_exact", int NUM_dm2=60, int NUM_amp=60,
                      double DM2_LO=-2, double DM2_HI=2, double AMP_LO=-4, double AMP_HI=0,
                      TString tag="")
{
  TH1D hd("hd","",NUM_dm2,DM2_LO,DM2_HI), ha("ha","",NUM_amp,AMP_LO,AMP_HI);
  TChain ch("tree"); ch.Add(indir+"/*.root");
  int idm2,iamp; double cls_data; std::vector<double> *cls_sm=0;
  ch.SetBranchAddress("idm2",&idm2); ch.SetBranchAddress("iamp",&iamp);
  ch.SetBranchAddress("cls_data",&cls_data); ch.SetBranchAddress("cls_sm",&cls_sm);
  int nent=ch.GetEntries(); ch.GetEntry(0); int NSM=cls_sm->size();
  printf(" entries %d, SM toys %d\n", nent, NSM);

  std::vector<std::vector<double>> CD(NUM_dm2, std::vector<double>(NUM_amp,-1));
  std::vector<std::vector<std::vector<double>>> CS(NUM_dm2,
      std::vector<std::vector<double>>(NUM_amp, std::vector<double>()));
  for(int i=0;i<nent;i++){ ch.GetEntry(i);
    CD[idm2][iamp]=cls_data; CS[idm2][iamp]=*cls_sm; }

  auto cross=[&](double c0,double c1,int j)->double{
    double x0=ha.GetBinCenter(j), x1=ha.GetBinCenter(j+1);
    return pow(10, x0+(x1-x0)*(0.05-c0)/(c1-c0)); };

  TGraph *g_data=new TGraph(), *g_med=new TGraph();
  TGraph *g1L=new TGraph(), *g1R=new TGraph(), *g2L=new TGraph(), *g2R=new TGraph();
  const double s1=0.6827, s2=0.9545;

  for(int i=0;i<NUM_dm2;i++){
    double y=pow(10, hd.GetBinCenter(i+1));
    for(int j=1;j<NUM_amp;j++){
      if(CD[i][j-1]>0.05 && CD[i][j]<=0.05){ g_data->SetPoint(g_data->GetN(), cross(CD[i][j-1],CD[i][j],j), y); break; } }
    std::vector<double> lim;
    for(int m=0;m<NSM;m++)
      for(int j=1;j<NUM_amp;j++){
        if(CS[i][j-1].empty()||CS[i][j].empty()) continue;
        double c0=CS[i][j-1][m], c1=CS[i][j][m];
        if(c0>0.05 && c1<=0.05){ lim.push_back(cross(c0,c1,j)); break; } }
    if(lim.size()<10) continue;
    std::sort(lim.begin(), lim.end());
    int n=lim.size();
    auto q=[&](double qq)->double{ for(int k=0;k<n;k++) if(k+1>=n*qq) return lim[k]; return lim[n-1]; };
    g_med->SetPoint(g_med->GetN(), q(0.5), y);
    g1L->SetPoint(g1L->GetN(), q(0.5-s1/2), y);  g1R->SetPoint(g1R->GetN(), q(0.5+s1/2), y);
    g2L->SetPoint(g2L->GetN(), q(0.5-s2/2), y);  g2R->SetPoint(g2R->GetN(), q(0.5+s2/2), y);
  }
  printf(" data pts %d, median pts %d\n", g_data->GetN(), g_med->GetN());

  TFile fr("/home/jjo/Work/MicroBooNE/natureReproducer/data/BNBNuMI_3p1_data_median_brazil_asimov_Apeparance.root");
  TGraph *rd=(TGraph*)fr.Get("gh_meas"), *rm=(TGraph*)fr.Get("gh_median_val");
  auto evalx=[](TGraph*g,double y)->double{ if(!g) return -1;
    for(int k=0;k<g->GetN()-1;k++){ double x0,y0,x1,y1; g->GetPoint(k,x0,y0); g->GetPoint(k+1,x1,y1);
      if((y0-y)*(y1-y)<=0&&y0!=y1) return exp(log(x0)+(log(x1)-log(x0))*(log(y/y0))/(log(y1/y0))); } return -1.0; };
  printf("\n   dm2     med_ours   med_pub  ratio |  data_ours data_pub  ratio\n");
  for(int i=0;i<=16;i++){ double y=pow(10,-0.3+1.3*i/16.0);
    double a=evalx(g_med,y),b=evalx(rm,y),c=evalx(g_data,y),d=evalx(rd,y);
    printf("%7.2f  %9.3e %9.3e %5.2f | %9.3e %9.3e %5.2f\n",y,a,b,(a>0&&b>0)?a/b:-1,c,d,(c>0&&d>0)?c/d:-1); }
  std::vector<double> rm_, rd_;
  for(int i=0;i<=200;i++){ double y=pow(10,-1.5+3.4*i/200.0);
    double a=evalx(g_med,y),b=evalx(rm,y),c=evalx(g_data,y),d=evalx(rd,y);
    if(a>0&&b>0) rm_.push_back(a/b); if(c>0&&d>0) rd_.push_back(c/d); }
  std::sort(rm_.begin(),rm_.end()); std::sort(rd_.begin(),rd_.end());
  printf("\n OVERALL  data: N=%d median %.3f [%.3f,%.3f] | median-sens: N=%d median %.3f [%.3f,%.3f]\n",
    (int)rd_.size(), rd_[rd_.size()/2], rd_[(int)(0.16*rd_.size())], rd_[(int)(0.84*rd_.size())],
    (int)rm_.size(), rm_[rm_.size()/2], rm_[(int)(0.16*rm_.size())], rm_[(int)(0.84*rm_.size())]);

  TFile fo("contour_3p1_rowwise"+tag+".root","recreate");
  g_data->Write("g_data_thiswork"); g_med->Write("g_sens_median_thiswork");
  g1L->Write("g_1s_LL"); g1R->Write("g_1s_RR"); g2L->Write("g_2s_LL"); g2R->Write("g_2s_RR");
  if(rd) rd->Write("ref_data"); if(rm) rm->Write("ref_median");
  fo.Close();
}
