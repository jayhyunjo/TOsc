// compare_rows.C
// One table: for every row test on disk, the data limit, median, band edges and
// width, each as a ratio to the published curves. Configurations compared:
//   final2  = data-profiled fixed hypothesis (current default)
//   fcboth  = per-toy profiling
//   asimov  = Asimov-profiled fixed hypothesis
void compare_rows(){
  TString B="/home/jjo/Work/MicroBooNE/ToAll_basic_nospeedup/";
  TFile fr("/home/jjo/Work/MicroBooNE/natureReproducer/data/BNBNuMI_3p1_data_median_brazil_asimov_Apeparance.root");
  TGraph *pmeas=(TGraph*)fr.Get("gh_meas"), *pmed=(TGraph*)fr.Get("gh_median_val");
  TGraph *pl=(TGraph*)fr.Get("gh_median_m1sigma"), *pr=(TGraph*)fr.Get("gh_median_p1sigma");
  auto ev=[](TGraph*g,double y)->double{ if(!g) return -1;
    for(int k=0;k<g->GetN()-1;k++){ double x0,y0,x1,y1; g->GetPoint(k,x0,y0); g->GetPoint(k+1,x1,y1);
      if((y0-y)*(y1-y)<=0&&y0!=y1) return exp(log(x0)+(log(x1)-log(x0))*(log(y/y0))/(log(y1/y0))); } return -1.0; };

  auto one=[&](TString dir,int ROW,const char*lab){
    TChain ch("tree"); ch.Add(B+dir+"/*.root");
    if(ch.GetEntries()==0) return;
    int idm2,iamp; double cls_data; std::vector<double>*cs=0;
    ch.SetBranchAddress("idm2",&idm2); ch.SetBranchAddress("iamp",&iamp);
    ch.SetBranchAddress("cls_data",&cls_data); ch.SetBranchAddress("cls_sm",&cs);
    std::vector<double> CD(60,-1); std::vector<std::vector<double>> C(60);
    for(int i=0;i<ch.GetEntries();i++){ ch.GetEntry(i); if(idm2==ROW){ CD[iamp]=cls_data; C[iamp]=*cs; } }
    int NSM=0; for(auto&v:C) if(!v.empty()){NSM=v.size();break;}
    if(!NSM) return;
    double dm2=pow(10,-2+(ROW+0.5)*4./80);
    auto cross=[&](double c0,double c1,int j){ double x0=-4+(j-0.5)*4./60,x1=-4+(j+0.5)*4./60;
      return x0+(x1-x0)*(0.05-c0)/(c1-c0); };
    double xd=-1;
    for(int j=1;j<60;j++){ if(CD[j-1]<0||CD[j]<0) continue;
      if(CD[j-1]>0.05&&CD[j]<=0.05){ xd=pow(10,cross(CD[j-1],CD[j],j)); break; } }
    std::vector<double> L;
    for(int m=0;m<NSM;m++) for(int j=1;j<60;j++){ if(C[j-1].empty()||C[j].empty()) continue;
      double c0=C[j-1][m],c1=C[j][m];
      if(c0>0.05&&c1<=0.05){ L.push_back(cross(c0,c1,j)); break; } }
    if((int)L.size()<50) return;
    std::sort(L.begin(),L.end()); int n=L.size();
    double om=pow(10,L[n/2]), ol=pow(10,L[(int)(.1587*n)]), orr=pow(10,L[(int)(.8413*n)]);
    double pdd=ev(pmeas,dm2), pmm=ev(pmed,dm2), pll=ev(pl,dm2), prr=ev(pr,dm2);
    printf("  %-8s dm2=%7.3f | data %5.2f | med %5.2f | left %5.2f | right %5.2f | width %.3f vs %.3f\n",
      lab, dm2, (xd>0&&pdd>0)?xd/pdd:-1, om/pmm, ol/pll, orr/prr, log10(orr/ol), log10(prr/pll));
  };
  printf("All ratios are ours/published. Target: everything 1.00 within 2-3 percent.\n\n");
  int rows[5]={30,40,47,61,72};
  printf("data-profiled fixed hypothesis (current default, full grid):\n");
  for(int r : rows) one("output_3p1_final2",r,"final2");
  printf("\nper-toy profiling (-fc 1):\n");
  for(int r : rows) one(Form("test_fcboth_r%d",r),r,"fcboth");
  printf("\nAsimov-profiled fixed hypothesis (-prof 1):\n");
  for(int r : rows) one(Form("test_asimov_r%d",r),r,"asimov");
}
