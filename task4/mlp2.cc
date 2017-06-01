{
  // Load library
  
  if (!gROOT->GetClass("TMultiLayerPerceptron")) {
    gSystem->Load("libMLP");
  }
  
  TPad* mlpa_canvas = new TCanvas("mlpa_canvas","",1000,700);
  mlpa_canvas->Divide(2,2);

  // Create tree with train data
  TFile *fin = new TFile("tree.root");
  TTree *tin = NULL;
  fin->GetObject("tree", tin);
  
  TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron(
    "@elayer:16:pid",
		tin,
		"(Entry$%10) && pid<2",
		"!(Entry$%10) && pid<2"
		);
 
  mlpa_canvas->cd(1);
  mlp->Train(400, "text,graph,update=10,current");
  mlp->Export("mlp2_func");

  //gROOT->LoadMacro("mlp2_func.cxx");
  //gROOT->LoadMacro("mlp_func.cxx");
  
  //  TPad* mlpa_canvas = draw_canvas("mlpa_canvas","",1000,700);
  TMLPAnalyzer ana(mlp);
  ana.GatherInformations();
  ana.CheckNetwork();

  // shows the network structure
  mlpa_canvas->cd(2);
  mlp->Draw();
  
  //mlpa_canvas->cd(2);
  //ana.DrawDInputs();
  
  mlpa_canvas->cd(3);
  //ana.DrawNetwork(0,"type==1","type==0");

  int nbin = 200;
  int thres = nbin/2+1;
  TH1F *bg = new TH1F("bgh", "NN output", nbin, -.5, 1.5);
  TH1F *sig = new TH1F("sigh", "Signal", nbin, -.5, 1.5);
  bg->SetDirectory(0);
  sig->SetDirectory(0);
  int pid;
  float elayer[10];
  float escale[10];
  for (int i = 0; i < 10; ++i) {
    escale[i] = 1.0*(gRandom->Uniform() - 0.5) + 1.0;
  }
  tin->SetBranchAddress("pid", &pid);
  tin->SetBranchAddress("elayer", &elayer);
  int nsig = 0, nbg = 0;
  Double_t params[10];
  for (int i = 0, N = tin->GetEntries(); i < N; ++i) {
    tin->GetEntry(i);
    if (pid >= 2) {
      continue;
    }
    TH1F *h;
    if (pid > 0) {
      h = sig;
      nsig += 1;
    } else {
      h = bg;
      nbg += 1;
    }
    for (int j = 0; j < 10; ++j) {
      params[j] = elayer[j]*escale[j];
    }
    h->Fill(mlp->Evaluate(0, params));
  }
  bg->SetLineColor(kBlue);
  bg->SetFillStyle(3008);   bg->SetFillColor(kBlue);
  sig->SetLineColor(kRed);
  sig->SetFillStyle(3003); sig->SetFillColor(kRed);
  bg->SetStats(0);
  sig->SetStats(0);
  bg->Draw();
  sig->Draw("same");
  TLegend *legend = new TLegend(.75, .80, .95, .95);
  double eff = sig->Integral(thres,nbin+2)/sig->Integral(0,nbin+2)*100.0;
  double err =  bg->Integral(thres,nbin+2)/ bg->Integral(0,nbin+2)*100.0;
  legend->AddEntry(sig, Form("Signal (%4.1f%)",eff));
  legend->AddEntry(bg, Form("Background (%4.1f%)",err));
  legend->Draw();
  
  double *x = new double[nbin], *y = new double[nbin];
  int csig = 0, cbg = 0;
  for (int i = 0; i < nbin; ++i) {
    csig += sig->GetBinContent(i);
    cbg += bg->GetBinContent(i);
    x[i] = double(csig)/nsig;
    y[i] = double(cbg)/nbg;
  }
  
  mlpa_canvas->cd(4);
  TGraph *graph = new TGraph(nbin,x,y);
  graph->Draw("AL*");
}
