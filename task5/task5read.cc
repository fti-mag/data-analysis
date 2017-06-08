#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

using namespace TMVA;

void task5read() {

  gSystem->Load("libTMVA");

  Reader* reader = new Reader("Color");

  TChain *tin = new TChain("tree");
  tin->Add("tree.root");

  int pid; double esum, elayer[10];
  tin->GetBranch("pid")->SetAddress(&pid);
  tin->GetBranch("esum")->SetAddress(&esum);
  tin->GetBranch("elayer")->SetAddress(elayer);

  float el[10];
  for( int i=0; i<10; i++ ) 
    reader->AddVariable(Form("elayer[%d]",i),&(el[i]));
  
  const char *method = "BDT";
  std::string path(method);
  path = "weights/myTMVA_" + path + ".weights.xml";
  reader->BookMVA(method, path.c_str());

  double l = -0.2, r = 0.4;
  int nbin = 100;
  TH1D *h0 = new TH1D("h0","h0",nbin,l,r);
  TH1D *h1 = new TH1D("h1","h1",nbin,l,r);

  for( int i=0; i<tin->GetEntries(); i++ ) {
    tin->GetEvent(i);
    if ( pid != 0 && pid != 1 ) continue;
    for( int j=0; j<10; j++ ) el[j]=elayer[j];
    double out = reader->EvaluateMVA(method);
    if( pid==0 ) h0->Fill(out);
    else h1->Fill(out);
  }

  TCanvas *c = new TCanvas("c","",900,600);
  c->Divide(2,1);

  h0->SetLineColor(kRed);
  h0->SetFillStyle(3003);
  h0->SetFillColor(kRed);

  h1->SetLineColor(kBlue);
  h1->SetFillColor(kBlue);
  h1->SetFillStyle(3008);

  c->cd(1);
  h0->Draw();
  h1->Draw("same");

  double *x = new double[nbin], *y = new double[nbin];
  int csig = 0, cbg = 0;
  int nsig = h0->GetEntries(), nbg = h1->GetEntries();
  for (int i = 0; i < nbin; ++i) {
    csig += h0->GetBinContent(i);
    cbg += h1->GetBinContent(i);
    x[i] = double(csig)/nsig;
    y[i] = double(cbg)/nbg;
  }
  
  c->cd(2);
  TGraph *graph = new TGraph(nbin,x,y);
  graph->Draw("AL*");
}
