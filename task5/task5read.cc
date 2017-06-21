#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

using namespace TMVA;

char digit[] = "0123456789";

void book(Reader *reader, const char *method) {
  std::string path(method);
  path = "weights/myTMVA_" + path + ".weights.xml";
  reader->BookMVA(method, path.c_str());
}

void task5read() {

  gSystem->Load("libTMVA");

  Reader* reader = new Reader("Color");

  TChain *tin = new TChain("tree");
  tin->Add("tree.root");

  int pid; float esum, elayer[10];
  tin->GetBranch("pid")->SetAddress(&pid);
  tin->GetBranch("esum")->SetAddress(&esum);
  tin->GetBranch("elayer")->SetAddress(elayer);

  float el[10];
  for( int i=0; i<10; i++ ) 
    reader->AddVariable(Form("elayer[%d]",i),&(el[i]));
  
  TCanvas *c = new TCanvas("c","",900,600);
  c->Divide(2,2);
  
  int nbin = 100;
  double b[6] = {-3.0, 3.0, -0.5, 1.5, -0.5, 0.5};
  TH1D *hist[6];
  for (int i = 0; i < 3; ++i) {
    std::string name;
    name = std::string("sig") + digit[i];
    hist[2*i] = new TH1D(name.c_str(),name.c_str(),nbin,b[2*i],b[2*i + 1]);
    name = std::string("bg") + digit[i];
    hist[2*i + 1] = new TH1D(name.c_str(),name.c_str(),nbin,b[2*i],b[2*i + 1]);
  }
  
  const char *methods[] = {"Likelihood", "PDERS", "BDT"};
  for (int type = 0; type < 3; ++type) {
    const char *method = methods[type];
    book(reader, method);

    for( int i=0; i<tin->GetEntries(); i++ ) {
      tin->GetEvent(i);
      if ( pid != 0 && pid != 1 ) continue;
      for( int j=0; j<10; j++ ) el[j]=elayer[j];
      double out = reader->EvaluateMVA(method);
      if( pid==0 ) hist[2*type]->Fill(out);
      else hist[2*type + 1]->Fill(out);
      //printf("%f\n", out);
    }

    hist[2*type]->SetLineColor(kRed);
    hist[2*type]->SetFillStyle(3003);
    hist[2*type]->SetFillColor(kRed);

    hist[2*type + 1]->SetLineColor(kBlue);
    hist[2*type + 1]->SetFillColor(kBlue);
    hist[2*type + 1]->SetFillStyle(3008);

    c->cd(type + 1);
    hist[2*type]->Draw();
    hist[2*type + 1]->Draw("same");
  }

  c->cd(4);
  TGraph *graph[3];
  for (int type = 0; type < 3; ++type) {
    double *x = new double[nbin], *y = new double[nbin];
    int csig = 0, cbg = 0;
    int nsig = hist[2*type]->GetEntries();
    int nbg = hist[2*type + 1]->GetEntries();
    for (int i = 0; i < nbin; ++i) {
      csig += hist[2*type]->GetBinContent(i);
      cbg += hist[2*type + 1]->GetBinContent(i);
      x[i] = double(csig)/nsig;
      y[i] = double(cbg)/nbg;
    }
    
    graph[type] = new TGraph(nbin,x,y);
    int color[] = {kRed, kBlue, kGreen};
    graph[type]->SetLineColor(color[type]);
    graph[type]->SetMarkerColor(color[type]);
    if(type == 0) {
      graph[type]->Draw("AL*");
    } else {
      graph[type]->Draw("LP*");
    }
    delete[] x;
    delete[] y;
  }
}
