#include <math>

const Int_t nbin = 100;
const Double_t xmin=0, xmax=100, xstep=(xmax-xmin)/nbin;

Double_t gs(Double_t x, Double_t x0, Double_t sig) {
  return 1/sqrt(2*TMath::Pi())/sig*exp(-pow((x-x0),2)/(2*sig*sig));
}

Double_t lgs(Double_t x, Double_t x0, Double_t sig, Double_t eta) {
  
  if( fabs(eta)<0.001 ) {
    return gs(x,x0,sig);
  }

  const double l4 = log(4.0);
  
  double x1 = 1-(x-x0)*eta/sig;
  if( x1<=0 ) return 0;

  double s0 = log(eta*sqrt(l4)+sqrt(1+l4*eta*eta))/sqrt(l4);
  double s02 = s0*s0;

  double res = (1/sqrt(2*TMath::Pi())/sig)*(eta/s0)*exp(-0.5*(pow(log(x1),2)/s02+s02));

  // Normalization
  if( eta>0 ) {
    double norm = TMath::Freq(log(1-(0.0-x0)*eta/sig)/s0-s0);
    res /= norm;
  }

  return res;

}

Double_t func_mu(Double_t *x, Double_t *par) {
  return par[0]*lgs(x[0],par[1],par[2],par[3]);
  //return par[0]*gs(x[0],par[1],par[2]);
}

Double_t func_e(Double_t *x, Double_t *par) {
  return par[0]*lgs(x[0],par[1],par[2],par[3]);
  //return par[0]*gs(x[0],par[1],par[2]);
}

Double_t func(Double_t *x, Double_t *par) {
  return func_e(x, par) + func_mu(x, par + 4);
}

void fit_mu() {
  TFile *file = new TFile("tree.root");
  TTree *tree = NULL;
  file->GetObject("tree", tree);

  // Fill histogram
  TH1D *h1 = new TH1D("h1","mu-",nbin,xmin,xmax);
  tree->Draw("esum>>h1","pid==2");

  TF1 *f1 = new TF1("f1",func_mu,xmin,xmax,4);
  //f1->SetParameters(3333.0,150.0,50.0,-0.2);
  f1->SetParameters(3333.0,35.0,5.0,-0.2);
  //TF1 *f1 = new TF1("f1",func_mu,xmin,xmax,3);
  //f1->SetParameters(1000.0,55.0,1.0);

  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->Divide(1,1);
  c1->cd(1);
  
  h1->Fit(f1,"L");
}
