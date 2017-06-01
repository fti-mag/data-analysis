#include <math>

const Int_t nbin = 300;
const Double_t xmin=0, xmax=300, xstep=(xmax-xmin)/nbin;
Double_t bin[nbin];

//---------------------------------------------

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
}

Double_t func(Double_t *x, Double_t *par) {
   return func_mu(x, par) + func_mu(x, par + 4);
}

/*
// chi-square
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t x[10];

  Double_t lfunc = 0;
  for( int i=0; i<nbin; i++ ) {
    if( bin[i]>0 ) {
      x[0] = (i+0.5)*xstep;
      double fval = func(x,par);
      double diff = bin[i] - fval;
      lfunc += diff*diff;
    }
  }
  //lfunc += par[0];

  f = lfunc;
}
*/

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t x[10];

  Double_t lfunc = 0;
  for( int i=0; i<nbin; i++ ) {
    if( bin[i]>0 ) {
      x[0] = (i+0.5)*xstep;
      double fval = func(x,par);
      if( fval>0 ) {
	      lfunc -= bin[i]*log(fval);
      }
    }
  }
  lfunc += par[0] + par[4];

  f = lfunc;
}

void fit_mu_minuit() {
  TFile *file = new TFile("tree.root");
  TTree *tree = NULL;
  file->GetObject("tree", tree);

  // Fill histogram
  TH1D *h1 = new TH1D("h1","e- and mu-",nbin,xmin,xmax);
  tree->Draw("esum>>h1","pid==1||pid==2");

  TCanvas *c1 = new TCanvas("c1","",800,800);
  c1->Divide(1,2);

  // Draw data histogram
  c1->cd(1);
  h1->Draw();

  // Get data
  for( int i=0; i<nbin; i++ ) bin[i] = h1->GetBinContent(i+1);

  // Initialize TMinuit with a maximum of 8 params

  const int npar = 8;
  TMinuit *gMinuit = new TMinuit(npar);  
  gMinuit->SetFCN(fcn);
 
  Double_t arglist[10];
  Int_t ierflg = 0;
 
  arglist[0] = 0.5;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
 
  // Set starting values and step sizes for parameters		
  const Double_t vstart[npar] = {3333.0,120.0,80.0,0.2, 3333.0,35.0,5.0,-0.2};
  const Double_t step[npar] = {10,1,0.1,0.01, 10,1,0.1,0.01};
  const char *name[npar] = {"Ne", "Ee", "Se", "Ae", "Nmu", "Emu", "Smu", "Amu"};
  for (int i = 0; i < npar; ++i) {
    gMinuit->mnparm(i, name[i], vstart[i], step[i], 0,0,ierflg);
  }

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
 
  // Get the results
  Double_t pval[npar], err[npar];
  for( int i=0; i<npar; i++ ) {
    gMinuit->GetParameter(i, pval[i], err[i]);
  }

  // Plot fit function
  TF1 *f = new TF1("f",func,xmin,xmax,npar);
  f->SetParameters(pval);
  f->SetNpx(200);
  f->Draw("same");

  // Plot contour plot between 2 and 3 parameters

  c1->cd(2);
  TGraph *gr1 = (TGraph*)gMinuit->Contour(40,0,4);
  gr1->Draw("AL");

}
