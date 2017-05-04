double sqr(double x) {
  return x*x;
}

double sqr12(double x) {
  return sqr(x - 0.5);
}

double func(double x) {
  return exp(-x)*sqr12(x);
}

const double E = 2.7182818284590452353602874713527;
const double norm = 5.0/4.0 - 13.0/(4.0*E);
const double x0 = (46.0 - 17.0*E)/(13.0 - 5.0*E);
const double sig0 = sqrt((201.0 - 74.0*E)/(13.0 - 5.0*E) - x0*x0);

double rand_unif() {
  return gRandom->Uniform();
}

double rand_sqr() {
  double x = 0.25*(gRandom->Uniform() - 0.5);
  return (1 - 2*(x < 0.0))*pow(fabs(x), 1.0/3.0) + 0.5;
}

//#define DISTORTION

double getRandom() {
  double x, y;
  do { 
    x = rand_sqr();
    y = rand_unif();
  } while (exp(-x) < y);
#ifndef DISTORTION
  return x;
#else
  if (gRandom->Uniform() < 0.01) {
    return x + 10*gRandom->Gaus(0,1);
  } else {
    return x;
  }
#endif
}

void task1(){

  std::cout << "x0: " << x0 << std::endl;
  std::cout << "sig0: " << sig0 << std::endl;

  // 1. Fill histogram with random numbers distributed according to
  //    predefined distribution  
  // 2. Calculate median and mean
  // 3. Show distribution of average of 2, 5 and 100 numbers

  // Set simple presentation style

  gROOT->SetStyle("Plain");

  // Constants
  
  const int nbin = 50;
  const int ndemo = 1000;
  const int nevent = 100;
  const int nexperiment = 1000;

  // Book histograms

  TH1D *h1 = new TH1D("h1","Initial distribution",nbin,0,1);  
  TH1D *h2 = new TH1D("h2","Distribution of median",nbin,0,1);
  TH1D *h3 = new TH1D("h3","Distribution of mean",nbin,0,1);
  TH1D *h4 = new TH1D("h4","Distribution of sum of 2 numbers",nbin,-3,3);
  TH1D *h5 = new TH1D("h5","Distribution of sum of 5 numbers",nbin,-3,3);
  TH1D *h6 = new TH1D("h6","Distribution of sum of 100 numbers",nbin,-3,3);
  TF1 *fgauss = new TF1("fgauss","gausn",-3,3);
  fgauss->SetParameters(ndemo*6/nbin,0,1);
  TF1 *ffunc = new TF1("ffunc","[0]*func(x)",0,1);
  ffunc->SetParameter(0, ndemo/nbin/norm);

  // Fill initial distribution

  for( int i=0; i<ndemo; i++ ) {
    double x = getRandom();
    h1->Fill(x);
  }

  // Fill distributions of mean
  double *data = new double[nevent];
  int *indices = new int[nevent];
  for( int iex=0; iex<nexperiment; iex++ ) {
    double mean = 0;
    double median = 0;
    for( int i=0; i<nevent; i++ ) {
      double x = getRandom();
      mean += x;
      data[i] = x;
    }
    mean /= nevent;
    TMath::Sort(nevent, data, indices);
    median = 0.5*(data[indices[(nevent - 1)/2]] + data[indices[nevent/2]]);
    h3->Fill(mean);
    h2->Fill(median);
  }
  delete[] data;
  delete[] indices;

  // Fill distributions of median

  // Fill demo distributions
  double sum;
  int ns[] = {2, 5, 100};
  TH1D *hs[] = {h4, h5, h6};
  for(int i = 0; i < ndemo; i++) {
    int j = 0;
    sum = 0.0;
    for (int k = 0; k < sizeof(ns)/sizeof(ns[0]); ++k) {
      for(j = 0; j < ns[k]; j++) {
        sum += (getRandom() - x0);
      }
      hs[k]->Fill(sum/(sqrt(ns[k])*sig0));
    }
  }

  // Draw histograms

  TCanvas *c1 = new TCanvas("c1","Example #1",800,600);

  c1->Divide(3,2);

  c1->cd(1);
  h1->Draw();
  ffunc->Draw("same");

  c1->cd(2);
  h2->Draw();

  c1->cd(3);
  h3->Draw();

  c1->cd(4);
  h4->Draw();
  fgauss->Draw("same");

  c1->cd(5);
  h5->Draw();
  fgauss->Draw("same");

  c1->cd(6);
  h6->Draw();
  fgauss->Draw("same");
}
  
