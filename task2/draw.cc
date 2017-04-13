
void draw(){
  // Set simple presentation style
  gROOT->SetStyle("Plain");
  
  TFile *file = new TFile("tree.root");
  TTree *tree = NULL;
  file->GetObject("tree", tree);

  TCanvas *c = new TCanvas("c","Canvas",800,600);

  c->Divide(3,3);

  c->cd(1);
  tree->Draw("esum","pid==0");
  
  c->cd(2);
  tree->Draw("esum","pid==1");
  
  c->cd(3);
  tree->Draw("esum","pid==2");
  
  c->cd(4);
  tree->Draw("elayer:Iteration$","pid==0","box");
  
  c->cd(5);
  tree->Draw("elayer:Iteration$","pid==1","box");
  
  c->cd(6);
  tree->Draw("elayer:Iteration$","pid==2","box");
  
  c->cd(7);
  tree->Draw("elayer:Iteration$","pid==0","prof");
  
  c->cd(8);
  tree->Draw("elayer:Iteration$","pid==1","prof");
  
  c->cd(9);
  tree->Draw("elayer:Iteration$","pid==2","prof");
}
