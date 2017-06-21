#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

void task5() {

  gSystem->Load("libTMVA");

  TMVA::Tools::Instance();  

  TFile* fout = TFile::Open("tmva_out.root","RECREATE");

  TMVA::Factory *factory = new TMVA::Factory(
    "myTMVA",
    fout, 
    "!V:!Silent"
    //"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification"
    );

  TChain *tin = new TChain("tree");
  tin->AddFile("tree.root");

  factory->SetInputTrees(tin, "pid==0", "pid==1");

  factory->AddVariable("el0 := elayer[0]",'F');
  factory->AddVariable("el1 := elayer[1]",'F');
  factory->AddVariable("el2 := elayer[2]",'F');
  factory->AddVariable("el3 := elayer[3]",'F');
  factory->AddVariable("el4 := elayer[4]",'F');
  factory->AddVariable("el5 := elayer[5]",'F');
  factory->AddVariable("el6 := elayer[6]",'F');
  factory->AddVariable("el7 := elayer[7]",'F');
  factory->AddVariable("el8 := elayer[8]",'F');
  factory->AddVariable("el9 := elayer[9]",'F');

  factory->PrepareTrainingAndTestTree("","");

  //factory->BookMethod(TMVA::Types::kLikelihood, "Likelihood", "H:!V:VarTransform=N,D:PDFInterpol=Spline2:TransformOutput=True");
  //factory->BookMethod(TMVA::Types::kPDERS, "PDERS", "H:!V:VolumeRangeMode=RMS:KernelEstimator=Gauss");
  factory->BookMethod(TMVA::Types::kBDT, "BDT", "H:!V:NTrees=800:MaxDepth=5:SeparationType=GiniIndex:BoostType=AdaBoost");

  factory->TrainAllMethods();

  factory->TestAllMethods();

  factory->EvaluateAllMethods();

}
