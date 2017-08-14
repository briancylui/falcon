// ------------------------------------------------------------------------
// FALCON TEST FILE
// ------------------------------------------------------------------------
#include <iostream>
#include <cmath>
#include <cassert>

#include "TKDTreeBinning.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH2Poly.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TMap.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "FalconTester.h"

// For TMVARegression with MLP
#include <cstdlib>
#include <map>
#include <string>

#include "TString.h"
#include "TObjString.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"

// #include<omp.h>

using namespace std;

void  FalconTester::Add(string filename)
{
  inputFiles.push_back(filename);
}

// Match Jets with Genjets based on min(dR)
void  FalconTester::DeltaRMatch(double dRcut)
{
  TChain chain("Delphes");
  for(size_t c=0; c < inputFiles.size(); c++)
    {
      cout << "Add: " << inputFiles[c] << endl;
      chain.Add(inputFiles[c].c_str());
    }
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  
  for(Int_t entry = 0; entry < numberOfEntries; entry++) //Entry loop
    {
      treeReader->ReadEntry(entry);
      int genjet_entries=branchGenJet->GetEntries();

      double mindR = 999;
      for(int gj=0; gj<genjet_entries; gj++) //GenJet loop
        {
	  Jet* genjet = (Jet*) branchGenJet->At(gj);
	  int jet_entries=branchJet->GetEntries();

	  Jets.push_back(JetMapItem());
	  if ( Jets.size() % 10000 == 0 )
	    cout << "\tjet count = " << Jets.size() << endl;
	  
	  JetMapItem& item = Jets.back(); // NB: get a reference not a copy
	  item.genjet.SetPtEtaPhiM(genjet->PT,
				   genjet->Eta,
				   genjet->Phi,
				   genjet->Mass);
	  
	  for(int j=0; j<jet_entries; j++) //Jet loop
            {
	      Jet *jet = (Jet*) branchJet->At(j);
	      TLorentzVector jet_v;
	      jet_v.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
	      double deltar = item.genjet.DeltaR(jet_v);
	      if( deltar < dRcut ) item.jet.push_back(jet_v);
	      if ( deltar < mindR ) mindR = deltar;
            }
	  drmin->Fill(mindR); 
        }
    }
}

void FalconTester::Match(std::string lookfilename,
			 double dRcut)
{
  assert(inputFiles.size()>0);
    
  cout << "1. match genjets to jets..." << endl;
  DeltaRMatch(dRcut);
  
  int totaljets = Jets.size();
  
  std::cout << "\tjet count = " << totaljets << std::endl;

  cout << "2. create lookup table" << endl;
  
  TFile rfile(lookfilename.c_str(), "RECREATE");
  rfile.cd();
  TTree* tree = new TTree("Falcon", "Lookup Table");
  int matched;
  double genPt, genEta, genPhi;
  double Pt, Eta, Phi, Mass;
  tree->Branch("matched", &matched, "matched/I");
  tree->Branch("genPt",   &genPt,   "genPt/D");
  tree->Branch("genEta",  &genEta,  "genEta/D");
  tree->Branch("genPhi",  &genPhi,  "genPhi/D");
  
  tree->Branch("Pt",      &Pt,      "Pt/D");
  tree->Branch("Eta",     &Eta,     "Eta/D");
  tree->Branch("Phi",     &Phi,     "Phi/D");
  tree->Branch("Mass",    &Mass,    "Mass/D");

  for (int entry=0; entry < totaljets; entry++)//loop over entries
    {
      if ( Jets.size() % 10000 == 0 )
	cout << "\tjet count = " << Jets.size() << endl;

      // get references to matched jets
      TLorentzVector&      genjet = Jets[entry].genjet;
      vector<TLorentzVector>& jet = Jets[entry].jet;

      genPt  = genjet.Pt();
      genEta = genjet.Eta();
      genPhi = genjet.Phi();
      
      if ( jet.size() > 0 )
	{
	  matched = 1;
	  Pt  = jet[0].Pt();
	  Eta = jet[0].Eta();
	  Phi = jet[0].Phi();
	  Mass= jet[0].M();
	}
      else
	{
	  matched = 0;
	  Pt  = 0;
	  Eta = 0;
	  Phi = 0;
	  Mass= 0;
	}
	
      rfile.cd();
      tree->Fill();
    }
  cout << "3. write out lookup table" << endl;
  tree->Write();
  rfile.Close();
}


void FalconTester::Build(std::string filename)
{
  cout << "1. open file " << filename << endl;

  TFile rfile(filename.c_str());
  rfile.cd();
  TTree* tree = (TTree*)rfile.Get("Falcon");
  int matched;
  double genPt, genEta, genPhi;
  double Pt, Eta, Phi, Mass;
  tree->SetBranchAddress("matched", &matched);
  tree->SetBranchAddress("genPt",   &genPt);
  tree->SetBranchAddress("genEta",  &genEta);
  tree->SetBranchAddress("genPhi",  &genPhi);
  
  tree->SetBranchAddress("Pt",      &Pt);
  tree->SetBranchAddress("Eta",     &Eta);
  tree->SetBranchAddress("Phi",     &Phi);
  tree->SetBranchAddress("Mass",    &Mass);
  
  int totaljets = tree->GetEntries();
  
  const UInt_t DATASZ  = totaljets;
  const UInt_t DATADIM = 3;
  const UInt_t NBINS   = totaljets;

  cout << "2. allocate ";
  Double_t* smp = new Double_t[DATASZ * DATADIM];
  if(smp==nullptr)
    {
      std::cerr<<"Can not allocate memory: "<<DATASZ * DATADIM<<" doubles\n";
      exit(0);
    }
  cout << sizeof(Double_t)*float(DATASZ*DATADIM)/1000000
       << " Mbytes of memory "
       << endl;

  cout << "3. fill memory with gen level data" << endl;

  for (int entry=0; entry < totaljets; entry++)//loop over entries
    {
      tree->GetEntry(entry);

      smp[DATASZ*0+entry] = genPt;
      smp[DATASZ*1+entry] = genEta;
      smp[DATASZ*2+entry] = genPhi;
    }

  cout << "4. bin gen level data...please be patient!" << endl;

//   TStopwatch swatch;
//   swatch.Start();
  kdt = new TKDTreeBinning(DATASZ, DATADIM, smp, NBINS);
//   cout << "\treal time: " << swatch.RealTime() << " seconds" << endl;
  
//  #pragma omp parallel for default(shared)
  for (int entry=0; entry < totaljets; entry++)//loop over entries
    {
      tree->GetEntry(entry);
//       cout << "\tjet count = " << entry << endl;
      if ( entry % 20000 == 0 ) cout << "\tjet count = " << entry << endl;
      
      double point[3];
      point[0] = genPt;
      point[1] = genEta;
      point[2] = genPhi;
      int bin  = kdt->FindBin(point);
      RecoJet jet;
      
      if ( matched )
	{
	  jet.PT   = Pt;
	  jet.Eta  = Eta;
	  jet.Phi  = Phi;
	  jet.Mass = Mass;
	}
      else
	{
	  jet.PT   = 0;
	  jet.Eta  = 0;
	  jet.Phi  = 0;
	  jet.Mass = 0;
	}	
      table[bin] = jet; 
    }
  std::cout << "\tjet count = " << totaljets << std::endl;  
  cout << "\tdone!" << endl;
  rfile.Close();
}

RecoJet FalconTester::MapJet(double pt, double eta, double phi)
{
  double point[3] = {pt, eta, phi};
  int index = kdt->FindBin(point);
  if ( table.find(index) != table.end() )
    return table[index];
  else
    return RecoJet();
}

void FalconTester::Learn(std::string filename)
{
  // First Build.                                                               
  Build(filename);

  cout << "1. open file " << filename << endl;

  TFile rfile(filename.c_str());
  rfile.cd();
  TTree* tree = (TTree*)rfile.Get("Falcon");
  int matched;
  double genPt, genEta, genPhi;
  double Pt, Eta, Phi, Mass;
  tree->SetBranchAddress("matched", &matched);
  tree->SetBranchAddress("genPt",   &genPt);
  tree->SetBranchAddress("genEta",  &genEta);
  tree->SetBranchAddress("genPhi",  &genPhi);

  tree->SetBranchAddress("Pt",      &Pt);
  tree->SetBranchAddress("Eta",     &Eta);
  tree->SetBranchAddress("Phi",     &Phi);
  tree->SetBranchAddress("Mass",    &Mass);


  int totaljets = tree->GetEntries();

  const UInt_t DATASZ  = totaljets;
  const UInt_t DATADIM = 3;
  const UInt_t NBINS   = totaljets;

  // Use MLP : 3 inputs, 4 outputs, 3-layer neural network.                     
  // First hidden layer: 10 nodes; second hidden layer: 5 nodes.                
  cout << "2. Use MLP to train a 3-layer neural network for multi-target 
regression." << endl;

  // This loads the library.                                                    
  TMVA::Tools::Instance();

  cout << endl << "==> Start TMVARegression" << endl;

  // Create a new root output file.                                             
  TString outfileName("TMVAReg.root");
  TFile* outputFile = TFile::Open(outfileName, "RECREATE");

  // Create the factory object.                                                 
  // The first argument is the base of the name of all the weightfiles 
  // in the directory weight, while the second argument is the output 
  // file for the training results.
  // All TMVA output can be suppressed by removing the "!" (not) 
  // in front of the "Silent" argument in the option string.                                       
  TMVA::Factory *factory = new TMVA::Factory("TMVARegression", outputFile, 
					     "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression");

  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

  // Define the input varialbes that shall be used for the MVA training         
  dataloader->AddVariable("genPt", "genPt", "units", 'F');
  dataloader->AddVariable("genEta", "genEta", "units", 'F');
  dataloader->AddVariable("genPhi", "genPhi", "units", 'F');

  // Add the variables carrying the regression targets                          
  dataloader->AddTarget("Pt");
  dataloader->AddTarget("Eta");
  dataloader->AddTarget("Phi");
  dataloader->AddTarget("Mass");

  // Read training and test data                                                
  // Load the signal and background event samples from ROOT trees               
  // Already opened the input file above as a TFile*: rfile.                    
  cout << "--- TMVARegression: Using input file: " << rfile->GetName() << endl;

  // Register the regression tree                                               
  // Done above: TTree* tree.                                                   

  // Add an arbitrary number of regression trees                                
  Double_t regWeight = 1.0;
  dataloader->AddRegressionTree(tree, regWeight);

  // Apply additional cuts on the signal and background samples (can be different)                                                                             
  TCut mycut = "matched==1";

  dataloader->PrepareTrainingAndTestTree(mycut, "nTrain_Regression=1000:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V");

  // Book MVA methods                                                           
  //                                                                            
  // Please lookup the various method configuration options in the corresponding cxx files, eg:                                                                
  // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html                                                                               
  // it is possible to preset ranges in the option string in which the cut optimisation should be done:                                                        
  // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable                                                                       

  factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:VarTransfo\
rm=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMeth\
od=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests\
=15:!UseRegulator" );

  // ------------------------------------------------------------------
  // Now you can tell the factory to train, test, and evaluate the MVAs         
  // Train MVAs using the set of training events                                
  factory->TrainAllMethods();
  // Evaluate all MVAs using the set of test events                             
  factory->TestAllMethods();
  // Evaluate and compare performance of all configured MVAs                    
  factory->EvaluateAllMethods();
  // --------------------------------------------------------------             
  // Save the output                                                            
  outputFile->Close();
  cout << "==> Wrote root file: " << outputFile->GetName() << endl;
  cout << "==> TMVARegression is done!" << endl;
  delete factory;
  delete dataloader;

  // Launch the GUI for the root macros                                         
  if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName );

  std::cout << "\tjet count = " << totaljets << std::endl;
  cout << "\tdone!" << endl;
  rfile.Close();
}

void FalconTester::LearnJet(double pt, double eta, double phi)
{
  double point[3] = {pt, eta, phi};
  int index = kdt->FindBin(point);
  if ( tablelearnt.find(index) != tablelearnt.end() )
    return tablelearnt[index];
  else
    return RecoJet();
}
