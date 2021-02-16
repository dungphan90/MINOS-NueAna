#include "NueAna/NueExpBuilder.h"
#include "NueAna/NueStandard.h"
#include "TRandom3.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "TMath.h"

NueExpBuilder::NueExpBuilder(){ 
  outSet = false;
  fReweight = false;
  fBaseline = 735;
  fDeltaMSquare = 0.0027;
  fUe3Square = 0.025;
  fTheta23 = TMath::Pi()/4;

  fSeed = -10;
  cutSet = "Default";
}

void NueExpBuilder::AddFiles(std::string infiles)
{
   files.push_back(infiles);
}

void NueExpBuilder::SetOutputFile(std::string file)
{
   outSet = true;
   outFile = file;
}   

void NueExpBuilder::SetDeltaMSquare(double dm2) { fDeltaMSquare = dm2; fReweight = true;}
void NueExpBuilder::SetUe3Square(double dUe32) { fUe3Square = dUe32; fReweight = true;}
void NueExpBuilder::SetTheta23(double t23) { fTheta23 = t23; fReweight = true;}
void NueExpBuilder::SetBaseline(double bl) { fBaseline = bl; fReweight = true;}

void NueExpBuilder::SetMaxProb(double max) {fMaxProb = max;}
double NueExpBuilder::GetMaxProb() {return fMaxProb;}

void NueExpBuilder::SetMeanNumberOfEvents(double num) {fMeanNumber = num;}
void NueExpBuilder::SetNumberOfEvents(int num) {fNumberOfEvents = num; }
                                                                                
void NueExpBuilder::GenerateExperiment(int nExp)
{
  NueRecord *nr = new NueRecord();
  NuePOT *np = new NuePOT();
  NuePOT *total = new NuePOT();

  TChain *chain = new TChain("ana_nue");

  for(unsigned int i = 0; i < files.size(); i++)
  {
      chain->Add(files[i].c_str());
  }

  chain->SetBranchAddress("NueRecord",&nr);
  chain->SetBranchStatus("fHeader*",1);
  chain->SetBranchStatus("fHeader.fVld*",1);
  chain->SetBranchStatus("fluxweight*", 1);
  chain->SetBranchStatus("mctrue*", 1);

  chain->SetBranchStatus("dtree*", 0);
  chain->SetBranchStatus("mcnnv*", 0);
  chain->SetBranchStatus("timing*", 0);
  chain->SetBranchStatus("cdi*", 0);
  chain->SetBranchStatus("mri*", 0);
  chain->SetBranchStatus("mri*", 0);
  chain->SetBranchStatus("treepid*", 0);
  chain->SetBranchStatus("highhit*", 0);
  chain->SetBranchStatus("mstvar*", 0);
  chain->SetBranchStatus("fracvars*", 0);
  chain->SetBranchStatus("shield*", 0);
  chain->SetBranchStatus("angcluster*", 0);
  chain->SetBranchStatus("shwfit*", 0);
  chain->SetBranchStatus("hitcalc*", 0);
  chain->SetBranchStatus("anainfo*", 0);
  chain->SetBranchStatus("mda*", 0);
  chain->SetBranchStatus("srshower*", 0);
  chain->SetBranchStatus("srtrack*", 0);
  
  std::vector<double> weights;
  double max = 0;

  Int_t n = chain->GetEntries();

  for(int i = 0; i < n; i++){
   chain->GetEvent(i);
   if(i%10000 == 0) std::cout<<"On call "<<i<<"/"<<n<<std::endl;
   double weight =1.0;

   if(EvaluateCuts(nr)){
     weight *= NueStandard::GetRPWBeamWeight(nr);

     if(fReweight && nr->mctrue.nuFlavor > -100)
     {
        int nuFlavor = nr->mctrue.nuFlavor;
        int  nonOsc = nr->mctrue.nonOscNuFlavor;
        float energy = nr->mctrue.nuEnergy;

        Float_t newWeight = NueConvention::Oscillate(nuFlavor, nonOsc, energy,                                735, fDeltaMSquare, fTheta23, fUe3Square);

        nr->mctrue.Ue3Squared = fUe3Square;
        nr->mctrue.DeltamSquared23 = fDeltaMSquare;
        nr->mctrue.Theta23  = fTheta23;
        nr->mctrue.fOscProb = newWeight;
        weight *= newWeight;
      }else{ weight *= nr->mctrue.fOscProb; }
   }else{weight = 0.0;}

   weights.push_back(weight);
   if(weight > fMaxProb) 
      std::cout<<"MaxWeight too small "<<weight<<"  "<<fMaxProb<<std::endl;
   if(weight > max) max = weight;
  }
  std::cout<<"Maximum weight: "<<max<<std::endl;


  TChain *schain = new TChain("ana_nue");
  schain->SetBranchAddress("NueRecord",&nr);
                                                                        
  TChain *pchain = new TChain("pottree");
  pchain->SetBranchAddress("NuePOT", &np);

  TRandom3 fRand;

  if(fSeed > -1) fRand.SetSeed(fSeed);

 
  for(unsigned int i = 0; i < files.size(); i++)
  {
      pchain->Add(files[i].c_str());
      schain->Add(files[i].c_str());
  }


  for(int i = 0; i < nExp; i++){

    //check if the target number was given
    fNumberOfEvents = fRand.Poisson(fMeanNumber);

    char suffix[10];
    sprintf(suffix, "Exp%d", i);

    std::string outName = outFile + std::string(suffix) + ".root";
                                                                                
    TFile *save = new TFile(outName.c_str(),"RECREATE");
    TTree *tree = new TTree("ana_nue","ana_nue");
    TTree::SetBranchStyle(1);
    TBranch* br = tree->Branch("NueRecord","NueRecord", &nr );
    br->SetAutoDelete(kFALSE);
                                                                        
    int count = 0;
    int calls = 0;
    std::vector<int> pos;

    while( count < fNumberOfEvents){
      int entry = (int) fRand.Uniform(n);
      double weight = weights[entry];
      calls++;

      if(calls%200 == 0) std::cout<<"On call "<<calls<<"/"<<count<<std::endl;

      double prob = fRand.Uniform(fMaxProb);
      if(prob <= weight){
         chain->GetEntry(entry);
         tree->Fill();
         count++;
      }
    }

    std::cout<<"Finished "<<suffix<<" with "<<count <<" events after "<<calls<<std::endl;
                                                                        
    TTree *ptree = new TTree("pottree","pottree");
    TBranch* br2 = ptree->Branch("NuePOT","NuePOT", &total );
    br2->SetAutoDelete(kFALSE);
    pchain->GetEntry(0);  
    total->beamtype = np->beamtype;
    total->pot = fTargetPOT;
    ptree->Fill();
                                                                        
    save->cd();
    tree->Write();
    ptree->Write();
    save->Close();
    
    
//    delete tree;
//    delete ptree;
  }
}

void NueExpBuilder::SetCuts(std::string type, int level)
{
   cutSet = type;
   cutLevel = level;   
}

bool NueExpBuilder::EvaluateCuts(NueRecord *nr)
{
   if(cutSet == "Default") return true;

   if(cutSet == "Fiducial") 
     return NueStandard::IsInFid(nr);

   if(cutSet == "Standard"){
      bool ret;
      switch(cutLevel){
        case 0: return NueStandard::PassesSelection(nr, Selection::kDataQual);
        case 1: return NueStandard::PassesSelection(nr, Selection::kFid);
        case 2: 
             ret = NueStandard::PassesSelection(nr, Selection::kFid);
             ret = ret &&  NueStandard::PassesNonHEPreSelection(nr);
             return ret;
        case 3: return NueStandard::PassesSelection(nr, Selection::kPre);
      }
   }
   if(cutSet == "MRE" || cutSet == "MRCC"){
      bool dq = NueStandard::PassesSelection(nr, Selection::kDataQual);
      bool mrefid = NueStandard::PassesMREFiducial(nr);
      bool mrepre = NueStandard::PassesMREPreSelection(nr) && mrefid && dq;
      bool ret;
      switch(cutLevel){
        case 0: return mrefid && dq;
        case 1: return mrepre;
        case 2: 
          return NueStandard::PassesSelection(nr, Selection::kFid) && mrepre;
        case 3:
             ret = NueStandard::PassesSelection(nr, Selection::kFid);
             ret = ret &&  NueStandard::PassesNonHEPreSelection(nr);
             return ret && mrepre;
        case 4: 
         return NueStandard::PassesSelection(nr, Selection::kPre) && mrepre;
      }
   }
   if(cutSet == "Presel")
      return NueStandard::PassesSelection(nr, Selection::kPre);

   if(cutSet == "PreselNoHE")
      return  NueStandard::PassesNonHEPreSelection(nr);

   if(cutSet == "Systematic"){
      // O - standard systematic envelope
      // 1 - standard systematic w/out HE cut
      if(!NueStandard::IsInFid(nr)) return false;

      switch(cutLevel){
        case 0: return NueStandard::PassesSysPreSelection(nr);
        case 1: return NueStandard::PassesSysPreSelectionNoHE(nr);
      }
   }


   if(cutSet == "CC")
      return NueStandard::PassesSelection(nr, Selection::kCC);




  std::cout<<"Invalid Cut Level for "<<cutSet<<": "<<cutLevel<<std::endl;
  return false;
}
