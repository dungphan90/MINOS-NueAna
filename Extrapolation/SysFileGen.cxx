
// THis will be the actual engine that handles a Full Extrapolation
#include <vector>
#include "NueData.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueAna/Extrapolation/SysFileGen.h"
#include "NueAna/NueStandard.h"

#include "TGraphAsymmErrors.h"

using namespace std;

SysFileGen::SysFileGen()
{
   outFileName = "DefaultOut.root";
   fTargetPOT = 3.25e20;
   fNearPOT = 0.0;
   fFarPOT = 0.0;
}

void SysFileGen::SetNueRecoBins(int nx,Double_t lx,Double_t ux)
{
  fNXBins = nx;
  fXBins = new Double_t[fNXBins+1];
  Float_t bwidth = (ux-lx)/float(fNXBins);
  for(int i=0;i<fNXBins+1;i++) fXBins[i] = lx + float(i)*bwidth;
}
                                                                                
void SysFileGen::SetTrueBins(int ny,Double_t ly,Double_t uy)
{
  fNYBins = ny;
  fYBins = new Double_t[fNYBins+1];
  Float_t bwidth = (uy-ly)/float(fNYBins);
  for(int i=0;i<fNYBins+1;i++) fYBins[i] = ly + float(i)*bwidth;
}
                                                                                
void SysFileGen::Initialize()
{
   this->Init();
}

void SysFileGen::SetOscParams(Double_t theta12,Double_t theta23,Double_t theta13,
                                 Double_t deltam12,Double_t deltam23,Double_t deltaCP,
                                 Int_t massH) 
{   
  fTheta12 = theta12; fTheta23 = theta23; fTheta13 = theta13;
  fDeltaMSq12 = deltam12; fDeltaMSq23 = deltam23;
  fDeltaCP = deltaCP; fMassHierarchy = massH; 
}

void SysFileGen::RunSystematicStudy(Selection::Selection_t sel)
{
   //For a given selection (Presel, ANN, MCNN)
   // there will be a set of "systematics"
   // normal, parameter + 1 sigma, par - 1 sigma
   // For each of these build the necessary histograms then write
   // them to file

   //Because this is much more arbitrary fSystematics is a vector of strings
   // These strings should describe the systematic under consideration

   for(unsigned int j = 0; j < fSystematics.size(); j++){
        std::cout<<"systematic "<<j<<endl;
        fCurrentSys = fSystematics[j];
        PrepareHistograms(sel);
        ResetHistograms();  //Clear out whatever happened last time
        FillSysHistograms();
        WriteToFile(sel);
   }
}

void SysFileGen::PrepareHistograms(Selection::Selection_t /*sel*/)
{
   //This is a function entirely up the the systematic so you need to fill this yourself 
   // build whatever histograms you're heart desires - 
   //   should not impact any of the default data structures
 
 
   return;
}

void SysFileGen::FillSysHistograms()
{
  std::map<Background::Background_t,SysHists*>::iterator mapBeg = fHistMap.begin();
  std::map<Background::Background_t,SysHists*>::iterator mapEnd = fHistMap.end();

  //Now looping over each background assign in the appropriate histograms

  while(mapBeg != mapEnd){
    //     SysHists* local = mapBeg->second;

     // Now here is one of the tricky parts
     // If you are going to overwrite local make sure you are using the same bins
     //*************************************************


     // Do whatever you have to do
 

     //***************************************************


     mapBeg++;
  }    
}

void SysFileGen::SetOutputFile(string name)
{
   outFileName = name;
}


void SysFileGen::WriteToFile(Selection::Selection_t sel)
{
  static TFile* file = 0;  
  static TTree* tree = 0;
  static char selection[256];
  static char sysname[256];

  if(file == 0){
     file = new TFile(outFileName.c_str(),"RECREATE");
     tree = new TTree("energytree","energytree");
     tree->Branch("Selection",selection,"Selection/C");
     tree->Branch("nearPOT",&fNearPOT,"nearPOT/D");
     tree->Branch("farPOT",&fFarPOT,"farPOT/D");
     tree->Branch("SysName", sysname,"SysName/C");
     tree->Branch("Theta12",&fTheta12,"Theta12/D");
     tree->Branch("Theta13",&fTheta12,"Theta13/D");
     tree->Branch("Theta23",&fTheta12,"Theta23/D");
     tree->Branch("DeltaMSq23",&fDeltaMSq23,"DeltaMSq23/D");
     tree->Branch("DeltaMSq12",&fDeltaMSq12,"DeltaMSq12/D");
     tree->Branch("DeltaCP",&fDeltaCP,"DeltaCP/D");
     tree->Branch("MassHierarchy",&fMassHierarchy,"MassHierarchy/I");
  }
  file->cd();
                                                                                
  sprintf(selection,"%s",Selection::AsString(sel));
  sprintf(sysname, "%s", fCurrentSys.c_str());

  std::map<Background::Background_t,SysHists*>::iterator FNbeg = fHistMap.begin();
  std::map<Background::Background_t,SysHists*>::iterator FNend = fHistMap.end();

  while(FNbeg!=FNend) {
    string fnh_name = FNbeg->second->fDirectory->GetName();
    fnh_name += "_" + fCurrentSys + "_" + string(Selection::AsString(sel));

    TDirectory *filedir = file->mkdir(fnh_name.c_str());
    filedir->cd();
    TList *list = FNbeg->second->fDirectory->GetList();
    TIter iter(list->MakeIterator());
    TObject *ob = 0;
    while((ob = iter())) ob->Write();
    file->cd();
    FNbeg++;
  }

  file->cd();

  tree->Fill();

  if(fCurrentSys == fSystematics[fSystematics.size()-1])
     tree->Write();
}

void SysFileGen::InitializeHistograms()
{

  vector<Background::Background_t> bgs;
  bgs.push_back(Background::kNC);
  bgs.push_back(Background::kNuMuCC);
  bgs.push_back(Background::kNuTauCC);
  bgs.push_back(Background::kNueCC);
  bgs.push_back(Background::kBNueCC);
  bgs.push_back(Background::kSelCC);

  for(unsigned int i = 0; i < bgs.size(); i++){
    gDirectory->cd("/");
    string fnh_name = string(Background::AsString(bgs[i]));
    SysHists *fnh = new SysHists(fnh_name.c_str());

    InitializeSysHists(fnh);
 
    (fHistMap)[bgs[i]] = fnh;
  }
}

//Assumes you have already built a SysHists
void SysFileGen::InitializeSysHists(SysHists* one)
{
   gDirectory->cd("/"); 
  
   int nRecoBins = fNXBins; 
   Double_t* RecoBins = fXBins; 
 
   one->fDirectory->cd(); 
   one->fND_RecoEnergy = new TH1D("ND_RecoEnergy","ND Reco Energy",nRecoBins, RecoBins);
   one->fND_RecoEnergy->Sumw2(); 
   one->fFD_RecoEnergy = new TH1D("FD_RecoEnergy","FD Reco Energy",nRecoBins, RecoBins);
   one->fFD_RecoEnergy->Sumw2(); 

   one->fND_RecoVsTrue = new TH2D("ND_RecoVsTrue", "ND Reco v True E", nRecoBins, RecoBins,
                                       fNYBins, fYBins); 
   one->fND_RecoVsTrue->Sumw2();  
 
   one->fFD_RecoVsTrue = new TH2D("FD_RecoVsTrue", "FD Reco v True E", nRecoBins, RecoBins,
                                       fNYBins, fYBins); 
   one->fFD_RecoVsTrue->Sumw2(); 
} 


void SysFileGen::ResetHistograms()
{

  if(fHistMap.size() == 0) InitializeHistograms();

   //Really I should change this to iterate over the fHistMap

  std::map<Background::Background_t, SysHists*>::iterator first = fHistMap.begin();
  std::map<Background::Background_t, SysHists*>::iterator last  = fHistMap.end();

  while(first != last)
  {
    //Background::Background_t bg = first->first;
    SysHists* fnh = first->second;

    fnh->fDirectory->cd();
    fnh->fND_RecoEnergy->Reset("ICE");
    fnh->fFD_RecoEnergy->Reset("ICE");
    fnh->fND_RecoVsTrue->Reset("ICE");
    fnh->fFD_RecoVsTrue->Reset("ICE");

    first++;
  }
}
