
// THis will be the actual engine that handles a Full Extrapolation
#include <vector>
#include "NueData.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueAna/Extrapolation/NueExtrapolationJB.h"
#include "NueAna/NueStandard.h"

#include "TGraphAsymmErrors.h"

using namespace std;

NueExtrapolationJB::NueExtrapolationJB():
 fNZBins(0), fZBins(0)
{
   fCurrentSys = 0;
   outFileName = "DefaultOut.root";
   fTargetPOT = 3.25e20;
   fUseMCAsData = true;
   fUseMCAsCCData = true;
   fNearCCPOT = 0.0;
   fFarCCPOT = 0.0;

   fNearDataPOT = 0.0;
   fNearCCDataPOT = 0.0;

   fExtrapMethod = 1;
}

void NueExtrapolationJB::SetNueRecoBins(int nx,Double_t lx,Double_t ux)
{
  fNXBins = nx;
  fXBins = new Double_t[fNXBins+1];
  Float_t bwidth = (ux-lx)/float(fNXBins);
  for(int i=0;i<fNXBins+1;i++) fXBins[i] = lx + float(i)*bwidth;
}
                                                                                
void NueExtrapolationJB::SetTrueBins(int ny,Double_t ly,Double_t uy)
{
  fNYBins = ny;
  fYBins = new Double_t[fNYBins+1];
  Float_t bwidth = (uy-ly)/float(fNYBins);
  for(int i=0;i<fNYBins+1;i++) fYBins[i] = ly + float(i)*bwidth;
}
                                                                                
void NueExtrapolationJB::SetCCRecoBins(int nz,Double_t lz,Double_t uz)
{
  fNZBins = nz;
  fZBins = new Double_t[fNZBins+1];
  Float_t bwidth = (uz-lz)/float(fNZBins);
  for(int i=0;i<fNZBins+1;i++) fZBins[i] = lz + float(i)*bwidth;
}

void NueExtrapolationJB::Initialize()
{
   this->Init();
}

bool NueExtrapolationJB::LoadFiles()
{
   std::map<TChain*, int>::iterator mapBeg =  fChainDataMap.begin();
   std::map<TChain*, int>::iterator mapEnd =  fChainDataMap.end();

   int chainCount = 0;
   int events = 0;

   while(mapBeg != mapEnd)
   {
      if(!mapBeg->first) cout<<"Error invalid chain"<<endl;
      TChain *chain = mapBeg->first;
      if(strstr(chain->GetName(), "ana_nue")!=0)  this->SetUpNueAnaChain(chain);
      if(strstr(chain->GetName(), "NueMini")!=0)  this->SetUpNueMiniChain(chain);

      int nEntries = chain->GetEntries();
      for(int i=0;i<nEntries;i++){
        if(i%10000==0) std::cout << "Processed "
                           << 100*Float_t(i)/Float_t(nEntries)
                           << "% of Set " << chainCount << std::endl;
         chain->GetEntry(i);
//         if(i > 10000) break;// && 
//            (fData[mapBeg->second]->GetDetector() != 2
//            || !fData[mapBeg->second]->IsNueData())) break; 
 
         if(strstr(chain->GetName(), "ana_nue")!=0){
            NueConvention::NueEnergyCorrection(fRecord);
            fData[mapBeg->second]->AddEvent(fRecord);
         }
         if(strstr(chain->GetName(),"NueMini")!=0)
            fData[mapBeg->second]->AddEvent(fMini);

      }
      mapBeg++;
      events += chain->GetEntries();
      chainCount++;
   }

   for(unsigned int i = 0; i <  fData.size(); i++) {
      if(fData[i]->GetDetector() == Detector::kFar) continue;
      if(fData[i]->IsData() || fUseMCAsData){
         if(fData[i]->IsNueData()) fNearDataPOT   = fData[i]->GetPOT();
         else                      fNearCCDataPOT = fData[i]->GetPOT();
      }
      if(!fUseMCAsData){ 
         if(fData[i]->IsNueData()) fNearDataPOT   = 1e19;
         else                      fNearCCDataPOT = 1e19;
      } 
   }

//   MSG("NueExtrapolationJB", Msg::kInfo) << "Loaded Files "
//      << chainCount<<" Chains yieled "<<events<<" events."<<endl;
   std::cout<<"Finished Loading all files"<<std::endl;

   TFile Input("NueAna/data/xsec_minos_modbyrs4_v3_5_0_mk.root", "READ");

   string names[5] = {"tot", "qe", "res", "dis", "coh"};

   for(int i = 0; i < 5; i++){
     string id = "h_numu_cc_" + names[i];
     Input.GetObject(id.c_str(), fNuMuCCXSec[i]);
     id = "h_nutau_cc_" + names[i];
     Input.GetObject(id.c_str(), fNuTauCCXSec[i]);
     id = "h_nue_cc_" + names[i];
     Input.GetObject(id.c_str(),fNueCCXSec[i]);

     id = "h_numubar_cc_" + names[i];
     Input.GetObject(id.c_str(), fNuMuBarCCXSec[i]);
     id = "h_nutaubar_cc_" + names[i];
     Input.GetObject(id.c_str(), fNuTauBarCCXSec[i]);
     id = "h_nuebar_cc_" + names[i];
     Input.GetObject(id.c_str(),fNueBarCCXSec[i]);

     fNueBarCCXSec[i]->SetDirectory(0);
     fNueCCXSec[i]->SetDirectory(0);
     fNuMuBarCCXSec[i]->SetDirectory(0);
     fNuMuCCXSec[i]->SetDirectory(0);
     fNuTauBarCCXSec[i]->SetDirectory(0);
     fNuTauCCXSec[i]->SetDirectory(0);
   }

   for(int k = 0; k < fNuMuCCXSec[0]->GetNbinsX()+1; k++){
      fXSecPos[k] = fNuMuCCXSec[0]->GetBinLowEdge(k);
   }

   return true;
}

void NueExtrapolationJB::AddChain(TChain *chain, double POT, bool isNue)
{
   if(strstr(chain->GetName(), "ana_nue")!=0)  this->SetUpNueAnaChain(chain);
   if(strstr(chain->GetName(), "NueMini")!=0)  this->SetUpNueMiniChain(chain);

   chain->GetEntry(0);

   ReleaseType::Release_t rel = ReleaseType::kUnknown;
   Detector::Detector_t det = Detector::kUnknown;
   BeamType::BeamType_t beam = BeamType::kUnknown;
 
   if(strstr(chain->GetName(),"ana_nue") != 0){
      rel = fRecord->GetHeader().GetRelease();
      det = fRecord->GetHeader().GetVldContext().GetDetector();
      beam = fRecord->GetHeader().GetBeamType();
   }
   if(strstr(chain->GetName(),"NueMini") !=0){
      rel = fMini->fRelease;
      det = fMini->fDet;
      beam = fMini->fBeam;
   }

   int pos = -1;

   for(unsigned int i = 0; i <  fData.size(); i++) {
      if(fData[i]->GetBeamType() == beam && fData[i]->GetRelease() == rel
           && fData[i]->GetDetector() == det && fData[i]->IsNueData() == isNue)
      {     pos = i; break;  }
   }
   if(pos < 0){
     NueData* dataSet = new NueData(beam, det, rel, isNue);
     fData.push_back(dataSet);
     pos = fData.size() - 1;
   }

   if(det == Detector::kNear && isNue) fNearPOT += POT;
   if(det == Detector::kFar && isNue) fFarPOT += POT;
   if(det == Detector::kNear && !isNue) fNearCCPOT += POT;
   if(det == Detector::kFar && !isNue) fFarCCPOT += POT;   
   
   fData[pos]->AddPOT(POT);

   fChainDataMap[chain] = pos;
}

void NueExtrapolationJB::PrepareExtrapHistograms(Selection::Selection_t sel)
{
   ResetExtrapHistograms();

   Double_t totWeight = 1.0;
   Selection::Selection_t tSel = sel;

   for(unsigned int i = 0; i < fData.size(); i++){
     if(fData[i]->IsData() ) continue;   //Don't need to remake Data Hists
     Detector::Detector_t det = fData[i]->GetDetector();

     NueHeader nh;
     fData[i]->SetupNueHeader(nh); 
     NueRecord* nr = new NueRecord(nh);

     bool isNue = fData[i]->IsNueData();
     tSel = sel;
     if(!isNue)  tSel = Selection::kCC;     

     for(int j = 0; j < fData[i]->NumberOfEntries(); j++)
     {
       fData[i]->FillRecord(nr, j);

       Background::Background_t bg =
          Background::TranslateFromNueClass(nr->mctrue.fNueClass);

       Double_t OscSysVal = fCurrentSys->GetSysValue(Systematic::kOscProb);

       if(fExtrapMethod != 1){         
         //If we are going for speed we still oscillate the NC since it doesn't matter
         // we need the full BNue and NuMu without oscillations to get true vs reco spectrum
         // for the nue and the nutau the oscillations only show up in the selected CC spectrum
         //  or in the true to reco which we should oscillate since it changes the map slightly
         if(bg == Background::kNuMuCC || bg == Background::kBNueCC
            || tSel == Selection::kCC)
         fCurrentSys->SetSysValue(Systematic::kOscProb, 0);
       }                                                                                      
       totWeight = fCurrentSys->UpdateRecord(nr, tSel);

       if(fExtrapMethod != 1){
         // reset osc state
         if(bg == Background::kNuMuCC || bg == Background::kBNueCC
            || tSel == Selection::kCC)
         fCurrentSys->SetSysValue(Systematic::kOscProb, OscSysVal);
       }
       
       double trueEnergy = nr->mctrue.nuEnergy;
       double recoEnergy = this->GetNueEnergy(nr, tSel);

       if(!isNue){
         bg = Background::kSelCC;

         if(nr->mctrue.fNueClass == 1 && det == Detector::kFar){
           fHistMap[bg]->fFD_numu_TrueEnergyBase->Fill(trueEnergy,totWeight);
         }
       }
      
       bool PassCuts = NueStandard::PassesPreSelection(nr) && NueStandard::PassesPIDSelection(nr, tSel);
       if(!isNue) PassCuts = NueStandard::PassesCCSelection(nr);

       if(det == Detector::kFar){
         //nutau and nue osc need the full fiducial voume true to reco spectrum
         //the others want true to reco after they have been selected
         if(bg == Background::kNuTauCC || bg == Background::kNueCC ){
           fHistMap[bg]->fFD_RecoVsTrue->Fill(recoEnergy,trueEnergy,totWeight);
         }else{
           // If we are using the accurate method we don't need to do anything special
           // to make things fast though we separately need to track the numu's that survive
           // and the numu's that appear from bnue -> numu

           if(PassCuts && fExtrapMethod == 1)  
                fHistMap[bg]->fFD_RecoVsTrue->Fill(recoEnergy,trueEnergy,totWeight);
           if(PassCuts && fExtrapMethod != 1)
           {
             if(bg == Background::kNuMuCC)
             {
               bool wasNumu = TMath::Abs(nr->mctrue.nonOscNuFlavor) == 14;
               if(wasNumu)  fHistMap[bg]->fFD_RecoVsTrue->Fill(recoEnergy,trueEnergy,totWeight);
               else fHistMap[bg]->fFD_RecoVsTrue_BNue->Fill(recoEnergy,trueEnergy,totWeight);
             }else
               fHistMap[bg]->fFD_RecoVsTrue->Fill(recoEnergy,trueEnergy,totWeight);
           }
         }
       }
       
       if(PassCuts)
       {
          if(det == Detector::kNear){ 
             fHistMap[bg]->fND_RecoEnergy->Fill(recoEnergy,totWeight);
             fHistMap[bg]->fND_TrueEnergy->Fill(trueEnergy,totWeight);
          }
          if(det == Detector::kFar){
             fHistMap[bg]->fFD_RecoEnergy->Fill(recoEnergy,totWeight);
             fHistMap[bg]->fFD_TrueEnergy->Fill(trueEnergy,totWeight);
          }
       }
     }
     delete nr;
   }

   return;
}

void NueExtrapolationJB::MakeDataHistograms(Selection::Selection_t sel)
{
   // So this isn't correct at least for Near Det
   // There I will be given Histograms for Data
     Double_t totWeight = 1.0;
     Selection::Selection_t tSel = sel;

     NueSystematic *nueSys = new NueSystematic("Standard");
     nueSys->AddSystematic(Systematic::kOscProb,Systematic::GetDefaultValue(Systematic::kOscProb));
     nueSys->SetOscParams(0.554,0.7854,0.,8.2e-5,2.4e-3,0,1);
     nueSys->AddSystematic(Systematic::kSKZP,0);

     for(unsigned int i = 0; i < fData.size(); i++){
       if(fData[i]->IsData()) continue;   //Don't need to deal with Data
  
       Detector::Detector_t det = fData[i]->GetDetector();
//       if(det != Detector::kNear) continue;

       bool isNue = fData[i]->IsNueData();
       tSel = sel;
       if(!isNue)  tSel = Selection::kCC;

       NueHeader nh;
       fData[i]->SetupNueHeader(nh);
       NueRecord* nr = new NueRecord(nh);

       for(int j = 0; j < fData[i]->NumberOfEntries(); j++)
       {
          fData[i]->FillRecord(nr, j);                                                                                                      
          totWeight = nueSys->UpdateRecord(nr, tSel);
          double recoEnergy = this->GetNueEnergy(nr, tSel);           

          Background::Background_t bg =
             Background::TranslateFromNueClass(nr->mctrue.fNueClass);
          if(!isNue)  bg = Background::kSelCC;

          bool PassCuts = NueStandard::PassesPreSelection(nr) && NueStandard::PassesPIDSelection(nr, tSel);
          if(!isNue) PassCuts = NueStandard::PassesCCSelection(nr);
                                                                                                         
          if(PassCuts){
           if(det == Detector::kNear)
              fDataMap[bg]->fND_RecoEnergy->Fill(recoEnergy,totWeight);  
           if(det == Detector::kFar)
              fDataMap[bg]->fFD_RecoEnergy->Fill(recoEnergy,totWeight);
          }
       }
       delete nr;
     }
   
}

void  NueExtrapolationJB::LoadDataHistograms(string /*fname*/, Selection::Selection_t sel)
{
  if(fUseMCAsData)  MakeDataHistograms(sel);
  else{
    TFile *fCC = new TFile("DataHists/NearDetCCData.root", "READ");
    TH1D* numuCC = 0;

    fCC->GetObject("ccdata", numuCC);
    numuCC->SetDirectory(0);

    

    fDataMap[Background::kSelCC]->fND_RecoEnergy = new TH1D("ND_Data_RecoEnergy_SelCC", "ND_Data_RecoEnergy_SelCC", fNZBins, fZBins);
    fDataMap[Background::kSelCC]->fND_RecoEnergy->SetDirectory(0);
 
    int k = 0;
    for(int i = 0; i < fNZBins; i++)
    {
       double pos = fDataMap[Background::kSelCC]->fND_RecoEnergy->GetBinCenter(i);
       while(pos > numuCC->GetBinCenter(k)) k++;
       
       if(pos == numuCC->GetBinCenter(k)){
         fDataMap[Background::kSelCC]->fND_RecoEnergy->SetBinContent(i, numuCC->GetBinContent(k));
         fDataMap[Background::kSelCC]->fND_RecoEnergy->SetBinError(i, numuCC->GetBinError(k));
       }       
    }

    string method = fSeparation;
    TFile *f1 = 0;
    TH1D* nc = 0;
    TH1D* cc = 0;
    TH1D* bnue = 0;

    nc = new TH1D("ND_RecoEnergy","ND Reco Energy",fNXBins, fXBins);
    nc->Sumw2();
    cc = new TH1D("FD_RecoEnergy","FD Reco Energy",fNXBins, fXBins);
    cc->Sumw2();
    bnue = new TH1D("FD_RecoEnergy_BNue","FD Reco Energy",fNXBins, fXBins);
    bnue->Sumw2();
                                                                                                  
    nc->SetDirectory(0);
    cc->SetDirectory(0);
    bnue->SetDirectory(0);

    f1 = new TFile(inDataFileName.c_str());

    if(f1 == 0){
       cout<<"Error loading files"<<endl;
       return;
    }

    TGraphAsymmErrors* INcc;
    TGraphAsymmErrors* INnc;
    TH1D* INbnue;
                                                                                
    string idstring = "TG_NC_" + string(Selection::AsString(sel));
    f1->GetObject(idstring.c_str(), INnc);
    idstring = "TG_CC_" + string(Selection::AsString(sel));
    f1->GetObject(idstring.c_str(), INcc);
    idstring = "MC_BNue_" + string(Selection::AsString(sel));
    f1->GetObject(idstring.c_str(), INbnue);
    double x,y;
    if(INcc == 0 || INnc == 0 || INbnue == 0){
         cout<<"Error loading histograms"<<endl;
    }

    for(int i = 0; i < INnc->GetN(); i++)
    {
       INnc->GetPoint(i,x,y);
       nc->Fill(x,y);
       nc->SetBinError(nc->FindBin(x), INnc->GetErrorY(i));
       INcc->GetPoint(i,x,y);
       cc->Fill(x,y);
       cc->SetBinError(cc->FindBin(x), INcc->GetErrorY(i));
    }
                                                                                
    for(int i = 1; i < 20; i++){
       double pos = INbnue->GetBinCenter(i);
       double val = INbnue->GetBinContent(i);
                                                                                
       bnue->Fill(pos, val);
       bnue->SetBinError(bnue->FindBin(pos), INbnue->GetBinError(i));
    }
                                                                                
    if(nc == 0 || cc == 0 || bnue == 0) cout<<"Failure to load files "<<nc<<"  "<<cc<<"  "<<bnue<<endl;
    nc->SetDirectory(0);
    cc->SetDirectory(0);
    bnue->SetDirectory(0);

    fDataMap[Background::kNuMuCC]->fND_RecoEnergy = (TH1D*) cc->Clone("ND_Data_RecoEnergy_NuMuCC");
    fDataMap[Background::kNuMuCC]->fND_RecoEnergy->SetDirectory(0);

    fDataMap[Background::kNC]->fND_RecoEnergy = (TH1D*) nc->Clone("ND_Data_RecoEnergy_NC");
    fDataMap[Background::kNC]->fND_RecoEnergy->SetDirectory(0);

    fDataMap[Background::kBNueCC]->fND_RecoEnergy = (TH1D*) bnue->Clone("ND_Data_RecoEnergy_BNueCC");
    fDataMap[Background::kBNueCC]->fND_RecoEnergy->SetDirectory(0);

     delete nc;
     delete cc;
     delete bnue;
     delete f1;
     delete fCC;
     delete numuCC;
  }

   TFile* EffInput = new TFile("NueAna/data/tau.root", "READ");;

   string name="tau_eff_"+ string(Selection::AsString(sel));  

   EffInput->GetObject(name.c_str(), fTauEff);
   if(fTauEff == 0) cout<<"error loading efficiency "<<name<<endl; 
   fTauEff->SetDirectory(0);

   delete EffInput;

   EffInput = new TFile("NueAna/data/MREEff.root", "READ");;
                                                                                
   name="eff_"+ string(Selection::AsString(sel));
                                                                                
   EffInput->GetObject(name.c_str(), fMREEff);
   if(fMREEff == 0) cout<<"error loading efficiency "<<name<<endl;
   fMREEff->SetDirectory(0);
                                                                                
   delete EffInput;

   cout<<"Finished Loading the Data"<<endl;
}

void NueExtrapolationJB::MakePrediction()
{
  std::map<Background::Background_t, TH1D*>::iterator mapBeg =  fPredMap.begin();
  std::map<Background::Background_t, TH1D*>::iterator mapEnd =  fPredMap.end();
                                                                                                                                 
  while(mapBeg != mapEnd){
     if(mapBeg->second){
            delete (mapBeg->second);
            mapBeg->second = NULL;
     }
     mapBeg->second = (TH1D*) GetPrediction(mapBeg->first);
     mapBeg++;
  }    
}


TH1* NueExtrapolationJB::GetPrediction(Background::Background_t bg)
{
  TH1D* extrapHist = NULL;
  TString name = "FarPrediction_" + string(Background::AsString(bg));
                                                                                
  //NuMu and NC get a straight F/N extrapolation
  if(bg == Background::kNuMuCC || bg == Background::kNC){
    extrapHist = (TH1D*) fDataMap[bg]->fND_RecoEnergy->Clone(name);
    
    TH1D* fHelperRatio = NULL;
    fHelperRatio = (TH1D*) this->GetFNRatio(bg);

    for(int i=0;i<extrapHist->GetNbinsX()+1;i++){
      Double_t ratio = fHelperRatio->GetBinContent(i);
      extrapHist->SetBinContent(i,extrapHist->GetBinContent(i)*ratio);
      extrapHist->SetBinError(i,extrapHist->GetBinError(i)*ratio);
      if(i > 10) i += 1000;
    }

    delete fHelperRatio;
    if(fNearDataPOT == 0) cout<<"No pot"<<endl;
    extrapHist->Scale(fTargetPOT/fNearDataPOT);
    return extrapHist;
  }

  if(bg == Background::kBNueCC)  // Just return the FarMC
  {
     if(fExtrapMethod == 1){
       extrapHist = (TH1D*) fHistMap[bg]->fFD_RecoEnergy->Clone(name);
     }else{
       TH2D* farHist = (TH2D*) fHistMap[bg]->fFD_RecoVsTrue;
       extrapHist = (TH1D*) CreateOscHist(farHist, 12, 12);
     }                                                                                                               
     extrapHist->Scale(fTargetPOT/fFarPOT);
     return extrapHist;
  }

   //If nue, nutau use more complicated approach:
  if(bg ==Background::kNuTauCC ||  bg ==Background::kNueCC){
     extrapHist = (TH1D*) fHistMap[bg]->fFD_RecoEnergy->Clone(name);
     extrapHist->Reset("ICE");
                                                                                                      
     TH1D* trueHist   = (TH1D*) fHistMap[bg]->fFD_TrueEnergy->Clone("TrueE");
     trueHist->Reset("ICE");
                                                                                                        
     TH2D* TrueToReco = (TH2D*) fHistMap[bg]->fFD_RecoVsTrue->Clone("TtoR");
     for(int j = 0; j < TrueToReco->GetNbinsY()+1; j++){
       double total = 0;
       for(int i = 0; i < TrueToReco->GetNbinsX()+1; i++){
          total += TrueToReco->GetBinContent(i,j);
       }
       if(total > 0){
         for(int i = 0; i < TrueToReco->GetNbinsX()+1; i++){
          TrueToReco->SetBinContent(i,j, TrueToReco->GetBinContent(i,j)/total);
         }
       }
     }
                                                                                                        
     if(fExtrapMethod == 1) this->BuildAppTrueHistExact(bg, trueHist);
     else                   this->BuildAppTrueHistFast(bg, trueHist);

     TH1D* eff = fTauEff;
     if(bg ==Background::kNueCC) eff = fMREEff;
                                                                                                        
     int nBin = eff->GetNbinsX();
     double* pos = new double[nBin];
     double* Eff = new double[nBin];
     for(int k = 0; k < nBin; k++){
        pos[k] = eff->GetBinCenter(k);
        Eff[k] = eff->GetBinContent(k);
     }
                                                                                                        
     for(int k = 0; k < TrueToReco->GetNbinsY(); k++)  //Bins in TrueGeV
     {
        double trueEnt = trueHist->GetBinContent(k);
        for(int m = 0; m < TrueToReco->GetNbinsX()+1; m++){  //bins in Reco E NueGeV
           double recoE = TrueToReco->GetXaxis()->GetBinCenter(m);
           if(recoE > 9){ m += 1000; continue; }
           double TtoR = TrueToReco->GetBinContent(m,k);
           extrapHist->Fill(recoE, trueEnt*TtoR*Eff[m]);
           if(recoE != pos[m]) cout<<"Skew on efficiency "<<m<<"  "<<pos[m]<<"  "<<recoE<<endl;
        }
     } //End of True To Reco

     extrapHist->Scale(fTargetPOT/fFarCCPOT);
     delete [] pos;
     delete [] Eff;
     delete TrueToReco;
     delete trueHist;                                                                                                     
     return extrapHist;
  }
     
  return NULL;
  
}

void NueExtrapolationJB::AddNueSystematic(NueSystematic *nueSys)
{
  fSystematics.push_back(nueSys);
}

void NueExtrapolationJB::RunSystematicStudy(Selection::Selection_t sel)
{
  vector<Selection::Selection_t> temp;
  temp.push_back(sel);
 
  return RunSystematicStudy(temp);
}


void NueExtrapolationJB::RunSystematicStudy(vector<Selection::Selection_t> &sel)
{
   string fname;
   InitializeExtrapHistograms();
   InitializePredictionHistograms();
   LoadFiles();

   bool isNeugen = false;
   for(unsigned int j = 0; j < fSystematics.size(); j++){
     Systematic::Systematic_t sys = Systematic::kKNO;
     if(fSystematics[j]->GetSysValue(sys) != Systematic::GetDefaultValue(sys))
        isNeugen = true;
     sys = Systematic::kMA_QE;
     if(fSystematics[j]->GetSysValue(sys) != Systematic::GetDefaultValue(sys))
        isNeugen = true;
     sys = Systematic::kMA_RES;
     if(fSystematics[j]->GetSysValue(sys) != Systematic::GetDefaultValue(sys))
        isNeugen = true;
   }

   if(isNeugen)  InitializeNeugen();

   for(unsigned int i = 0; i < sel.size(); i++){
     LoadDataHistograms(fname, sel[i]);

     for(unsigned int j = 0; j < fSystematics.size(); j++){
        fCurrentSys = fSystematics[j];
        std::cout<<"systematic "<<j<<endl;
        PrepareExtrapHistograms(sel[i]);
        MakePrediction();
        WriteToFile(sel[i]);
     }
   }
}

void NueExtrapolationJB::RunExtrapStudy(Selection::Selection_t sel)
{
   string fname;
   InitializeExtrapHistograms();
   InitializePredictionHistograms();
   LoadFiles();
   
   LoadDataHistograms(fname, sel);

   for(unsigned int j = 0; j < fSystematics.size(); j++){
      fCurrentSys = fSystematics[j];
      if(fExtrapMethod != 3 || j == 0) PrepareExtrapHistograms(sel);
      if(fExtrapMethod == 3 && j == 0){
        for(unsigned int k = 0; k < fData.size(); k++)   fData[k]->Clear();
      }

      MakePrediction();

      double th13, delta, dm23, dum;
      int mH;

      fCurrentSys->GetOscParams(dum,dum,th13,dum,dm23,delta,mH);
      std::map<Background::Background_t, TH1D*>::iterator mapBeg =  fPredMap.begin();
      std::map<Background::Background_t, TH1D*>::iterator mapEnd =  fPredMap.end();

      cout<<th13<<"  "<<delta/3.1415926<<"  "<<dm23<<"  "<<mH<<"  ";

      double nc, numu, nue, bnue, tau;
      nc = numu= nue = bnue = tau = 0.0;
      double total = 0.0;
      while(mapBeg != mapEnd){
        if(mapBeg->second){
          double sum = mapBeg->second->GetSum();
          if(mapBeg->first == Background::kNC) nc = sum;
          if(mapBeg->first == Background::kNuMuCC) numu = sum;
          if(mapBeg->first == Background::kNuTauCC) tau = sum;
          if(mapBeg->first == Background::kNueCC) nue = sum;
          if(mapBeg->first == Background::kBNueCC) bnue = sum;
          total += sum;
        }
        mapBeg++;
      }
     cout<<nc<<"  "<<numu<<"  "<<bnue<<"  "<<tau<<"  "<<nue<<"  "<<total<<endl;
   }
}


void NueExtrapolationJB::SetOutputFile(string name)
{
   outFileName = name;
}


void NueExtrapolationJB::WriteToFile(Selection::Selection_t sel)
{
  static TFile* file = 0;  
  static TTree* tree = 0;
  static char selection[256];

  if(file == 0){
     file = new TFile(outFileName.c_str(),"RECREATE");
     tree = new TTree("energytree","energytree");
     tree->Branch("Selection",selection,"Selection/C");
     tree->Branch("nearPOT",&fNearPOT,"nearPOT/D");
     tree->Branch("farPOT",&fFarPOT,"farPOT/D");
     tree->Branch("nearCCPOT",&fNearCCPOT,"nearCCPOT/D");
     tree->Branch("farCCPOT",&fFarCCPOT,"farCCPOT/D");
  }
  file->cd();
                                                                                
  sprintf(selection,"%s",Selection::AsString(sel));

  std::map<Background::Background_t,FNHistsM2*>::iterator FNbeg = fHistMap.begin();
  std::map<Background::Background_t,FNHistsM2*>::iterator FNend = fHistMap.end();

  while(FNbeg!=FNend) {
    string fnh_name = FNbeg->second->fDirectory->GetName();
    fnh_name += "_" + string(fCurrentSys->GetName()) + "_" + string(Selection::AsString(sel));

    TDirectory *filedir = file->mkdir(fnh_name.c_str());
    filedir->cd();
    TList *list = FNbeg->second->fDirectory->GetList();
    TIter iter(list->MakeIterator());
    TObject *ob = 0;
    while((ob = iter())) ob->Write();
    file->cd();
    FNbeg++;
  }

  //Now the predictions
  string fnh_name = "Predictions_"+ string(fCurrentSys->GetName()) + "_" + string(Selection::AsString(sel));

  TDirectory* filedir = file->mkdir(fnh_name.c_str());
  filedir->cd();
  std::map<Background::Background_t, TH1D*>::iterator mapBeg =  fPredMap.begin();
  std::map<Background::Background_t, TH1D*>::iterator mapEnd =  fPredMap.end();

  std::map<Background::Background_t, FNHistsM2*>::iterator dataBeg =  fDataMap.begin();
  std::map<Background::Background_t, FNHistsM2*>::iterator dataEnd =  fDataMap.end();
                                                                                
                                                                                
  while(mapBeg != mapEnd){
     if(mapBeg->second) mapBeg->second->Write();
     mapBeg++;
     dataBeg->second->fND_RecoEnergy->Write();
     dataBeg->second->fFD_RecoEnergy->Write();
     dataBeg++;
  }
  file->cd();

  fCurrentSys->MakeBranches(tree);                                                                                
  tree->Fill();

  if(fCurrentSys == fSystematics[fSystematics.size()-1])
     tree->Write();
}

void NueExtrapolationJB::InitializeExtrapHistograms()
{
  Int_t max_bg_index = 0;

  int nRecoBins;
  Double_t* RecoBins;

  while(strcmp(Background::
               AsString(Background::EBackground(max_bg_index)),
               "?Unknown?")!=0) {
    gDirectory->cd("/");
    string fnh_name = string(Background::
                             AsString(Background::EBackground(max_bg_index)));
    FNHistsM2 *fnh = new FNHistsM2(fnh_name.c_str());
 
    nRecoBins = fNXBins;
    RecoBins = fXBins;

    if(Background::EBackground(max_bg_index) == Background::kSelCC){
       nRecoBins = fNZBins;
       RecoBins  = fZBins;
    }

    fnh->fDirectory->cd();
    fnh->fND_RecoEnergy = new TH1D("ND_RecoEnergy","ND Reco Energy",nRecoBins, RecoBins);
    fnh->fND_RecoEnergy->Sumw2();
    fnh->fFD_RecoEnergy = new TH1D("FD_RecoEnergy","FD Reco Energy",nRecoBins, RecoBins);
    fnh->fFD_RecoEnergy->Sumw2();
    fnh->fND_TrueEnergy = new TH1D("ND_TrueEnergy","ND True Energy",fNYBins,fYBins);
    fnh->fND_TrueEnergy->Sumw2();
    fnh->fFD_TrueEnergy = new TH1D("FD_TrueEnergy","FD True Energy",fNYBins,fYBins);
    fnh->fFD_TrueEnergy->Sumw2();
                                                                                
    fnh->fFD_RecoVsTrue = new TH2D("FD_RecoVsTrue", "FD Reco v True E", nRecoBins, RecoBins,
                                       fNYBins, fYBins);
    fnh->fFD_RecoVsTrue->Sumw2();
 
    if(Background::EBackground(max_bg_index) == Background::kSelCC){
      fnh->fFD_numu_TrueEnergyBase = new TH1D("fFD_numu_TrueEnergyBase","FD True Numu in Fid vs. True Energy",fNYBins,fYBins);
      fnh->fFD_numu_TrueEnergyBase->Sumw2();
    }

    if(fExtrapMethod != 1 && Background::EBackground(max_bg_index) == Background::kNuMuCC){
      fnh->fFD_RecoVsTrue_BNue = new TH2D("FD_RecoVSTrue_BNue","FD BNue to Numu Reco vs. True Energy",nRecoBins, RecoBins,fNYBins,fYBins);
      fnh->fFD_RecoVsTrue_BNue->Sumw2();
    }

                                                                                
    (fHistMap)[Background::EBackground(max_bg_index)] = fnh;
    max_bg_index++;
  }
}

void NueExtrapolationJB::InitializePredictionHistograms()
{
  vector<Background::Background_t> bgs;
  bgs.push_back(Background::kNC);
  bgs.push_back(Background::kNuMuCC);
  bgs.push_back(Background::kNuTauCC);
  bgs.push_back(Background::kNueCC);
  bgs.push_back(Background::kBNueCC);
  bgs.push_back(Background::kSelCC);
                                                                                                                                 
  int nRecoBins;
  Double_t* RecoBins;
                                                                                                                                 
  for(unsigned int i = 0; i < bgs.size(); i++)
  {
    string fnh_name = "Data Sample " + string(Background::AsString(bgs[i]));

    nRecoBins = fNXBins;
    RecoBins = fXBins;
     
    TH1D* temp = NULL; //new TH1D(fnh_name.c_str(),"Reco Energy",nRecoBins, RecoBins);
//    temp->Sumw2();

    (fPredMap)[bgs[i]] = temp;

    gDirectory->cd("/");

    FNHistsM2 *fnh = new FNHistsM2(fnh_name.c_str());
                                                                                
    if(bgs[i] == Background::kSelCC){
       nRecoBins = fNZBins;
       RecoBins  = fZBins;
    }

    TString name = "ND_Data_RecoEnergy_" + string(Background::AsString(bgs[i]));

    fnh->fND_RecoEnergy = new TH1D(name,"ND Reco Energy",nRecoBins, RecoBins);
    fnh->fND_RecoEnergy->Sumw2();

    name = "FD_Data_RecoEnergy_" + string(Background::AsString(bgs[i]));
    fnh->fFD_RecoEnergy = new TH1D(name,"FD Reco Energy",nRecoBins, RecoBins);
    fnh->fFD_RecoEnergy->Sumw2();

    fDataMap[bgs[i]] = fnh;
  }
}


void NueExtrapolationJB::ResetExtrapHistograms()
{

  if(fHistMap.size() == 0) InitializeExtrapHistograms();

   //Really I should change this to iterate over the fHistMap

  Int_t max_bg_index = 0;
  while(strcmp(Background::
               AsString(Background::EBackground(max_bg_index)),
               "?Unknown?")!=0) {

    Background::Background_t bg = Background::EBackground(max_bg_index);
    FNHistsM2* fnh = fHistMap[bg];

    fnh->fDirectory->cd();
    fnh->fND_RecoEnergy->Reset("ICE");
    fnh->fFD_RecoEnergy->Reset("ICE");
    fnh->fND_TrueEnergy->Reset("ICE");
    fnh->fFD_TrueEnergy->Reset("ICE");
    fnh->fFD_RecoVsTrue->Reset("ICE");

    if(bg == Background::kSelCC) fnh->fFD_numu_TrueEnergyBase->Reset("ICE");
    if(bg == Background::kNuMuCC && fExtrapMethod != 1) fnh->fFD_RecoVsTrue_BNue->Reset("ICE");                                                                                                        
    max_bg_index++;
  }

}

TH1 *NueExtrapolationJB::GetFNRatio(Background::Background_t bg)
{
  Double_t NearPOT = 0.0;
  Double_t FarPOT = 0.0;

  for(unsigned int i = 0; i < fData.size(); i++){
     if(fData[i]->IsData() ) continue;  
     bool isNue = fData[i]->IsNueData();

     if(bg == Background::kSelCC && isNue) continue;
     else if(!isNue) continue;

     Detector::Detector_t det = fData[i]->GetDetector();
         
     if(det == Detector::kNear) NearPOT = fData[i]->GetPOT();
     if(det == Detector::kFar) FarPOT = fData[i]->GetPOT();
  }

  TH1 *helperHistFD = (TH1*) fHistMap[bg]->fFD_RecoEnergy;
  TH1* farRecoHist = 0;
                                                                                                        
  if(bg == Background::kNuMuCC && fExtrapMethod != 1)
  {
     TH2D* farHist = (TH2D*) fHistMap[bg]->fFD_RecoVsTrue;
     farRecoHist = this->CreateOscHist(farHist, 14, 14);

     farHist = (TH2D*) fHistMap[bg]->fFD_RecoVsTrue_BNue;
     TH1* bnue = this->CreateOscHist(farHist, 12, 14);
     farRecoHist->Add(bnue);
     delete bnue;
                                                                                                        
     helperHistFD = farRecoHist;
  }
                                                                                                        
  TH1 *helperHistND = (TH1*) fHistMap[bg]->fND_RecoEnergy;
  TH1D* fHelperRatio = (TH1D*) helperHistFD->Clone("helperRatio");
  fHelperRatio->SetDirectory(0);
  fHelperRatio->Divide(helperHistND);
  if(farRecoHist != 0) delete farRecoHist;

  if(FarPOT!=0) fHelperRatio->Scale(NearPOT/FarPOT);
  return (TH1*) fHelperRatio;
}

TH1 *NueExtrapolationJB::GetCCRatio(string hist, Background::Background_t /*bg*/)
{
  TH1D *helperHistFD = 0;
  TH1D *helperHistND = 0;

  Double_t NearPOT = 1.0;
  Double_t FarPOT = 1.0;
                                                                                                         
  if(hist == "PE"){  // P/E 
    helperHistFD = (TH1D*) fHistMap[Background::kSelCC]->fFD_numu_TrueEnergyBase;
    helperHistND = (TH1D*) fHistMap[Background::kSelCC]->fFD_TrueEnergy;
  }
 
  if(hist == "NN"){   // NearData / NearMC
    helperHistFD = (TH1D*) fDataMap[Background::kSelCC]->fND_RecoEnergy;
    helperHistND = (TH1D*) fHistMap[Background::kSelCC]->fND_RecoEnergy;
                                                                                                         
    for(unsigned int i = 0; i < fData.size(); i++){
      if(fData[i]->GetDetector() != Detector::kNear) continue;
      if(fData[i]->IsNueData()) continue;
                                                                                                         
      if(fData[i]->IsData() ) continue;   //Don't need to remake Data Hists
                                                                                                         
      FarPOT  = fNearCCDataPOT;
      if(!fData[i]->IsData())                NearPOT = fData[i]->GetPOT();
    }
  }

  TH1D* fHelperRatio = (TH1D*) helperHistFD->Clone("helperRatio");
  fHelperRatio->SetDirectory(0);
  fHelperRatio->Divide(helperHistND);

  if(FarPOT!=0) fHelperRatio->Scale(NearPOT/FarPOT);
  return (TH1*) fHelperRatio;
}

void NueExtrapolationJB::InitializeNeugen()
{
  for(unsigned int i = 0; i < fData.size(); i++){
    if(fData[i]->IsData() ) continue;   //Don't need to remake Data Hists
  
    NueHeader nh;
    fData[i]->SetupNueHeader(nh);
    NueRecord* nr = new NueRecord(nh);

    Systematic::Systematic_t sys = Systematic::kMA_QE;
    double defV = Systematic::GetDefaultValue(sys);

    NueSystematic nuesys("temp");

    for(int j = 0; j < fData[i]->NumberOfEntries(); j++)
    {
      fData[i]->FillRecord(nr, j);
  
      double xsec = nuesys.DoNeugenCalc(nr, sys, defV);
      fData[i]->SetNeugenStdXsec(xsec, j);
    }
    delete nr;
  }
} 


TH1* NueExtrapolationJB::CreateOscHist(TH2* hist, int nonOsc, int nuFlavor)
{
 //create a NueRecord for later use
  NueRecord* nr = 0;
  Background::Background_t bg = Background::TranslateFromMC(1, nuFlavor, nonOsc);

  static int TrueEndBin[10] = {fNYBins,fNYBins,fNYBins,fNYBins,fNYBins, 
			       fNYBins,fNYBins,fNYBins,fNYBins,fNYBins};
  static int RecoEndBin[10] = {10,10,10,10,10, 10,10,10,10,10};

  for(unsigned int i = 0; i < fData.size(); i++){
     if(fData[i]->IsData() ) continue;
     bool isNue = fData[i]->IsNueData();
     if(!isNue) continue;
     Detector::Detector_t det = fData[i]->GetDetector();
     if(det == Detector::kFar){
       NueHeader nh;
       fData[i]->SetupNueHeader(nh);
       nr = new NueRecord(nh);
     }
  }
 
  TH2D* farHist = (TH2D*) hist;
  TH1D* extrapHist = (TH1D*) fHistMap[bg]->fFD_RecoEnergy->Clone("farreco");
  extrapHist->Reset();                                              

  int tempTrue = 0;
  int pos = bg;
  if(nonOsc == 12 && nuFlavor == 14) pos = 8;
 
  for(int i = 0; i < farHist->GetNbinsY()+1; i++){
    double E = farHist->GetYaxis()->GetBinCenter(i);
    double oscWeight = 1.0;
    nr->mctrue.nuEnergy = E;
    
    nr->mctrue.nonOscNuFlavor = nonOsc;
    nr->mctrue.nuFlavor = nuFlavor; 
    oscWeight *= fCurrentSys->GetAppearanceWeight(nr, bg);

    if(i > TrueEndBin[pos]) i = farHist->GetNbinsY()+1;

    for(int j = 0; j < farHist->GetNbinsX()+1; j++){
       double events = farHist->GetBinContent(j,i);
       events *= oscWeight;
       double recoE = farHist->GetXaxis()->GetBinCenter(j);
       extrapHist->Fill(recoE,events);

       if(events > 0) {tempTrue = i;}
       if(j > RecoEndBin[pos]) j = farHist->GetNbinsX()+1;
    }
  }

  if(TrueEndBin[pos] != tempTrue){
     TrueEndBin[pos] = tempTrue; 
//     cout<<"Setting max pos to "<<tempTrue<<" for "<<pos<<endl;
  }

  if(nr != 0) delete nr;
  return (TH1*) extrapHist;
}

void NueExtrapolationJB::BuildAppTrueHistExact(Background::Background_t bg, TH1D* trueHist)
{
  TH1D* NearONear  = (TH1D*) this->GetCCRatio("NN", Background::kSelCC);
  TH1D* POverE     = (TH1D*) this->GetCCRatio("PE", Background::kSelCC);

  trueHist->Reset();
  for(unsigned int i = 0; i < fData.size(); i++){
    //Only Need the CC Far MC spectrum
    if(fData[i]->IsData() ) continue;   //Don't need to remake Data Hists
    if(fData[i]->GetDetector() != Detector::kFar) continue;
    if(fData[i]->IsNueData()) continue;   //Only loop over the Far CC events

    NueHeader nh;
    fData[i]->SetupNueHeader(nh);
    NueRecord* nr = new NueRecord(nh);

    for(int j = 0; j < fData[i]->NumberOfEntries(); j++)
    {
       fData[i]->FillRecord(nr, j);
       double totWeight = 1.0;

       int resCode = nr->mctrue.resonanceCode;
       if(resCode < 1001 || resCode > 1004)  resCode = 0;

       resCode = 0;  //just using the totals for now

       //Reweight immediately for the x-sec
       TH1F *denom =  fNuMuCCXSec[resCode%1000];
       TH1F *num =  fNuTauCCXSec[resCode%1000];

       if(bg == Background::kNueCC)  num =  fNueCCXSec[resCode%1000];

       if(nr->mctrue.nonOscNuFlavor < 0){
          denom =  fNuMuBarCCXSec[resCode%1000];
          num = fNuTauBarCCXSec[resCode%1000];
          if(bg == Background::kNueCC) num =  fNueBarCCXSec[resCode%1000];
       }

       double trueEnergy = nr->mctrue.nuEnergy;

       for(int k = 0; k < num->GetNbinsX(); k++){
          double start = fXSecPos[k];
          double end   = fXSecPos[k+1];
          if(trueEnergy > start && trueEnergy < end)
          {
              float low = denom->GetBinContent(k);
              float high = num->GetBinContent(k);
              if(low > 0) totWeight *= high/low;
              k = num->GetNbinsX();
          }
       }

       // below threshold just move on
       if(totWeight == 0) continue;
       totWeight *= fCurrentSys->UpdateRecord(nr, Selection::kCC, bg);
       if(!NueStandard::PassesCCSelection(nr)) continue;  //Only things passing CC sel
       double recoEnergy = this->GetNueEnergy(nr, Selection::kCC);
       for(int k = 1; k < NearONear->GetNbinsX(); k++){
          if(recoEnergy > fXBins[k-1] && recoEnergy < fXBins[k]){
              totWeight *= NearONear->GetBinContent(k);  k += 1000; 
          }
       }

       for(int k = 1; k < POverE->GetNbinsX(); k++)  //Bins in TrueGeV
       {
          double start = fYBins[k-1];
          double end   = fYBins[k];

          if(trueEnergy > start && trueEnergy < end)
          {
               totWeight *= POverE->GetBinContent(k);
               trueHist->Fill(trueEnergy, totWeight);
               k = POverE->GetNbinsX();
          }
       }
    }// End of loop over entries in chain
    delete nr;
  }// End of loop over all data

  delete NearONear;
  delete POverE;
}

void NueExtrapolationJB::BuildAppTrueHistFast(Background::Background_t bg, TH1D* trueHist)
{
   TH1D* NearONear  = (TH1D*) this->GetCCRatio("NN", Background::kSelCC);
   TH1D* POverE     = (TH1D*) this->GetCCRatio("PE", Background::kSelCC);

   trueHist->Reset();

   TH2D* RecoToTrue = (TH2D*) fHistMap[Background::kSelCC]->fFD_RecoVsTrue->Clone("RtoT");

   NueRecord* nr = 0;
   for(unsigned int i = 0; i < fData.size(); i++){
     //Only Need the CC Far MC spectrum
     if(fData[i]->IsData() ) continue;   //Don't need to remake Data Hists
     if(fData[i]->GetDetector() != Detector::kFar) continue;
     if(fData[i]->IsNueData()) continue;   //Only loop over the Far CC events
                                                                                                        
     NueHeader nh;
     fData[i]->SetupNueHeader(nh);
     nr = new NueRecord(nh);
   }

   for(int i = 0; i < RecoToTrue->GetNbinsX()+1; i++){
      double nOverNWeight = NearONear->GetBinContent(i);
      if(nOverNWeight < 0) cout<<"negative weight"<<endl;
      for(int j = 0; j < RecoToTrue->GetNbinsY()+1; j++){
         double trueE = RecoToTrue->GetYaxis()->GetBinCenter(j);
         double entries = RecoToTrue->GetBinContent(i,j)*nOverNWeight;
         if(entries < 0) cout<<"entries less than zero"<<i<<"  "<<j<<entries<<"  "<<nOverNWeight<<endl;
         trueHist->Fill(trueE,entries);
      }
   }

   Double_t totWeight = 1.0;

   //Reweight immediately for the x-sec
   TH1F *denom =  fNuMuCCXSec[0];
   TH1F *num =  fNuTauCCXSec[0];
   if(bg == Background::kNueCC)  num =  fNueCCXSec[0];

   int xSecInd = 0;

   for(int i = 0; i < trueHist->GetNbinsX()+1; i++){
     double trueE = trueHist->GetBinCenter(i);
     if(trueE < 0) continue;
     totWeight = 1.0;

     //First do the xsec reweighting (which may not be synchronized)
     if(trueE > fXSecPos[fNuMuCCXSec[0]->GetNbinsX()]) continue;
     while(trueE > fXSecPos[xSecInd+1]) xSecInd++;

     double start = fXSecPos[xSecInd];
     double end   = fXSecPos[xSecInd+1];
     if(trueE > start && trueE < end)
     {
        float low = denom->GetBinContent(xSecInd);
        float high = num->GetBinContent(xSecInd);
        if(low > 0) totWeight *= high/low;
     }
     //Next add in the oscillation and tau prod weight
     nr->mctrue.nuEnergy = trueE;
     totWeight *= fCurrentSys->GetAppearanceWeight(nr, bg);
     //Finally add the POverE weight
     totWeight *= POverE->GetBinContent(i);
     trueHist->SetBinContent(i, trueHist->GetBinContent(i)*totWeight);
   }// End of loop over all data

                                                                                                       
  delete NearONear;
  delete POverE;
  delete RecoToTrue;
  delete nr;

}
