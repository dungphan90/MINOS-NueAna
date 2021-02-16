#include "NueAna/NueUtilities.h"  //Header contains functions that are commonly needed
                                  //in macros to process AnaNue files.  For example:
                                  //   chaining files and getting a POT counting
                                  //Contains AnaNueProcessor which tries to encapsulate
                                  //these functions
                                  //
                                  //Header also happens to contains other headers that
                                  //are needed for many of the common NueAna functions

#include "NueAna/NueAnaTools/NueConvention.h" //Header contains function to recalculate E-scale

#include "NueAna/MultiBinAna/SysFileGen2D.h"  //Makes LisaFile Histograms


using namespace std;  //so you don't have to put std:: in front of things like cout

void FillHistos(bool isFDTauNue, string nameForStudy, string sysNames, string filelist, SysFileGen2D *sysFile_Fid, SysFileGen2D *sysFile_Pre);

void fillVector(vector<string> & v_string, string file);



void make_sys_file(
                    string nameForStudy, string sysNamesFile, string fileListsFileND,
                    string fileListsFileFD_f210, string fileListsFileFD_f213, string fileListsFileFD_f214
                   )
{
 vector<string> sysNames;
 fillVector(sysNames, sysNamesFile);
 if(sysNames.size() < 1)
 {
   cout << "\n\nERROR: Need at least one systematic name"<<endl;
   return; 
 }

 vector<string> fileListsND;
 fillVector(fileListsND, fileListsFileND);
 if(fileListsND.size() != sysNames.size())
 {
   cout << "\n\nERROR: fileListsND.size() != sysNames.size()"<<endl;
   return; 
 }

 vector<string> fileListsFD_f210;
 fillVector(fileListsFD_f210, fileListsFileFD_f210);
 if(fileListsFD_f210.size() != sysNames.size())
 {
   cout << "\n\nERROR: fileListsFD_f210.size() != sysNames.size()"<<endl;
   return; 
 }

 vector<string> fileListsFD_f213;
 fillVector(fileListsFD_f213, fileListsFileFD_f213);
 if(fileListsFD_f213.size() != sysNames.size())
 {
   cout << "\n\nERROR: fileListsFD_f213.size() != sysNames.size()"<<endl;
   return; 
 }

 vector<string> fileListsFD_f214;
 fillVector(fileListsFD_f214, fileListsFileFD_f214);
 if(fileListsFD_f214.size() != sysNames.size())
 {
   cout << "\n\nERROR: fileListsFD_f214.size() != sysNames.size()"<<endl;
   return; 
 }

 //Create SysFile
 SysFileGen2D *sysFile_Fid = new SysFileGen2D;
  sysFile_Fid->SetCutLevel(Selection::kFid);
 SysFileGen2D *sysFile_Pre = new SysFileGen2D;
  sysFile_Pre->SetCutLevel(Selection::kPre);
 
 //Add your different sample names:
 for(unsigned int i=0; i<sysNames.size(); ++i)
 {
  sysFile_Fid->AddSystematic(sysNames[i]);
  sysFile_Pre->AddSystematic(sysNames[i]);
 }
 sysFile_Fid->AddPID("ANN11");
 sysFile_Fid->AddPID("ANN14");
 sysFile_Fid->AddPID("PAR");
 sysFile_Fid->AddPID("LEM");

 sysFile_Pre->AddPID("ANN11");
 sysFile_Pre->AddPID("ANN14");
 sysFile_Pre->AddPID("PAR");
 sysFile_Pre->AddPID("LEM");


 //Set oscillation parameters
 NueStandard::SetDefaultOscParam();
 

 for(unsigned int sample=0; sample<sysNames.size(); ++sample)
 {
  cout << "\n\nOn sample: " << sysNames[sample] << endl;
  sysFile_Fid->fCurrentSysName = sysFile_Fid->fSysNames[sample];
  sysFile_Pre->fCurrentSysName = sysFile_Pre->fSysNames[sample];

  sysFile_Fid->ResetHistograms();
  sysFile_Pre->ResetHistograms();
 
  FillHistos(false, nameForStudy, sysNames[sample], fileListsND[sample], sysFile_Fid, sysFile_Pre);
  FillHistos(false, nameForStudy, sysNames[sample], fileListsFD_f210[sample], sysFile_Fid, sysFile_Pre);
  FillHistos(true, nameForStudy, sysNames[sample], fileListsFD_f213[sample], sysFile_Fid, sysFile_Pre);
  FillHistos(true, nameForStudy, sysNames[sample], fileListsFD_f214[sample], sysFile_Fid, sysFile_Pre);

  sysFile_Fid->WriteToFile();
  sysFile_Pre->WriteToFile();

 }//loop over special samples

}//end make_sys_file


void fillVector(vector<string> & v_string, string file)
{
 ifstream is_file(file.c_str());
 if( is_file )
 {
   while( is_file.good() && !is_file.eof() )
   {
     std::string line;
     getline(is_file,line,'\n');
     if ( line.size()==0 ) continue;
     v_string.push_back(line);
   }
 }
}



void FillHistos(bool isFDTauNue, string nameForStudy, string sysNames, string filelist, SysFileGen2D *sysFile_Fid, SysFileGen2D *sysFile_Pre)
{
  if( filelist == "NULL")
  {
    //Intentionally didn't supply files
    //Do nothing
    return;
  }
  
  NueUtilities::AnaNueProcessor anp(filelist);  //Create chain of ana_nue trees, pottrees in the files
  while( anp.NextEntry() ) //loop over events in ana_nue chain
  {
   //Event information is stored in pointer to a NueRecord that is accessible through anp.nueRecord
   //Some event information can be accessed with anp.isMC(), anp.isFar(), anp.isNC(), etc.   
    anp.PrintProgress(1); //cout progress after increments of 1% of the entries have been analyzed
    
    NueConvention::NueEnergyCorrection(anp.nueRecord);
 
    //Set event weight 
    double eventWeight = 1.0;
   
    //Apply MC corrections for flux, etc.
    eventWeight *= NueStandard::GetMCWeights(anp.nueRecord);

    //If your using event weights to evaluate systematics then insert your code here:
      //eventWeight *= blah
    //If your readjusting values:
      //adjust them here
    //Suggestion: key off of the value of sysNames
   
    //normalize POT
      //This normalization is only appropriate if Run 1, 2, and 3 are in separate samples
    eventWeight *= anp.Normalization_RunSeparatedSample();

    //Apply oscillations for Far
    //This method of applying oscillations depends on how you do add your samples together
    //This method works if you keep the f210*, f213*, and f214* files separate
    double oscWeight = 1.0;
    oscWeight *= anp.GetOscWeight_f210f213f214Separate();


    //Apply Fiducial Cuts (which includes Data Quality)
    if(!NueStandard::PassesFarDataTiming(anp.nueRecord)) continue;  //event is within spill window
    if(!NueStandard::PassesSelection(anp.nueRecord, Selection::kFid)) continue;


    Background::Background_t bkg = Background::kUnknown;
    if(anp.isSignal() && isFDTauNue) bkg = Background::kNueCC;
    else if(anp.isTau() && isFDTauNue) bkg = Background::kNuTauCC;
    else if(anp.isCC() && !isFDTauNue) bkg = Background::kNuMuCC;
    else if(anp.isNC() && !isFDTauNue) bkg = Background::kNC;
    else if(anp.isBeamNue() && !isFDTauNue) bkg = Background::kBNueCC;

    //Fill histograms that pass fiducial cuts here
    if(anp.isNear() && bkg != Background::kUnknown)
    {
     sysFile_Fid->fSysHistsMap[bkg]->ND_TrueVsReco->Fill(sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),anp.nueRecord->mctrue.nuEnergy,eventWeight);

     sysFile_Fid->fSysHistsMap[bkg]->ND_RecoVsPID["ANN11"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN2PE_DAIKON04),sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Fid->fSysHistsMap[bkg]->ND_RecoVsPID["ANN14"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN14_DAIKON04),sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Fid->fSysHistsMap[bkg]->ND_RecoVsPID["PAR"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kParticlePID),sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Fid->fSysHistsMap[bkg]->ND_RecoVsPID["LEM"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kMCNN),sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
    }
    else if(anp.isFar() && bkg != Background::kUnknown)
    {
     sysFile_Fid->fSysHistsMap[bkg]->FD_TrueVsReco->Fill(sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),anp.nueRecord->mctrue.nuEnergy,eventWeight);
 
     sysFile_Fid->fSysHistsMap[bkg]->FD_RecoVsPID["ANN11"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN2PE_DAIKON04),sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Fid->fSysHistsMap[bkg]->FD_RecoVsPID["ANN14"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN14_DAIKON04),sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Fid->fSysHistsMap[bkg]->FD_RecoVsPID["PAR"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kParticlePID),sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Fid->fSysHistsMap[bkg]->FD_RecoVsPID["LEM"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kMCNN),sysFile_Fid->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
    }

    //Apply Preselection Cuts (which includes Fiducial and DQ cuts)
    if(!NueStandard::PassesSelection(anp.nueRecord, Selection::kPre)) continue;

    //Fill histograms that pass preselection cuts here
    if(anp.isNear() && bkg != Background::kUnknown)
    {
     sysFile_Pre->fSysHistsMap[bkg]->ND_TrueVsReco->Fill(sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),anp.nueRecord->mctrue.nuEnergy,eventWeight);

     sysFile_Pre->fSysHistsMap[bkg]->ND_RecoVsPID["ANN11"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN2PE_DAIKON04),sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Pre->fSysHistsMap[bkg]->ND_RecoVsPID["ANN14"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN14_DAIKON04),sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Pre->fSysHistsMap[bkg]->ND_RecoVsPID["PAR"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kParticlePID),sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Pre->fSysHistsMap[bkg]->ND_RecoVsPID["LEM"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kMCNN),sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
    }
    else if(anp.isFar() && bkg != Background::kUnknown)
    {
     sysFile_Pre->fSysHistsMap[bkg]->FD_TrueVsReco->Fill(sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),anp.nueRecord->mctrue.nuEnergy,eventWeight);

     sysFile_Pre->fSysHistsMap[bkg]->FD_RecoVsPID["ANN11"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN2PE_DAIKON04),sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Pre->fSysHistsMap[bkg]->FD_RecoVsPID["ANN14"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN14_DAIKON04),sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Pre->fSysHistsMap[bkg]->FD_RecoVsPID["PAR"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kParticlePID),sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
     sysFile_Pre->fSysHistsMap[bkg]->FD_RecoVsPID["LEM"]->Fill(
       NueStandard::GetPIDValue(anp.nueRecord, Selection::kMCNN),sysFile_Pre->GetNueRecoEnergy(anp.nueRecord),oscWeight*eventWeight);
    }
  }//Loop over entries

  static bool firstTime = true;
  if(firstTime)
  {
   //Do this one time
   sysFile_Fid->SetNearPOT(NueStandard::kNormalizedNearPOT*1e12);
   sysFile_Pre->SetNearPOT(NueStandard::kNormalizedNearPOT*1e12);
   string run ="";

   if(anp.didRun1() && !anp.didRun2() && !anp.didRun3())
   {
    sysFile_Pre->SetFarPOT(NueStandard::kNormalizedFarPOT_Run1*1e12);
    sysFile_Fid->SetFarPOT(NueStandard::kNormalizedFarPOT_Run1*1e12);
    run = "1";
   }
   else if(!anp.didRun1() && anp.didRun2() && !anp.didRun3())
   { 
    sysFile_Pre->SetFarPOT(NueStandard::kNormalizedFarPOT_Run2*1e12);
    sysFile_Fid->SetFarPOT(NueStandard::kNormalizedFarPOT_Run2*1e12);
    run = "2";
   }
   else if(!anp.didRun1() && !anp.didRun2() && anp.didRun3())
   {
    sysFile_Pre->SetFarPOT(NueStandard::kNormalizedFarPOT_Run3*1e12);
    sysFile_Fid->SetFarPOT(NueStandard::kNormalizedFarPOT_Run3*1e12);
    run = "3";
   }

   string outfiletag_Fid = "SysFile2D_" + nameForStudy + "_Fid_" + run + ".root";
   string outfiletag_Pre = "SysFile2D_" + nameForStudy + "_Pre_" + run + ".root";
   sysFile_Fid->SetOutputFile(outfiletag_Fid);
   sysFile_Pre->SetOutputFile(outfiletag_Pre);
   
   firstTime = false;
  }//end one time config

}
