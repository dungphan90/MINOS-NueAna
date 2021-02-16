////////////////////////////////////////////////////////////////////////
// $Id: NueUtilities.cxx,v 1.17 2015/12/22 05:36:18 wingmc Exp $
//
// The NueUtilities namespace is intended to contain common functions
//  and utilities that analyzers frequently find themselves implementing 
//  For example:
//     a)Chaining files together
//     b)Getting the total POT for the chained files
//     c)Looping over events in a chain
//     d)Making histograms
//
//     Examples a, b, and c are encapsulated by the AnaNueProcessor class
//     Example d is encapsulated by the Histograms class

//  Note: There is an overhead associated with copying a AnaNueProcessor object
//  Hence, it is best to pass it by reference to avoid making a copy
//
//  Note: You cannot copying a Histogram object, the copy and assignment
//  constructors are private.  If you have to pass it, then pass by reference
//
// Author: Gregory Pawloski 
// Created: August 28, 2009
////////////////////////////////////////////////////////////////////////

#include "NueAna/NueUtilities.h"

#include "MessageService/MsgService.h"

CVSID("");

 TChain* NueUtilities::CreateChain(std::string treeName,std::string fileListName)
 {
  //This function can take:
  //1) a single root file which must end with ".root"
  //2) a collection of root files using a wild card
  //   for example: "MyAnaNueFiles/AnaNue*"
  //3) a ASCII file with a list of AnaNue file names
  
  TChain *chain= new TChain(treeName.c_str());
  
  if(fileListName.find("*")!=string::npos)
  {
   //Filenames using a wildcard:
   chain->Add(fileListName.c_str());
  }
  else if( ".root" == fileListName.substr(fileListName.length()-5,5) )
  {
   //Filenames using a wildcard:
   chain->Add(fileListName.c_str());
  }
  else
  {  
   ifstream fileList(fileListName.c_str());
   if( fileList )
   {
    //Must be an ASCII file with a list of AnaNue file names
    while( fileList.good() && !fileList.eof() )
    {
     std::string line;
     getline(fileList,line,'\n');
     if ( line.size()==0 ) continue;
     chain->Add(line.c_str());
    }
   }
  }
  return chain;
 }


 Double_t NueUtilities::GetChainPOT(TChain *chain)
 {
  NuePOT *nuepot = NULL;
  chain->SetBranchAddress("NuePOT",&nuepot);
 
  Double_t pot = 0;
  for( int i = 0; i<chain->GetEntries(); ++i )
  {
    chain->GetEntry(i);
    pot += nuepot->pot;
  }
  chain->ResetBranchAddresses();

  return(pot);
 }

 Double_t NueUtilities::GetWeightedChainPOT(TChain *chain)
 {
  NuePOT *nuepot = NULL;
  chain->SetBranchAddress("NuePOT",&nuepot);
 
  Double_t pot = 0;
  for( int i = 0; i<chain->GetEntries(); ++i )
  {
    chain->GetEntry(i);
    switch(nuepot->beamtype){
      case BeamType::kL010z185i_i124: pot += 1.24932*nuepot->pot; break;
      case BeamType::kL010z185i_i191: pot += 1.62047*nuepot->pot; break;
      case BeamType::kL010z185i_i213: pot += 1.19021*nuepot->pot; break;
      case BeamType::kL010z185i_i224: pot += 1.69376*nuepot->pot; break;
      case BeamType::kL010z185i_i232: pot += 1.28621*nuepot->pot; break;
      case BeamType::kL010z185i_i243: pot += 1.26755*nuepot->pot; break;
      case BeamType::kL010z185i_i257: pot += 1.07360*nuepot->pot; break;
      case BeamType::kL010z185i_i282: pot += 1.11812*nuepot->pot; break;
      case BeamType::kL010z185i_i303: pot += 1.06092*nuepot->pot; break;
      case BeamType::kL010z185i_i324: pot += 2.60657*nuepot->pot; break;
      default:                        pot += nuepot->pot;         break;
    }
  }
  chain->ResetBranchAddresses();

  return(pot);
 }


 NueUtilities::AnaNueProcessor::AnaNueProcessor()
 {
   filelist = "";
   ana_nue_chain = NULL;
   pottree_chain = NULL;

   totalPOTNotFilled=true;

   nueRecord = NULL;
   numEntries = 0;
   currentEntry = -1;
   lastfraction = -1;

   ProcessedMC = false;
   ProcessedData = false;
   ProcessedFar = false;
   ProcessedNear = false;
   ProcessedNC = false;
   ProcessedCC = false;
   ProcessedBeamNue = false;
   ProcessedTau = false;
   ProcessedSignal = false;
   ProcessedRun1 = false;
   ProcessedRun2 = false;
   ProcessedRun3 = false;

   ProcessedRun11 = false;
   ProcessedRun12 = false;
   ProcessedRun13 = false;

 }

 NueUtilities::AnaNueProcessor::AnaNueProcessor(std::string fileListName)
 {
   filelist = fileListName;
   ana_nue_chain= NueUtilities::CreateChain("ana_nue",fileListName);
   pottree_chain= NueUtilities::CreateChain("pottree",fileListName);

   totalPOTNotFilled=true;

   nueRecord = NULL;
   ana_nue_chain->SetBranchAddress("NueRecord",&nueRecord);
   
   numEntries = ana_nue_chain->GetEntries();
   currentEntry = -1;
   lastfraction = -1;

   ProcessedMC = false;
   ProcessedData = false;
   ProcessedFar = false;
   ProcessedNear = false;
   ProcessedNC = false;
   ProcessedCC = false;
   ProcessedBeamNue = false;
   ProcessedTau = false;
   ProcessedSignal = false;
   ProcessedRun1 = false;
   ProcessedRun2 = false;
   ProcessedRun3 = false;

   ProcessedRun11 = false;
   ProcessedRun12 = false;
   ProcessedRun13 = false;

 }


 //Copy constructor
 NueUtilities::AnaNueProcessor::AnaNueProcessor(const NueUtilities::AnaNueProcessor &original)
 {
   *this = original;
 }

 NueUtilities::AnaNueProcessor & NueUtilities::AnaNueProcessor::operator=(const NueUtilities::AnaNueProcessor &original)   
 {  
  // check for "self assignment" and do nothing in that case
  if (this != &original)
  {
   filelist = original.GetFileListName();

   ana_nue_chain= NueUtilities::CreateChain("ana_nue",filelist);
   pottree_chain= NueUtilities::CreateChain("pottree",filelist);

   totalPOTNotFilled=true;

   nueRecord = NULL;
   ana_nue_chain->SetBranchAddress("NueRecord",&nueRecord);

   ProcessedMC = original.didMC();
   ProcessedData = original.didData();
   ProcessedFar = original.didFar();
   ProcessedNear = original.didNear();
   ProcessedNC = original.didNC();
   ProcessedCC = original.didCC();
   ProcessedBeamNue = original.didBeamNue();
   ProcessedTau = original.didTau();
   ProcessedSignal = original.didSignal();
   ProcessedRun1 = original.didRun1();
   ProcessedRun2 = original.didRun2();
   ProcessedRun3 = original.didRun3();

   ProcessedRun11 = original.didRun11();
   ProcessedRun12 = original.didRun12();
   ProcessedRun13 = original.didRun13();

   
   numEntries = ana_nue_chain->GetEntries();
   lastfraction = original.GetLastFraction();

   this->SetEntry(original.CurrentEntry());
  }// if not itself
  return *this;

 }
 


 NueUtilities::AnaNueProcessor::~AnaNueProcessor()
 {
   delete ana_nue_chain;
   delete pottree_chain;

 }

 Double_t NueUtilities::AnaNueProcessor::GetTotalPOT()
 {
  if(totalPOTNotFilled)
  {
   totalPOT = NueUtilities::GetChainPOT(pottree_chain);
   totalPOTNotFilled = false;
  }
  return (totalPOT);
  
 }

 Bool_t NueUtilities::AnaNueProcessor::NextEntry()
 {
  ++currentEntry;
  return (SetEntry(currentEntry));
 }

 Bool_t NueUtilities::AnaNueProcessor::SetEntry(Long64_t entry)
 {
  currentEntry = entry;
  if(currentEntry <0 || currentEntry >= numEntries)
  {
    //nueRecord= NULL; //just don't do anything
    return(false);
  }
  ana_nue_chain->GetEntry(currentEntry);
   ProcessedMC |= isMC();
   ProcessedData |= isData();
   ProcessedFar |= isFar();
   ProcessedNear |= isNear();
   ProcessedNC |= isNC();
   ProcessedCC |= isCC();
   ProcessedBeamNue |= isBeamNue();
   ProcessedTau |= isTau();
   ProcessedSignal |= isSignal();
   ProcessedRun1 |= isRun1();
   ProcessedRun2 |= isRun2();
   ProcessedRun3 |= isRun3();

   ProcessedRun11 |= isRun11();
   ProcessedRun12 |= isRun12();
   ProcessedRun13 |= isRun13();

  //Could of had this function return the result of GetEntry(currentEntry),
  //but that would just cause loops to end and people would probably think that
  //they just hit the end of the entries. However, at this point in the code if
  //GetEntry(currentEntry) returns a nonpositve value, then something wrong happend
  //and I want the program to crash

  return(true);
  
 }

 void NueUtilities::AnaNueProcessor::PrintProgress(unsigned int percentInterval)
 {
   Int_t fraction = (100*currentEntry)/numEntries;
   
   if ( (fraction-lastfraction) >= (Int_t)percentInterval )
   {
    std::cout<<"\n"<<fraction<<"% done"<<std::endl;
    lastfraction = fraction;
   }

   if( currentEntry == (numEntries-1) ) 
   {
    std::cout<<"\nOn last entry"<<std::endl;
    lastfraction = fraction;
   }
 }

 void NueUtilities::AnaNueProcessor::PrintFileNames()
 {
  Long64_t num = 0;
  std::cout << "\nProcessing the following files:";
  ifstream fileList(filelist.c_str());
  if( fileList )
  {
   while( fileList.good() && !fileList.eof() )
   {
     std::string line;
     getline(fileList,line,'\n');
     if ( line.size()==0 ) continue;
     std::cout << "\n " << line;
     ++num;
   }
  }
  std::cout << "\nTotal number of files = " << num << endl;
 }


 Bool_t NueUtilities::AnaNueProcessor::isMC()
 {
  return(nueRecord->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC);
 }

 Bool_t NueUtilities::AnaNueProcessor::isData()
 {
  return(nueRecord->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kData);
 }

 Bool_t NueUtilities::AnaNueProcessor::isFar()
 {
  return(nueRecord->GetHeader().GetVldContext().GetDetector() == Detector::kFar);
 }

 Bool_t NueUtilities::AnaNueProcessor::isNear()
 {
  return(nueRecord->GetHeader().GetVldContext().GetDetector() == Detector::kNear);
 }

 Bool_t NueUtilities::AnaNueProcessor::isNC()
 {
  return(isMC() && nueRecord->mctrue.interactionType==0);
 }

 Bool_t NueUtilities::AnaNueProcessor::isCC()
 {
  return(isMC() && !isNC() && abs(nueRecord->mctrue.nuFlavor)==14);
 }

 Bool_t NueUtilities::AnaNueProcessor::isBeamNue()
 {
  return(isMC() && !isNC() && abs(nueRecord->mctrue.nuFlavor)==12 && abs(nueRecord->mctrue.nonOscNuFlavor)==12);
 }

 Bool_t NueUtilities::AnaNueProcessor::isTau()
 {
  return(isMC() && !isNC() && abs(nueRecord->mctrue.nuFlavor)==16);
 }

 Bool_t NueUtilities::AnaNueProcessor::isSignal()
 {
  return( isMC() && !isNC() && abs(nueRecord->mctrue.nuFlavor)==12 && abs(nueRecord->mctrue.nonOscNuFlavor)==14 );
 }

 Bool_t NueUtilities::AnaNueProcessor::isNu()
 {
  return( isMC() && nueRecord->mctrue.nuFlavor > 0 );
 }

 Bool_t NueUtilities::AnaNueProcessor::isNuBar()
 {
  return( isMC() && nueRecord->mctrue.nuFlavor < 0 && nueRecord->mctrue.nuFlavor > -9000 );
 }

 Bool_t NueUtilities::AnaNueProcessor::isQuasiElastic()
 {
  return( isCC() && nueRecord->mctrue.resonanceCode==1001);
 }



//FHC Run Calls
 Bool_t NueUtilities::AnaNueProcessor::isRun1()
 {
  return( NueStandard::IsRun1(nueRecord) );
 }

 Bool_t NueUtilities::AnaNueProcessor::isRun2()
 {
  return( NueStandard::IsRun2(nueRecord) );
 }

 Bool_t NueUtilities::AnaNueProcessor::isRun3()
 {
  return( NueStandard::IsRun3(nueRecord) );
 }

//RHC Run Calls
Bool_t NueUtilities::AnaNueProcessor::isRun4RHC()
{
  return( NueStandard::IsRun4RHC(nueRecord) );
}

Bool_t NueUtilities::AnaNueProcessor::isRun7RHC()
{
  return( NueStandard::IsRun7RHC(nueRecord) );
}
Bool_t NueUtilities::AnaNueProcessor::isRun9RHC()
{
  return( NueStandard::IsRun9RHC(nueRecord) );
}

//FHC MINOS+
 Bool_t NueUtilities::AnaNueProcessor::isRun11()
 {
  return( NueStandard::IsRun11(nueRecord) );
 }

 Bool_t NueUtilities::AnaNueProcessor::isRun12()
 {
  return( NueStandard::IsRun12(nueRecord) );
 }

 Bool_t NueUtilities::AnaNueProcessor::isRun13()
 {
  return( NueStandard::IsRun13(nueRecord) );
 }



//FHC POT Weighting
 Double_t NueUtilities::AnaNueProcessor::Normalization_RunSeparatedSample()
 {
  //normalize POT assuming sample only contains events from one run period
   Double_t eventWeight =1;
   if(isNear())
   {
     //Normalize Near Data and MC
      eventWeight = NueStandard::kNormalizedNearPOT/GetTotalPOT();
   }
   else if(isMC())
   {
     //Only Normalize Far MC
     //This normalization is only appropriate if Run 1, 2, and 3 are in separate samples
     if( NueStandard::IsRun1(nueRecord) )      eventWeight = NueStandard::kNormalizedFarPOT_Run1/GetTotalPOT();
     else if( NueStandard::IsRun2(nueRecord) ) eventWeight = NueStandard::kNormalizedFarPOT_Run2/GetTotalPOT();
     else if( NueStandard::IsRun3(nueRecord) ) eventWeight = NueStandard::kNormalizedFarPOT_Run3/GetTotalPOT();
     //Adam is going to cheat and add MINOS+ functions to this nonsense
     else if( NueStandard::IsRun11(nueRecord) ) eventWeight = NueStandard::kNormalizedFarPOT_Run11/GetTotalPOT();
     else if( NueStandard::IsRun12(nueRecord) ) eventWeight = NueStandard::kNormalizedFarPOT_Run12/GetTotalPOT();
     else if( NueStandard::IsRun13(nueRecord) ) eventWeight = NueStandard::kNormalizedFarPOT_Run13/GetTotalPOT();
   }
    
   return(eventWeight);
 }
 
 Double_t NueUtilities::AnaNueProcessor::Normalization_MixedRunSample()
 {
  //normalize POT assuming sample only contains mixture of events from different run periods
   Double_t eventWeight = 1;
   if(isNear())
   {
     //Normalize Near Data and MC
      eventWeight = NueStandard::kNormalizedNearPOT/GetTotalPOT();
   }
   else if(isMC())
   {
     //Only Normalize Far MC
     //This normalization is only appropriate if runs mixed -- and is thoroughly inappropriate for MINOS+
     eventWeight = NueStandard::kNormalizedFarPOT/GetTotalPOT();
   }
    
   return(eventWeight);
 }


//RHC POT Weighting
Double_t NueUtilities::AnaNueProcessor::Normalization_RunSeparatedSampleRHC()
{
  //normalize POT assuming sample contains events from only one run period
  Double_t eventWeight = 1;
  if(isNear())
    {
      eventWeight = NueStandard::kNormalizedNearPOT/GetTotalPOT();
    }
  else if(isMC())
    {
      if( NueStandard::IsRun4RHC(nueRecord) ) eventWeight = NueStandard::kNormalizedFarPOT_Run4RHC/GetTotalPOT();
      else if( NueStandard::IsRun7RHC(nueRecord) ) eventWeight = NueStandard::kNormalizedFarPOT_Run7RHC/GetTotalPOT();
    }

  return(eventWeight);
}

Double_t NueUtilities::AnaNueProcessor::Normalization_MixedRunSampleRHC()
{
  //normalize POT assuming mixture of events from different run periods
  Double_t eventWeight = 1;
  if(isNear())
    {
      eventWeight = NueStandard::kNormalizedNearPOT/GetTotalPOT();
    }
  else if(isMC())
    {
      eventWeight = NueStandard::kNormalizedFarPOT_RHC/GetTotalPOT();
    }

  return(eventWeight);
}

//end RHC POT Weighting



Double_t NueUtilities::AnaNueProcessor::GetNDDataWeightedPOT(Bool_t mixed, unsigned int percentInterval){
  Double_t chainpot = 0;
  lastfraction = 0;
  SetEntry(0);
  Bool_t pass = true;
  if( ( !isRun3() && !mixed && !isRun4RHC() && !isRun7RHC() && !isRun9RHC() ) || isMC() ){
    chainpot = GetTotalPOT();
    pass = false;
    cout << "===========> Not Run3 Data! Getting Standard POT <===========" << endl;
  }
  else cout << "===========> Beginning POT weight loop <===========" << endl;
  SetEntry(-1);
  while( NextEntry() && pass ) //loop over events in ana_nue chain
  {
    PrintProgress(percentInterval);

    if(NueStandard::PassesPOTStandards(nueRecord)){
      if(nueRecord->GetHeader().GetEventNo() == 0 || nueRecord->GetHeader().GetEvents()==0)
      {
        double dataweight = NueStandard::GetNDDataWeights(nueRecord);
        chainpot += dataweight*nueRecord->bmon.bI;
        SetEntry(CurrentEntry()+TMath::Max(0,nueRecord->GetHeader().GetEvents()-1));
      }
    }
  }
  lastfraction = 0;
  SetEntry(-1);
  return chainpot;
}

 Double_t NueUtilities::AnaNueProcessor::GetOscWeight_f210f213f214Separate()
 {
  //This method of applying oscillations depends on how you do add your samples together
  //This method works if you keep the f210*, f213*, and f214* files separate
   Double_t eventWeight =1;

   if( isMC() && isFar() && !isNC() )
   {
    eventWeight = NueStandard::GetOscWeight(nueRecord->mctrue.nuFlavor, nueRecord->mctrue.nonOscNuFlavor, nueRecord->mctrue.nuEnergy);
   }

   return(eventWeight);
 }


 Double_t NueUtilities::AnaNueProcessor::GetRecoE()
 {
   return(nueRecord->srevent.phNueGeV);
 }


 Double_t NueUtilities::AnaNueProcessor::GetTrueE()
 {
   return(nueRecord->mctrue.nuEnergy);
 }


//------------------------------------------------------------------------------
//
//class: Histograms
//
//A class that is used to handle histograms which includes:
//  Opening a ROOT file to store them in
//  Initializing whatever histograms that you are interested in
//  Filling the histograms
//  Writing the histograms
//  Closing the file 
//
//This is an attempt to avoid rewriting lines of code for multiple histograms
//that are very similar.  This hopefully avoids copy and paste errors in that case
//The exact histograms are defined in another class that inherits from this one
//
//Note don't make copies of an object of this class (why would you want to)
//Only pass by reference if you find that you have too
//------------------------------------------------------------------------------


TH1* NueUtilities::Histograms::GetHistogram(TString setName,TString subSetName,TString name)
{
 for( unsigned int dimensions=1; dimensions <= 3; ++dimensions)
 {
  if( _histos.find(dimensions) != _histos.end() )
  {
   if( _histos[dimensions].find(setName) != _histos[dimensions].end() )
   {
    if( _histos[dimensions][setName].find(subSetName) != _histos[dimensions][setName].end() )
    {
     if( _histos[dimensions][setName][subSetName].find(name) != _histos[dimensions][setName][subSetName].end() )
     {
       return( _histos[dimensions][setName][subSetName][name] );
     }//Does name exist
    }//Does subSetName exist
   }//Does setName exist
  }//Does Dimension exist
 }

 cout << "\n\n ERROR: Cannot find histogram for:"<< endl;
 cout << " "<< setName << " " << subSetName << " " << name << endl << endl;
 exit(1);
 return(NULL);
}


////////////////////////////
//Methods to Make Histograms
////////////////////////////


//Make 1D Histogram
void NueUtilities::Histograms::MakeHisto(TString setName, TString subSetName, TString name,
                                         int numBinsX, double lowValueX, double highValueX)
{
 MakeHisto(setName,subSetName,name,
           numBinsX, lowValueX, highValueX,
           0,0,0,
           0,0,0,
           1);
}

//Make 2D Histogram
void NueUtilities::Histograms::MakeHisto(TString setName, TString subSetName, TString name,
                                         int numBinsX, double lowValueX, double highValueX,
                                         int numBinsY, double lowValueY, double highValueY)
{
 MakeHisto(setName,subSetName,name,
           numBinsX, lowValueX, highValueX,
           numBinsY, lowValueY, highValueY,
           0,0,0,
           2);
}

//Make 3D Histogram
void NueUtilities::Histograms::MakeHisto(TString setName, TString subSetName, TString name,
                                         int numBinsX, double lowValueX, double highValueX,
                                         int numBinsY, double lowValueY, double highValueY,
                                         int numBinsZ, double lowValueZ, double highValueZ)
{
 MakeHisto(setName,subSetName,name,
           numBinsX, lowValueX, highValueX,
           numBinsY, lowValueY, highValueY,
           numBinsZ, lowValueZ, highValueZ,
           3);
}

//Make Histograms called by protected MakeHisto member functions
void NueUtilities::Histograms::MakeHisto(TString setName,TString subSetName,TString name,
                                         int numBinsX, double lowValueX, double highValueX,
                                         int numBinsY, double lowValueY, double highValueY,
                                         int numBinsZ, double lowValueZ, double highValueZ,
                                         int dimensions)
{
 if( _histos[dimensions][setName][subSetName].find(name) != _histos[dimensions][setName][subSetName].end() )
 {
  cout << "\n\nError overwriting " << "th"<< dimensions << "d_"+setName+"_"+subSetName+"_"+name << endl;
  delete (_histos[dimensions][setName][subSetName][name]);
  exit(1);
 }


  TString currentDir = gDirectory->GetPath();
  _file->cd();

  if( 1 == dimensions)
  {
           _histos[dimensions][setName][subSetName][name] = new TH1D("th1d_"+setName+"_"+subSetName+"_"+name,"",
                                                            numBinsX,lowValueX,highValueX);
           ((TH1D*)_histos[dimensions][setName][subSetName][name])->Sumw2();  //Always do a Sumw2() to get errors right
  }
  else if( 2 == dimensions)
  {
           _histos[dimensions][setName][subSetName][name] = new TH2D("th2d_"+setName+"_"+subSetName+"_"+name,"",
                                                            numBinsX,lowValueX,highValueX, numBinsY,lowValueY,highValueY);
           ((TH2D*)_histos[dimensions][setName][subSetName][name])->Sumw2();  //Always do a Sumw2() to get errors right
  }
  else if( 3 == dimensions)
  {
           _histos[dimensions][setName][subSetName][name] = new TH3D("th3d_"+setName+"_"+subSetName+"_"+name,"",
                                                            numBinsX,lowValueX,highValueX, numBinsY,lowValueY,highValueY, numBinsZ,lowValueZ,highValueZ);
           ((TH3D*)_histos[dimensions][setName][subSetName][name])->Sumw2();  //Always do a Sumw2() to get errors right
  }
  else
  {
           cout << "\n\nERROR: Don't know how to make histogram of " << dimensions << " dimensions\n"<<endl;
           gROOT->cd(currentDir);
           exit(1);
  }

  gROOT->cd(currentDir);

  if( !DoesSetExist(setName) ) _setNames.push_back(setName);

}


//Make Histograms by Retrieving Histograms in File
//Function assumes Set/Subset/histo file structure
void NueUtilities::Histograms::AutoRetrieveHistograms()
{
  TString currentDir = gDirectory->GetPath();
  _file->cd();

  // loop over all keys in this directory
  TIter nextSetKey( _file->GetListOfKeys() );
  TKey *setKey;
  while((setKey = (TKey*)nextSetKey()))
  {
    TObject *set = setKey->ReadObj();
    if( set->IsA()->InheritsFrom( "TDirectory" ) )
    {
      TString setName = set->GetName();
      _file->cd(setName);
      TIter nextSubsetKey( gDirectory->GetListOfKeys() );
      TKey *subsetKey;
      while((subsetKey = (TKey*)nextSubsetKey()))
      {
        TObject *subset = subsetKey->ReadObj();
        if( subset->IsA()->InheritsFrom( "TDirectory" ) )
        {
          TString subsetName = subset->GetName();
          _file->cd(setName+"/"+subsetName);
          TIter nextHistoKey( gDirectory->GetListOfKeys() );
          TKey *histoKey;
          while((histoKey = (TKey*)nextHistoKey()))
          {
            TObject *histo = histoKey->ReadObj();
            if( histo->IsA()->InheritsFrom( "TH1" ) )
            {
              TString histoName = histo->GetName();
              int dimensions = 1;
              if( histo->IsA()->InheritsFrom( "TH3" ) ) dimensions = 3;
              else if( histo->IsA()->InheritsFrom( "TH2" ) ) dimensions = 2;
              RetrieveHistogram(setName,subsetName,histoName,dimensions);
            }//If found directory subset
          } // while subset keys
          _file->cd(setName);  //return to set dir
        }//If found directory subset
      } // while subset keys
     _file->cd();  //return to top dir
    }//If found directory set
  } // while set keys

  //return to original directory
  gROOT->cd(currentDir);
}


void NueUtilities::Histograms::RetrieveHistogram(TString setName, TString subSetName, TString histoName, unsigned int dimensions)
{
 //remove "thnd_" from the histoname
 TString name = histoName(5,histoName.Length()-5);
 if( _histos[dimensions][setName][subSetName].find(name) != _histos[dimensions][setName][subSetName].end() )
 {
  cout << "\n\nError overwriting " << "th"<< dimensions << "d_"+setName+"_"+subSetName+"_"+name << endl;
  cout << "How does that happen when reading a file??????"<< endl;
  exit(1);
 }

 _histos[dimensions][setName][subSetName][name] = (TH1*) _file->Get(setName+"/"+subSetName+"/"+histoName);
 if( !DoesSetExist(setName) ) _setNames.push_back(setName);
}


////////////////////////////
//Methods to Fill Histograms
////////////////////////////


//Fill 1D Histogram
void NueUtilities::Histograms::FillHisto(TString setName,TString subSetName,TString name,
                                         double valX, double eventWeight)
{
 FillHisto(setName,subSetName,name,valX, 0, 0, eventWeight, 1);
}

//Fill 2D Histogram
void NueUtilities::Histograms::FillHisto(TString setName,TString subSetName,TString name,
                                         double valX, double valY, double eventWeight)
{
 FillHisto(setName,subSetName,name,valX, valY, 0, eventWeight, 2);
}

//Fill 3D Histogram
void NueUtilities::Histograms::FillHisto(TString setName,TString subSetName,TString name,
                                         double valX, double valY, double valZ, double eventWeight)
{
 FillHisto(setName,subSetName,name,valX, valY, valZ, eventWeight, 3);
}

//Fill Histograms called by protected FillHisto member functions
void NueUtilities::Histograms::FillHisto(TString setName,TString subSetName,TString name,
                                         double valX, double valY, double valZ,
                                         double eventWeight, int dimensions)
{
  if( _histos.find(dimensions) == _histos.end() )
  {
   cout << "\n\n ERROR: Cannot find set of histograms for TH"<< dimensions<< "D"<< endl;
   cout << " Call to fill "<< setName << " " << subSetName << " " << name << endl << endl;
   exit(1);
  }


  if( _histos[dimensions].find(setName) == _histos[dimensions].end() )
  {
   cout << "\n\n ERROR: Cannot find SetName of TH"<< dimensions<< "D histograms for:"<< endl;
   cout << " Fill to "<< setName << " " << subSetName << " " << name << endl << endl;
   exit(1);
  }


  if( _histos[dimensions][setName].find(subSetName) == _histos[dimensions][setName].end() )
  {
   cout << "\n\n ERROR: Cannot find SubSetName of TH"<< dimensions<< "D histograms for:"<< endl;
   cout << " Fill to "<< setName << " " << subSetName << " " << name << endl << endl;
   exit(1);
  }

  if( _histos[dimensions][setName][subSetName].find(name) == _histos[dimensions][setName][subSetName].end() )
  {
   cout << "\n\n ERROR: Cannot find HistoName of TH"<< dimensions<< "D histograms for:"<< endl;
   cout << " Fill to "<< setName << " " << subSetName << " " << name << endl << endl;
   exit(1);
 }

  
  if( 1 == dimensions)
  {
           ((TH1D*)_histos[dimensions][setName][subSetName][name])->Fill(valX, eventWeight);
  }
  else if( 2 == dimensions)
  {
           ((TH2D*)_histos[dimensions][setName][subSetName][name])->Fill(valX, valY, eventWeight);
  }
  else if( 3 == dimensions)
  {
           ((TH3D*)_histos[dimensions][setName][subSetName][name])->Fill(valX, valY, valZ, eventWeight);
  }
  else
  {
           cout << "\n\nERROR: Don't know how to fill histogram of " << dimensions << " dimensions\n"<<endl;
           exit(1);
  }
}





////////////////////////////
//Methods to Write Histograms
////////////////////////////


//Write TH1D histograms here
void NueUtilities::Histograms::WriteHisto(TString setName,TString subSetName,TString name)
{
 bool histoNotFound = true;

 for( unsigned int dimensions=1; dimensions <= 3; ++dimensions)
 {
  if( _histos.find(dimensions) != _histos.end() )
  {
   if( _histos[dimensions].find(setName) != _histos[dimensions].end() )
   {
    if( _histos[dimensions][setName].find(subSetName) != _histos[dimensions][setName].end() )
    {
     if( _histos[dimensions][setName][subSetName].find(name) != _histos[dimensions][setName][subSetName].end() )
     {
      histoNotFound = false;

      TString currentDir = gDirectory->GetPath();
      _file->cd();

      if( !_file->Get( setName ) )
      {
       gDirectory->mkdir(setName);
      }
      gDirectory->cd(setName);

      if( !_file->Get( setName+"/"+subSetName)) 
      {
       gDirectory->mkdir(subSetName);
      }
      gDirectory->cd(subSetName);

      _histos[dimensions][setName][subSetName][name]->SetName( Form("th%dd_",dimensions)+name );
      _histos[dimensions][setName][subSetName][name]->Write();

      gROOT->cd(currentDir);
       
     }//Does name exist
    }//Does subSetName exist
   }//Does setName exist
  }//Dimension
 }

  if( histoNotFound )
  {
   cout << "\n\n ERROR: Cannot write histogram for:"<< endl;
   cout << " "<< setName << " " << subSetName << " " << name << endl << endl;
   exit(1);
  }
}





/////////////////////////////////
//Wrappers to Histograms Methods
/////////////////////////////////


void NueUtilities::Histograms::HistoTask(TString action, TString setName, TString subSetName, TString name,
                           int numBinsX, double lowValueX, double highValueX,
                           double fillValX, double eventWeight)
{
 HistoTask(action, setName, subSetName, name,
           numBinsX, lowValueX, highValueX,
           0,0,0,
           0,0,0,
           1, fillValX, 0, 0, eventWeight);
}

void NueUtilities::Histograms::HistoTask(TString action, TString setName, TString subSetName, TString name,
                           int numBinsX, double lowValueX, double highValueX,
                           int numBinsY, double lowValueY, double highValueY,
                           double fillValX, double fillValY, double eventWeight)
{
 HistoTask(action, setName, subSetName, name,
           numBinsX, lowValueX, highValueX,
           numBinsY, lowValueY, highValueY,
           0,0,0,
           2, fillValX, fillValY, 0, eventWeight);
}

void NueUtilities::Histograms::HistoTask(TString action, TString setName, TString subSetName, TString name,
                           int numBinsX, double lowValueX, double highValueX,
                           int numBinsY, double lowValueY, double highValueY,
                           int numBinsZ, double lowValueZ, double highValueZ,
                           double fillValX, double fillValY, double fillValZ, double eventWeight)
{
 HistoTask(action, setName, subSetName, name,
           numBinsX, lowValueX, highValueX,
           numBinsY, lowValueY, highValueY,
           numBinsZ, lowValueZ, highValueZ,
           3, fillValX, fillValY, fillValZ, eventWeight);
}

void NueUtilities::Histograms::HistoTask(TString action, TString setName, TString subSetName, TString name,
                           int numBinsX, double lowValueX, double highValueX,
                           int numBinsY, double lowValueY, double highValueY,
                           int numBinsZ, double lowValueZ, double highValueZ,
                           int dimensions, double fillValX, double fillValY, double fillValZ, double eventWeight)
{

 if( "make" == action)
 {
  MakeHisto(setName, subSetName, name,
            numBinsX, lowValueX, highValueX,
            numBinsY, lowValueY, highValueY,
            numBinsZ, lowValueZ, highValueZ,
            dimensions);
 }
 else if( "fill" == action)
 {
  FillHisto(setName, subSetName, name,
            fillValX, fillValY, fillValZ,
            eventWeight, dimensions);
 }
 else if( "write" == action)
 {
  WriteHisto(setName, subSetName, name);
 }
 else
 {
  cout << "\n\nERROR: Unknown HistoTask = " << action << "\n"<<endl;
  exit(1);
 }

}





////////////////////////////////
//Constructor, Destructors, Etc
////////////////////////////////


NueUtilities::Histograms::Histograms()
{
 _file=NULL;
}

NueUtilities::Histograms::~Histograms()
{
 CloseFile();
}

void NueUtilities::Histograms::CloseFile()
{
 if(_file != NULL)
 {
  //This should delete the histogram pointers
  _file->Close();
  delete _file;
 }
 _file=NULL;
 _histos.clear(); //shouldn't cause memory leak if _file->Close() worked
 _setNames.clear();
}

bool NueUtilities::Histograms::OpenFile(TString outputFile, TString option)
{
 CloseFile();
 _file = new TFile(outputFile, option);

 if(_file != NULL && _file->IsOpen())
 {
 
   //If reading the file, then get the histograms in the file
   if(
        (option.Data()[0] == 'r' || option.Data()[0] == 'R')
      &&(option.Data()[1] == 'e' || option.Data()[1] == 'E')
      &&(option.Data()[2] == 'a' || option.Data()[2] == 'A')
      &&(option.Data()[3] == 'd' || option.Data()[3] == 'D')
     )
   {
    AutoRetrieveHistograms();
   }
  
   return(true);
 }

 return( false );
}

bool NueUtilities::Histograms::DoesSetExist(TString setName) const
{
  return (_setNames.end() != find (_setNames.begin(), _setNames.end(), setName) );
}

unsigned int NueUtilities::Histograms::NumSets() const
{
 return (_setNames.size());
}

TString NueUtilities::Histograms::GetSetName(unsigned int setIndex) const
{
 if( setIndex < _setNames.size() )
 {
   return(_setNames[setIndex]);
 }
 return("");
}



