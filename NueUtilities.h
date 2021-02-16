////////////////////////////////////////////////////////////////////////
// $Id: NueUtilities.h,v 1.10 2015/12/22 05:36:18 wingmc Exp $
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
#ifndef  NUEUTILITIES
#define NUEUTILITIES

//------------------------------------------------------------------------------
//include C++ STL headers
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>  //Going to store histograms in map 
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//include ROOT related headers
#include "TChain.h"
#include "TROOT.h"  //Access to gROOT object
#include "TFile.h"  //Going to store output in a ROOT file 
#include "TKey.h"
#include "TH1D.h"   //Can store 1D histogram in file
#include "TH2D.h"   //Can store 2D histogram in file
#include "TH3D.h"   //Can store 3D histogram in file
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//include minossoft related headers
#include "NueAna/NuePOT.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NueStandard.h"
//------------------------------------------------------------------------------


namespace NueUtilities
{

 TChain* CreateChain(std::string treeName,std::string fileListName);

 Double_t GetChainPOT(TChain *chain);

 Double_t GetWeightedChainPOT(TChain *chain);

 class AnaNueProcessor
 {

  public:
   AnaNueProcessor();
   AnaNueProcessor(std::string fileListName);
   AnaNueProcessor(const AnaNueProcessor &original);
   AnaNueProcessor & operator=(const AnaNueProcessor &original);
   ~AnaNueProcessor();

   //Methods related to total sample
   Double_t GetTotalPOT();

   //Did the AnaNueProcessor run over any of these type of events
   inline Bool_t didMC() const { return(ProcessedMC); };
   inline Bool_t didData() const { return(ProcessedData); };
   inline Bool_t didFar() const { return(ProcessedFar); };
   inline Bool_t didNear() const { return(ProcessedNear); };
   inline Bool_t didNC() const { return(ProcessedNC); };
   inline Bool_t didCC() const { return(ProcessedCC); };
   inline Bool_t didBeamNue() const { return(ProcessedBeamNue); };
   inline Bool_t didTau() const { return(ProcessedTau); };
   inline Bool_t didSignal() const { return(ProcessedSignal); };
   inline Bool_t didRun1() const { return(ProcessedRun1); };
   inline Bool_t didRun2() const { return(ProcessedRun2); };
   inline Bool_t didRun3() const { return(ProcessedRun3); };

   inline Bool_t didRun11() const { return(ProcessedRun11); };
   inline Bool_t didRun12() const { return(ProcessedRun12); };
   inline Bool_t didRun13() const { return(ProcessedRun13); };


   //Methods related event processing
   Bool_t NextEntry();
   Bool_t SetEntry(Long64_t entry);

   Long64_t GetEntries() const {return(numEntries);};
   Long64_t CurrentEntry() const {return(currentEntry);};

   void PrintProgress(unsigned int percentInterval);
   std::string GetFileListName() const {return(filelist);};
   Int_t GetLastFraction() const {return(lastfraction);};

   inline Int_t GetNumberOfFiles(){ return( ana_nue_chain->GetNtrees() ); };
   void PrintFileNames();

   //Methods related to current event
   Bool_t isMC();
   Bool_t isData();

   Bool_t isFar();
   Bool_t isNear();

   Bool_t isNC();
   Bool_t isCC();
   Bool_t isBeamNue();
   Bool_t isTau();
   Bool_t isSignal();

   Bool_t isNu();
   Bool_t isNuBar();

   Bool_t isQuasiElastic();

   Bool_t isRun1();
   Bool_t isRun2();
   Bool_t isRun3();
   Bool_t isRun4RHC();
   Bool_t isRun7RHC();
   Bool_t isRun9RHC();

   Bool_t isRun11();
   Bool_t isRun12();
   Bool_t isRun13();
   
   Double_t Normalization_RunSeparatedSample();
   Double_t Normalization_MixedRunSample();
   Double_t Normalization_RunSeparatedSampleRHC();
   Double_t Normalization_MixedRunSampleRHC();


   Double_t GetNDDataWeightedPOT(Bool_t mixed = false, unsigned int percentInterval = 5);

   Double_t GetOscWeight_f210f213f214Separate();

   Double_t GetRecoE();
   Double_t GetTrueE();

   NueRecord *nueRecord;

  private:
   std::string filelist;
   TChain *ana_nue_chain;
   TChain *pottree_chain;
   Double_t totalPOT;
   bool totalPOTNotFilled;
   
   Long64_t numEntries;
   Long64_t currentEntry;
   Int_t lastfraction;

   Bool_t ProcessedMC;
   Bool_t ProcessedData;
   Bool_t ProcessedFar;
   Bool_t ProcessedNear;
   Bool_t ProcessedNC;
   Bool_t ProcessedCC;
   Bool_t ProcessedBeamNue;
   Bool_t ProcessedTau;
   Bool_t ProcessedSignal;
   Bool_t ProcessedRun1;
   Bool_t ProcessedRun2;
   Bool_t ProcessedRun3;

   Bool_t ProcessedRun11;
   Bool_t ProcessedRun12;
   Bool_t ProcessedRun13;

 };//End Class AnaNueProcessor





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

class Histograms
{
 public:
  Histograms();
  ~Histograms();

  bool DoesSetExist(TString setName) const;
  unsigned int NumSets() const;
  TString GetSetName(unsigned int setIndex) const;
  TH1* GetHistogram(TString set,TString subset,TString name);

 protected:

  bool OpenFile(TString outputFile, TString option="READ");
  void CloseFile();

  //Make 1D Histogram
  void MakeHisto(TString setName, TString subSetName, TString name,
                 int numBinsX, double lowValueX, double highValueX);
  //Make 2D Histogram
  void MakeHisto(TString setName, TString subSetName, TString name,
                 int numBinsX, double lowValueX, double highValueX,
                 int numBinsY, double lowValueY, double highValueY);
  //Make 3D Histogram
  void MakeHisto(TString setName,TString subSetName,TString name,
                 int numBinsX, double lowValueX, double highValueX,
                 int numBinsY, double lowValueY, double highValueY,
                 int numBinsZ, double lowValueZ, double highValueZ);

  //Make Histograms called by protected MakeHisto member functions
  void MakeHisto(TString setName,TString subSetName,TString name,
                 int numBinsX, double lowValueX, double highValueX,
                 int numBinsY, double lowValueY, double highValueY,
                 int numBinsZ, double lowValueZ, double highValueZ,
                 int dimensions);

  //Make Histograms by Retrieving Histograms in File
  void AutoRetrieveHistograms();
  void RetrieveHistogram(TString setName, TString subsetName, TString histoName, unsigned int dimensions);


  //Fill 1D Histogram
  void FillHisto(TString setName,TString subSetName,TString name,
                 double valX, double eventWeight);

  //Fill 2D Histogram
  void FillHisto(TString setName,TString subSetName,TString name,
                 double valX, double valY, double eventWeight);

  //Fill 3D Histogram
  void FillHisto(TString setName,TString subSetName,TString name,
                 double valX, double valY, double valZ, double eventWeight);

  //Fill Histograms called by protected FillHisto member functions
  void FillHisto(TString setName,TString subSetName,TString name,
                 double valX, double valY, double valZ,
                 double eventWeight, int dimensions);


  //Write Histogram to file
  void WriteHisto(TString setName,TString subSetName,TString name);



  //Wrapper to above Histogram methods

  void HistoTask(TString action, TString setName, TString subSetName, TString name,
                 int numBinsX, double lowValueX, double highValueX,
                 double fillValX, double eventWeight);
  void HistoTask(TString action, TString setName, TString subSetName, TString name,
                 int numBinsX, double lowValueX, double highValueX,
                 int numBinsY, double lowValueY, double highValueY,
                 double fillValX, double fillValY, double eventWeight);
  void HistoTask(TString action, TString setName, TString subSetName, TString name,
                 int numBinsX, double lowValueX, double highValueX,
                 int numBinsY, double lowValueY, double highValueY,
                 int numBinsZ, double lowValueZ, double highValueZ,
                 double fillValX, double fillValY, double fillValZ, double eventWeight);
  void HistoTask(TString action, TString setName, TString subSetName, TString name,
                 int numBinsX, double lowValueX, double highValueX,
                 int numBinsY, double lowValueY, double highValueY,
                 int numBinsZ, double lowValueZ, double highValueZ,
                 int dimensions, double fillValX, double fillValY, double fillValZ, double eventWeight);


 private:
  //Don't make a copy of this class! Just pass by reference if needed!
  //Here the assignment and copy constructors are made private to prevent you from doing that
  Histograms(const Histograms&){};
  void operator=(const Histograms&){};
  



  //data members
  std::map< unsigned int, std::map< TString, std::map< TString, std::map< TString,TH1*> > > > _histos;  //Organized collection of TH1D, TH2D, and TH3D pointers
  std::vector< TString> _setNames;
  TFile* _file;

};//End Class Histograms









} //End namespace NueUtilities

#endif


