////////////////////////////////////////////////////////////////////////
//
// Module to put MCNN info into tree (from MCNN file format)
//
// template from boehm@physics.harvard.edu
// filled in by rbpatter@caltech.edu
////////////////////////////////////////////////////////////////////////
#ifndef MERGEMCNNMOD_H
#define MERGEMCNNMOD_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <vector>
#include "TTree.h"
#include "MCNNAnalysis/NNReader.h"
#include "MCNNAnalysis/NueAnaReader.h"
#include "MCNNAnalysis/MCNNFiller.h"
#include "NueAna/NueAnalysisCuts.h"

using namespace std;


class TH1F;
class TFile;
class NueRecord;


class MCNNMergeModule : public JobCModule
{
public:
  MCNNMergeModule();
  ~MCNNMergeModule();

public:
  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);
  void EndJob();
  void BeginJob();

  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);

private:

  std::string kInputMCNNFile;
  Bool_t kOkay;
  Int_t mcnn_entry;

  Int_t kNumBestMatches;
  Double_t kYCut;
  std::string kPDFFile;

  // Pedro's wrappers, which do the MCNN calculations, etc.
  NNReader *_nnreader;
  NueAnaReader *_anareader;
  MCNNFiller *_filler;


};
#endif 
////////////////////////////////////////////////////////////////////////
