////////////////////////////////////////////////////////////////////////
//
// Module to put MCNN info into tree (from MCNN file format)
//
// template from boehm@physics.harvard.edu
// filled in by rbpatter@caltech.edu
////////////////////////////////////////////////////////////////////////
#include "TFile.h"
#include "Conventions/Detector.h"
#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include <fstream>
#include "AnalysisNtuples/ANtpDefaultValue.h"                                   
#include "TClass.h"
#include <string>
#include "MCNNMergeModule.h"
#include "NueAna/NueAnalysisCuts.h"
#include "NueAna/NueAnaTools/NueConvention.h"

const int CALend=110;
const int SM1end=230;
const int SM2end=484;

JOBMODULE(MCNNMergeModule, "MCNNMergeModule",
          "Import MCNN info into AnaNue tree");
CVSID("$Id: MCNNMergeModule.cxx,v 1.4 2009/05/24 22:19:53 rbpatter Exp $");
//......................................................................

MCNNMergeModule::MCNNMergeModule():
  kInputMCNNFile(""), kOkay(1), mcnn_entry(0)
{}

//......................................................................

MCNNMergeModule::~MCNNMergeModule()
{}

//......................................................................
void MCNNMergeModule::BeginJob()
{

  kOkay = kTRUE;
  mcnn_entry = 0;

  if (kInputMCNNFile != "") {

    std::cout << "MCNNMergeModule: Opening MCNN file for merging" << std::endl;
    std::cout << "--> FILE = " << kInputMCNNFile.c_str() << std::endl;
    _nnreader = new NNReader("nn",true);
    Int_t count = _nnreader->Add(kInputMCNNFile.c_str());
    if (!count||!_nnreader->GetEntries()) {
      std::cout << "--> No valid MCNN info found in file (or file not found)!" << std::endl;
      kOkay = kFALSE;
    }

    // initialize Pedro's NueAnaReader object
    _anareader = new NueAnaReader();
    _anareader->SetMaxNumBestMatches(kNumBestMatches);

    // initialize MCNNFiller object ...
    // ... these should become config options
    _filler = new MCNNFiller();
    _filler->SetBestMatchesToFill(kNumBestMatches);
    _filler->SetMCNNpidFile(kPDFFile.c_str());
    _filler->Setycut(kYCut);

    std::cout << "MCNNMergeModule -- starting up with this configuration:" << std::endl;
    std::cout << "  NumBestMatches = " << kNumBestMatches << std::endl;
    std::cout << "  YCut           = " << kYCut << std::endl;
    std::cout << "  PDFFile        = " << kPDFFile.c_str() << std::endl;

  }
  else {
    std::cout << "MCNNMergeModule: No MCNN file specified, so no merging will be done." << std::endl;
  }

}


JobCResult MCNNMergeModule::Reco(MomNavigator* mom)
{
  // okay?
  if (!kOkay) {
    std::cout << "MCNNMergeModule has a fatal error." << std::endl;
    return JobCResult::kFailed;
  }

  //get all NueRecords from mom 
  //may have more than one per go since mom reads in a snarl's worth of data
  //so, this is a little more complicated than just asking for a NueRecord
  TObject *obj=0;
  
  vector<NueRecord *> records;

  TIter objiter = mom->FragmentIter();
  while((obj=objiter.Next())){
    NueRecord *nr = dynamic_cast<NueRecord *>(obj);
    if(nr){
      MSG("MCNNMergeModule",Msg::kDebug)<<"Found a NueRecord in MOM"<<endl;
    }
    else{
      MSG("MCNNMergeModule",Msg::kDebug)<<"Didn't find a NueRecord in MOM"<<endl;
      continue;
    }

    // do it...
    if (kInputMCNNFile != "") {

      // Set the NueRecord object in _anareader object
      _anareader->SetNueRecord(nr);
      _nnreader->GetEntry(mcnn_entry++);

      // check that things are lined up
      MSG("MCNNMergeModule",Msg::kDebug) << "--------------" << std::endl;
      MSG("MCNNMergeModule",Msg::kDebug) << _anareader->nuerecord->GetHeader().GetRun() << " " << _nnreader->run << std::endl;
      MSG("MCNNMergeModule",Msg::kDebug) << _anareader->nuerecord->GetHeader().GetSnarl() << " " << _nnreader->snarl << std::endl;
      MSG("MCNNMergeModule",Msg::kDebug) << _anareader->nuerecord->GetHeader().GetEventNo() << " " << _nnreader->evt << std::endl;

      if ((_anareader->nuerecord->GetHeader().GetRun()     != _nnreader->run  )||
	  (_anareader->nuerecord->GetHeader().GetSnarl()   != _nnreader->snarl)||
	  (_anareader->nuerecord->GetHeader().GetEventNo() != _nnreader->evt  )) {
	MSG("MCNNMergeModule",Msg::kError) << "MCNNMergeModule: ERROR -- trees not in sync!" << std::endl;
	MSG("MCNNMergeModule",Msg::kError) << _anareader->nuerecord->GetHeader().GetRun() << " " << _nnreader->run << std::endl;
	MSG("MCNNMergeModule",Msg::kError) << _anareader->nuerecord->GetHeader().GetSnarl() << " " << _nnreader->snarl << std::endl;
	MSG("MCNNMergeModule",Msg::kError) << _anareader->nuerecord->GetHeader().GetEventNo() << " " << _nnreader->evt << std::endl;
	return JobCResult::kFailed;
      }
      else {
	// MCNN needs a calibrated energy in the NueRecord...
	NueConvention::NueEnergyCorrection(_anareader->nuerecord);
	// Fill the MCNN info!
	_filler->FillMCNN(_anareader,_nnreader);
      }
	  

    }
    
  }

  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}


////////////////////////////////////////////////////////////////////////
void MCNNMergeModule::EndJob()
{
  // Deleting these caused seg fault.  Whatever.  Just let the job end, and there'll be no problem.
  //delete _nnreader;
  //delete _anareader;
}

const Registry& MCNNMergeModule::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("MCNNMergeModule",Msg::kDebug)<<"In MCNNMergeModule::DefaultConfig"<<endl;

  static Registry r;
 
  // Set values in configuration
  r.UnLockValues();
  r.Set("InputMCNNFile", "");
  r.Set("NumBestMatches", 50);
  r.Set("YCut", 0.9);
  r.Set("PDFFile", "${SRT_PRIVATE_CONTEXT}/MCNNAnalysis/macros/files/MCNNpdf_2D_0.5GeVbins_MCNNv2_3.0pe_3.0pelib.root");
  r.LockValues();
                                                                                
  return r;
}

void MCNNMergeModule::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("MCNNMergeModule",Msg::kDebug)<<"In MCNNMergeModule::Config"<<endl;
                                                                                 
  const char* tmps;
  Int_t tmpi;
  Double_t tmpd;
  if (r.Get("InputMCNNFile", tmps)) kInputMCNNFile = tmps;
  if (r.Get("PDFFile", tmps)) kPDFFile = tmps;
  if (r.Get("NumBestMatches", tmpi)) kNumBestMatches = tmpi;
  if (r.Get("YCut", tmpd)) kYCut = tmpd;

}



