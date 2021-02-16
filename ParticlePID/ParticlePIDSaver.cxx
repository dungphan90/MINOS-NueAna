////////////////////////////////////////////////////////////////////////
/// $Id: ParticlePIDSaver.cxx,v 1.6 2009/08/13 11:29:10 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "NueAna/ParticlePID/ParticlePIDSaver.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro


JOBMODULE(ParticlePIDSaver, "ParticlePIDSaver",
          "copies PRecords from MOM into NueRecord");
CVSID("$Id: ParticlePIDSaver.cxx,v 1.6 2009/08/13 11:29:10 scavan Exp $");




//......................................................................

ParticlePIDSaver::ParticlePIDSaver() 
{

  ///
  /// (Document me!)
  ///
}
//......................................................................

ParticlePIDSaver::~ParticlePIDSaver()
{
  ///
  /// (Document me!)
  ///

}

//......................................................................

void ParticlePIDSaver::BeginJob()
{
  ///
  /// (Document me!)
  ///
}

//......................................................................

void ParticlePIDSaver::EndJob()
{
  ///
  /// (Document me!)
  ///

  if(kPOTTreeName=="")return; //don't want to copy pottree...


   TFile *fpf = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->FindObject(kPOTTreeName.c_str()));
   if(fpf){
	printf("got outfile\n");
   }


   TFile *fpi = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->FindObject(kPOTTreeNameIn.c_str()));
   if(fpi){
        printf("got infile\n");
        TTree *pt = (TTree*)fpi->Get("pottree");	   
	if(pt)
	{
		printf("got pottree\n");
		fpf->cd();
		pt->Write();
	}
   }




}

//......................................................................

JobCResult ParticlePIDSaver::Reco(MomNavigator* mom)
{
  ///
  /// (Document me!)
  ///


  std::vector<TObject* > nrv = ( mom->GetFragmentList("NueRecord"));
  std::vector<TObject* > hft =( mom->GetFragmentList("PRecord","Normal"));
  std::vector<TObject* > hftm =( mom->GetFragmentList("PRecord","MRCC"));


for(unsigned int s =0;s<hft.size();s++)
{

	PRecord *pr=0;
	pr=dynamic_cast<PRecord *>(hft[s]);

	if(!pr)continue;
	
	NueRecord *nr=0;
	if(s<nrv.size())
		nr=dynamic_cast<NueRecord *>(nrv[s]);


	//make sure we are synced...
	int match=0;
	if(nr)match=isMatched(nr,pr);
	
	if(!match)
	{
		for(int m=0;m<(int)nrv.size();m++)
		{
			nr=dynamic_cast<NueRecord *>(nrv[m]);
			int mm = isMatched(nr,pr);
			if(mm){
				match=1;
				break;
			}
		}
	}
  
	if(!match)continue;//unable to finding matching nue record... probably not passing nue quality cuts


/*
	if(!match)
	{
		printf("unable to match in ParticlePIDSaver!\n");
		exit(1);
	}
*/	
	
	(nr->precord) = *pr;

	//is there an mrcc record?
	for(unsigned int j=0;j<hftm.size();j++)
	{
	        PRecord *pr=0;
        	pr=dynamic_cast<PRecord *>(hftm[j]);
		if(!pr)continue;

	        //make sure we are synced...
       		int match=1;
        	match = match && (nr->GetHeader().GetRun() == pr->GetHeader().GetRun());
        	match = match && (nr->GetHeader().GetSnarl() == pr->GetHeader().GetSnarl());
        	match = match && (nr->GetHeader().GetEventNo() == pr->GetHeader().GetEvent());

		if(match)
		{

        		(nr->precordMRCC) = *pr;

		}

	}

}







  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................

const Registry& ParticlePIDSaver::DefaultConfig() const
{
  ///
  /// Supply the default configuration for the module
  ///
  static Registry r; // Default configuration for module


  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());




  // Set values in configuration
  r.UnLockValues();
  r.Set("POTTreeName","");
  r.Set("POTTreeNameIn","");

  r.LockValues();


  return r;
}

//......................................................................

void ParticlePIDSaver::Config(const Registry& r)
{
  ///
  /// Configure the module given the Registry r
  ///
  
    const char* tmps;
    if(r.Get("POTTreeName",tmps))kPOTTreeName=tmps;
    if(r.Get("POTTreeNameIn",tmps))kPOTTreeNameIn=tmps;


  
}

//......................................................................

void ParticlePIDSaver::Reset()
{
  ///
  /// (Document me!)
  ///
}

////////////////////////////////////////////////////////////////////////

		


int ParticlePIDSaver::isMatched(NueRecord *nr, PRecord *pr)
{
	int match = 1;
        match = match && (nr->GetHeader().GetRun() == pr->GetHeader().GetRun());
        match = match && (nr->GetHeader().GetSnarl() == pr->GetHeader().GetSnarl());
        match = match && (nr->GetHeader().GetEventNo() == pr->GetHeader().GetEvent());

	return match;
}


