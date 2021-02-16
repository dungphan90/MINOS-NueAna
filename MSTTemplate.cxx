////////////////////////////////////////////////////////////////////////
// $Id: MSTTemplate.cxx,v 1.4 2008/11/19 18:22:51 rhatcher Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TMath.h"
#include "Conventions/Detector.h"
#include "NueAna/MSTTemplate.h"
#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "HistMan/HistMan.h"

const int CALend=110;
const int SM1end=230;
const int SM2end=484;
const int NXBINS=31;
const int NYBINS=10;
const int NZBINS=31;
const int NYMBINS=21;


JOBMODULE(MSTTemplate, "MSTTemplate",
          "Makes hists to compare MST variables between detectors and truth classes");
CVSID("$Id: MSTTemplate.cxx,v 1.4 2008/11/19 18:22:51 rhatcher Exp $");
//......................................................................

MSTTemplate::MSTTemplate():
  counter(0),
  fname(),
  fout(0),
  nlambdanele(0),
  nlambdanother(0),
  mipdistele(0),  
  mipdistother(0)
{
  //x corresponds to the number of edges in the event
  //it goes from 0-30, then all events with 60<=#edges<300
  //are grouped together in one bin
  xbins=new double[NXBINS];
  for(int i=0;i<NXBINS-1;i++){
    xbins[i]=1.*i;
  }
  xbins[NXBINS-1]=300;
  
  //y corresponds to the values of the weights
  //it goes from 0-40, in steps of 5, 40<weights<1000 are
  //grouped together in one bin
  ybins=new double[NYBINS];
  for(int i=0;i<NYBINS-1;i++){
    ybins[i]=i*5.;
  }
  ybins[NYBINS-1]=1000.;
  
  //ymbins corresponds to the mip distributions, not great yet
  ymbins=new double[NYMBINS];
  for(int i=0;i<NYMBINS-1;i++){
    ymbins[i]=i*5;
  }
  ymbins[NYMBINS-1]=150;
  
  //z also corresponds to the number of edges in the event
  //it goes from 0-30, then all events with 30<#edges<100
  //are grouped together in one bin
  zbins=new double[NZBINS];
  for(int i=0;i<NZBINS-1;i++){
    zbins[i]=i;
  }
  zbins[NZBINS-1]=100;

}

//......................................................................

MSTTemplate::~MSTTemplate()
{
  delete [] xbins;

}

//......................................................................
void MSTTemplate::BeginJob()
{

  fout=new TFile(fname.c_str(),"RECREATE");
  
   //create the template hisograms
   nlambdanele=new TProfile2D("nlambdanele","nlambdanele",
			      NXBINS-1,xbins,NYBINS-1,ybins);
   nlambdanother=new TProfile2D("nlambdanother","nlambdanother",
				NXBINS-1,xbins,NYBINS-1,ybins);
   
   mipdistele=new TProfile2D("mipdistele","mipdistele",
			     NXBINS-1,xbins,NYMBINS-1,ymbins);
   mipdistother=new TProfile2D("mipdistother","mipdistother",
			       NXBINS-1,xbins,NYMBINS-1,ymbins);

}


JobCResult MSTTemplate::Ana(const MomNavigator* mom)
{
  if(counter%1000==0){
    cout<<"on entry "<<counter<<endl;
  }
  counter++;

   //get all NueRecords from mom 
   //may have more than one per go since mom reads in a snarl's worth of data
   //so, this is a little more complicated than just asking for a NueRecord
   TObject *obj=0;
   TIter objiter = mom->FragmentIter();
   while((obj=objiter.Next())){
      NueRecord *nr = static_cast<NueRecord *>(obj);
      if(nr){
	 MSG("MSTTemplate",Msg::kDebug)<<"Found a NueRecord in MOM"<<endl;
      }
      else{
	 MSG("MSTTemplate",Msg::kError)<<"Didn't find a NueRecord in MOM"<<endl;
	 continue;
      }
      MSG("MSTTemplate",Msg::kDebug)<<"Found a NueRecord "<<nr<<endl;
      
      //figure out whether event is signal or background
      bool signal = false;
      bool bg = false;
      if(abs(nr->mctrue.nuFlavor)==12&&nr->mctrue.interactionType==1&&
	 nr->mctrue.resonanceCode==1001){
	signal=true;
      }
      if(nr->mctrue.interactionType==0){
	bg=true;
      }

      //if it's neither signal or bg, go on
      if(!(signal||bg)){
	continue;
      }

      //need at least 1 sucessfully reconstructed event
      if(nr->GetHeader().GetEventNo()<0){
	continue;
      }
      //only look at events for which mst actually gets calculated
      if(nr->GetHeader().GetTrackLength()>25){
	continue;
      }
      //only look at contained events (vertex contained, end doesn't go past 
      //calorimeter or sm gap
      if(nr->GetHeader().GetVldContext().GetDetector()==Detector::kFar){
	if(sqrt(pow(nr->srevent.vtxX,2)+pow(nr->srevent.vtxY,2))>3.87){
	  continue;
	}
	if(nr->srevent.vtxZ<1||(nr->srevent.vtxZ>14&&nr->srevent.vtxZ<17)||
	   nr->srevent.vtxZ>20){
	  continue;
	}
	if((nr->srevent.begPlane<SM1end&&nr->srevent.endPlane>SM1end-2)||
	   (nr->srevent.begPlane>SM1end&&nr->srevent.endPlane>SM2end-2)){
	  continue;
	}	  
      }
      else if(nr->GetHeader().GetVldContext().GetDetector()==Detector::kNear){
	if(sqrt(pow(nr->srevent.vtxX-1.5,2)+pow(nr->srevent.vtxY,2))>1.0){
	  continue;
	}
	if(nr->srevent.vtxZ<4||nr->srevent.vtxZ>6.5){
	  continue;
	}
	if(nr->srevent.endPlane>CALend-2){
	  continue;
	}	  
      }
      
    
      //make some temporary histograms
      TH1F ewhist("ewhist","ewhist",NYBINS-1,ybins);
      TH1F emhist("emhist","emhist",NYMBINS-1,ymbins);
      TH1F owhist("owhist","owhist",NYBINS-1,ybins);
      TH1F omhist("omhist","omhist",NYMBINS-1,ymbins);

      //fill the temp hists.
      for(int i=0;i<nr->mstvars.enn1;i++){
	 if(nr->mstvars.eallw1[i]!=0){
	    ewhist.Fill(nr->mstvars.eallw1[i]);
	 }
	 if(nr->mstvars.eallm1[i]!=0){	 
	    emhist.Fill(nr->mstvars.eallm1[i]);     
	 }
      }
      for(int i=0;i<nr->mstvars.onn1;i++){
	 if(nr->mstvars.oallw1[i]!=0){
	    owhist.Fill(nr->mstvars.oallw1[i]);
	 }
	 if(nr->mstvars.oallm1[i]!=0){
	    omhist.Fill(nr->mstvars.oallm1[i]);     
	 }
      }
      

      if(signal){
	 //fill signal profile hist.
	 for(int j=1;j<=ewhist.GetNbinsX()+1;j++){
	    nlambdanele->Fill(nr->mstvars.enn1,
			      ewhist.GetBinCenter(j),
			      ewhist.GetBinContent(j));
	 }
	 for(int j=1;j<=owhist.GetNbinsX()+1;j++){
	    nlambdanele->Fill(nr->mstvars.onn1,
			      owhist.GetBinCenter(j),
			      owhist.GetBinContent(j));
	 }
	 for(int j=1;j<=emhist.GetNbinsX()+1;j++){
	    mipdistele->Fill(nr->mstvars.enn1,
			     emhist.GetBinCenter(j),
			     emhist.GetBinContent(j));
	 }
	 for(int j=1;j<=omhist.GetNbinsX()+1;j++){
	    mipdistele->Fill(nr->mstvars.onn1,
			     omhist.GetBinCenter(j),
			     omhist.GetBinContent(j));
	 }

      }
      else if(bg){
	 //fill bg profile hist.
	 for(int j=1;j<=ewhist.GetNbinsX()+1;j++){
	    nlambdanother->Fill(nr->mstvars.enn1,
				ewhist.GetBinCenter(j),
				ewhist.GetBinContent(j));
	 }
	 for(int j=1;j<=owhist.GetNbinsX()+1;j++){
	    nlambdanother->Fill(nr->mstvars.onn1,
				owhist.GetBinCenter(j),
				owhist.GetBinContent(j));
	 }
	 for(int j=1;j<=emhist.GetNbinsX()+1;j++){
	    mipdistother->Fill(nr->mstvars.enn1,
			       emhist.GetBinCenter(j),
			       emhist.GetBinContent(j));
	 }
	 for(int j=1;j<=omhist.GetNbinsX()+1;j++){
	    mipdistother->Fill(nr->mstvars.onn1,
			       omhist.GetBinCenter(j),
			       omhist.GetBinContent(j));
	 }
      }

      //reset the temp histograms
      ewhist.Reset();
      owhist.Reset();
      emhist.Reset();
      omhist.Reset();
   }
  
   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
void MSTTemplate::EndJob()
{
  if(fout!=0){
    if(fout->IsOpen()){
      fout->Write();
      fout->Close();
    }
  }
}

//......................................................................

const Registry& MSTTemplate::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("MSTTemplate",Msg::kDebug)<<"In MSTTemplate::DefaultConfig"<<endl;

  static Registry r; // Default configuration for module

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  // Set values in configuration
  r.UnLockValues();
  r.Set("MSTTmpltFile","templates.root");
  r.LockValues();

  return r;
}

//......................................................................

void MSTTemplate::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("MSTTemplate",Msg::kDebug)<<"In MSTTemplate::Config"<<endl;

  const char* tmps;

  if(r.Get("MSTTmpltFile",tmps)) { fname = (string)(tmps);}
}
