////////////////////////////////////////////////////////////////////////
// $Id: FillPIDFromText.cxx,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <cassert>
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "CandNtupleSR/NtpSRRecord.h"
#include "NueAna/NuePID.h"
#include "NueAna/FillPIDFromText.h"

JOBMODULE(FillPIDFromText, "FillPIDFromText",
          "fills pid info from a text file");
CVSID("$Id:");
//......................................................................

TxtEntry::TxtEntry():
   run(0),
   subrun(0),
   snarl(0),
   event(0)
{}

TxtEntry::TxtEntry(int r, int sbr, int s, int e):
   run(r),
   subrun(sbr),
   snarl(s),
   event(e)
{}

TxtEntry::TxtEntry(int r, int sbr, int s, int e, int f, int rsn):
   run(r),
   subrun(sbr),
   snarl(s),
   event(e),
   flav(f),
   res(rsn)
{}

TxtEntry::TxtEntry(int r, int sbr, int s, int e, float l):
   run(r),
   subrun(sbr),
   snarl(s),
   event(e),
   like(l)
{}


TxtEntry::~TxtEntry()
{}
FillPIDFromText::FillPIDFromText():
   counter(0),
   fTextFile(),
   fDecider(NuePIDHeader::kUnknown),
   elist(),
   kSelRes(-1),
   kSelFlav(-1)

{}

//......................................................................

FillPIDFromText::~FillPIDFromText()
{}

//......................................................................

void FillPIDFromText::BeginJob()
{
   ReadTextFile();
}

//......................................................................

JobCResult FillPIDFromText::Reco(MomNavigator* mom)
{
   if(counter%1000==0){
      MSG("FillPIDFromText",Msg::kInfo)<<"On entry "<<counter<<endl;
   }
   counter++;

   //make sure the elist set was filled
   if(elist.size()==0){
      MSG("FillPIDFromText",Msg::kError)<<"No entries in elist, abort"<<endl;
      assert(0);
   }

   //read in a NtpSRRecord
   NtpSRRecord *sr = static_cast<NtpSRRecord *>(mom->GetFragment("NtpSRRecord"));
   if(!sr){
      MSG("FillPIDFromText",Msg::kError)<<"Couldn't get a NtpSRRecord from mom"<<endl;
      return JobCResult::kFailed;
   }
   //get vldcontext
   VldContext vc=sr->GetHeader().GetVldContext();

   Int_t dec=-1;
   Float_t likelihood=1;

   int evtn=sr->evthdr.nevent;
   if(evtn==0){
      //do something clever, 
      //ok, not very clever, just do it once
      //make a pid header
      NuePIDHeader h(vc);
      h.SetSnarl(sr->GetHeader().GetSnarl());
      h.SetRun(sr->GetHeader().GetRun());
      h.SetSubRun(sr->GetHeader().GetSubRun());
      h.SetEventNo(-1);
      h.SetEvents(evtn);
      h.SetDecider(fDecider);

      //make a pid object
      NuePID *pid=new NuePID(h);


      //is this event in the set?
      TxtEntry test(h.GetRun(),h.GetSubRun(),h.GetSnarl(),0);
      std::set<TxtEntry>::iterator t=elist.find(test);
      if(t==elist.end()){
	 //not in our list
	 pid->IsNue=-1;
	 pid->likelihood=1;
      }
      else{
	 //is in our list

	if(fDecider==6) {
	  Int_t decflav=-1;
	  Int_t decres=-1;
	  if(kSelFlav==-1){decflav=1;}
	  else if((*t).flav==kSelFlav&&kSelFlav>-1){decflav =1;}
	  if(kSelRes==-1){decres=1;}	  
	  else if((*t).res==kSelRes&&kSelRes>-1){decres = 1;}

	  if(decflav==-1||decres==-1){dec=-1;}
	  else {dec=1;}
	}
	if(fDecider==7) {likelihood=(*t).like; dec=1;}
	MSG("FillPIDFromText",Msg::kDebug) << (*t).run 
			 << (*t).subrun << " " << (*t).snarl << " " 
			 << (*t).event << " " << (*t).flav << " " << (*t).res 
			 << " " << dec << endl;
	 pid->IsNue=dec;
	 pid->likelihood=likelihood;
      }
      MSG("FillPIDFromText",Msg::kDebug)<<"PID decision: "<<pid->IsNue<<endl;
      //give pid object to mom to write to file
      mom->AdoptFragment(pid);
      return JobCResult::kPassed;
   }

   //loop over events is NtpSRRecord
   for(int i=0;i<evtn;i++){
      //make a pid header
      NuePIDHeader h(vc);
      h.SetSnarl(sr->GetHeader().GetSnarl());
      h.SetRun(sr->GetHeader().GetRun());
      h.SetSubRun(sr->GetHeader().GetSubRun());
      h.SetEventNo(i);
      h.SetEvents(evtn);
      h.SetDecider(fDecider);

      //make a pid object
      NuePID *pid=new NuePID(h);

      //is this event in the set?
      TxtEntry test(h.GetRun(),h.GetSubRun(),h.GetSnarl(),h.GetEventNo());
      std::set<TxtEntry>::iterator t=elist.find(test);
      if(t==elist.end()){
	 //not in our list
	 pid->IsNue=-1;
	 pid->likelihood=1;
      }
      else{
	if(fDecider==6) {
	  Int_t decflav=-1;
	  Int_t decres=-1;
	  if(kSelFlav==-1){decflav=1;}
	  else if((*t).flav==kSelFlav&&kSelFlav>-1){decflav =1;}
	  if(kSelRes==-1){decres=1;}	  
	  else if((*t).res==kSelRes&&kSelRes>-1){decres = 1;}

	  if(decflav==-1||decres==-1){dec=-1;}
	  else {dec=1;}
	}
	if(fDecider==7) {likelihood=(*t).like; dec=1;}
	MSG("FillPIDFromText",Msg::kDebug)<< (*t).run 
	     << (*t).subrun << " " << (*t).snarl << " " 
	     << (*t).event << " " << (*t).flav << " " << (*t).res 
	     << " " << dec << endl;
	 //is in our list
	 pid->IsNue=dec;
	 pid->likelihood=likelihood;
      }
      MSG("FillPIDFromText",Msg::kDebug)<<"PID decision: "<<pid->IsNue<<endl;
      //give pid object to mom to write to file}
      mom->AdoptFragment(pid);
   }

  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................

const Registry& FillPIDFromText::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
  static Registry r; // Default configuration for module

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  // Set values in configuration
  r.UnLockValues();
  r.Set("TextFile","pid.txt");
  r.Set("Decider",0);
  r.Set("SelRes",-1);
  r.Set("SelFlav",-1);

  r.LockValues();

  return r;
}

//......................................................................

void FillPIDFromText::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  const char* tmps;

  if (r.Get("TextFile",tmps)) { fTextFile = (string)(tmps); }
  int d;
  if (r.Get("Decider",d)) { fDecider = (NuePIDHeader::Decider_t)(d); }
  if (r.Get("SelRes",d)) { kSelRes = d; }
  if (r.Get("SelFlav",d)) { kSelFlav = d; }
}

//......................................................................

void FillPIDFromText::ReadTextFile()
{

   int nread=0;
   std::ifstream in(fTextFile.c_str());
   if(!in){
      MSG("FilPIDFromText",Msg::kError)<<"Could not open "<<fTextFile<<endl;
      return;
   }

   int run, subrun, snarl, event, flav, res;
   float like;
   in>>run;
   while(!in.eof()){
      nread++;
      in>>subrun>>snarl>>event;
      if(fDecider==6) {in >> flav >> res; 
	MSG("FillPIDFromText",Msg::kDebug)<< run << " " 
             << subrun << " " << snarl << " " 
	     << event << " " << flav << " " << res << endl;}
      if(fDecider==7) in >> like;

      TxtEntry *a;
      if (fDecider==6) {	TxtEntry t(run,subrun,snarl,event, flav,res);a=&t;}
      else if (fDecider==7) {	TxtEntry t(run,subrun,snarl,event, like);a=&t;}
      else{                     TxtEntry t(run,subrun,snarl,event);a=&t;}

      elist.insert(*a);
      in>>run;
   }

   MSG("FillPIDFromText",Msg::kDebug)<<"Read in  "<<nread<<" entries "<<endl;
   MSG("FillPIDFromText",Msg::kDebug)<<"elist has  "<<elist.size()<<" entries "<<endl;
}

////////////////////////////////////////////////////////////////////////
