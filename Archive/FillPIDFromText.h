////////////////////////////////////////////////////////////////////////
// $Id: FillPIDFromText.h,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef FILLPIDFROMTEXT_H
#define FILLPIDFROMTEXT_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <set>
#include "NueAna/NuePIDHeader.h"

class TxtEntry;


class FillPIDFromText : public JobCModule
{
public:
  FillPIDFromText();
  ~FillPIDFromText();

public:
  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);
  virtual void BeginJob();
  void ReadTextFile();

  // Module configuration
  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);


private:
  int counter;
  string fTextFile;
  NuePIDHeader::Decider_t fDecider;
  std::set<TxtEntry> elist;
  // Module member data
  int kSelRes;
  int kSelFlav;

 

};

class TxtEntry
{
public:
   TxtEntry();
   TxtEntry(int r, int sbr, int s, int e);
   TxtEntry(int r, int sbr, int s, int e, int f, int rsn);
   TxtEntry(int r, int sbr, int s, int e, float l);
   ~TxtEntry();
   
   bool operator < (const TxtEntry &e2) const {
      if(run<e2.run) return true;
      if(run>e2.run) return false;
      if(subrun<e2.subrun) return true;
      if(subrun>e2.subrun) return false;
      if(snarl<e2.snarl) return true;
      if(snarl>e2.snarl) return false;      
      if(event<e2.event) return true;
      if(event>e2.event) return false;
      return false;
   }

   bool operator == (const TxtEntry &e2) const {
      if(run==e2.run&&subrun==e2.subrun&&snarl==e2.snarl&&event==e2.event) return true;
      return false;
   }

   int run;
   int subrun;
   int snarl;
   int event;

   int flav;
   int res;
   float like;
};

#endif // FILLPIDFROMTEXT_H
////////////////////////////////////////////////////////////////////////
