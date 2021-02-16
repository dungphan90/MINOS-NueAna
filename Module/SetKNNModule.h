#ifndef SETKNNMODULE_H
#define SETKNNMODULE_H

//
// $Id: SetKNNModule.h,v 1.2 2008/10/25 22:56:02 boehm Exp $
//

// MINOS
#include "JobControl/JobCModule.h"
#include "Registry/Registry.h"

// Local
#include "PhysicsNtuple/Handle.h"
#include "PhysicsNtuple/Store/Interface.h"

class MomNavigator;

class StorekNNData: virtual public Anp::Base
{
   public:
      StorekNNData();
      virtual ~StorekNNData();

      bool SetValidity(VldContext &vld);
      const VldContext & GetValidity();
      bool Add(int i, const std::string &key, float data);      
      bool Get(int i, const std::string &key, float &data);
      bool SetPrefix(std::string in);
      bool Clear();

   private:
      VldContext fValidity;
      std::string fPrefix;
      std::map<std::string, float> fData;
};

class SetKNNModule : public JobCModule
{
public:
   
   SetKNNModule();
   virtual ~SetKNNModule();
   
   JobCResult Reco(MomNavigator *mom);

   void Config(const Registry &reg); 

   void BeginJob();

private:

   void Analyze(int i, TObject *object, bool print = false);
   bool AnalyzeRecord(NtpStRecord* ntprec);
   
private:
   Anp::Interface *fInterface;

   unsigned int fNPass;
   unsigned int fNFail;

   bool fFastMode;
   bool fPrintEvent;
   bool fPrintTrack;
   bool fStrip;
   std::string fPath;

   bool fIsMRCC;

   Registry fConfig;
};

#endif

