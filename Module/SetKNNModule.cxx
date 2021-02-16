// $Id: SetKNNModule.cxx,v 1.5 2011/02/17 21:16:42 whitehd Exp $

// C/C++
#include <cassert>
#include <iomanip>

// ROOT
#include "TClonesArray.h"

// MINOS
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "JobControl/JobCModuleRegistry.h"
#include "JobControl/JobCommand.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "StandardNtuple/NtpStRecord.h"

// Local
#include "PhysicsNtuple/Default.h"
#include "PhysicsNtuple/Factory.h"
#include "PhysicsNtuple/Store/FillBasic.h"
#include "PhysicsNtuple/Store/Interface.h"

// Local
#include "SetKNNModule.h"

CVSID("$Id: SetKNNModule.cxx,v 1.5 2011/02/17 21:16:42 whitehd Exp $");

JOBMODULE(SetKNNModule, "SetKNNModule", "SetKNNModule");

//---------------------------------------------------------------------------------------------
StorekNNData::StorekNNData()
{
}

//---------------------------------------------------------------------------------------------
StorekNNData::~StorekNNData()
{
}

//---------------------------------------------------------------------------------------------
bool StorekNNData::Add(int i, const std::string &key, const float data)
{
   char temp[3];
   sprintf(temp, "%02d", i);
   const std::string full_key = fPrefix + "_" + key + "_" + string(temp);
   
//   if(true) cout<<"Storing: "<<full_key<<"  "<<data<<endl;

   return fData.insert(map<string, float>::value_type(full_key, data)).second;
}

//---------------------------------------------------------------------------------------------
bool StorekNNData::SetValidity(VldContext &vld)
{
   fValidity = vld;
   return true;
}

//---------------------------------------------------------------------------------------------
const VldContext& StorekNNData::GetValidity()
{ 
   return fValidity;
} 

//---------------------------------------------------------------------------------------------
bool StorekNNData::SetPrefix(std::string in)
{
   fPrefix = in;
   return true;
} 

//---------------------------------------------------------------------------------------------
bool StorekNNData::Get(int i, const std::string &key, float &data)
{
   char temp[3];
   sprintf(temp, "%02d", i);
   const std::string full_key = fPrefix + "_" + key + "_" + string(temp);

   map<string, float>::const_iterator fit = fData.find(full_key);
   if(fit != fData.end())
   {
      data = fit -> second;
      return true;
   }

   return false;
}

//---------------------------------------------------------------------------------------------
bool StorekNNData::Clear()
{
//   cout<<"flushing the data!"<<endl;
   fData.clear();
   return true;
}  

//---------------------------------------------------------------------------------------------
SetKNNModule::SetKNNModule()
   :fInterface(new Anp::Interface()),
    fNPass(0),
    fNFail(0),
    fPrintEvent(false),
    fPrintTrack(false),
    fStrip(false),
    fIsMRCC(false),
    fConfig(false)
{

//   fPath = "junk/afs/fnal.gov/files/home/room1/rustem/data/muon-knn/";
   fPath = "/minos/app/nue/Releases/Griffin/Data/kNN_cedar/";
   fFastMode = false;

   Anp::Factory<StorekNNData>::
      Instance().Hold("kNNData", Anp::Handle<StorekNNData>(new StorekNNData));
}

//---------------------------------------------------------------------------------------------
SetKNNModule::~SetKNNModule() 
{
   //
   // Delete interface object
   //
   delete fInterface;
   fInterface = 0;

   //
   // Destroy StorekNNData if it exists
   //
   Anp::Factory<StorekNNData>::Instance().Remove("kNNData");

   MSG("SetKNNModule", Msg::kDebug) 
      << endl
      << "**************************************************" << std::endl
      << "    SetKNNModule" << std::endl
      << "      Number of passed records " << fNPass << std::endl
      << "      Number of failed records " << fNFail << std::endl
      << "**************************************************" << std::endl;
}

//---------------------------------------------------------------------------------------------
JobCResult SetKNNModule::Reco(MomNavigator *mom)
{  
  bool foundST=false;
  std::vector<TObject*> records = mom->GetFragmentList("NtpStRecord");

  if(fFastMode){
    std::vector<TObject*> nrRecord =  mom->GetFragmentList("NueRecord");
    if(nrRecord.size() == 0) return  JobCResult::kFailed;
  }

  bool pass = true;

  static bool first = true;

  Anp::Handle<StorekNNData> data = Anp::Factory<StorekNNData>::Instance().Get("kNNData");

  if(!data.valid())
  {
      cerr << " SetKNNModule::Analyze - Handle<StorekNNData> is invalid" << endl;
      return JobCResult::kFailed;
  }
  data->Clear();

  for(unsigned int i=0;i<records.size() && pass; i++)
  {
      NtpStRecord *str=dynamic_cast<NtpStRecord *>(records[i]);
      if(str){  foundST=true;  }
      else continue;
   
      if(first){
          Registry reg(false);
          string pub = getenv("SRT_PUBLIC_CONTEXT");
          string config = pub + "/PhysicsNtuple/Config/Config2008Real.txt";
          reg.Set("InterfaceConfigPath", config.c_str());
 
          string file = fPath;


          int release = str->GetRelease();
          bool diakon00 = ( ReleaseType::IsDaikon(release) && 0==ReleaseType::GetMCSubVersion(release) );

          if(str->GetHeader().GetVldContext().GetDetector() == Detector::kNear)
          {
           if(diakon00)
           {
             file += "knn.physics.near.daikon_00.cedar_phy.L010z185i.root";
           }
           else
           {
             file += "knn.physics.near.daikon_04.cedar_phy_bhcurv.L010z185i.root";
           }
          }


          if(str->GetHeader().GetVldContext().GetDetector() == Detector::kFar)
          {
           if(diakon00)
           {
             file += "knn.physics.far.daikon_00.cedar_phy.L010z185i.root";
           }
           else
           {
             file += "knn.physics.far.daikon_04.cedar_phy_bhcurv.L010z185i.root";
           }
          }
 
            reg.Set("FillkNNFilePath", file.c_str());

          fInterface->Config(reg);
          first = false;
       }    


      if(records.size() > 1 && !fIsMRCC){
         fIsMRCC = true;
         MSG("SetKNNModule", Msg::kInfo)<<"Multiple NtpStRecords Detected - Shifting to MRCC Mode"<<endl;
      }

      pass = pass && AnalyzeRecord(str);
  }

  if(!pass) return JobCResult::kFailed;
  return JobCResult::kAOK;
}
   
bool SetKNNModule::AnalyzeRecord(NtpStRecord* ntprec)
{
   Anp::Handle<StorekNNData> data = Anp::Factory<StorekNNData>::Instance().Get("kNNData");
   if(!data.valid())
   {
      cerr << " SetKNNModule::Analyze - Handle<StorekNNData> is invalid" << endl;
      return false;
   }

   std::string name = ntprec->GetName();
   data->SetPrefix("SNTP");
   if(fIsMRCC && name == "Primary") data->SetPrefix("OldSNTP");

   VldContext vld = ntprec->GetHeader().GetVldContext();
   data->SetValidity(vld);

   // Fill interface with data for a new snarl
   if(fInterface -> FillSnarl(ntprec))
   {
      ++fNPass;
   }
   else
   {
      ++fNFail;
      return false; 
   } 

   TClonesArray *event_array = ntprec -> evt;
   TClonesArray *track_array = ntprec -> trk;
   if(!event_array || !track_array)
   {           
      MSG("SetKNNModule", Msg::kWarning) << "Invalid TClonesArray object(s)" << endl;
      return false;
   }

   // Iterate over all events and access analysis variables
   for(int ievent = 0; ievent < event_array -> GetEntries(); ++ievent)
   {
      NtpSREvent *event = dynamic_cast<NtpSREvent *>(event_array -> At(ievent));
      if(!event)
      {
	 MSG("SetKNNModule", Msg::kError) << "NtpSREvent dynamic_cast failed" << endl;
	 continue;
      }

      SetKNNModule::Analyze(ievent, event, fPrintEvent);
      assert(event -> index == ievent && "mismatched index");
   }

   return true;
}

//---------------------------------------------------------------------------------------------
void SetKNNModule::Config(const Registry& reg)
{
   MSG("SetKNNModule", Msg::kVerbose) << "SetKNNModule::Config()..." << std::endl;

   //
   // Store Registry copy
   //
   fConfig.UnLockValues();
   fConfig.UnLockKeys();
   fConfig.Merge(reg);
   fConfig.LockKeys();
   fConfig.LockValues();

   const char* tmps;
   if(reg.Get("kNNPath", tmps)) fPath = tmps;
   
   int itmp;
   if(reg.Get("ForceMRCC", itmp)){
     if(itmp == 1){
        fIsMRCC = true;
         MSG("SetKNNModule", Msg::kInfo)<<"Forced to Shift to MRCC Mode"<<endl;
     }
   }

   if(reg.Get("RunFast", itmp)){
     if(itmp == 1){
        fFastMode = true;
         MSG("SetKNNModule", Msg::kInfo)<<"Shifting to Fast Mode"<<endl;
     }
   }


   Anp::Read(reg, "SetKNNModulePrintEvent", fPrintEvent);
   Anp::Read(reg, "SetKNNModulePrintTrack", fPrintTrack);
   Anp::Read(reg, "SetKNNModuleStrip", fStrip);
}

//---------------------------------------------------------------------------------------------
void SetKNNModule::BeginJob()
{
   MSG("SetKNNModule", Msg::kVerbose) << "SetKNNModule::BeginJob()..." << std::endl;

   //
   // Configure interface object at the beginning of analysis job
   //
//   fInterface -> Config(fConfig);
}

//------------------------------------------------------------------------------------------
void SetKNNModule::Analyze(int i, TObject *object, const bool print)
{
   //
   // Extract variables from PhysicsNtuple interface
   //
   const float numubar = fInterface -> GetVar("numubar", object);
   const float rel_ang = fInterface -> GetVar("rel_ang", object);
   const float knn_pid = fInterface -> GetVar("knn_pid", object);
   const float knn_01  = fInterface -> GetVar("knn_01",  object);
   const float knn_10  = fInterface -> GetVar("knn_10",  object);
   const float knn_20  = fInterface -> GetVar("knn_20",  object);
   const float knn_40  = fInterface -> GetVar("knn_40",  object);

   if(print)
   {
      cout << " For event: "<<i
           << "   anti neutrino selection = " << numubar << endl
	   << "   relative angle = " << rel_ang << endl
	   << "   knn muon pid = " << knn_pid << endl
	   << "   number of scintillator planes = " << knn_01 << endl
	   << "   mean track signal = " << knn_10 << endl
	   << "   low mean signal over high mean signal = " << knn_20 << endl
	   << "   track mean signal over track window mean signal = " << knn_40 << endl;
   }
   
   Anp::Handle<StorekNNData> data = Anp::Factory<StorekNNData>::Instance().Get("kNNData");

   if(!data.valid())
   {
      cerr << " SetKNNModule::Analyze - Handle<StorekNNData> is invalid" << endl;
      return;
   }

   data -> Add(i, "numubar", numubar);
   data -> Add(i, "rel_ang", rel_ang);
   data -> Add(i, "knn_pid", knn_pid);
   data -> Add(i, "knn_01", knn_01);
   data -> Add(i, "knn_10", knn_10);
   data -> Add(i, "knn_20", knn_20);
   data -> Add(i, "knn_40", knn_40);
}
