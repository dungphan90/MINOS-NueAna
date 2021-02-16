#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "TH1F.h"
#include "TFile.h"


#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro

#include "TClass.h"
#include "TDataType.h"
#include "TDataMember.h"
#include "TRealData.h"
#include "NueAna/ParticlePID/ParticleFinder/Systematic/SystematicGains.h"
#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"

#include "TClonesArray.h"

#include "CandNtupleSR/NtpSRStrip.h"  
#include "CandNtupleSR/NtpSREvent.h" 
#include "CandNtupleSR/NtpSRShower.h" 
#include "CandNtupleSR/NtpSRTrack.h" 

#include "TObject.h"


JOBMODULE(SystematicGains, "SystematicGains",
          "Adjust strip energies for a given gain adjustment");
          
CVSID("$Id: SystematicGains.cxx,v 1.1 2009/06/19 14:32:40 scavan Exp $");


int SystematicGains::printit=-1;

//......................................................................

SystematicGains::SystematicGains()
{
        if(printit<0)
        {
                if(MsgService::Instance()->IsActive("SystematicGains",Msg::kDebug))printit=1;
                else printit=0;
        }
	allshift=0;
	rndshift=0;
}

//......................................................................

SystematicGains::~SystematicGains()
{}

//......................................................................
void SystematicGains::BeginJob()
{

}


void SystematicGains::EndJob()
{
}

JobCResult SystematicGains::Reco(MomNavigator* mom)
{
   bool foundST=false;
                                                                               
//   VldContext vc;
   //first ask mom for a NtpStRecord

//iterate over all NtpStRecords in MOM

std::vector<TObject*> records = mom->GetFragmentList("NtpStRecord");


for(unsigned int i=0;i<records.size();i++)
{

   NtpStRecord *str=dynamic_cast<NtpStRecord *>(records[i]); //(mom->GetFragment("NtpStRecord","Primary"));
   if(str){
     foundST=true;
  //   vc=str->GetHeader().GetVldContext();
   }else{
        continue;
   }
   //cout<<"Running over "<<str->GetName()<<endl;

        NtpSREventSummary * evthdr = dynamic_cast<NtpSREventSummary *>(&str->evthdr);

        for(int ievt=0;ievt<evthdr->nevent;ievt++)
        {
        



   			NtpSREvent* event = dynamic_cast<NtpSREvent*>(str->evt->At(ievt));  
   			int nstrips=event->nstrip;

			for(int m=0;m<nstrips;m++)
            {                               
            	NtpSRStrip* strip = (NtpSRStrip*)(str->stp->At(event->stp[m]));
				if(!strip)
                {
                	printf("CANT FIND STRIP AT IDX %d\n",event->stp[m]);
                	continue;
                }
                                
				        
		double adjfrac0=0;
		double adjfrac1=0;
		if(strip->ph0.pe)
		{
			adjfrac0 = 1. + allshift ;
			if(rndshift)adjfrac0 += rnd.Gaus(0,rndshift);		                        
                	strip->ph0.pe*=adjfrac0;
				strip->ph0.siglin*=adjfrac0;
         	        strip->ph0.sigcor*=adjfrac0;
				strip->ph0.raw*=adjfrac0;
                }

		if(strip->ph1.pe)
		{
			adjfrac1 = 1. - allshift ;
			if(rndshift)adjfrac1 += rnd.Gaus(0,rndshift);
                	strip->ph1.pe*=adjfrac1;
				strip->ph1.siglin*=adjfrac1;
                	strip->ph1.sigcor*=adjfrac1;
				strip->ph1.raw*=adjfrac1;
                }
                if(printit)printf("adjusting size 0 by %f and side 1 by %f\n", adjfrac0, adjfrac1);
                                        
                                        
           }
           

        }
}
   
   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}




const Registry& SystematicGains::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("SystematicGains",Msg::kDebug)<<"In Trimmer::DefaultConfig"<<endl;

  
  static Registry r;
  std::string name = this->JobCModule::GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  r.UnLockValues();
  r.Set("TotalShift",(double)0.);
  r.Set("RandomShift",(double)0.);          

  r.LockValues();             
               
                                                                              
  return r;
}

void SystematicGains::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("SystematicGains",Msg::kDebug)<<"In Trimmer::Config"<<endl;
                                                                               
  bool islocked = this->GetConfig().ValuesLocked();
  if (islocked) this->GetConfig().UnLockValues();
 
  double valint;
  if(r.Get("TotalShift", valint)) 
  {
    allshift=valint;
  }

  if(r.Get("RandomShift", valint)) 
  {
    rndshift=valint;
  }

  if (islocked) this->GetConfig().LockValues();

  cout<< "SystematicGains   TotalShift "<<allshift<<" RandomShift "<<rndshift<<endl;

}





