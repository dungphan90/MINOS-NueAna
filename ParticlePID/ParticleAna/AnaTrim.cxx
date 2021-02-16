////////////////////////////////////////////////////////////////////////
/// $Id: AnaTrim.cxx,v 1.2 2009/06/23 22:42:24 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "NueAna/ParticlePID/ParticleAna/AnaTrim.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro

#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "NueAna/ParticlePID/ParticleAna/PRecordAna.h"




#include <math.h>


JOBMODULE(AnaTrim, "AnaTrim",
          "does the ana trim for ParticleFinder");
CVSID("$Id: AnaTrim.cxx,v 1.2 2009/06/23 22:42:24 scavan Exp $");
//......................................................................

AnaTrim::AnaTrim() 
{




	
  ///
  /// (Document me!)
  ///
}
//......................................................................

AnaTrim::~AnaTrim()
{
  ///
  /// (Document me!)
  ///


}

//......................................................................

void AnaTrim::BeginJob()
{
  ///
  /// (Document me!)
  ///
}

//......................................................................

void AnaTrim::EndJob()
{
  ///
  /// (Document me!)
  ///

  
}

//......................................................................

JobCResult AnaTrim::Reco(MomNavigator* mom)
{
  ///
  /// (Document me!)
  ///

  std::vector<TObject* > hft =( mom->GetFragmentList("PRecord"));




	for(unsigned int s =0;s<hft.size();s++)
	{
		PRecord *n=0;
		n=dynamic_cast<PRecord *>(hft[s]);	

//		int infid=1;
//	
//		infid = infid && (n->event.vtx_z>0.47692 && n->event.vtx_z<14.6 || n->event.vtx_z>16.0 && n->event.vtx_z <29.8);
//		infid = infid && (sqrt(n->event.vtx_u*n->event.vtx_u+n->event.vtx_v*n->event.vtx_v)<4);
//		infid = infid && (sqrt(n->event.vtx_u*n->event.vtx_u+n->event.vtx_v*n->event.vtx_v)>0.5);
	
	
	
		int passcuts=1;
		
//		passcuts=passcuts && infid;
//		passcuts=passcuts && n->particles.ntot>1;
//		passcuts=passcuts && n->particles.totvise > 25.;
//		passcuts=passcuts && n->particles.totvise < 25.*8.;
			
	
		if(!passcuts)continue;
			n->SetName("TrimmedPA");
		
	//	n->mctrue.trainweight *=3.5/6.5/10;
	
	

	}

		


  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................

const Registry& AnaTrim::DefaultConfig() const
{
  ///
  /// Supply the default configuration for the module
  ///
  static Registry r; // Default configuration for module

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());
 // r.Set("POTTreeFileName","pottree.root");
  // Set values in configuration
  r.UnLockValues();
  r.LockValues();


  return r;
}

//......................................................................

void AnaTrim::Config(const Registry& /*r*/)
{
  ///
  /// Configure the module given the Registry r
  ///
  
 
   // const char* tmps;
    //if(r.Get("POTTreeFileName",tmps)){kPOTTreeName=tmps;}
  
}

//......................................................................

void AnaTrim::Reset()
{
  ///
  /// (Document me!)
  ///
}

////////////////////////////////////////////////////////////////////////

