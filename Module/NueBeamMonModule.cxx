#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "NueAna/Module/NueBeamMonModule.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NuePOT.h"

#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h"
#include "Validity/VldContext.h"
#include "Conventions/Detector.h"

#include "BeamDataUtil/BDSpillAccessor.h"
#include "BeamDataUtil/BeamMonSpill.h"
#include "BeamDataUtil/BMSpillAna.h"
#include "SpillTiming/SpillTimeFinder.h"

#include "DataUtil/MCInfo.h"
#include "DcsUser/CoilTools.h"
#include "DataUtil/GetTempTags.h"

#include "NueAna/NueStandard.h"
#include "TSystem.h"

#include <fstream>
#include <iomanip>
#include <vector>

JOBMODULE(NueBeamMonModule, "NueBeamMon",
          "Fill the beam monitoring data in the NueAna ntuples");
CVSID("$Id: NueBeamMonModule.cxx,v 1.16 2008/11/19 18:22:51 rhatcher Exp $");


NueBeamMonModule::NueBeamMonModule()
  :fBMSpillAna(),firstevent(true), zbeamtype(0), lastrun(-1)
{
  beamtype = BeamType::kUnknown;
  nuepot = new NuePOT();
  fBMSpillAna.UseDatabaseCuts();
}

NueBeamMonModule::~NueBeamMonModule()
{
  if(nuepot){
    delete nuepot;
    nuepot=0;
  }
}

const Registry& NueBeamMonModule::DefaultConfig() const
{
    static Registry r;
    r.UnLockValues();
    r.Set("ZarkoBeamType",0);
    r.Set("POTTreeFileName","pottree.root");
    r.Set("BeamString", "");
    r.LockValues();

    return r;
}

void NueBeamMonModule::Config(const Registry& r)
{
    fBMSpillAna.Config(r);
    const char* tmps;
    if(r.Get("POTTreeFileName",tmps)){kPOTTreeName=tmps;}
//    int imps;
//    if(r.Get("ZarkoBeamType",imps)) {beamtype=imps;}

}

JobCResult NueBeamMonModule::Reco(MomNavigator* mom) 
{

    // Instantiate a BDSpillAccessor object, necessary to access the
    // BeamMonSpill database.
    BDSpillAccessor& sa = BDSpillAccessor::Get();

    // Same to access SpillTimeND table 
    SpillTimeFinder &stf = SpillTimeFinder::Instance();    

    // Get the nue records from mom -- basically copied from the
    // Trimmer class
    TObject *obj=0;
    vector<NueRecord *> records;

    int nnuerec=0;
    bool contains_mc=false;
    TIter objiter = mom->FragmentIter();
    while(( obj=objiter.Next() )){
        NueRecord *nr = dynamic_cast<NueRecord *>(obj);

        if (!nr){
            MSG("NueBeamMon",Msg::kDebug)<<"Didn't find a NueRecord in MOM"<<endl;
            continue;
        }        
       MSG("NueBeamMon",Msg::kDebug)<<"Found a NueRecord in MOM"
                                      <<" Snarl "<<nr->GetHeader().GetSnarl()
                                      <<" Event "<<nr->GetHeader().GetEventNo()<<endl;      
        
        VldContext evt_vldc = nr->GetHeader().GetVldContext();
        ReleaseType::Release_t rel = nr->GetHeader().GetRelease();
        beamtype = nr->GetHeader().GetBeamType();

	if(firstevent){
	  firstevent=false;
	  nuepot->beamtype=beamtype;
          zbeamtype = BeamType::ToZarko(beamtype);
	}

	if(nr->GetHeader().GetRun()!=lastrun){
          nuepot->nruns++;
          lastrun = nr->GetHeader().GetRun();
          if(evt_vldc.GetSimFlag() != SimFlag::kData &&
               evt_vldc.GetDetector()==Detector::kFar){
               double temp = MCInfo::GetMCPoT(Detector::kFar, beamtype, rel);
               nuepot->pot += temp;
               nuepot->pot_nocut += temp; 
          }
        }
        if(nnuerec==0){
          nuepot->nsnarls++;
        }
                                                                                
        // Only fill if this is data
        if (evt_vldc.GetSimFlag() != SimFlag::kData){
            contains_mc=true;
            if(nnuerec==0&&evt_vldc.GetDetector()==Detector::kNear){
              nuepot->pot+=MCInfo::GetMCPoT(Detector::kNear, beamtype, rel);
            }
            nnuerec++;
            continue;
        }

        // Get the closest spill out of the database
        const BeamMonSpill* bms = sa.LoadSpill(evt_vldc.GetTimeStamp());
        VldTimeStamp bms_vts;
        if (!bms) {
            MSG("NueBeamMon",Msg::kError) << "No BeamMonSpill found for " << evt_vldc << endl;
            bms_vts = VldTimeStamp::GetEOT();
        }
        else {
            bms_vts = bms->SpillTime();
        }

        fBMSpillAna.SetSpill(*bms);
        // Get the time of the event, the value in the header is
        // accurate enough for this purpose
        fBMSpillAna.SetSnarlTime(evt_vldc.GetTimeStamp());

        // First reset the values of the BeamMon branch
        nr->bmon.Reset();
        // Fill in the variables
        if (fBMSpillAna.SelectSpill())
            nr->bmon.goodBeamMon=1;
        else
            nr->bmon.goodBeamMon=0;

        nr->bmon.tortgt = bms->fTortgt;
        nr->bmon.trtgtd = bms->fTrtgtd;
        nr->bmon.tor101 = bms->fTor101;
        nr->bmon.tr101d = bms->fTr101d;

//        cout<<fBMSpillAna.SelectSpill()<<"  "<<"   "<<bms->BeamType()<<endl;
        nr->bmon.bI = nr->bmon.GetPot();
        
        double xbpm,ybpm,xrms,yrms;
        bms->BpmAtTarget(xbpm,ybpm,xrms,yrms);
        nr->bmon.hpos2 = xbpm;
        nr->bmon.vpos2 = ybpm;

        for (Int_t i=0;i<6;++i){
            nr->bmon.batchposx[i]=bms->fTargBpmX[i];
            nr->bmon.batchposy[i]=bms->fTargBpmY[i];
            nr->bmon.batchint[i]=bms->fBpmInt[i];
        }

        nr->bmon.hpos1 = bms->fTargProfX;
        nr->bmon.vpos1 = bms->fTargProfY;
        
        nr->bmon.hbw = bms->fProfWidX;
        nr->bmon.vbw = bms->fProfWidY;

        nr->bmon.hornI = bms->fHornCur;

        nr->bmon.bmst_vts = bms->SpillTime();
        
        VldTimeStamp vtsdif = nr->bmon.bmst_vts-evt_vldc.GetTimeStamp();
        nr->bmon.dt_bmst=vtsdif.GetSeconds();

        nr->bmon.stnd_vts=stf.GetTimeOfNearestSpill(evt_vldc);
        nr->bmon.dt_stnd=stf.GetTimeToNearestSpill(evt_vldc);

        // count the pots for the good spills, only take the value
        // from the first event in the snarl
        bool goodCoil = CoilTools::IsOK(evt_vldc) && !CoilTools::IsReverse(evt_vldc);    
        goodCoil = goodCoil || (evt_vldc.GetDetector()==Detector::kFar);

        nr->bmon.goodDataQual = (int) NueStandard::PassesPOTStandards(nr);  
        // includes the goodBeamMon, SpillType, CoilCuts and FarDet special cuts 

        if (nnuerec==0 && nr->bmon.goodDataQual == 1){
//            fTotPot+=nr->bmon.bI;
	    nuepot->pot+=nr->bmon.bI;
        }
        if(nnuerec==0) nuepot->pot_nocut +=nr->bmon.bI;        
        ++nnuerec;
    }
    if (!contains_mc && nnuerec==0)
        MSG("NueBeamMon",Msg::kWarning)<<"No NueRecord found in MOM"<<endl;
    
    return JobCResult::kPassed;
}

void NueBeamMonModule::EndJob(){
   MSG("NueBeamMon",Msg::kInfo)<<"Number of POT in this job: "<<nuepot->pot<<endl;

   TDirectory *savedir = gDirectory;

   TFile *fpf = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->FindObject(kPOTTreeName.c_str()));
   if(fpf){
     fpf->cd();
     TTree *pottree = new TTree("pottree","pottree");
     pottree->Branch("NuePOT",&nuepot);
     pottree->Fill();
     pottree->Write();
     savedir->cd();
   }
   else{
     MSG("NueBeamMon",Msg::kError)<<"Could not find the file to write the pottree to, there will be no pottree"<<endl;
   }

}
