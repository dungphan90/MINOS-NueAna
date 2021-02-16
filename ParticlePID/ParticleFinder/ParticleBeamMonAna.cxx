#include "NueAna/ParticlePID/ParticleFinder/ParticleBeamMonAna.h"

#include "Conventions/BeamType.h"
#include "Conventions/Detector.h"
#include "Conventions/ReleaseType.h"

#include "BeamDataUtil/BDSpillAccessor.h"
#include "BeamDataUtil/BeamMonSpill.h"
#include "BeamDataUtil/BMSpillAna.h"
#include "SpillTiming/SpillTimeFinder.h"

#include "MessageService/MsgService.h"
#include "DcsUser/CoilTools.h"
#include "TSystem.h"

#include "DataUtil/GetTempTags.h"


CVSID("$Id: ParticleBeamMonAna.cxx,v 1.3 2013/07/31 04:42:21 rhatcher Exp $");



ParticleBeamMonAna::ParticleBeamMonAna()
{
	fname="";
}

ParticleBeamMonAna::~ParticleBeamMonAna()
{}


void ParticleBeamMonAna::ana(ParticleObjectHolder * poh, ParticleBeamMon * e)
{

    BDSpillAccessor& sa = BDSpillAccessor::Get();
    SpillTimeFinder &stf = SpillTimeFinder::Instance();


    
    BeamType::BeamType_t beamtype;
    int zbeamtype; 
    
    BMSpillAna fBMSpillAna;




	VldContext evt_vldc = poh->GetHeader().GetVldContext();



    std::string file = DataUtil::GetTempTagString(poh, "file");
    std::string filename = gSystem->BaseName(file.c_str());



    if(evt_vldc.GetSimFlag() == SimFlag::kData)
      beamtype = DetermineBeamType(evt_vldc);
    else
      beamtype = DetermineBeamType(filename);

    //last ditch effort...
    	if(beamtype==BeamType::kUnknown)
    	{
       		beamtype = DetermineBeamType(fname);   	
    	}
    
    e->beamtype=beamtype;


    // Only fill if this is data
    if (evt_vldc.GetSimFlag() != SimFlag::kData)
        return;
   








    const BeamMonSpill* bms = sa.LoadSpill(evt_vldc.GetTimeStamp());
    VldTimeStamp bms_vts;
    if (!bms) {
    	MSG("NueParticleBeamMon",Msg::kError) << "No ParticleBeamMonSpill found for " << evt_vldc << endl;
    	bms_vts = VldTimeStamp::GetEOT();
    }else{
        bms_vts = bms->SpillTime();
    }
  
 //   beamtype = bms->BeamType();
    
//    e->beamtype=beamtype;
    
    zbeamtype = BeamType::ToZarko(beamtype);
  
  
        
    fBMSpillAna.SetSpill(*bms);
    // Get the time of the event, the value in the header is
    // accurate enough for this purpose
    fBMSpillAna.SetSnarlTime(evt_vldc.GetTimeStamp());

        e->ResetAll();                     // Fill in the variables
        if (fBMSpillAna.SelectSpill())
            e->goodParticleBeamMon=1;
        else
            e->goodParticleBeamMon=0;

        e->tortgt = bms->fTortgt;
        e->trtgtd = bms->fTrtgtd;
        e->tor101 = bms->fTor101;
        e->tr101d = bms->fTr101d;

//        cout<<fBMSpillAna.SelectSpill()<<"  "<<"   "<<bms->BeamType()<<endl;
        e->bI = e->GetPot();

        double xbpm,ybpm,xrms,yrms;
        bms->BpmAtTarget(xbpm,ybpm,xrms,yrms);
        e->hpos2 = xbpm;
        e->vpos2 = ybpm;

        for (Int_t i=0;i<6;++i){
            e->batchposx[i]=bms->fTargBpmX[i];
            e->batchposy[i]=bms->fTargBpmY[i];
            e->batchint[i]=bms->fBpmInt[i];
        }

        e->hpos1 = bms->fTargProfX;
        e->vpos1 = bms->fTargProfY;

        e->hbw = bms->fProfWidX;
        e->vbw = bms->fProfWidY;

        e->hornI = bms->fHornCur;

        e->bmst_vts = bms->SpillTime();

        VldTimeStamp vtsdif = e->bmst_vts-evt_vldc.GetTimeStamp();
        e->dt_bmst=vtsdif.GetSeconds();

        e->stnd_vts=stf.GetTimeOfNearestSpill(evt_vldc);
        e->dt_stnd=stf.GetTimeToNearestSpill(evt_vldc);

        // count the pots for the good spills, only take the value
        // from the first event in the snarl
        bool goodCoil = CoilTools::IsOK(evt_vldc) && !CoilTools::IsReverse(evt_vldc);
        goodCoil = goodCoil || (evt_vldc.GetDetector()==Detector::kFar);

     
	
}


BeamType::BeamType_t ParticleBeamMonAna::DetermineBeamType(VldContext vc){

/*    BDSpillAccessor& sa = BDSpillAccessor::Get();
    VldContext evt_vldc = nr->GetHeader().GetVldContext();

    const BeamMonSpill* bms = sa.LoadSpill(vc.GetTimeStamp());
    VldTimeStamp bms_vts;
    if (!bms) {
       MSG("NueBeamMon",Msg::kError) << "No BeamMonSpill found for " << evt_vldc << endl;
       bms_vts = VldTimeStamp::GetEOT();
    }
    else {
        bms_vts = bms->SpillTime();
    }
    return bms->BeamType();
*/
    int time = vc.GetTimeStamp().GetSec();
    
    //printf("time %d\n",time);

    BeamType::BeamType_t beam = BeamType::kUnknown;

    if(time >= 1107216000 && time < 1109539850) beam = BeamType::kL100z200i;
    if(time >= 1109540615 && time < 1109899325) beam = BeamType::kL250z200i;
    if(time >= 1109899938 && time < 1110239564) beam = BeamType::kL100z200i;
    if(time >= 1110323763 && time < 1111622400) beam = BeamType::kL000z200i;
    if(time >= 1114892377 && time < 1115927583) beam = BeamType::kL100z200i;
    if(time >= 1115937438 && time < 1116604821) beam = BeamType::kL250z200i;
    if(time >= 1116618256 && time < 1122659668) beam = BeamType::kL010z185i;
    if(time >= 1122659886 && time < 1122922688) beam = BeamType::kL010z170i;
    if(time >= 1122922890 && time < 1123112674) beam = BeamType::kL010z200i;
    if(time >= 1123112803 && time < 1139605423) beam = BeamType::kL010z185i;
    if(time >= 1139605543 && time < 1140022084) beam = BeamType::kL010z000i;
    if(time >= 1140026702 && time < 1140908579) beam = BeamType::kL010z185i;
    // End of Run 1
    if(time >= 1149180600 && time < 1150047780) beam = BeamType::kL150z200i;
    if(time >= 1150047780 && time < 1151690460) beam = BeamType::kL250z200i;
    if(time >= 1153956600 && time < 1155510000) beam = BeamType::kL250z200i;
    if(time >= 1158004800 && time < 1158019870) beam = BeamType::kL010z200i;
    if(time >= 1158019870 && time < 1161892800) beam = BeamType::kL010z185i;
    if(time >= 1161892800 && time < 1184351737) beam = BeamType::kL010z185i;
    if(time >= 1184351737 && time < 1184708040) beam = BeamType::kL010z000i;
    //End of Run 2

    if(time >= 1184800000 ) beam = BeamType::kL010z185i;

    return beam;

}

BeamType::BeamType_t ParticleBeamMonAna::DetermineBeamType(std::string file)
{

//printf("looking in %s\n",file.c_str());
   BeamType::BeamType_t beam = BeamType::kUnknown;
   if(file.find("L010185")!=string::npos){ beam=BeamType::kL010z185i; }
   if(file.find("L100200")!=string::npos){ beam=BeamType::kL100z200i; }
   if(file.find("L250200")!=string::npos){ beam=BeamType::kL250z200i; }
   if(file.find("L150200")!=string::npos){ beam=BeamType::kL150z200i; }
   if(file.find("L010200")!=string::npos){ beam=BeamType::kL010z200i; }
   if(file.find("L010170")!=string::npos){ beam=BeamType::kL010z170i; }
   if(file.find("L010000")!=string::npos){ beam=BeamType::kL010z000i; }

   return beam;
}



