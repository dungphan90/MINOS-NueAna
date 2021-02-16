#include<iostream>
#include "NueAna/ParticlePID/ParticleFinder/ParticleBeamMon.h"

#include "DataUtil/GetTempTags.h"

ClassImp(ParticleBeamMon)

ParticleBeamMon::ParticleBeamMon(const RecCandHeader& hdr) : RecRecordImp<RecCandHeader>(hdr)
{
	ResetAll();	
}


ParticleBeamMon::ParticleBeamMon() : RecRecordImp<RecCandHeader>()
{
	ResetAll();
}

ParticleBeamMon::~ParticleBeamMon()
{}

void ParticleBeamMon::ResetAll()
{
    goodParticleBeamMon=-1;
    bI=-9999.;
    tortgt=-9999.;
    trtgtd=-9999.;
    tor101=-9999.;
    tr101d=-9999.;
    hpos2=-9999.;
    vpos2=-9999.;
    for (Int_t i=0;i<6;++i){
        batchposx[i]=-9999.;
        batchposy[i]=-9999.;
        batchint[i]=-9999.;
    }
    hpos1=-9999.;
    vpos1=-9999.;
    hbw=-9999.;
    vbw=-9999.;
    htan=-9999.;
    vtan=-9999.;
    hornI=-9999.;
    nuTarZ=-9999.;
    time=0; 
    bmst_vts=0;
    stnd_time=0;
    stnd_vts=0;
    dt_bmst=-9999;
    dt_stnd=-9999;
    goodDataQual = -9999;
    beamtype=BeamType::kUnknown;
}

void ParticleBeamMon::PrintInfo(Option_t */*t*/) const
{
    std::cout <<" goodParticleBeamMon " << goodParticleBeamMon << std::endl
              <<" bI "<<bI<<std::endl
              <<" tortgt "<<tortgt<<std::endl
              <<" tortgt "<<trtgtd<<std::endl
              <<" tor101 "<<tor101<<std::endl
              <<" tr101d "<<tr101d<<std::endl
              <<" hbw "<<hbw<<std::endl
              <<" vbw "<<vbw<<std::endl
              <<" hpos1 "<<hpos1<<std::endl
              <<" vpos1 "<<vpos1<<std::endl
              <<" hpos2 "<<hpos2<<std::endl
              <<" vpos2 "<<vpos2<<std::endl
              <<" htan "<<htan<<std::endl
              <<" vtan "<<vtan<<std::endl
              <<" hornI "<<hornI<<std::endl
              <<" nuTarZ "<<nuTarZ<<std::endl
              <<" time "<<time<<std::endl
              <<" goodDataQual"<<goodDataQual<<std::endl;
}


double ParticleBeamMon::GetPot()
{
    double pot = tortgt;
    if (pot==0) pot = trtgtd;
    if (pot==0) pot = tor101;
    if (pot==0) pot = tr101d;
    return pot;
}

