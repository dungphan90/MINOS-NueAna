#include "NueAna/ParticlePID/ParticleAna/PRecord.h"

ClassImp(PRecord)

PRecord::PRecord() : RecRecordImp<RecCandHeader>()
{


}

PRecord::PRecord(const RecCandHeader& hdr) : RecRecordImp<RecCandHeader>(hdr)
{


}


PRecord::~PRecord()
{




}

void PRecord::Reset()
{
	particles.Clear();
        event.Clear();
        mctrue.Clear();
        truthcompare.Clear();
        mrccinfo.Clear();
}



