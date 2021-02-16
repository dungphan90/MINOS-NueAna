/**
 * \class BagTreeAna
 *
 * \ingroup NueAna
 *
 * \version $Revision: 1.1 $
 *
 * \date $Date: 2007/10/22 15:59:27 $
 *
 * Created on: Thu Mar 16 20:10:10 2006
 *
 * $Id: BagTreeAna.cxx,v 1.1 2007/10/22 15:59:27 boehm Exp $
 *
 */

#include "NueAna/BagTreeAna.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueAna/NueRecord.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "MessageService/MsgService.h"

CVSID("$Id: BagTreeAna.cxx,v 1.1 2007/10/22 15:59:27 boehm Exp $");

DecisionTreeReader BagTreeAna::heBag;

ClassImp(BagTreeAna)

BagTreeAna::BagTreeAna(NueRecord &nr, BagTree &tp):
    nueRec(nr),
    fbt(tp)
{
   static bool first = true;
   if(first){
//     if(!heBag.CreateFromFile("/afs/fnal.gov/files/home/room3/boehm/NueAnaRuns/TreeReader/Test/BagTree_l5_n500_c5.txt"));

     first = false;
   }



}

BagTreeAna::~BagTreeAna()
{}

void BagTreeAna::Analyze()
{
    NueRecord *nr=&nueRec;

   Double_t AvgIDU,AvgIDV,AvgProbEMU,AvgProbEMV;
   double E2to1U, E2to1V; //, ubeam, vbeam;                                                                               
   double fract_road;
   double fv_3_planes, fv_4_planes, fv_5_planes, fv_6_planes;
   double fv_4_count, fv_8_count, fv_10_count, fv_12_count;
   double abtrackphmean;
   double uv_molrad_sum;
   double uvbeam;
   double uvbeamsum;
   double SSVar1, SSVar2, SSVar3;
   double mstVar1;
   double SSVar4, SSVar5;
   double trkrat, shwrat;
   double ob1, eb1;
                                                                                
   double  input[30];
                                                                                

    int evtplanes = nr->srevent.planes;
    int shwplanes = nr->srshower.planes;
    int trklikeplanes = nr->srtrack.trklikePlanes;
                                                                                
    trkrat = shwrat = 0;
    if(evtplanes > 0){
             trkrat = double(trklikeplanes)/double(evtplanes);
             shwrat = double(shwplanes)/double(evtplanes);
    }
                                                                               
    double t1 = nr->shwfit.UBeamLike;
    double t2 = nr->shwfit.VBeamLike;
    if(t1 < -1000) uvbeam = uvbeamsum = t1;
    else if(t2 < -1000) uvbeam = uvbeamsum = t2;
    else{   uvbeam = t1*t2;   uvbeamsum= t1 + t2;    }
                                                                                
    t2 = nr->shwfit.u_molrad_peak;
    t2 = nr->shwfit.v_molrad_peak;
                                                                               
    if(t1 < -1000) uv_molrad_sum = t1;
    else if(t2 < -1000) uv_molrad_sum = t2;
    else{   uv_molrad_sum = t1 + t2;          }
                                                                                
    fract_road = nr->fracvars.fract_road;
    fv_3_planes = nr->fracvars.fract_3_planes;
    fv_4_planes = nr->fracvars.fract_4_planes;
    fv_5_planes = nr->fracvars.fract_5_planes;
    fv_6_planes = nr->fracvars.fract_6_planes;
    fv_4_count = nr->fracvars.fract_4_counters;
    fv_8_count = nr->fracvars.fract_8_counters;
    fv_10_count = nr->fracvars.fract_10_counters;
    fv_12_count = nr->fracvars.fract_12_counters;
            abtrackphmean = nr->anainfo.abtrackphmean;
                                                                                
            AvgIDU = nr->subshowervars.PHAvgIDU;
            AvgIDV = nr->subshowervars.PHAvgIDV;
            AvgProbEMU = nr->subshowervars.PHAvgProbEMU;
            AvgProbEMV = nr->subshowervars.PHAvgProbEMV;
            E2to1U = nr->subshowervars.E2to1U;
            E2to1V = nr->subshowervars.E2to1V;
            ob1 = nr->mstvars.ob1;
            eb1 = nr->mstvars.eb1;
      SSVar1 = SSVar2 = SSVar3 = SSVar4 = SSVar5 = -9999;
            mstVar1 = -9999;
                                                                                
            if(AvgIDU > -1000 &&  AvgIDV > -1000) SSVar1 = AvgIDU + AvgIDV;
            if(AvgProbEMU > -1000 && AvgProbEMV > -1000){
              SSVar2 = AvgProbEMU*AvgProbEMU + AvgProbEMV*AvgProbEMV;
              SSVar3 = AvgProbEMU + AvgProbEMV;
            }
            if(E2to1U > -1000 && E2to1V > -1000){
              SSVar4 = E2to1U + E2to1V;
              SSVar5 = E2to1V*E2to1V + E2to1U*E2to1U;
            }
                                                                                
     if(ob1 > -1000 && eb1 > -1000)
               mstVar1 = ob1*ob1 + eb1*eb1;
                                                                                
                                                                                
     input[0] =   nr->shwfit.uv_molrad_peak;
     input[1] =   nr->shwfit.uv_kurt;
     input[2] =  nr->fracvars.fract_road;
     input[3] =  nr->fracvars.fract_3_planes;
     input[4] =    fv_4_planes;
     input[5] =    fv_5_planes;
     input[6] =    fv_6_planes;
     input[7] =    fv_4_count;
     input[8] =    fv_8_count;
     input[9] =    fv_10_count;
     input[10] =   fv_12_count;
     input[11] =    nr->shwfit.uv_skew;
     input[12] =  abtrackphmean;
     input[13] =   uv_molrad_sum;
     input[14] =    uvbeam;
     input[15] =  uvbeamsum;
     input[16] =    nr->shwfit.par_b;
     input[17] =   nr->shwfit.uv_rms;
     input[18] =   SSVar1;
     input[19] =   SSVar2;
     input[20] =   SSVar3;
     input[21] =  SSVar4;
     input[22] =   SSVar5;
     input[23] =   mstVar1;
     input[24] =   trkrat;
     input[25] =   shwrat;
     input[26] =    ob1;
     input[27] =    eb1;
 
//     fbt.hePID = heBag.Evaluate(input,28);
                            
}
  
void BagTreeAna::Reset()
{}
