/**
 *
 * $Id: TimingVars.cxx,v 1.3 2007/03/15 17:00:56 boehm Exp $
 *
 * \class TimingVars
 *
 * \package NueAna
 *
 * \brief Hold variables related to the TimingVars package
 **/
#include <iostream>
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "NueAna/TimingVars.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "MessageService/MsgService.h"

CVSID("$Id: TimingVars.cxx,v 1.3 2007/03/15 17:00:56 boehm Exp $");

ClassImp(TimingVars)

TimingVars::TimingVars():

   ShwMaxBin(ANtpDefaultValue::kInt),
   TrkMaxBin(ANtpDefaultValue::kInt),
   TotalMaxBin(ANtpDefaultValue::kInt),
   ShwMaxBinCont(ANtpDefaultValue::kFloat),
   TrkMaxBinCont(ANtpDefaultValue::kFloat),
   TotalMaxBinCont(ANtpDefaultValue::kFloat),
   ShwSecondMaxBin(ANtpDefaultValue::kFloat),
   TrkSecondMaxBin(ANtpDefaultValue::kFloat),
   TotalSecondMaxBin(ANtpDefaultValue::kFloat),
   ShwEventPulseHeight(ANtpDefaultValue::kFloat),
   TrkEventPulseHeight(ANtpDefaultValue::kFloat),
   TotalEventPulseHeight(ANtpDefaultValue::kFloat),
   ShwPercentofTotal(ANtpDefaultValue::kFloat),
   TrkPercentofTotal(ANtpDefaultValue::kFloat),
   TotalPercentofTotal(ANtpDefaultValue::kFloat),
   ShwbinsPrior(ANtpDefaultValue::kInt),
   TrkbinsPrior(ANtpDefaultValue::kInt),
   TotalbinsPrior(ANtpDefaultValue::kInt),
   ShwBinsOver40(ANtpDefaultValue::kInt),
   TrkBinsOver40(ANtpDefaultValue::kInt),
   TotalBinsOver40(ANtpDefaultValue::kInt)
//   lenepl(0),
//   tenestu(0),
//   tenestv(0),
//   efit(0),
//   hfit(0),
//   info1(0)
{}

TimingVars::TimingVars(const TimingVars &s):
  TObject(),

   ShwMaxBin(s.ShwMaxBin),
   TrkMaxBin(s.TrkMaxBin),
   TotalMaxBin(s.TotalMaxBin),
   ShwMaxBinCont(s.ShwMaxBinCont),
   TrkMaxBinCont(s.TrkMaxBinCont),
   TotalMaxBinCont(s.TotalMaxBinCont),
   ShwSecondMaxBin(s.ShwSecondMaxBin),
   TrkSecondMaxBin(s.TrkSecondMaxBin),
   TotalSecondMaxBin(s.TotalSecondMaxBin),
   ShwEventPulseHeight(s.ShwEventPulseHeight),
   TrkEventPulseHeight(s.TrkEventPulseHeight),
   TotalEventPulseHeight(s.TotalEventPulseHeight),
   ShwPercentofTotal(s.ShwPercentofTotal),
   TrkPercentofTotal(s.TrkPercentofTotal),
   TotalPercentofTotal(s.TotalPercentofTotal),
   ShwbinsPrior(s.ShwbinsPrior),
   TrkbinsPrior(s.TrkbinsPrior),
   TotalbinsPrior(s.TotalbinsPrior),
   ShwBinsOver40(s.ShwBinsOver40),
   TrkBinsOver40(s.TrkBinsOver40),
   TotalBinsOver40(s.TotalBinsOver40)

{
/*
  if(s.lenepl!=0){
    lenepl=new TH1F(*(s.lenepl));
  }
  else{
    lenepl=0;
  }
  if(s.tenestu!=0){
    tenestu=new TH1F(*(s.tenestu));
  }
  else{
    tenestu=0;
  }
  if(s.tenestv!=0){
    tenestv=new TH1F(*(s.tenestv));
  }
  else{
    tenestv=0;
  }
  if(s.efit!=0){
    efit=new TF1(*(s.efit));
  }
  else{
    efit=0;
  }
  if(s.hfit!=0){
    hfit=new TF1(*(s.hfit));
  }
  else{
    hfit=0;
  }
  if(s.info1!=0){
    info1=new TPaveText(*(s.info1));
  }
  else{
    info1=0;
  }
*/
}


TimingVars::~TimingVars()
{
/*
   if(lenepl!=0){
      delete lenepl;
      lenepl=0;
   }
   if(tenestu!=0){
      delete tenestu;
      tenestu=0;
   }
   if(tenestv!=0){
      delete tenestv;
      tenestv=0;
   }

   if(efit!=0){
      delete efit;
      efit=0;
   }
   if(hfit!=0){
      delete hfit;
      hfit=0;
   }
   if(info1!=0){
      delete info1;
      info1=0;
   }
*/
}

void TimingVars::Reset()
{
   ShwMaxBin = ANtpDefaultValue::kInt;
   TrkMaxBin = ANtpDefaultValue::kInt;
   TotalMaxBin = ANtpDefaultValue::kInt;
   ShwMaxBinCont = ANtpDefaultValue::kFloat;
   TrkMaxBinCont = ANtpDefaultValue::kFloat;
   TotalMaxBinCont = ANtpDefaultValue::kFloat;
   ShwSecondMaxBin = ANtpDefaultValue::kFloat;
   TrkSecondMaxBin = ANtpDefaultValue::kFloat;
   TotalSecondMaxBin = ANtpDefaultValue::kFloat;
   ShwEventPulseHeight = ANtpDefaultValue::kFloat;
   TrkEventPulseHeight = ANtpDefaultValue::kFloat;
   TotalEventPulseHeight = ANtpDefaultValue::kFloat;
   ShwPercentofTotal = ANtpDefaultValue::kFloat;
   TrkPercentofTotal = ANtpDefaultValue::kFloat;
   TotalPercentofTotal = ANtpDefaultValue::kFloat;
   ShwbinsPrior = ANtpDefaultValue::kInt;
   TrkbinsPrior = ANtpDefaultValue::kInt;
   TotalbinsPrior = ANtpDefaultValue::kInt;
   ShwBinsOver40 = ANtpDefaultValue::kInt;
   TrkBinsOver40 = ANtpDefaultValue::kInt;
   TotalBinsOver40 = ANtpDefaultValue::kInt;

}

