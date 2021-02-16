#include <iostream>

using namespace std;

#include "TClonesArray.h"
#include "NueAna/MCNNBestMatch.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "MCNtuple/NtpMCStdHep.h"

ClassImp(MCNNBestMatch)

std::ostream &operator << (std::ostream& os, const MCNNBestMatch& mcnnbmatch);
  
MCNNBestMatch::MCNNBestMatch() :
  stdhep(0)
{
  Reset();
}

MCNNBestMatch::MCNNBestMatch(const MCNNBestMatch *mbm):
  run(mbm->run),
  snarl(mbm->snarl),
  interactionType(mbm->interactionType),
  nuFlavor(mbm->nuFlavor),
  true_enu(mbm->true_enu),
  true_y(mbm->true_y),
  resonanceCode(mbm->resonanceCode),
  dlnL(mbm->dlnL),
  fracQmatched(mbm->fracQmatched),
  stdhep(mbm->stdhep)
{
}

MCNNBestMatch::~MCNNBestMatch()
{
}

void MCNNBestMatch::Clear(Option_t* /* option */)
{
  if(stdhep) {stdhep->Clear("C");}
}

void MCNNBestMatch::Reset()
{
  run = ANtpDefVal::kInt;
  snarl = ANtpDefVal::kInt;
  interactionType = ANtpDefVal::kInt;
  nuFlavor = ANtpDefVal::kInt;
  true_enu = ANtpDefVal::kFloat;
  true_y = ANtpDefVal::kFloat;
  resonanceCode = ANtpDefVal::kInt;
  dlnL = ANtpDefVal::kFloat;
  fracQmatched = ANtpDefVal::kFloat;
  stdhepsize = 0;
  //stdhep = new TClonesArray("NtpMCStdHep");//<-- TClonesArray of NtpMCStdHep
}

void MCNNBestMatch::Print(Option_t * /*option*/) const
{
  std::cout << "Run,snarl: " << run << "," << snarl << std::endl;
  std::cout << "nuFlavor,interactionType: " << nuFlavor << "," << interactionType << std::endl;
  std::cout << "true_enu,y: " << true_enu << "," << true_y << std::endl;
  std::cout << "resonanceCode: " << resonanceCode << std::endl;
  std::cout << "dlnL,fracQmatched: " << dlnL << "," << fracQmatched << std::endl;
  std::cout << "stdhepsize: " << stdhepsize << std::endl;
}
