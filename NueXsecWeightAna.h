#ifndef NUEXSECWEIGHTANA_H
#define NUEXSECWEIGHTANA_H

#include "NueAna/NueAnaBase.h"
#include "NueAna/NueXsecWeight.h"

class Registry;
class MCReweight;

class NueXsecWeightAna: public NueAnaBase
{

 public:
  NueXsecWeightAna(NueXsecWeight &fw);
  virtual ~NueXsecWeightAna();
  void SetMCReweight(MCReweight *m) {mcr=m;}
  void SetRegistry(Registry *r) {rwtreg=r;}
  void Analyze(int evt, RecRecordImp<RecCandHeader> *st);

 private:
  NueXsecWeight &fNueXsecWeight;
  Registry *rwtreg;
  MCReweight *mcr;

};

#endif// NUEXSECWEIGHTANA_H
