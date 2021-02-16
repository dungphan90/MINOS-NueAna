#ifndef FRACVARANA_H
#define FRACVARANA_H

#include "NueAna/NueAnaBase.h"
#include "NueAna/FracVar.h"
#include "TMultiLayerPerceptron.h"

class TPad;
class TNtuple;

class NtpSRRecord;
class NtpSREvent;
class FracVarAna : public NueAnaBase
{

 public:

  FracVarAna(FracVar &fv);
  virtual ~FracVarAna();

  void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
  bool LoadLargestTrackFromEvent(NtpSREvent* event,
				 RecRecordImp<RecCandHeader>* record,
				 int& trkidx);
  bool LoadShowerAtTrackVertex(NtpSREvent* event,
			       RecRecordImp<RecCandHeader>* record,
			       int  trkidx,
			       int& shwidx);
  bool LoadLargestShowerFromEvent(NtpSREvent* event,
				  RecRecordImp<RecCandHeader>* record,
				  int& shwidx);
  void SetDisplay(Int_t fd){fDisplay = fd;};
  void Draw(TPad *pad);
  void Print(TPad *pad);

 private:
  FracVar &fFracVar;
  Int_t fDisplay;
  TNtuple *display;
  static TMultiLayerPerceptron *fneuralNet;
};

#endif// FRACVARANA_H
