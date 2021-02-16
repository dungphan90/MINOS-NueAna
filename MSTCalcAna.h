#ifndef MSTCALCANA_H
#define MSTCALCANA_H

#include <vector>
#include <set>
#include "NueAna/NueAnaBase.h"
#include "NueAna/MSTCalc.h"
#include "NueAna/NueAnaTools/DCGraph.h"
#include "NueAna/NueAnaTools/DCHit.h"

class TH3F;
class TH1;
class TProfile;
class TProfile2D;
class TFile;
class UberRecord;

using std::vector;
using std::set;

class MSTCalcAna: public NueAnaBase
{

public:
   MSTCalcAna(MSTCalc &mst);
   virtual ~MSTCalcAna();

   //   void Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord *mc=0, NtpTHRecord *th=0);
   void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);


   //so I can still use this on CalDet data
   void Analyze(UberRecord *ur);

   void ComputeGraphQuantities();
   void SetSigTemplate(TProfile2D *e) {etemplate=e;}
   void SetBGTemplate(TProfile2D *b) {btemplate=b;}

   void SetSigMIPTemplate(TProfile2D *e) {emtemplate=e;}
   void SetBGMIPTemplate(TProfile2D *b) {bmtemplate=b;}
   void SetMSTParams(float minsc, float maxm, 
		     float minfsc, float maxmlowz, float scm){
     minsig=minsc;
     maxmetric=maxm;
     minfarsigcor=minfsc;
     maxmetriclowz=maxmlowz;
     SetParams(scm);
   }

   void FindAllGraphs(set<DCHit> &evenhistset,set <DCHit> &oddhitset);
   
   //   void FillHitSets(int evtn, NtpSRRecord *srobj, 
   //		    set<DCHit> &evenhitset, set<DCHit> &oddhitset);

   void FillHitSets(int evtn, RecRecordImp<RecCandHeader> *srobj, 
		    set<DCHit> &evenhitset, set<DCHit> &oddhitset);


   //so I can still use on caldet data
   void FillHitSets(UberRecord *ur,
		    set<DCHit> &evenhitset, set<DCHit> &oddhitset);


   void DrawMaxGraph(int view);
   DCGraph<DCHit> *GetMaxGraph(int view);
   TH1D *FindProjection(TProfile2D *tmplt, int nedges, const char* name);
   double ComputeChi2(TH1 *h1, TH1 *h2, 
		      double &chi2, int &ndof);
   double ComputeLLike(TH1 *data, TH1 *tmpl, int &ndof);
   float FindAlpha(TH1* data, TH1 *ehist, TH1 *bhist);
   void FindLikeAlpha(TH1* data, TH1 *ehist, TH1 *bhist,
		      double &alpha, double &beta);

   void GraphLoop(DCGraph<DCHit> *g,int view);

   void Print() const;
   void Reset();
private:
   MSTCalc &fMSTCalc;
   float minsig;
   float maxmetric;
   float minfarsigcor;
   float maxmetriclowz;
     
   vector<DCGraph<DCHit> *>even; //!
   vector<DCGraph<DCHit> *>odd; //!

   TProfile2D *etemplate; //!
   TProfile2D *btemplate; //!
   TProfile2D *emtemplate; //!
   TProfile2D *bmtemplate; //!

};
#endif//MSTCALCANA_H
