#ifndef NUEGUI_H
#define NUEGUI_H
#include <map>
#include <TGClient.h>
#include <TGButton.h>
//#include <TGGroupFrame.h>
#include "NueAna/NueAnaTools/Selection.h"
#include "NueAna/Extrapolation/Systematic.h"
#include "NueAna/Extrapolation/Background.h"
#include "NueAna/Extrapolation/Comparator.h"

class NueGroupFrame;

class NueGui : public TGMainFrame 
{

private:

  NueGroupFrame *fGroupFrameSel,*fGroupFrameSys,*fGroupFrameBg;
  TGTextButton    *fButton0,*fButton1,*fButton2,*fButton3;
  std::map<Selection::Selection_t,TGCheckButton*> fSelRBut;
  std::map<Systematic::Systematic_t,TGCheckButton*> fSysRBut;
  std::map<Background::Background_t,TGCheckButton*> fBgRBut;  
  TGCheckButton   *fSumBut;
  TGCheckButton   *fPrintBut;
  TGCheckButton   *fSplitBut;
  TGCheckButton   *fOscBut;
  TGLayoutHints   *fLayout;
  
  Comparator *fComp;
  
  Bool_t fDoTotal;
  Bool_t fDoPrint;
  Bool_t fSplitBNue;
  Bool_t fTheta23;
  Selection::Selection_t fCurSel;
  Background::Background_t fCurBg;
  Systematic::Systematic_t fCurSys;
  
  void StartComparator();
  
 public:
  
  NueGui(const TGWindow *p, UInt_t w, UInt_t h);
  ~NueGui();
  Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

  Bool_t DoTotal() { return fDoTotal;}
  Selection::Selection_t GetCurSel() {return fCurSel;}
  Systematic::Systematic_t GetCurSys() {return fCurSys;}
  Background::Background_t GetCurBg() {return fCurBg;}

  void UnCheckAllBut(Background::Background_t);
  void UnCheckAllBut(Systematic::Systematic_t);
  void UnCheckAllBut(Selection::Selection_t);

  ClassDef(NueGui,0)

};


class NueGroupFrame : public TGGroupFrame
{

 private:
  
 public:
  NueGroupFrame(NueGui* p, const char*);
  ~NueGroupFrame() {}
  Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

  ClassDef(NueGroupFrame,0)
};


#endif  //NUEGUI_H
