/**
 *
 * $Id: ParticleDisplay.h,v 1.2 2009/06/23 18:19:15 rhatcher Exp $
 *
 * \class ParticleDisplay
 *
 * \package NueAna
 *
 * \a Job Module which makes a display integrated with Midad for Nue analysis.
 *
 * Contact: scavan@fnal.gov
 *
 * Created on: Mon Nov  4 11:53:18 2002
 *
 */

#ifndef PARTICLEDISPLAY_H
#define PARTICLEDISPLAY_H

#include "TLine.h"
#include "TArrow.h"

#include <JobControl/JobCModule.h>

#include "NueAna/FracVar.h"
#include "NueAna/FracVarAna.h"
#include "NueAna/Shwfit.h"
#include "NueAna/ShwfitAna.h"
#include "NueAna/HitCalc.h"
#include "NueAna/HitCalcAna.h"
#include "NueAna/AngCluster.h"
#include "NueAna/AngClusterAna.h"
#include "NueAna/AngClusterFit.h"
#include "NueAna/AngClusterFitAna.h"
#include "NueAna/NueAnalysisCuts.h"
#include "NueAna/StdHepInfo.h"
#include "NueAna/StdHepInfoAna.h"
#include "Conventions/ReleaseType.h"
#include <Riostream.h>
#include <vector>
#include <map>


#include "NueAna/ParticlePID/ParticleFinder/Display/HoughView.h"
#include "NueAna/ParticlePID/ParticleFinder/Display/ChainView.h"
#include "NueAna/ParticlePID/ParticleFinder/Display/HitView.h"
#include "NueAna/ParticlePID/ParticleFinder/Display/FitView.h"
#include "NueAna/ParticlePID/ParticleFinder/Display/ViewParticle3D.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleAna/PRecord.h"


class TCanvas;
class TNtuple;
class TH2D;
class THStack;
class TPad;
class TText;
class TLatex;
class TLegend;
class TMarker;
class TButton;
class TBox;
class TGraph;
class TMultiGraph;
class TPolyLine;
class MomNavigator;
class SteelOutline;
class NueRecord;
class NtpSRRecord;
class NtpMCRecord;
class NtpTHRecord;
class NtpSREvent;
class NtpMCTruth;
class NtpTHEvent;
class NtpStRecord;

class GuiTextButton;
class GuiCheckButton;
class GuiTextEntry;
class GuiBox;
class GuiLabel;
class GuiTextView;
class SelectPad;

const int PHSIZE = 5500;

class ParticleDisplay : public JobCModule
{

public:

  ParticleDisplay();
  virtual ~ParticleDisplay();
  

  
  // Build the display here.
  virtual void BeginRun();
  
  // Update our display here
  JobCResult Ana(const MomNavigator *mom);
  
  // Module configuration
  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);
 


 
 private:
  






 // ANtpAnalysisInfo objects
//   ANtpAnalysisInfoNue anainfo;
//   ANtpAnalysisInfoAna anaia;
  //////////////////////////////////////////////////////////////////

  //buttons in the left box of the main frame
  GuiTextButton* fNextEventbut;
  GuiTextButton* fPrevEventbut;
  Int_t clickbutton;
  
  GuiLabel* fEventNo;
  GuiTextEntry* fEventEntry;

  GuiTextButton* fCuts; //preselection

  GuiTextButton* fIO;   //I/O

  GuiTextButton* fNextSelEvt;
  GuiTextButton* fPrevSelEvt;

  GuiTextButton* fMCTruth;
  GuiTextButton* fFixMCInfo;

  GuiTextButton* fIntRecoCalc1;
  GuiTextButton* fIntRecoCalc2;
  ///////////////////////////////////////////

  Int_t RunNo;
  Int_t SubRunNo;
  Int_t fSlice;
  Int_t fNumSlices;
  Int_t fEvent;
  Int_t fEvent_old;
  Int_t fNumEvents;
  Int_t fSnarl;
  //used in 'fast scan' style
  Int_t RunNo_old;
  Int_t SubRunNo_old;

  Detector::Detector_t fDetectorType;
  SimFlag::SimFlag_t fSimFlag;
  NueAnalysisCuts fCut;


  ParticleObjectHolder  *poh;
  NtpStRecord * st;

    

  bool foundST;
  bool foundPOH;

 

  void BuildDisplay();
  void UpdateDisplay(bool, bool);

  Int_t GetEvent(Int_t evt);


  void DrawInteractionDiagram(Int_t);
  void SetUpStripButtons();
  void IntRecoCalc();

  void NextEvent();
  void PrevEvent();
  void GotoEvent();

  int imctruth;
  vector<TLine *> paru;
  vector<TLine *> parv;
  void showmctruth();

  
  int ifixmcinfo;  
  void fixmcinfo();
  
  void GenerateOverlayMultiGraphs();
  void UpdateOverlayGraphs(int evt);
  void UpdateEventOverlayColors(int evt, TMultiGraph *mgr);
  void UpdateSliceOverlayColors(int evt, TMultiGraph *mgr);
   void  FillEventGraph(NtpSREvent* Event, ParticleObjectHolder *str,
          TMultiGraph* uzGraph, TMultiGraph *vzGraph, TString name = "NULL");
   void SetFrame(TMultiGraph* mg);
 
  void GhostGraph(TGraph *gr);
  void ColorGraph(TGraph *gr, int color, int style = 8, Double_t size = 0.8);
 
                                                  
  bool showStpTrueMu;
  bool showStpTrueShw;
  bool showStpScaled;
  bool showStpRecoTrk;
  bool showStpAll;          

  bool showNewEvent;
  bool showOrigEvent;


  Int_t kDPlaneCut;
  Int_t kLoPhNStripCut;
  Int_t kLoPhNPlaneCut;
  Double_t kPhStripCut;
  Double_t kPhPlaneCut;
  Double_t kCPhPlaneCut;


  Double_t kPIDCut;  
  Int_t kScanMode; //0: default  1: click 'next event' to record
  Int_t kTestMode; //0: default not a test; 1: test mode i.e. no truth info 
  Int_t kHideRunSnarl; //0: default, don't hide 1: hide
  Int_t kDrawClu;  //0: don't; 1:do
  Int_t kIntReco;  //0: no, 1 yes

  bool foundmeu;
  double SIGCORRMEU;

  //layouts below the main display
  //              fHBox
  //fVBox1 | fVBox2 | fVBox3 | fHBox2 
  //                  fHBox3         
  //		      fHBox4         
  //                  ......         
  //                  fHBox8         
  //                  fHBox9         
  GuiBox* fHBox;
  GuiBox* fVBox1;
  GuiBox* fVBox2;
  GuiBox* fVBox3;
  GuiBox* fHBox1; //not used
  GuiBox* fHBox2;
  GuiBox* fHBox3;
  GuiBox* fHBox4;
  GuiBox* fHBox5;
  GuiBox* fHBox6;
  GuiBox* fHBox7;
  GuiBox* fHBox8;
  GuiBox* fHBox9;

  //voting stuff
  GuiTextButton* evtp1;
  GuiTextButton* evtp2;
  GuiTextButton* evtp3;
  GuiTextButton* evtp4;
  GuiTextButton* evtp5;
  GuiTextButton* evtp6;
  GuiTextButton* topo1;
  GuiTextButton* topo2;
  GuiTextButton* topo3;
  GuiTextButton* topo4;
  GuiTextButton* vote;
  GuiLabel* fEncoded;
  GuiTextEntry* fComment;
  GuiLabel* fComLab;
  GuiTextView* fRecoInfo;
  Int_t ievtp;
  Int_t itopo;



  Int_t iEvtp[7];
  Int_t iTopo[5];

  Int_t selecevtp;
  Int_t selectopo;


  vector<Int_t> runno;
  vector<Int_t>::iterator runnoitr;
  vector<Int_t> subrunno;
  vector<Int_t>::iterator subrunnoitr;
  vector<Int_t> snarlno;
  vector<Int_t>::iterator snarlnoitr;
  vector<Int_t> eventno;
  vector<Int_t>::iterator eventnoitr;
  vector<Int_t> eventtype;
  vector<Int_t>::iterator typeitr;
  vector<Int_t> eventtopo;
  vector<Int_t>::iterator topoitr;
  vector<string> eventcomm;
  vector<string>::iterator commitr;

  GuiLabel* fFileR;
  GuiLabel* fTYPETOPO;

  GuiLabel* fFileW;

  //duplications
  GuiTextButton* fNextEventbut1;
  GuiTextButton* fPrevEventbut1;

  void SetMode();



  void updateselec();

  void OpenFileRead();

  void NextSelEvt();
  void PrevSelEvt();

  ReleaseType::Release_t fRel;

  double ph0[PHSIZE];
  double ph1[PHSIZE];


	ChainView chainview;
	HitView hitview;
	FitView fitview;
	ViewParticle3D view3d;
	HoughView houghview;


 std::vector<TObject* > hft;


};                              // end of class ParticleDisplay

#endif  // PARTICLEDISPLAY_H
