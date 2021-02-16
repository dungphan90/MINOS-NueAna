/**
 *
 * $Id: NueDisplayModule.h,v 1.61 2009/06/30 23:53:42 pawloski Exp $
 *
 * \class NueDisplayModule
 *
 * \package NueAna
 *
 * \a Job Module which makes a display integrated with Midad for Nue analysis.
 *
 * Contact: vahle@hep.ucl.ac.uk  tjyang@stanford.edu
 *
 * Created on: Mon Nov  4 11:53:18 2002
 *
 */

#ifndef NUEDISPLAYMODULE_H
#define NUEDISPLAYMODULE_H

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
class NtpMRRecord;
class NtpMREvent;
class GuiTextButton;
class GuiCheckButton;
class GuiTextEntry;
class GuiBox;
class GuiLabel;
class GuiTextView;
class SelectPad;

const int PHSIZE = 5500;

class NueDisplayModule : public JobCModule
{

public:

  NueDisplayModule();
  virtual ~NueDisplayModule();
  
  typedef std::deque<Float_t> DeqFloat_t; 
  typedef std::deque< std::deque <Int_t> > DeqDeqInt_t;
  
  // Build the display here.
  virtual void BeginRun();
  
  // Update our display here
  JobCResult Ana(const MomNavigator *mom);
  
  // Module configuration
  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);
 
  //For canvas9
  void AdjustOverlayRange(int range); 
  void MRShowMRStrips(int val);
  void MRShowEvent(int val); 
 
 private:
  
  void DeleteOldCrap();

  //each of the following canvases will appear one tab in the main frame
  //main canvas
  TCanvas* fCanvas0;           // We will get this from gMint.
  TPad* fUVview, *fButtonPad, *fHistPad;
  TPad* fInfo0;
  TPad *fStdHepCan;
  TH2D *fUZview, *fVZview;
  TH2D *fTrkUZ, *fTrkVZ;  //track strips
  TH2D *fShwUZ, *fShwVZ;  //shower strips
  TH2D *fSlcUZ, *fSlcVZ;  //slice information
  TNtuple* ftrkshw;
  SteelOutline* fSteelOutline;   //detector outline
  TList * VS; //veto shield
  TPolyLine *pu1_outline; TPolyLine *fu1_outline; 
  TPolyLine *pv1_outline; TPolyLine *fv1_outline;
  TPolyLine *pu2_outline; TPolyLine *fu2_outline; 
  TPolyLine *pv2_outline; TPolyLine *fv2_outline;
  TPolyLine * coil;
  StdHepInfo stdhepinfo;
  StdHepInfoAna shia;
  TButton* fZoom;

  //magnified images canvas
  TCanvas* fCanvas1;
  TPad* fBC, *fHistcolz;
  TPad* fBL, *fHistlego;
  TPad* fStdLeg;
  TH2D *fUZcolz, *fVZcolz;
  TH2D *fUZlego, *fVZlego;
  TButton* cbut[5];
  TButton* lbut[5];
  Int_t ifixc, ifixl;
  //shamelessly stolen from Mad
  float evthighest_z;
  float evtlowest_z;
  float evthighest_t0;
  float evtlowest_t0;
  float evthighest_t1;
  float evtlowest_t1;

  //FracVar display canvas
  TCanvas* fCanvas2;
  TPad* fFracVar_plots, *fFracVar_info;
  TPad* fReco_plots;
  //shamelessly stolen from Mad
  float highest_z;
  float lowest_z;
  float highest_t0;
  float lowest_t0;
  float highest_t1;
  float lowest_t1;

  FracVar fracvars;
  FracVarAna fva;
  
  //Shwfit display canvas
  TCanvas* fCanvas3;
  TPad* fShwfit_plots, *fShwfit_plots_sub, *fShwfit_info;

  Shwfit shwfit;
  ShwfitAna sfa;
  
  HitCalc hitcalc;
  HitCalcAna hca;

  //AngCluster display canvas
  TCanvas* fCanvas4;
  TPad* fAngClusterFitAna_plots;
  
  AngCluster angcluster;
  AngClusterAna aca;
  
  AngClusterFit angclusterfit;
  AngClusterFitAna acfa;
  
  //SubShower canvas
  TCanvas *fCanvas5;
  TMultiGraph *ssGraphU;
  TMultiGraph *ssGraphV;  
  TLegend *cluLegU;
  TLegend *cluLegV;

  TCanvas *fCanvas6;
  TCanvas *fCanvas7;
  SelectPad *fSelectPad1;
  SelectPad *fSelectPad2;
  vector<TButton*> fStripButtonU;
  vector<TButton*> fStripButtonV;
  map<Int_t,Int_t> fStpIndexMapU;
  map<Int_t,Int_t> fStpIndexMapV;
  TArrow *fIntRecoLineU;
  TArrow *fIntRecoLineV;
  TH1F *fIntRecoHistU;
  TH1F *fIntRecoHistV;
  TH1F *fIntRecoHistZ;
  TH1F *fPredictedHistU;
  TH1F *fPredictedHistV;
  TH1F *fPredictedHistZ;
  TH2F *fUZPred;
  TH2F *fVZPred;
  TButton *fIntRecoDoSim;

  //timing etc
  TCanvas *fCanvas8;
  
  //Time Histograms
  THStack *TimeHst;
  TH1F *TimeHstTrk;
  TH1F *TimeHstShw;

  THStack *TimeHstUV;
  TH1F *TimeHstTrkU;
  TH1F *TimeHstTrkV;
  TH1F *TimeHstShwU;
  TH1F *TimeHstShwV;

  THStack *TimeHst2;
  TH1F *TimeHstTrk2;
  TH1F *TimeHstShw2;

  THStack *TimeHst2UV;
  TH1F *TimeHstTrk2U;
  TH1F *TimeHstTrk2V;
  TH1F *TimeHstShw2U;
  TH1F *TimeHstShw2V;

  TGraph* gr_dtds;
  TF1* tfit_dt_ds_pos;
  TF1* tfit_dt_ds_neg;

  TCanvas *fCanvas9;
  TMultiGraph *uzOverlayView;
  TMultiGraph *vzOverlayView;
  TMultiGraph *vzEventOverlay;
  TMultiGraph *uzEventOverlay;
  TMultiGraph *uzSliceOverlay;
  TMultiGraph *vzSliceOverlay;

  TLegend* leg;
  TBox* place[6];
  TButton* fPlusMinusOne;
  TButton* fPlusMinusTwo;
  TButton* fFullSnarl;
  int kShowRange;


  TCanvas* fMRCanvas1;
  TPad* mrInfoPanel;
  TPad* mrGraphPad;

  TMultiGraph * uzMREventOverlay;
  TMultiGraph * vzMREventOverlay;
  TMultiGraph *uzMREventOld;
  TMultiGraph *vzMREventOld;

  //Using MRu for upper part of MR Canvas
  TButton* fMRuShowAll;
  TButton* fMRuShowMR;
  TButton* fMRuShowOld;

  TLegend* fMRUpperLeg;
  TLegend* fMRLowerLeg;

  TButton* fMRdShowAll;
  TButton* fMRdShowTrueMu;
  TButton* fMRdShowTrueShw;
  TButton* fMRdShowScaled;
  TButton* fMRdShowReco;

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

  GuiTextButton* fPrint;
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

  NtpSRRecord *sr;
  NtpMCRecord *mc;
  NtpTHRecord *th;
  NueRecord *nr;
//  NuePID *pid;
  NtpStRecord *st;
  NtpMRRecord *mr;  
  NtpStRecord *stOld;
  NtpSREvent *eventOld;
  NtpMCTruth *mctruthOld;


  NtpSREvent* event;
  NtpMCTruth* mctruth;
  NtpTHEvent* thevent;
    

  bool foundvrmatch;
  bool foundpidmatch;
  bool foundSR;
  bool foundMC;
  bool foundTH;
  bool foundST;
  bool foundMR;
  bool foundSTOld;
 
  TLatex* info1;
  TLatex* info2;
  TLatex* info3;
  TLatex* info4;
  TLatex* info41;
  TLatex* info5;
  TLatex* info6;
  TLatex* info7;
  TLatex* info8;
  TLatex* info9;
  TLatex* info10;
  TLatex* info11;
  TLatex* info12;
  TLatex* info13;

  TLatex* mrInfo1;
  TLatex* mrInfo2;
  TLatex* mrInfo3;
  TLatex* mrInfo4;
  TLatex* mrInfo5;
  TLatex* mrInfo6;
  TLatex* mrInfo7;
  TLatex* mrInfo8;
  TLatex* mrInfo9;  


  //vertex positions
  TMarker* mcvtx_u;
  TMarker* mcvtx_v;
  TMarker* srvtx_u;
  TMarker* srvtx_v;
  TMarker* srvtx_xy;
  
  

  void BuildDisplay();
  void UpdateDisplay(bool, bool);
  void GetBasicInfo();
  Int_t GetEvent(Int_t evt);
  void Analyze(Int_t evt);
  void FillClusterGraphs();
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
  void plotmc();
  void delmc();
  
  int ifixmcinfo;  
  void fixmcinfo();
  
  void PrintPads();

  void GenerateOverlayMultiGraphs();
  void UpdateOverlayGraphs(int evt);
  void UpdateEventOverlayColors(int evt, TMultiGraph *mgr);
  void UpdateSliceOverlayColors(int evt, TMultiGraph *mgr);
  void GenerateMRMultiGraphs(int evtNum,
                        vector<int> &MREvent, vector<int> &OldEvent);
  
  void CreateMRMap(vector<int> &MRnewEvent, vector<int> &MRoldEvent);
  int UpdateMRGraphs(int evt);
  void  FillEventGraph(NtpSREvent* Event, NtpStRecord *str,
          TMultiGraph* uzGraph, TMultiGraph *vzGraph, TString name = "NULL");
  void ResetMRGraphs(); 

  void FillMREventGraph(NtpMREvent* mrevt, NtpStRecord *str,
           TMultiGraph* uzGraph, TMultiGraph *vzGraph, int type);
  void UpdateMREventOverlayColors(TMultiGraph *mgr);
  void UpdateMREventStripColors(TMultiGraph *mgr);
  void SetFrame(TMultiGraph* mg);
  void DrawLowerMRCanvas();
  void DrawUpperMRCanvas();

  void GhostGraph(TGraph *gr);
  void ColorGraph(TGraph *gr, int color, int style = 8, Float_t size = 0.8);
  void UpdateMRInformation(int evt,
                        vector<int> &MREvent, vector<int> &OldEvent);
 
 
                                                  
  bool showStpTrueMu;
  bool showStpTrueShw;
  bool showStpScaled;
  bool showStpRecoTrk;
  bool showStpAll;          

  bool showNewEvent;
  bool showMREvent;
  bool showOrigEvent;


  Int_t kDPlaneCut;
  Int_t kLoPhNStripCut;
  Int_t kLoPhNPlaneCut;
  Float_t kPhStripCut;
  Float_t kPhPlaneCut;
  Float_t kCPhPlaneCut;


  Float_t kPIDCut;  
  Int_t kScanMode; //0: default  1: click 'next event' to record
  Int_t kTestMode; //0: default not a test; 1: test mode i.e. no truth info 
  Int_t kHideRunSnarl; //0: default, don't hide 1: hide
  Int_t kDrawClu;  //0: don't; 1:do
  Int_t kIntReco;  //0: no, 1 yes

  bool foundmeu;
  float SIGMAPMEU;

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

  void vote10(){if (!hitlog){ievtp = 1; updateEncoded();}};
  void vote20(){if (!hitlog){ievtp = 2; updateEncoded();}};
  void vote30(){if (!hitlog){ievtp = 3; updateEncoded();}};
  void vote40(){if (!hitlog){ievtp = 4; updateEncoded();}};
  void vote50(){if (!hitlog){ievtp = 5; updateEncoded();}};
  void vote60(){if (!hitlog){ievtp = 6; updateEncoded();}};
  void vote01(){if (!hitlog){itopo = 1; updateEncoded();}};
  void vote02(){if (!hitlog){itopo = 2; updateEncoded();}};
  void vote03(){if (!hitlog){itopo = 3; updateEncoded();}};
  void vote04(){if (!hitlog){itopo = 4; updateEncoded();}};
  //you can't change your decision after you make it, otherwise you will skew up the record in the file unless a wise way is figured out
  void logvote();
  void makecomment();

  void updateEncoded();
  ofstream outfile;
  string filename;
  void OpenFileWrite();
  int iFileW;  //0: not open  1: opened

  int hitlog; //make a decision before leaving a comment
  int icomm;  //left a comment?

  //////////////////////////////////////////////////////////////

  //preslection
  int passfid;  //in the fiducial volume? 1:pass 0:fail
  int passtrk;  //no long track? 1:pass 0:fail
  int passtrklike;  //trklike track? 1:pass 0:fail
  int passshw;  //high energy shower? 1:pass 0:fail
  int passeng;  //event energy? 1:pass 0:fail

  int preselec;  //do a preselection? 1:yes 0:no
  bool PassCuts();
  void SetCuts();

  ///////////////////////////////////////////////////////////////

  //I/O
  Int_t iIO; //0:write 1:read
  GuiLabel* fFileLab;
  GuiTextEntry* fFileEnt;
  GuiLabel* fTypeLab;
  GuiLabel* fTopoLab;
  GuiCheckButton* Evtp1;
  GuiCheckButton* Evtp2;
  GuiCheckButton* Evtp3;
  GuiCheckButton* Evtp4;
  GuiCheckButton* Evtp5;
  GuiCheckButton* Evtp6;
  GuiCheckButton* Evtp7;
  GuiCheckButton* Topo1;
  GuiCheckButton* Topo2;
  GuiCheckButton* Topo3;
  GuiCheckButton* Topo4;
  GuiCheckButton* Topo5;

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

  //clearly not smart
  void selec10(){ if (!iEvtp[1]){iEvtp[1]=1;}else{iEvtp[1]=0;} updateselec();};
  void selec20(){ if (!iEvtp[2]){iEvtp[2]=1;}else{iEvtp[2]=0;} updateselec();};
  void selec30(){ if (!iEvtp[3]){iEvtp[3]=1;}else{iEvtp[3]=0;} updateselec();};
  void selec40(){ if (!iEvtp[4]){iEvtp[4]=1;}else{iEvtp[4]=0;} updateselec();};
  void selec50(){ if (!iEvtp[5]){iEvtp[5]=1;}else{iEvtp[5]=0;} updateselec();};
  void selec60(){ if (!iEvtp[6]){iEvtp[6]=1;}else{iEvtp[6]=0;} updateselec();};
  void selec70(){ if (!iEvtp[0]){iEvtp[0]=1;}else{iEvtp[0]=0;} updateselec();};
				         	      
  void selec01(){ if (!iTopo[1]){iTopo[1]=1;}else{iTopo[1]=0;} updateselec();};
  void selec02(){ if (!iTopo[2]){iTopo[2]=1;}else{iTopo[2]=0;} updateselec();};
  void selec03(){ if (!iTopo[3]){iTopo[3]=1;}else{iTopo[3]=0;} updateselec();};
  void selec04(){ if (!iTopo[4]){iTopo[4]=1;}else{iTopo[4]=0;} updateselec();};
  void selec05(){ if (!iTopo[0]){iTopo[0]=1;}else{iTopo[0]=0;} updateselec();};

  void updateselec();

  void OpenFileRead();

  void NextSelEvt();
  void PrevSelEvt();

  ReleaseType::Release_t fRel;

  float ph0[PHSIZE];
  float ph1[PHSIZE];

};                              // end of class NueDisplayModule

#endif  // NUEDISPLAYMODULE_H
