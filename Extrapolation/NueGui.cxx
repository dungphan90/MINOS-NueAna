#include <iostream>
#include "TColor.h"
#include "NueAna/Extrapolation/NueGui.h"

ClassImp(NueGui)

NueGui::NueGui(const TGWindow *p, UInt_t w, UInt_t h)
  : TGMainFrame(p, w, h)
{

  fDoTotal = 0;
  fDoPrint = 0;
  fSplitBNue = 0;
  fTheta23 = 0;
  fCurSel = Selection::kUnknown;
  fCurBg = Background::kUnknown;
  fCurSys = Systematic::kUnknown;
  fComp = NULL;

  this->SetBackgroundColor(TColor::Number2Pixel(5));

  fLayout = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY);
  fGroupFrameSel = new NueGroupFrame(this,"Selections ");
  fGroupFrameSel->SetBackgroundColor(TColor::Number2Pixel(5));
  fGroupFrameSys = new NueGroupFrame(this,"Systematics");
  fGroupFrameSys->SetBackgroundColor(TColor::Number2Pixel(5));
  fGroupFrameBg  = new NueGroupFrame(this,"Backgrounds");
  fGroupFrameBg->SetBackgroundColor(TColor::Number2Pixel(5));

  Int_t max_sel_index = 0;
  while(strcmp(Selection::AsString(Selection::
				   ESelection(max_sel_index)),"Unknown")!=0){
    Selection::Selection_t sel = Selection::ESelection(max_sel_index);    
    fSelRBut[sel] = new TGCheckButton(fGroupFrameSel,Selection::AsString(sel),150+max_sel_index);
    fSelRBut[sel]->SetBackgroundColor(TColor::Number2Pixel(5));
    fGroupFrameSel->AddFrame(fSelRBut[sel], new TGLayoutHints(kLHintsLeft | kLHintsTop));
    max_sel_index++;
  }
  this->AddFrame(fGroupFrameSel,fLayout);

  fSplitBut = new TGCheckButton(this, "Split BNue", 6);
  fSplitBut->SetBackgroundColor(TColor::Number2Pixel(5));
  this->AddFrame(fSplitBut, fLayout);

  fButton0 = new TGTextButton(this, "&Load", 10);
  fButton0->SetCommand("cout << \"Loading from \" << Selection::AsString(mainWin->GetCurSel()) << \" files\" << endl;");
  fButton0->SetBackgroundColor(TColor::Number2Pixel(5));
  this->AddFrame(fButton0, fLayout);

  Int_t max_bg_index = 0;
  while(strcmp(Background::AsString(Background::
				    EBackground(max_bg_index)),"?Unknown?")!=0){
    Background::Background_t bg = Background::EBackground(max_bg_index);    
    fBgRBut[bg] = new TGCheckButton(fGroupFrameBg,Background::AsString(bg),50+max_bg_index);
    fBgRBut[bg]->SetBackgroundColor(TColor::Number2Pixel(5));
    fGroupFrameBg->AddFrame(fBgRBut[bg],new TGLayoutHints(kLHintsLeft | kLHintsTop));   
    max_bg_index++;
  }
  this->AddFrame(fGroupFrameBg,fLayout);

  fSumBut = new TGCheckButton(this, "All Backgrounds", 4);
  fSumBut->SetBackgroundColor(TColor::Number2Pixel(5));
  this->AddFrame(fSumBut, fLayout);

  Int_t max_sys_index = 0;
  while(strcmp(Systematic::AsString(Systematic::
				    ESystematic(max_sys_index)),"Unknown")!=0) {
    Systematic::Systematic_t sys = Systematic::ESystematic(max_sys_index);    
    fSysRBut[sys] = new TGCheckButton(fGroupFrameSys,Systematic::AsString(sys),100+max_sys_index);
    fSysRBut[sys]->SetBackgroundColor(TColor::Number2Pixel(5));
    fGroupFrameSys->AddFrame(fSysRBut[sys],new TGLayoutHints(kLHintsLeft | kLHintsTop));
    max_sys_index++;
  }
  this->AddFrame(fGroupFrameSys,fLayout);

  fOscBut = new TGCheckButton(this, "Theta23", 7);
  fOscBut->SetBackgroundColor(TColor::Number2Pixel(5));
  this->AddFrame(fOscBut, fLayout);

  fButton1 = new TGTextButton(this, "&Draw Plots", 1);
  fButton1->SetCommand("cout << \"Drawing Plots: Systematic = \" << Systematic::AsString(mainWin->GetCurSys()) << \" Background = \"; if(mainWin->DoTotal()) cout << \"All\"; else cout << Background::AsString(mainWin->GetCurBg()); cout << endl;");
  fButton1->SetBackgroundColor(TColor::Number2Pixel(5));
  this->AddFrame(fButton1, fLayout);
  
  fButton2 = new TGTextButton(this, "&Draw Summary", 2);
  fButton2->SetCommand("cout << \"Drawing Summary\" << endl;");  
  fButton2->SetBackgroundColor(TColor::Number2Pixel(5));
  this->AddFrame(fButton2, fLayout);

  fPrintBut = new TGCheckButton(this, "Print?", 5);
  fPrintBut->SetBackgroundColor(TColor::Number2Pixel(5));
  this->AddFrame(fPrintBut, fLayout);

  fButton3 = new TGTextButton(this, "&Exit", 3);
  fButton3->SetCommand(".q" );
  fButton3->SetBackgroundColor(TColor::Number2Pixel(5));
  this->AddFrame(fButton3, fLayout);

  MapSubwindows();
  Layout();
  SetWindowName("F/N Extrapolation Display");
  SetIconName("F/N Extrapolation Display");
  MapWindow();

}

void NueGui::StartComparator() {

  if(fCurSel!=Selection::kANN30 && fCurSel!=Selection::kSSPID) {
    std::cout << "Selection method: " << Selection::AsString(fCurSel)
	      << " not currently available" << std::endl;
    return;
  }

  if(fComp) delete fComp;

  std::cout << "Processing... " << std::endl;
  
  fComp = new Comparator("RecoEnergy");
  fComp->AddBackground(Background::kNC);
  fComp->AddBackground(Background::kNuMuCC);
  if(fSplitBNue) {
    fComp->AddBackground(Background::kPiBNueCC);
    fComp->AddBackground(Background::kKaBNueCC);
  }
  else fComp->AddBackground(Background::kBNueCC);
  fComp->AddBackground(Background::kNuTauCC);

  Int_t max_sys_index = 0;
  while(strcmp(Systematic::AsString(Systematic::
				    ESystematic(max_sys_index)),"Unknown")!=0) {
    Systematic::Systematic_t sys = Systematic::ESystematic(max_sys_index);    
    //if(sys==Systematic::kShwDev) {max_sys_index++; continue;}
    string file = string(string("EnergySpectraHelperFiles/EnergySpectraHelper_") +
                         string(Selection::AsString(fCurSel)) +
                         string("_Alter") + string(Systematic::AsString(sys)) +
                         string(".root")).c_str();
                                                                                
    ifstream Test(file.c_str());
    if(!Test) std::cout<<"Unable to Find File: "<<file<<std::endl;
    else{
      Test.close();
      fComp->AddSysFile(sys,file);
    }
    max_sys_index++;
  }
  fComp->ComputeAll();
  fComp->DoPrint(fDoPrint);
  if(fTheta23) fComp->SetOscSysString("Theta23");
  else fComp->SetOscSysString("Delta23");
  std::cout << "Done!" << std::endl;
}

NueGui::~NueGui()
{
  delete fComp;
  delete fLayout;
  delete fGroupFrameSel;
  delete fGroupFrameSys;
  delete fGroupFrameBg;
  delete fButton0;
  delete fButton1;
  delete fButton2;
  delete fButton3;
  delete fPrintBut;
  delete fSumBut;
  delete fSplitBut;
  delete fOscBut;

  std::map<Selection::Selection_t,TGCheckButton*>::iterator beg_sel = fSelRBut.begin();
  while(beg_sel!=fSelRBut.end()){
    delete beg_sel->second;
    beg_sel++;
  }
  std::map<Background::Background_t,TGCheckButton*>::iterator beg_bg = fBgRBut.begin();
  while(beg_bg!=fBgRBut.end()){
    delete beg_bg->second;
    beg_bg++;
  }
  std::map<Systematic::Systematic_t,TGCheckButton*>::iterator beg_sys = fSysRBut.begin();
  while(beg_sys!=fSysRBut.end()){
    delete beg_sys->second;
    beg_sys++;
  }
}

Bool_t NueGui::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
  // Process events generated by the buttons in the frame.
  switch (GET_MSG(msg)) {
  case kC_COMMAND:    
    switch (GET_SUBMSG(msg)) {
    case kCM_BUTTON:
      if(parm1==1) {
	if(fDoTotal) fComp->DrawAll(int(fCurSys));
	else fComp->DrawAll(int(fCurSys),int(fCurBg));
      }
      else if(parm1==2) fComp->DrawSummary();
      else if(parm1==10) this->StartComparator();
      break;
    case kCM_CHECKBUTTON:
      if(parm1==4) {
	if(fDoTotal==false) fDoTotal = true;
	else fDoTotal = false;
      }
      else if(parm1==5) {
	if(fDoPrint==false) fDoPrint = true;
	else fDoPrint = false;
	if(fComp) fComp->DoPrint(fDoPrint);
      }
      else if(parm1==6) {
	if(fSplitBNue==false) fSplitBNue = true;
	else fSplitBNue = false;
      }
      else if(parm1==7) {
	if(fTheta23==false) fTheta23 = true;
	else fTheta23 = false;
	if(fTheta23) fComp->SetOscSysString("Theta23");
	else fComp->SetOscSysString("Delta23");
      }
      break;
    default:
      break;
    }
  default:
    break;
  }
  return kTRUE;
}

void NueGui::UnCheckAllBut(Background::Background_t bg)
{
  fCurBg = bg;
  std::map<Background::Background_t,TGCheckButton*>::iterator beg = fBgRBut.begin();
  while(beg!=fBgRBut.end()){
    if(beg->first!=bg) beg->second->SetState(kButtonUp);
    beg++;
  }
}

void NueGui::UnCheckAllBut(Selection::Selection_t sel)
{
  fCurSel = sel;
  std::map<Selection::Selection_t,TGCheckButton*>::iterator beg = fSelRBut.begin();
  while(beg!=fSelRBut.end()){
    if(beg->first!=sel) beg->second->SetState(kButtonUp);
    beg++;
  }
}

void NueGui::UnCheckAllBut(Systematic::Systematic_t sys)
{
  fCurSys = sys;
  std::map<Systematic::Systematic_t,TGCheckButton*>::iterator beg = fSysRBut.begin();
  while(beg!=fSysRBut.end()){
    if(beg->first!=sys) beg->second->SetState(kButtonUp);
    beg++;
  }
}

/////////////////////////////////////////////////////////////////////
ClassImp(NueGroupFrame)

NueGroupFrame::NueGroupFrame(NueGui *p, const char* title) 
  : TGGroupFrame(p,title)
{
}

Bool_t NueGroupFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
  switch (GET_MSG(msg)) {
  case kC_COMMAND:
    switch (GET_SUBMSG(msg)) {
    case kCM_CHECKBUTTON:
      if(parm1 < 100) {
	//this is background
	Background::Background_t bg = Background::EBackground(parm1-50);
	NueGui *ng = (NueGui*) this->GetParent();
	ng->UnCheckAllBut(bg);
      }
      else if (parm1 >= 100 && parm1<150) {
	//this is systematic
	Systematic::Systematic_t sys = Systematic::ESystematic(parm1-100);
	NueGui *ng = (NueGui*) this->GetParent();
	ng->UnCheckAllBut(sys);
      }
      if (parm1 >= 150 && parm1<200) {
	//this is selection
	Selection::Selection_t sel = Selection::ESelection(parm1-150);
	NueGui *ng = (NueGui*) this->GetParent();
	ng->UnCheckAllBut(sel);
      }
      break;
    default:
      break;
    }
  default:
    break;
  }
  return kTRUE;
}
