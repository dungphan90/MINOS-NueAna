#include "NueAna/Display/SelectPad.h"
#include "TButton.h"
#include "TList.h"
#include <iostream>
using namespace std;

ClassImp(SelectPad)

SelectPad::SelectPad()

    : TPad("sPad","sPad",0,0,1,1)
    , fHaveSelection(false)
    , fDrawingSelection(false)
    , fInside(false)
{
}


SelectPad::SelectPad(const char* name,  const char* title,
		     double xlow, double ylow, double xup, double yup)

    : TPad(name,title,xlow,ylow,xup,yup)
    , fHaveSelection(false)
    , fDrawingSelection(false)
    , fInside(false)
{
}

SelectPad::~SelectPad()
{
}

void SelectPad::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
   switch (event) {

   case kButton1Down:
       this->SetCursor(kHand);
       fInside = 1;
       xmin = px;
       ymin = py;
       break;

   case kButton1Motion:
     fInside=1;
     if (fDrawingSelection)
       gVirtualX->DrawBox(xmin,ymin,xmax,ymax,
			  TVirtualX::kHollow);
       xmax = px;
       ymax = py;
       fDrawingSelection = true;
       gVirtualX->DrawBox(xmin,ymin,xmax,ymax,
                          TVirtualX::kHollow);
       break;
       
   case kButton1Up: {
       this->SetCursor(kPointer);

       Int_t dtp = this->DistancetoPrimitive(xmax,ymax);

       if (!fDrawingSelection || dtp || !fInside) {
           cerr
               << "fDrawingSelection = " << fDrawingSelection
               << ", DistanceToPrimitive = " << dtp
               << ", fInside = " << fInside
               << endl;
	   fDrawingSelection = false;
	   fInside = false;
           this->ClearSelection();
           break;
       }

       fHaveSelection = true;

       this->ApplySelectionToButtons();
       this->Update();
       fInside = false;
       fDrawingSelection = false;
       this->ClearSelection();
       break;
   }
   case kMouseLeave:
       fInside = false;
       break;
   case kMouseEnter:
       fInside = true;
       break;
   default:
       this->TPad::ExecuteEvent(event,px,py);
       break;
   }
}

void SelectPad::ClearSelection(void)
{
    fHaveSelection = false;    
    xmin = 0;
    ymin = 0;
    xmax = 0;
    ymax = 0;
}

void SelectPad::ApplySelectionToButtons()
{
  if (!fHaveSelection) return;
  
  double x1 = this->AbsPixeltoX(xmin);
  double x2 = this->AbsPixeltoX(xmax);
  double y1 = this->AbsPixeltoY(ymin);
  double y2 = this->AbsPixeltoY(ymax);
  
  if (x1 > x2) { double tmp = x1; x1 = x2; x2 = tmp; }
  if (y1 > y2) { double tmp = y1; y1 = y2; y2 = tmp; }

  TList *theList = this->GetListOfPrimitives();
  TIterator *iter = theList->MakeIterator();
  TObject *ob = NULL;
  while((ob = iter->Next())){
    if(ob->InheritsFrom("TButton")) {
      TButton *but = (TButton*) ob;
      double bx1,bx2,by1,by2;
      but->GetPadPar(bx1,by1,bx2,by2);
      if (bx1 > bx2) { double tmp = bx1; bx1 = bx2; bx2 = tmp; }
      if (by1 > by2) { double tmp = by1; by1 = by2; by2 = tmp; }      
      if(((bx1>x1 && bx1<x2) || (bx2>x1 && bx2<x2)) &&
	 ((by1>y1 && by1<y2) || (by2>y1 && by2<y2))){
	//if(bx1>x1&&bx2<x2&&by1>y1&&by2<y2){
	but->ExecuteEvent(kButton1Down,this->XtoAbsPixel((bx1+bx2)/2.),
			  this->YtoAbsPixel((by1+by2)/2.));
	but->ExecuteEvent(kButton1Up,this->XtoAbsPixel((bx1+bx2)/2.),
			  this->YtoAbsPixel((by1+by2)/2.));
      }
    }
  }
}


