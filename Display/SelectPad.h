#ifndef NUE_SELECTPAD_H
#define NUE_SELECTPAD_H
#include "TPad.h"

class TObject;
class TGaxis;

class SelectPad : public TPad
{

 public:
  SelectPad();
  SelectPad(const char* name, const char* title,
	    double xlow=0, double ylow=0, 
	    double xup=1, double yup=1);
  virtual ~SelectPad();
  
 private:
  
  void ExecuteEvent(Int_t event, Int_t px, Int_t py);
  void ClearSelection(void);
  void ApplySelectionToButtons();
  
  Int_t xmin;
  Int_t xmax;
  Int_t ymin;
  Int_t ymax;
  bool fHaveSelection;
  bool fDrawingSelection;
  bool fInside;
  
  ClassDef(SelectPad,1)
};                              // end of class SelectPad


#endif  // NUE_SELECTPAD_H
