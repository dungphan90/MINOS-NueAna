#ifndef ANN_H
#define ANN_H

#include "TObject.h"

class Ann : public TObject
{

public:
   Ann();
   virtual ~Ann();

   void Reset();

   //ann variables
   Double_t pid;
   Double_t pid_30inp;
   Double_t pid_6inp;
   Double_t pid_11inp;
   Double_t pid_11inp_daikon04;
   Double_t pid_14inp_daikon04;

private:

   ClassDef(Ann,5)
};

#endif// ANN_H
