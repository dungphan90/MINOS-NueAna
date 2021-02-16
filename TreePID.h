/**
 * \class TreePID
 *
 * \ingroup NueAna
 *
 * \brief Calculate PID from Tree or Cut methods and perform classification.
 *
 * Contact: Mayly Sanchez (msanchez@physics.harvard.edu)
 *
 * \version $Revision: 1.4 $
 *
 * \date $Date: 2007/05/25 20:58:47 $
 *
 * Created on: Thu Mar 16 19:35:58 2006
 *
 * $Id: TreePID.h,v 1.4 2007/05/25 20:58:47 msanchez Exp $
 *
 */

#ifndef TREEPID_H
#define TREEPID_H

#include "TObject.h"

class TreePID : public TObject
{

public:

    TreePID();
    virtual ~TreePID();

    virtual void Draw(Option_t *option);
    virtual void Print(Option_t *option) const; 
    
    void Reset();

    //CutPID variables
    Int_t fPassBCuts;
    Int_t fCutPID;
    Int_t fCutPID1;
    Int_t fCutPID2;
    Int_t fCutPID3;
 
    Int_t fCutClass;
    Int_t fCutClass1;
    Int_t fCutClass2;
    Int_t fCutClass3;
    //TreePID variables

private:
    TreePID& operator=(const TreePID& rhs); //assignment
    
    ClassDef(TreePID,2);

};
// 
#endif  // TREEPID_H


