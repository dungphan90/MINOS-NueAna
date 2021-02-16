/**
 * \class BagTreeAna
 *
 * \ingroup NueAna
 *
 * \brief Calculate PID from Tree or Cut methods and perform classification.
 *
 * Author: Mayly Sanchez (msanchez@minos.phy.tufts.edu)
 *
 * \version $Revision: 1.1 $
 *
 * \date $Date: 2007/10/22 15:59:27 $
 *
 * Created on: Thu Mar 16 19:54:12 2006
 *
 * $Id: BagTreeAna.h,v 1.1 2007/10/22 15:59:27 boehm Exp $
 *
 */

#ifndef BAGTREEANA_H
#define BAGTREEANA_H

#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/BagTree.h"
#include "NueAna/NueAnaTools/DecisionTreeReader.h"

class NueRecord;
//class TreePID;

class BagTreeAna
{

public:

    BagTreeAna(NueRecord &nr, BagTree &tp);
    virtual ~BagTreeAna();
    
    void Analyze();
    void Reset();

private: 
    NueRecord &nueRec;
    BagTree &fbt;

    BagTreeAna(const BagTreeAna& rhs); //copy constructor
    BagTreeAna& operator=(const BagTreeAna& rhs); //assignment

    static DecisionTreeReader heBag;

    ClassDef(BagTreeAna,1)
};
// 
#endif  // TREEPIDANA_H

