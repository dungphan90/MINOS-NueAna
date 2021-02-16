/**
 * \class TreePIDAna
 *
 * \ingroup NueAna
 *
 * \brief Calculate PID from Tree or Cut methods and perform classification.
 *
 * Author: Mayly Sanchez (msanchez@minos.phy.tufts.edu)
 *
 * \version $Revision: 1.1 $
 *
 * \date $Date: 2006/03/20 06:59:13 $
 *
 * Created on: Thu Mar 16 19:54:12 2006
 *
 * $Id: TreePIDAna.h,v 1.1 2006/03/20 06:59:13 msanchez Exp $
 *
 */

#ifndef TREEPIDANA_H
#define TREEPIDANA_H

#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/TreePID.h"

class NueRecord;
//class TreePID;

class TreePIDAna
{

public:

    TreePIDAna(NueRecord &nr, TreePID &tp);
    virtual ~TreePIDAna();
    
    void Analyze();
    void Reset();

private: 
    NueRecord &nueRec;
    TreePID &fTreePID;

    TreePIDAna(const TreePIDAna& rhs); //copy constructor
    TreePIDAna& operator=(const TreePIDAna& rhs); //assignment

    ClassDef(TreePIDAna,1)
};
// 
#endif  // TREEPIDANA_H

