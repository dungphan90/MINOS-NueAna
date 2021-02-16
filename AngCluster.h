/// $Id: AngCluster.h,v 1.4 2007/05/18 21:01:09 danche Exp $
///
/// class AngCluster
///
/// NueAna package
///
/// Purpose: Return 3D Hit angular clustering variables.
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Thu Apr 28 2005


#ifndef ANGCLUSTER_H
#define ANGCLUSTER_H

#include "TObject.h"

class TH1F;
class TH2F;


class AngCluster : public TObject
{

public:

    // Con/de-structors:
    AngCluster();
    //   AngCluster(const AngCluster& rhs); // copy constructor
    ~AngCluster();

    virtual void Draw(Option_t *option);
    virtual void Print(Option_t *option) const;

    void Reset();

    //Clustering variables
    Float_t fACluRmsShwAxis;
    Float_t fACluRmsZAxis;
    Float_t fACluPrimEnergy;
    Float_t fACluPrimEnergyRatio;
    Float_t fACluShwDirX;
    Float_t fACluShwDirY;
    Float_t fACluShwDirZ;
    Float_t weightedPH0;
    Float_t weightedPH1;
    Float_t weightedPH2;
    Float_t weightedPH3;

private:

    // copy constructor, assignment:
    AngCluster& operator=(const AngCluster& rhs); // assignment

     ClassDef(AngCluster,2);

};                              // end of class AngCluster

#endif  // ANGCLUSTER_H
