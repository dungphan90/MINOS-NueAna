//////////////////////////////////////////////////////////////////////////////
// TrimEngine.C
//
// A macro for running bulk trimming jobs, first you organize the files
//   you want to trim, then setup the cuts desired.
// Finally the engine can then be called to run over the full files set
//
// Created:  J Boehm -- March 2007
//
// $Author: boehm $
//
// $Revision: 1.1 $
//
// $Id: TrimEngine.C,v 1.1 2007/04/03 17:13:35 boehm Exp $
//
// $Log: TrimEngine.C,v $
// Revision 1.1  2007/04/03 17:13:35  boehm
// Brand new trimmer which runs outside of job control.  This allows it to run about 5x faster than the old trimmer.  I've had to employ some tricks to make everything work together reasonably - so if you want to run this yourself you may want to contact me first - Josh
//
//
//
/////////////////////////////////////////////////////////////////////

#include "NueAna/macros/TrimParams.h"
#include <vector>
#include <string>

void TrimEngine(int i)
{
  gROOT->LoadMacro("NueAna/macros/LoadLibs.C");
  LoadLibs();

  AddFirstYearNearData();
  AddNearMCSet("Daikon/NearL010185/AnaNue-n130110");
  AddNearMCSet("Daikon/NearL010185/AnaNue-n130111");
  AddNearMCSet("Daikon/NearL010185/AnaNue-n130112");
  AddNearMCSet("Daikon/NearL010185/AnaNue-n130113");
  AddNearMCSet("Daikon/NearL010185/AnaNue-n130114");
  AddNearMCSet("Daikon/NearL010185/AnaNue-n130115");
  AddNearMCSet("Daikon/NearL010185/AnaNue-n130116");

  AddFarMCSet(0);
  AddFarMCSet(3);
  AddFarMCSet(4); 
                                                                                
  Trimmer trim;

  Registry r2;
                                                                                
  r2.Set("HiPlaneTrackCut",25); // max number of planes in track
  r2.Set("HiEnergyCut",150); // max total energy in Meu
  r2.Set("HiTrackLikeCut",16); // max number of tracklike planes
  r2.Set("LoEnergyShowerCut",1.00);
  r2.Set("FiducialCut", 1);

  for(int j = 0; j < fileSets[i].size(); j++){
     trim.AddFiles(fileSets[i][j]);
  }

  trim.Config(r2);
  trim.RunTrimmer();
}

