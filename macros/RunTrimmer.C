#include <vector>
#include <string>

void RunTrimmer()
{
  gROOT->LoadMacro("NueAna/macros/LoadLibs.C");
  LoadLibs();

  Trimmer trim;

  Registry r2;
                                                                                
  r2.Set("HiPlaneTrackCut",25); // max number of planes in track
  r2.Set("HiEnergyCut",150); // max total energy in Meu
  r2.Set("HiTrackLikeCut",16); // max number of tracklike planes
  r2.Set("LoEnergyShowerCut",1.00);
  r2.Set("FiducialCut", 1);

  string start = "/afs/fnal.gov/files/data/minos/d109/Dragon/";
  trim.AddFiles(start + "AnaNue-N00007786*.root");

  trim.Config(r2);
  trim.RunTrimmer();
}

