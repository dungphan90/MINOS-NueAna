// Macro for making AnaNue trees
//-----------------------------------------------------------------------------
// Can be run in several ways:
//----------------------------
// SNTP files (Birch/Carrot or Cedar/Daikon):
//      loon NueAna/macros/MakeAnaNueTree.C blah.sntp.root
// or
//      loon 'NueAna/macros/MakeAnaNueTree.C("blah.sntp.root")'
//
// MRNT files (Birch/Carrot single file era):
//      loon NueAna/macros/MakeAnaNueTree.C blah.mrnt.root
// or
//      loon 'NueAna/macros/MakeAnaNueTree.C("blah.mrnt.root")'
//
// SNTP + MRNT files (Cedar/Daikon split file era):
//      loon 'NueAna/macros/MakeAnaNueTree.C("blah.sntp.root","blah.mrnt.root")'
//
//------------------------------------------------------------------------------
#include "NueAna/NuePOT.h"

void RepairAnaNueMREData(string nue_file="",string old_file=""){

  ////////////////////////////////////////////////////////
  //load necessary libraries
  gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs();
  
  JobC j;

  ////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  //set the analysis path
  JobCPath& main_path = j.Path.Create("Main",
				      "SetKNNModule::Reco "
				      "FixModule::Reco "
				      "Output::Put");
  
  ////////////////////////////////////////////////////////
  //Set Input Streams
  j.Input.Set("Streams=NtpSt,NtpSR,NtpMC,NtpTH,NtpMR,NtpOld,snarltree,NueAna,ana_nue"); 

  j.Input.AddFile(nue_file.c_str(),"ana_nue");
  j.Input.DefineStream("NtpOld","NtpSt");
  j.Input.AddFile(old_file.c_str(),"NtpOld");


//  j.Msg.SetLevel("FixModule","Debug");

  ////////////////////////////////////////////////////////
  //configure NueModule

  Registry r(false);
  r.Set("kNNPath", "/afs/fnal.gov/files/home/room1/rustem/data/muon-knn/");
  r.Set("ForceMRCC", 1);
  r.Set("RunFast", 1);
  main_path.Mod("SetKNNModule").Config(r);
                                                                                                    
  //set the streams that will be written
  j.Path("Main").Mod("Output").Cmd("DefineStream ana_nue NueRecord");
  j.Path("Main").Mod("Output").Set("Streams=ana_nue");
  
  string filename;
  filename = nue_file;

  string base = filename.substr(filename.find_last_of("/")+1,filename.size()-5);
  
  string file = base + "_Fix.root";

  char setfilename[200];
  sprintf(setfilename,"FileName=%s", file.c_str());
  cout<<"Saving AnaNue Tree to "<<setfilename<<endl;
  j.Path("Main").Mod("Output").Set(setfilename);
  j.Path("Main").Mod("FixModule").Set(setfilename);
  char infilename[300];
  sprintf(infilename,"InFileName=%s", nue_file.c_str());
  cout<<infilename<<endl;
  j.Path("Main").Mod("FixModule").Set(infilename);
  
/*
  char pottreefn[1000];
  sprintf(pottreefn,"POTTreeFileName=%sAnaNue-%s-PECut-%f.root",path,file.c_str(),pecut);
  j.Path("Main").Mod("NueBeamMon").Set(pottreefn);
  */

  j.Input.Next();
  j.Path("Main").Run();
    
  ////////////////////////////////////////////////////////
  //report
  j.Path("Main").Report();
  
  ////////////////////////////////////////////////////////
}
