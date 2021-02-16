{ 
  gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs();
  gSystem->Load("libNueAnaExtrapolation.so");
  NueGui mainWin(gClient->GetRoot(), 150, 850);
}
