{
  gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs();
  gSystem->Load("libSensitivity.so");

  NueSensitivity  ns;
  ns.Run("InputMethod2.inp", "Output.root", 4.0);
}

