void RunLoonMacro(TString macro)
{
 gROOT->ProcessLine(".x NueAna/macros/LoadLibs.C"); 
 gROOT->ProcessLine(".x "+macro);
}
