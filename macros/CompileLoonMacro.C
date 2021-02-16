void CompileLoonMacro(TString LoadLibs, TString macro)
{
 gROOT->ProcessLine(".x "+LoadLibs); 
 gROOT->ProcessLine(".L "+macro+"++");
}
