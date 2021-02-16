#include <string>

void trim(string filename){

  gROOT->LoadMacro("$SRT_PRIVATE_CONTEXT/NueAna/ParticlePID/ParticleFinder/macros/LoadLibs.C"); LoadLibs();

	ParticleTrimmer a;
	
/*	
    JobCEnv& jce = JobCEnv::Instance();
    string filename = jce.GetFileName(0);
*/
	
	a.AddFiles(filename);
	
	a.RunTrimmer();

}
