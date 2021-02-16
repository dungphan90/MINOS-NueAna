{

JobC j;

////////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 
gSystem->Load("libHistMan.so");

////////////////////////////////////////////////////////
//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
//j.Msg.SetLevel("NueSensitivity","Debug");
//j.Msg.SetLevel("HistMan","Debug");

////////////////////////////////////////////////////////
//set the analysis path
j.Path.Create("Main",
	      "NueSensitivity::Ana ");
////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=NtpSR,NtpMC,NtpTH,pid_nue"); 

 char name[256];

 j.Input.AddFile("/data/MDC/farmock/sntp/R1.12/F21100001_0000.sntp.R1.12.root","NtpSR");
 j.Input.AddFile("/home/cbs/NueMDCStudy/PIDcutsmock/PIDNue-F21100001.root","pid_nue");
 j.Path("Main").Mod("NueSensitivity").Set("nNuMuFiles=5");
 j.Path("Main").Mod("NueSensitivity").Set("nNueFiles=5");
 j.Path("Main").Mod("NueSensitivity").Set("nNuTauFiles=5");

 /*
 //far:
 ifstream inNtpFar("fileLists/farNtpFiles_cuts.dat");
 while(inNtpFar>>name) j.Input.AddFile(name,"NtpSR,NtpMC,NtpTH");
 ifstream inPidFar("fileLists/farPidFiles_cuts.dat");
 while(inPidFar>>name) j.Input.AddFile(name,"pid_nue"); 
 */

 /*
 //near:
 ifstream inNtpNear("fileLists/nearNtpFiles_jbms.dat");
 while(inNtpNear>>name) j.Input.AddFile(name,"NtpSR,NtpMC,NtpTH");
 ifstream inPidNear("fileLists/nearPidFiles_jbms.dat");
 while(inPidNear>>name) j.Input.AddFile(name,"pid_nue");
 j.Path("Main").Mod("NueSensitivity").Set("nNearFiles=200");
 */

 /* 
 //near challenge:
 ifstream inNtpNear("fileLists/nearChallengeNtpFiles.dat");
 while(inNtpNear>>name) j.Input.AddFile(name,"NtpSR");
 ifstream inPidNear("fileLists/nearChallengePidFiles.dat");
 while(inPidNear>>name) j.Input.AddFile(name,"pid_nue");
 j.Path("Main").Mod("NueSensitivity").Set("nChallengeNearFiles=200");
 */

////////////////////////////////////////////////////////
//Run
j.Path("Main").Run();
//j.Path("Main").RunNpass(50);
//j.Path("Main").RunNin(50);

////////////////////////////////////////////////////////
//report
j.Path("Main").Report();

////////////////////////////////////////////////////////

}
