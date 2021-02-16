{

JobC j;

////////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 
gSystem->Load("libHistMan.so");

////////////////////////////////////////////////////////
//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
//j.Msg.SetLevel("NueReadTJPID","Debug");
//j.Msg.SetLevel("HistMan","Debug");

////////////////////////////////////////////////////////
//set the analysis path
j.Path.Create("Main",
	      "NueReadTJPID::Ana ");
////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=NtpSR,NtpMC,NtpTH,pid_nue"); 

////////////////////////////////////////////////////////
//Run
j.Path("Main").Run();
//j.Path("Main").RunNpass(50);
//j.Path("Main").RunNin(50);

////////////////////////////////////////////////////////
//report
j.Path("Main").Report();

////////////////////////////////////////////////////////
HistMan *hm = new HistMan("tjpid");
hm->WriteOut("TJpid_plots.root");


}
