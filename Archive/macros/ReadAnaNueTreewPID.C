{

JobC j;

////////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 
gSystem->Load("libHistMan.so");

////////////////////////////////////////////////////////
//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
//j.Msg.SetLevel("NueReadwPID","Debug");
j.Msg.SetLevel("NueReadwPID","Info");
j.Msg.SetLevel("HistMan","Info");

////////////////////////////////////////////////////////
//set the analysis path
j.Path.Create("Main",
	      "NueReadwPID::Ana ");
////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=ana_nue,pid_nue"); 

////////////////////////////////////////////////////////
//Run
j.Path("Main").Run();
//j.Path("Main").RunNpass(50);
//j.Path("Main").RunNin(50);

////////////////////////////////////////////////////////
//report
j.Path("Main").Report();

////////////////////////////////////////////////////////
HistMan *hm = new HistMan("nueana");
hm->WriteOut("nueana_plots.root");


}
