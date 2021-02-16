{

  JobC j;
  gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 
  
  j.Msg.SetLevel("Io","Warning");
  j.Msg.SetLevel("NueHandScan","Debug");

  j.Path.Create("Main",
		"NueHandScan::Ana ");  

  j.Input.Set("Streams=NtpSt,pid_nue");


  //Scan File Options:
  //--------------------
  //set the random seed. If this is zero, it will use clock time
  //random seed used is output in the Summary file, so can back-track:
  j.Path("Main").Mod("NueHandScan").Set("RandomSeed=0");
  //adds a text string to the output file names:
  j.Path("Main").Mod("NueHandScan").Set("FileTag=FD_R1.18_Run20");
  //all events which pass are candidates for being added to the random file
  //this selects the fraction of the total that passed to be added.
  //basically controls the output list size.
  j.Path("Main").Mod("NueHandScan").Set("FracPassed=0.8");
  //pre-scale factors to scale down contributions from particular event types:
  //default is 1 i.e. no scaling
  j.Path("Main").Mod("NueHandScan").Set("PreScaleFactorNue=1");
  j.Path("Main").Mod("NueHandScan").Set("PreScaleFactorNuMu=1");
  j.Path("Main").Mod("NueHandScan").Set("PreScaleFactorNuTau=1");
  j.Path("Main").Mod("NueHandScan").Set("PreScaleFactorBNue=1");
  j.Path("Main").Mod("NueHandScan").Set("PreScaleFactorNC=1");
  
  //Add files:
  j.Input.AddFile("/data/MDC/far/sntp/R1.18/f24100020_0000.sntp.R1_18.root","NtpSt");
  j.Input.AddFile("/data/MDC/far/sntp/R1.18/f24110020_0000.sntp.R1_18.root","NtpSt");
  j.Input.AddFile("/data/MDC/far/sntp/R1.18/f24130020_0000.sntp.R1_18.root","NtpSt");

  //Run the job!
  j.Path("Main").Run();
  j.Path("Main").Report();  

}
