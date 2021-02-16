//already have the compiled MiniPlotMaker.C+ loaded...
void near_mc_mrcc_hornON_ParticlePID(string path, int sys=0)
{

//	gSystem->CompileMacro("/data/minos/NueAnalysis/Code/MiniPlotMaker/MiniPlotMaker.C");

	MiniPlotMaker p;
	p.SetSystematic(sys);
	
	//add the necessary files
	// follow: p.AddInputFiles(file,group);






//	string path = "/minos/data/Dogwood1_Near/MRCC/out/minos/data/analysis/nue/2ndAnalysis/AnaNue_Files/BeforeLEM/Full/Untrimmed/Near/MC/L010185N/Standard/";

	p.AddInputFiles(path+"/An*.root");


	//set the output file
	p.SetOutputFile("out.root");
	
	
	//for testing.. set the maximum number to run over...
	//p.SetMaxEntries(100);
	
	//set the reco type  0:normal, 1:particlePID
	p.SetRecoType(1); 

	//set MC/data  0:data, 1:MC
	p.SetMC(1);

	//set MRCC  0:normal, 1:MRCC
	p.SetMRCC(1);

	//set detector 0:far, 1:near
	p.SetDetector(1);
	
	//set horn state 0:on, 1:off
	p.SetHorn(0);
			
	//normalize to what
	p.SetWantPOTs(1.e19);					
			
	//now make the plots
	p.MakePlots();

}
