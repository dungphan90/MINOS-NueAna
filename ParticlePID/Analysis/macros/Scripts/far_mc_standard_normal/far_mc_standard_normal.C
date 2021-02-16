//already have the compiled MiniPlotMaker.C+ loaded...
void far_mc_standard_normal(string path, int sys=0)
{

	//gSystem->CompileMacro("/data/minos/NueAnalysis/Code/MiniPlotMaker/MiniPlotMaker.C");

	MiniPlotMaker p;
	p.SetSystematic(sys);
	
	//add the necessary files
	// follow: p.AddInputFiles(file,group);


	p.AddInputFiles((path+"/AnaNue-f210*").c_str(),0);
	p.AddInputFiles((path+"/AnaNue-f213*").c_str(),1);
	p.AddInputFiles((path+"/AnaNue-f214*").c_str(),2);
	
	
	//set the output file
	p.SetOutputFile("out.root");
	
	
	//for testing.. set the maximum number to run over...
	//p.SetMaxEntries(100);
	
	//set the reco type  0:normal, 1:particlePID
	p.SetRecoType(0); 

	//set MC/data  0:data, 1:MC
	p.SetMC(1);

	//set MRCC  0:normal, 1:MRCC
	p.SetMRCC(0);

	//set detector 0:far, 1:near
	p.SetDetector(0);

			
	//normalize to what
	p.SetWantPOTs(7.e20);					
			
	//now make the plots
	p.MakePlots();

}
