//already have the compiled MiniPlotMaker.C+ loaded...
void near_data_mrcc_hornON_normal(string path, int sys=0)
{

//	gSystem->CompileMacro("/data/minos/NueAnalysis/Code/MiniPlotMaker/MiniPlotMaker.C");

	MiniPlotMaker p;
	p.SetSystematic(sys);
	
	//add the necessary files
	// follow: p.AddInputFiles(file,group);


//	p.AddInputFiles("/minos/data/mini_08_26/standard/Griffin_Production_v00/near/data/MRCC/official/*");


//	string path="/minos/data/Dogwood1_Sorted/Standard/out/minos/data/analysis/nue/2ndAnalysis/AnaNue_Files/BeforeLEM/Full/Untrimmed/Near/Data/L010185N/MRCC/AllRuns";

	//string path="/minos/data/Dogwood1_Near/Standard/out/minos/data/analysis/nue/2ndAnalysis/AnaNue_Files/BeforeLEM/Full/Untrimmed/Near/Data/L010185N/MRCC";

	p.AddInputFiles(path+"/Run1/An*.root");
	p.AddInputFiles(path+"/Run2/An*.root");
	p.AddInputFiles(path+"/Run3/An*.root");


	//set the output file
	p.SetOutputFile("out.root");
	
	
	//for testing.. set the maximum number to run over...
	//p.SetMaxEntries(100);
	
	//set the reco type  0:normal, 1:particlePID
	p.SetRecoType(0); 

	//set MC/data  0:data, 1:MC
	p.SetMC(0);

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
