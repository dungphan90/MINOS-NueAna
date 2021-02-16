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
j.Msg.SetLevel("CompareAll","Info");
j.Msg.SetLevel("HistMan","Info");

//j.Input.AddFile("/home/boehm/work/Data/June3Files/nue/*.root", "ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/nue/AnaNue-f24110001_0000.root", "ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/numu/AnaNue-f24100001_0000.root", "ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/nutau/AnaNue-f24130001_0000.root", "ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/near/AnaNue-n14301001_0000.root","ana_nue");
//j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/near/AnaNue-n14301001_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006842*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006843*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006844*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006845*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006846*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006847*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006848*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006849*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006850*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006852*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006853*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006854*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue-N00006855*.root","ana_nue");
 j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/2005-03/AnaNue*.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301001_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301002_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301003_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301004_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301005_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301006_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301007_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301008_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301009_0000.root","ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue-n14301010_0000.root","ana_nue");
 j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/newnear/AnaNue*.root","ana_nue");
// j.Input.AddFile("nf/*.root","ana_nue");
////////////////////////////////////////////////////////
//set the analysis path
j.Path.Create("Main",
	      "CompareMD::Ana ");

j.Path("Main").Mod("CompareMD").Set("NMCFiles=366");//have to set it by hand

////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=ana_nue"); 

////////////////////////////////////////////////////////
//Run
j.Path("Main").Run();
//j.Path("Main").RunNpass(50);
//j.Path("Main").RunNin(5);

////////////////////////////////////////////////////////
//report
j.Path("Main").Report();

////////////////////////////////////////////////////////

//Get the run number from the first Input file name
//doesn't have to be this way, I just thought I'd be tricky
//char *path = getenv("NUE_OUT_DIR");
//if(path=='\0'){
//   path="./";
//}
//JobCEnv& jce = JobCEnv::Instance();
//string filename = jce.GetFileName(0);
//string rnumstring = filename.substr(filename.find_last_of("_")-9,9);
//string subrunstring = filename.substr(filename.find_last_of("_")+1,4);
//char setfilename[5000];
//sprintf(setfilename,"%sallcomp-%s_%s.root",
//	path,rnumstring.c_str(),subrunstring.c_str());
HistMan *hm = new HistMan("allcomp");
hm->WriteOut("test.root");

}
