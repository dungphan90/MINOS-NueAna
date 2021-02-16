void TestShort()
{
	NNTrain n;
	char tmp[200];
	
	sprintf(tmp,"/minos/app/scavan/NueAnaTrain/filelists/stest/");
	char path[200];

	sprintf(path,"%s/AnaNue-f210*",tmp);
	n.SetInputFile(0,path);
        sprintf(path,"%s/AnaNue-f213*",tmp);
	n.SetInputFile(1,path);
        sprintf(path,"%s/AnaNue-f214*",tmp);
	n.SetInputFile(2,path);


	n.SetTrainFile("traintestshort.root");
	n.MakeTestTree();

}
