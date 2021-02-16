void MakeTrainFile(string inpath, string outfile)
{
	NNTrain n;

        char tmp[200];

        sprintf(tmp,inpath.c_str());//"/minos/app/scavan/NueAnaTrain/filelists/train/");

        char path[200];

        sprintf(path,"%s/AnaNue-f210*",tmp);
        n.SetInputFile(0,path);
        sprintf(path,"%s/AnaNue-f213*",tmp);
        n.SetInputFile(1,path);
        sprintf(path,"%s/AnaNue-f214*",tmp);
        n.SetInputFile(2,path);



	n.SetTrainFile(outfile);//"nout.root");
	n.MakeTrainTree(1);


}
