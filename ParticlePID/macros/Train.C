void Train(string format="14:7",int steps=20, int traintype=6)
{
	NNTrain n;


        char tmp[200];

        sprintf(tmp,"/minos/app/scavan/NueAnaTrain/filelists/train/");

        char path[200];
        
        sprintf(path,"%s/AnaNue-f210*",tmp);
        n.SetInputFile(0,path);
        sprintf(path,"%s/AnaNue-f213*",tmp);
        n.SetInputFile(1,path);
        sprintf(path,"%s/AnaNue-f214*",tmp);
        n.SetInputFile(2,path);
        
        



	n.SetTrainFile("nout.root");
	n.Train(steps,5,traintype,format);

}
