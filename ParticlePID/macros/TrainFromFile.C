void TrainFromFile(string file, string format="14:7", int sp=20, int typ=6)
{


	NNTrain n;
	n.SetTrainFile(file,1);

	n.Train(sp,5,typ,format);

}
