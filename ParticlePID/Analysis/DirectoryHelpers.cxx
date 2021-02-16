
#include "NueAna/ParticlePID/Analysis/DirectoryHelpers.h"

using namespace std;


void DirectoryHelpers::Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

	//printf("tokens: ");

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        string t = str.substr(lastPos, pos - lastPos);
        if(t!="")tokens.push_back(t);
        
      //  printf("-%s- ",t.c_str());
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
   // printf("\n");
}



TDirectory * DirectoryHelpers::GetDirectory(TDirectory *f, string dir, int create)
{
	if(!create)
	{ 
		TDirectory * r = f->GetDirectory(dir.c_str());
		if(r)return r;
		
		//maybe the string has repeating tokens?
		vector<string> paths;
		Tokenize(dir,paths,"/");
		
		r = f->GetDirectory("/");
		for(int i=0;i<(int)paths.size();i++)
		{
			r=r->GetDirectory(paths[i].c_str());
			if(!r)break;
		}
		return r;
		
	}
	
	//parse the dir
	
	vector<string> paths;
	Tokenize(dir,paths,"/");
	
	
	TDirectory *tmp=f;
	TDirectory *tmpNew=0;
	for(int i=0;i<(int)paths.size();i++)
	{
		tmpNew=tmp->GetDirectory(paths[i].c_str());
		if(!tmpNew)
		tmp->mkdir(paths[i].c_str());
		tmpNew=tmp->GetDirectory(paths[i].c_str());
		tmp=tmpNew;
	}	
	
	
	//clear directory?
	if(create==2)
	{
		tmp->Close();
		tmp->Write();
	}
	
	return tmp;	
}


TDirectory * DirectoryHelpers::GetDirectory(TFile *f, string dir, int create)
{
	TDirectory *ret = GetDirectory(f->GetDirectory("/"),dir,create);

	if(!ret)printf("failed to get directory %s\n",dir.c_str());
	
	return ret;
}




int DirectoryHelpers::MakeTrueDirectory(vector<string> paths, int pos)
{
	
	if(paths.size()>0)
	{
		string str ="/";
		if(pos > (int)paths.size())return 1;

		for(int i=0;i<pos;i++)str+=paths[i]+"/";
		
		int res = mkdir(str.c_str(),0755);

		return MakeTrueDirectory(paths,pos+1) || res;  //0 if OK
	
	}

	return 10;
}


int DirectoryHelpers::MakeTrueDirectory(string dir)
{
	printf("Making directory %s\n",dir.c_str());

	vector<string> paths;
	Tokenize(dir,paths,"/");
	
	return MakeTrueDirectory(paths,1);
}



