#include <string.h>
#include <iostream>
#include <vector>
#include <algorithm>
                                                                                
	
void CreateChain(char* filename, char *dirname, char *chain){
  
  //* CreateChain.C *//
  // Reads in all the ana_tree ntuples inside a user specified 
  // directory and adds the trees within to their corresponding chains.
  // writes chain to file named as input using all files in dirname
  //   usage:
  //          root> .L CreateChain.C
  //          root> CreateChain("filename","dirname")
  using namespace std;
                                                                                  
 vector <string> files;
    
	
  TChain *chain_ana = new TChain(chain);

  Int_t nfiles = 0;
  Char_t *file;
  Char_t path_file[100];

  TString me;
  
  //**Change this to reflect the ntuple's path**//
  //  Char_t dirname[]="./";
  
  void *dirhandle = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname));

  //Loop over the filenames inside the directory
  while (file = gSystem->GetDirEntry(dirhandle)) {
    if (strstr(file,".root") != 0) {
      strcpy(path_file,dirname);
      strcat(path_file,file);
      files.push_back(path_file);
      
      //Add trees to their correpsonding chains
//      chain_ana->Add(path_file,-1);  //Ana_Tree chain

      printf("file %d: %s\n",nfiles,path_file);
      nfiles++;
    }
  }

    //sort the string vector
  sort(files.begin(), files.end());
  
  for(int i = 0; i < files.size(); i++){
     cout<<files[i]<<endl;	  
     chain_ana->Add(files[i].c_str(),-1);          
  }
    if (nfiles == 0) {
    printf("No *.root files found in specified directory\n");
    break;
  }
    TFile *Target = TFile::Open(filename,"RECREATE");
    chain_ana->Write();

}

