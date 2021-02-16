#include <vector>
#include <string>

void TrimMRCC(TString option,string sfilelist,string path,bool separatebyRunPeriod)
{
  //option is 'Fid' or 'Pre'
  //sfilelist if a text file with the names of the files being merged into one trimmed file (for preselected trimmed files, which are made from fiducial trimmed files, there will only be one line in this file).  this text file will be created in the shell script which is submitted to condor.  there will be a separate shell script for each kind of file set (e.g. nd mc intensity 124, etc)
  //path is the directory where the output file will go
  
  gROOT->LoadMacro("$SRT_PRIVATE_CONTEXT/NueAna/macros/LoadLibs.C");
  LoadLibs();
  
  Trimmer trim;
  
  if(option=="Fid")//standard fiducial/data quality, mrcc fiducial and mrcc preselection
  {
    trim.SetCuts("MRCC",2);
  }
  else if(option=="Pre")
  {
    trim.SetCuts("MRCC",6);
  }
  else
  {
    cout<<"Error: choose 'Fid' or 'Pre'"<<endl;
    return;
  }
  
  string outfile;
  ifstream filelist(gSystem->ExpandPathName(sfilelist.c_str()));
  if(!filelist.good())
  {
    cout<<"problem reading filelist"<<endl;
    return;
  }
  
  string line;
  int nfiles=0;
  while(!filelist.eof())
  {
    getline(filelist,line);
    if(!line.empty())
    {
      if(nfiles==0)
      {
        outfile = line.substr(line.find_last_of("/")+1,line.find_last_of(".root")-line.find_last_of("/") - 5);
        if(option=="Fid")
        {
          outfile += "-Fid.root";
        }
        if(option=="Pre")
        {
          string tmp = line.substr(line.find_last_of(".root")-16+1,16);
          if(tmp=="_RunPeriod1.root" || tmp=="_RunPeriod2.root" || tmp=="_RunPeriod3.root")//need to keep this RunPeriod label in name
          {
            outfile = line.substr(line.find_last_of("/")+1,line.find("-Fid")-1-line.find_last_of("/"));
            outfile += "-Pre";
            outfile += tmp;
          }
          else//no RunPeriod label
          {
            outfile = line.substr(line.find_last_of("/")+1,line.find_last_of("-Fid.root")-line.find_last_of("/") - 9);
            outfile += "-Pre.root";
          }
          
        }
        path += "/";
        outfile.insert(0,path);
      }
      trim.AddFiles(line);
      nfiles++;
    }
  }
  
  filelist.close();
  
  trim.SetOutputFile(gSystem->ExpandPathName(outfile.c_str()));
  
  if(separatebyRunPeriod)
  {
    trim.SeparatebyRunPeriod();
  }
  
  trim.RunTrimmer();
  
  return;
}

