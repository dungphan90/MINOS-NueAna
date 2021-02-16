#ifndef DirectoryHelpers_h
#define DirectoryHelpers_h

#include "TDirectory.h"
#include "TFile.h"
#include <string>
#include <vector>

#include <sys/stat.h>
#include <sys/types.h>

using namespace std;


namespace DirectoryHelpers
{


void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters);


TDirectory * GetDirectory(TDirectory *f, string dir, int create=0);


TDirectory * GetDirectory(TFile *f, string dir, int create=0);



int MakeTrueDirectory(vector<string> paths, int pos=1);


int MakeTrueDirectory(string dir);

}

#endif

