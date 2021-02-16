////////////////////////////////////////////////////////////////////////
//
//  use this utility to grab NtpSt events from a list 
// 
// with Gather(string filelist, string runlist, string outfile, int doMRCC=0)
//  filelist is a file containing a list of /pnfs files to look in
// runlist is a list with each line containing "run snarl event"  so "40182 234 0"
// outfile is a root file to save the events to
// doMRCC=1 if you are gathering from mrcc files and want the NtpMR tree in addition to NtpSt
//
////////////////////////////////////////////////////////////////////////


#ifndef NTPSTEVENTGRABBER_H
#define NTPSTEVENTGRABBER_H

#include <vector>
#include "TChain.h"
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "StandardNtuple/NtpStRecord.h"
#include "MuonRemoval/NtpMRRecord.h"


using namespace std;



class NtpStEventGrabber 
{
public:
  NtpStEventGrabber();
  ~NtpStEventGrabber();

	void Gather(string filelist, string runlist, string outfile, int doMRCC=0);

	TFile *outfile;
	TTree *outtree;
	TTree *outtreeMRCC;

	int counter;


	NtpStRecord * ntpst;
	NtpMRRecord * ntpmr;

	vector<string> filelist;
	vector<int>filerunlist;

	vector<int>runlist;
	vector<int>snarllist;
	vector<int>eventlist;

        void LoadFileList(string filelist);
        void LoadRunList(string runlist);
	void DoGather();

	void SetFastBranchStatus(TChain *c);

	int doMRCC;

};
#endif // COMPAREALL_H
////////////////////////////////////////////////////////////////////////
