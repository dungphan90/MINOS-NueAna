///////////////////////////////////////////////
////
////            SystematicGains
////
////    This module will adjust strip energy to simulate changes in gains
////
////	This must be run before reconstruction and only using a reconstruction
////       that uses strip energies in NtpSt to do the reco
////
////    Steven Cavanaugh, May 28, 2009
////
///////////////////////////////////////////////

#ifndef SystematicGains_H
#define SystematicGains_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <vector>

#include "StandardNtuple/NtpStRecord.h"

#include "TRandom2.h"



using namespace std;


class NtpStRecord;



class SystematicGains : public JobCModule
{

public:

  SystematicGains();
  ~SystematicGains();

public:

	// Analysis and Reconstruction methods
	JobCResult Reco(MomNavigator* mom);
	void EndJob();
	void BeginJob();
	void Config(const Registry&);
	const Registry& DefaultConfig() const;


private:

	static int printit;
	double allshift;
	double rndshift;
	TRandom2 rnd;
};
#endif // COMPAREALL_H
////////////////////////////////////////////////////////////////////////





