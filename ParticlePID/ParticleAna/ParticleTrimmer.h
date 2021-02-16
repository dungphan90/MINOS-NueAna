#ifndef ParticleTrimmer_H
#define ParticleTrimmer_H

#include "NueAna/ParticlePID/ParticleFinder/POT.h"
#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "TChain.h"
#include "TFile.h"


#include <string>
using std::string;


class ParticleTrimmer
{
	public:
		ParticleTrimmer();
		~ParticleTrimmer();

		
		void AddFiles(string file);
		void RunTrimmer();

	private:
		TChain * chain_precord;
		TChain * chain_pot;
		
		PRecord * precord;
		POT * pot;
		
		POT * totalpot;
		
		void trimPOT();
		void trimPRecord();
		
		string outfile;
		
		TFile * out;
		
		int nearFileTypeCount;
		int farFileTypeCount[5];
		
};



#endif



