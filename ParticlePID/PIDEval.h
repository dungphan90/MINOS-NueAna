#ifndef PIDEVAL_H
#define PIDEVAL_H

#include "JobControl/JobCModule.h"
#include "NueAna/ParticlePID/ParticleAna/PRecord.h"

#include "NueAna/ParticlePID/ParticleAna/Particles.h"
#include "NueAna/ParticlePID/ParticleAna/Event.h"
#include "TMultiLayerPerceptron.h"
//#include "NueAna/ParticlePID/ParticleAna/OscCalc.h"
#include "NueAna/NueRecord.h"

#include <string>

using std::string;

class PIDEval : public JobCModule
{
	public:
	
		PIDEval();
		~PIDEval();
		static void LoadMLP();

		
		void ana(Particles * p, Event * e,NueRecord *nr);

	
  		// Handle job status changes
  		void BeginJob();
  		void EndJob();

  		// Analysis and Reconstruction methods
  		JobCResult Reco(MomNavigator* mom);

  		// Module configuration
  		const Registry& DefaultConfig() const;
  		void Config(const Registry& r);

  		// User interface methods
  		void Reset();



	private:
		static int first;
		static TMultiLayerPerceptron * mlp1;
		static TMultiLayerPerceptron * mlp2;
		static TMultiLayerPerceptron * mlp3;
		static TMultiLayerPerceptron * mlp4;
		static TMultiLayerPerceptron * mlp5;
		static TMultiLayerPerceptron * mlp6;
		string infile;
		string outfile;
		int cnt;
		int total;


	/*	OscCalc fOscCalc;
		float fBaseLine;
		double fDeltaMS12;
		double fTh12;
		double fTh23;
		double fDensity;
		void SetOscParamBase( float dm2, float ss13, float delta, int hierarchy);
		void recalc_oscprob(NueRecord *p);
*/

//		void adjRecoE(PRecord * p);

        string kPOTTreeName;
        string kPOTTreeNameIn;


};

#endif

