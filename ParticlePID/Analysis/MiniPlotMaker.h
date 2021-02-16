#ifndef MiniPlotMaker_h
#define MiniPlotMaker_h

#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "DataUtil/MCInfo.h"
#include "TObjArray.h"
#include "TChainElement.h"

#include "TH1D.h"

#include "TH2D.h"

#include <string>
#include <vector>

#include "NueAna/ParticlePID/Analysis/DirectoryHelpers.h"
#include "NueAna/ParticlePID/NueMiniPID.h"


#include "MCReweight/SKZPWeightCalculator.h"
#include "MCReweight/WeightCalculator.h"
#include "Registry/Registry.h"
#include "MCReweight/MCReweight.h"
#include "MCReweight/MCEventInfo.h"
#include "MCReweight/NeugenWeightCalculator.h"
#include "MCReweight/ReweightHelpers.h"



using namespace std;

#define MAXCUTS 40  //must be equal to number of cuts to use x 4 -- 4 for the mrcc qp cuts...

class MiniPlotMaker
{

	public:
		MiniPlotMaker();
		~MiniPlotMaker();
		
		void PrintConfigs();
		int AddInputFiles(string s,int group=0);

		void SetOutputFile(string s){outputFile=s;};
		
		void SetRecoType(int t){recoType=t;};
		void SetMC(int t){isMC=t;};
		void SetMRCC(int t){isMRCC=t;};
		void SetDetector(int t){detector=t;};
		void SetHorn(int t){horn=t;};
		void SetSystematic(int t);

		void SetMaxEntries(int i){maxEntries=i;};

		void MakePlots();
		
		void SetWantPOTs(double wp=1e19){wantPOTs=wp;};
		
		//should we overwrite the output file if it exists?
		void SetOverwrite(int i){overwrite=i;};

		//adjust the overflow bins for lisa histos
		void MoveOverFlow(TH2D *h);

	private:
		string outputFile;
		TFile * outfile;
		
		
		vector<string> GetFileList(TChain *c);
		string GetDirectory();

		TChain * infiles[3];
		TChain * infilesPOTTree[3];
		TChain * infilesForPOTCounting[3];
		int infileCount[3];
		int infilePOTTreeCount[3];
		int infileForPOTCountingCount[3];
		
		TDirectory * MakeDirectory(TFile *f);
		
		void CountPots(TChain *c, TChain *d, double rp[3]);
		
		void MakeHistos();
		void FillHistos(int cut, double weight, int rp);
		
		void SaveHistos(TDirectory *td);
		
		/////////////
		//job settings and helpers
		int recoType;
		int isMC;
		int isMRCC;
		int detector;
		int horn;
		int systematic;
		
		string GetRecoTypeString();
		string GetMCString();
		string GetMRCCString();
		string GetDetectorString();
		string GetHornString();
		string GetSystematic();
		
		string GetMRCCQPCut(int cut);
		/////////////
		
		int maxEntries;
		
		
		
		//////
		
		NueMiniPID *nm;
		
		void ProcessChainNormal(TChain * c, double chainweight[]);
		void ProcessChainNormalMRCC(TChain * c, double chainweight[]);
		void ProcessChainParticlePID(TChain * c, double chainweight[]);
		void ProcessChain2Pairs(TChain * c, double chainweight[]);

		//apply the systematic in var systematic to NueMiniPID nm
		void ApplySystematic();

		int GetRunPeriod(time_t ts);
		string GetRunName(int j);

		int GetNCuts();
		string GetCutName(int cut);
		void SetUpChain(TChain * c);
		
		double wantPOTs;
		string GetNueClassSuffix(int j);

		int overwrite;

		//histos
		// [run period] [cuts] [mctype]
		// run period 0-2 (run 1-3) 3 for all
		// mctype 0-4 (type) 5 all
		TH1D * hist_recoE[4][MAXCUTS][6];
		TH1D * hist_resCode[4][MAXCUTS][6];

		TH1D * hist_nuEnergy[4][MAXCUTS][6];
		TH1D * hist_ann14[4][MAXCUTS][6];
		TH1D * hist_ann11[4][MAXCUTS][6];
		TH1D * hist_ann11_firebird[4][MAXCUTS][6];

		TH1D * hist_pidA[4][MAXCUTS][6];
		TH1D * hist_pidB[4][MAXCUTS][6];
		TH1D * hist_pidC[4][MAXCUTS][6];
		TH1D * hist_pidD[4][MAXCUTS][6];
		TH1D * hist_pidE[4][MAXCUTS][6];
		TH1D * hist_pidF[4][MAXCUTS][6];

		TH2D * hist_pidA_recoE[4][MAXCUTS][6];
		TH2D * hist_pidB_recoE[4][MAXCUTS][6];
		TH2D * hist_pidC_recoE[4][MAXCUTS][6];
		TH2D * hist_pidD_recoE[4][MAXCUTS][6];
		TH2D * hist_pidE_recoE[4][MAXCUTS][6];
		TH2D * hist_pidF_recoE[4][MAXCUTS][6];
			

		//keep lisa histos seperate... because of the way to deal with the overflow bins
		TH2D * hist_lisa_trueE_recoE[4][MAXCUTS][6];
		TH2D * hist_lisa_recoE_pidF[4][MAXCUTS][6];
		TH2D * hist_lisa_recoE_ann11[4][MAXCUTS][6];
		TH2D * hist_lisa_recoE_ann11_firebird[4][MAXCUTS][6];
		TH2D * hist_lisa_recoE_ann14[4][MAXCUTS][6];



		TH1D * hist_ntot[4][MAXCUTS][6];
		TH1D * hist_pars[14][4][MAXCUTS][6];
		
		TH2D * hist_pids_PID_ann[4][6];
		
		
		string GetDirectoryString();
		
		static double potratios[3]; 
		static double potratios_mrcc[3];

		//for reweighting
		SKZPWeightCalculator *skzp;
	
		MCReweight *mcr;
		NeugenWeightCalculator *nwc;
		Registry *rwtconfig;	
		int lastReweightType;
		double GetReweight(int type);

		bool PassesNuePreselection();


		double inputPOT[3][3];//a record of the input files for the tree...

};

#endif

