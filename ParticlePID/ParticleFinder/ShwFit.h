#ifndef SHWFIT_H
#define SHWFIT_H

#include "TH1F.h"
#include "TF1.h"

#include <vector>

#include <map>

class ShwFit{
	public:
		ShwFit();
		~ShwFit();


		double par_a;
		double par_b;
		double par_e0;
	//	double par_zoffset;
		
		double par_a_err;
		double par_b_err;
		double par_e0_err;
	//	double par_zoffset_err;		
		
		double chisq;
		double ndf;
		
		double prob;
		
		double shwmax;
		double shwmaxplane;
		double conv;
	
		double pp_chisq;
		int pp_ndf;
		int pp_igood;
		double pp_p;

		double cmp_chisq;
		int cmp_ndf;
		double peakdiff;	
	
		void Reset();
		void Insert(double ph, double z);

		
		void Fit3d(int psf=0);


		
		//must insert in order from front to back!
		//must call reset before inserting!
		void Insert3d(double ph, double u, double v, double z);
		void SetDetector(int mydet);


		double pred_e_a;
		double pred_g_a;
		double pred_b;
		double pred_e0;
		double pred_e_chisq;
		double pred_e_ndf;
		double pred_g_chisq;
		double pred_g_ndf;	
		
		double zshift;
		
		double pre_over;
		double pre_under;
		double post_over;
		double post_under;	
		
//	private:
		static TH1F * lenepl;
		static TF1 * efit;
		Double_t GetMaximumX(TF1* efit, Double_t xmin=0, Double_t xmax=0);
		
		
		static TH1F * lenepl3d;

		
		std::vector<double> v_ph;
		std::vector<double> v_u;
		std::vector<double> v_v;
		std::vector<double> v_z;
		
		int detector;
		
		static std::map<double,int> * planemap;

		void BuildPlaneMap();	
		int GetPlaneFromZ(double z);
		
		int debug;



};

#endif

