#ifndef TH1DStat_h
#define TH1DStat_h


#include "TH1D.h"
#include "TArrayD.h"


class TH1DStat 
{
	public:
		TH1DStat(TH1D * t);
		~TH1DStat(){};
	
		void Fill(double v, double w, double s);
		TH1D * GetTH1D(int withstats=1);
		void Dump();
		
	private:
		TH1D * th1d;
		
		TArrayD td;
		TArrayD pw;

		bool wasErrorAdjusted;


};

#endif


