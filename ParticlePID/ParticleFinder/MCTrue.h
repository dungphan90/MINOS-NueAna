#ifndef MCTRUE_H
#define MCTRUE_H


#include "TObject.h"

#include "MCNtuple/NtpMCFluxInfo.h"

#include "MCNtuple/NtpMCStdHep.h"
#include <vector>

// ROOT forward declaration
class TObjArray;


class MCTrue : public TObject

{

	public:
	
		MCTrue();
		virtual ~MCTrue();
 		void Clear(Option_t* option = "");


		double vtx_u;
		double vtx_v;
		double vtx_z;
		double nuenergy;
		double oscprob;
		double totbeamweight;
		
		double osc_L;
		double osc_dm2;
		double osc_sinth23;
		double osc_sin2th13;
		double visible_energy;


		
		int type; // 0 nc 1 nue 2 cc 3 tau 4 bnue
		int inu;
		int iresonance;
		int iaction;
		int inunoosc;

  		std::vector<NtpMCStdHep>stdhep;   //-> array of NtpMCStdHep 
				
		NtpMCFluxInfo flux;



  
	private:
		void Init();


        ClassDef(MCTrue,1)


};

#endif

