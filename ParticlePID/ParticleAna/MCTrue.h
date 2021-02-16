#ifndef MCTRUE_H
#define MCTRUE_H

class MCTrue
{

	public:
	
		MCTrue();
		~MCTrue();

		double vtx_u;
		double vtx_v;
		double vtx_z;
			
		int inu;
		int iresonance;
		int iaction;
		int inunoosc;
		
		double nuenergy;
		double visenergy;
		
		double oscprob;
	

		
		double totbeamweight;
		
		double osc_L;
		double osc_dm2;
		double osc_sinth23;
		double osc_sin2th13;
		
		int type; // 0 nc 1 nue 2 cc 3 tau 4 bnue

		double trainweight; //oscprob * totbeamweight


};

#endif

