#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

class Particle
{
	public:
		Particle(){};
		~Particle(){};

		Int_t type; //charged pion, neutron, proton, em-like

		Double_t stripEnergy; //sum of strip energy

		Double_t fittedEnergy; //fitted energy based on particle type

		Int_t begPlane;
		Int_t endPlane;

		Int_t begStrip;
		Int_t endStrip;


		Double_t endz;
		Double_t endt;
		Double_t begz;
		Double_t begt;
		Double_t begu;
		Double_t begv;

		std::vector<int> strip;
		std::vector<int> plane;
		std::vector<double> energy;
		std::vector<double> z;
		std::vector<double> t;

};

#endif

