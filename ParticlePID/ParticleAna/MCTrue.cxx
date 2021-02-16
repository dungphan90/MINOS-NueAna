#include "NueAna/ParticlePID/ParticleAna/MCTrue.h"

MCTrue::MCTrue() : vtx_u(0), vtx_v(0), vtx_z(0),
					inu(0), iresonance(-1), iaction(-1), inunoosc(0),
					nuenergy(0), visenergy(0), 
					oscprob(1), totbeamweight(1),
					osc_L(0), osc_dm2(0), osc_sinth23(0), osc_sin2th13(0),
					type(-1), trainweight(0)
{}

MCTrue::~MCTrue()
{}

