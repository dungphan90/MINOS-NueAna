#ifndef SHWFIT_H
#define SHWFIT_H

#include "TObject.h"

class TH1F;
class TF1;
class TPad; 
class TPaveText;

class Shwfit : public TObject
{

public:
   Shwfit();
   Shwfit(const Shwfit &s);
   virtual ~Shwfit();

   virtual void Draw(Option_t *option);
   virtual void Print(Option_t *option) const;
   void Reset();
   void Clear(Option_t* option = "");
   void Draw(TPad *pad);

   //shwfit variables
   Float_t par_a;
   Float_t par_b;
   Float_t par_e0;
   Float_t chisq;
   Float_t shwmax;
   Int_t shwmaxplane;
   Int_t conv;
   Float_t chisq_ndf;
   Float_t e0_pe_ratio;
   Float_t caldet_comp;
   Float_t max_pe_plane;
   Float_t shwmaxplane_diff;
   Int_t hiPhStripCountM4;
   Int_t hiPhPlaneCountM4;
   Int_t hiPhStripCountM2;
   Int_t hiPhPlaneCountM2;
   Int_t hiPhStripCount;
   Int_t hiPhPlaneCount;
   Int_t hiPhStripCountP2;
   Int_t hiPhPlaneCountP2;
   Int_t hiPhStripCountP4;
   Int_t hiPhPlaneCountP4;

   //Dan's new variables
   Float_t Beta_Maxwell;
   Float_t Energy_Maxwell;
   Float_t chisq_Maxwell;
   Float_t ndf_Maxwell;
   Float_t Beta_Maxwell3;
   Float_t Energy_Maxwell3;
   Float_t chisq_Maxwell3;
   Float_t ndf_Maxwell3;
   Float_t trans_u_mean;
   Float_t trans_u_sigma;
   Float_t trans_u_chisq;
   Float_t trans_u_ndf;
   Float_t trans_v_mean;
   Float_t trans_v_sigma;
   Float_t trans_v_chisq;
   Float_t trans_v_ndf;
   Float_t E_ratio_half;
   Float_t n_ratio_half;
   Float_t E_ratio_2;
   Float_t n_ratio_2;
   Float_t pos_E_split;
   Float_t pos_n_split;


   // new contPlaneCount 
   Int_t contPlaneCount;
   Int_t contPlaneCount015;
   Int_t contPlaneCount030;
   Int_t contPlaneCount050;
   Int_t contPlaneCount075;
   Int_t contPlaneCount100;
   Int_t contPlaneCount200;

   //utrans parameters
   Float_t u_asym_peak;
   Float_t u_asym_vert;
   Float_t u_molrad_peak;
   Float_t u_molrad_vert;
   Float_t u_mean;
   Float_t u_rms;
   Float_t u_skew;
   Float_t u_kurt;


   //vtrans parameters
   Float_t v_asym_peak;
   Float_t v_asym_vert;
   Float_t v_molrad_peak;
   Float_t v_molrad_vert;
   Float_t v_mean;
   Float_t v_rms;
   Float_t v_skew;
   Float_t v_kurt;

   //uvtrans parameters
   Float_t uv_asym_peak;
   Float_t uv_asym_vert;
   Float_t uv_molrad_peak;
   Float_t uv_molrad_vert;
   Float_t uv_mean;
   Float_t uv_rms;
   Float_t uv_skew;
   Float_t uv_kurt;
   Float_t uv_ratio;

   Float_t vtxEnergy;
   Float_t energyPlane0;
   Float_t energyPlane1;
   Float_t energyPlane2;

   Float_t USlope;
   Float_t VSlope;
   Float_t UOffset;
   Float_t VOffset;
   Float_t UFitquality;
   Float_t VFitquality;
   Float_t UWDiff;
   Float_t VWDiff;
   Float_t USlopeMom;
   Float_t VSlopeMom;
   Float_t slopefix;

   Float_t UBeamLike;
   Float_t VBeamLike;
   Float_t UVSlope;

   Float_t ULongE;
   Float_t VLongE;
   Float_t LongE;

   Float_t complexity;
   Float_t wcomplexity;

   //9s cut 2pe
   //utrans parameters
   Float_t u_asym_peak_9s_2pe;
   Float_t u_asym_vert_9s_2pe;
   Float_t u_molrad_peak_9s_2pe;
   Float_t u_molrad_vert_9s_2pe;
   Float_t u_mean_9s_2pe;
   Float_t u_rms_9s_2pe;
   Float_t u_skew_9s_2pe;
   Float_t u_kurt_9s_2pe;

   //vtrans parameters
   Float_t v_asym_peak_9s_2pe;
   Float_t v_asym_vert_9s_2pe;
   Float_t v_molrad_peak_9s_2pe;
   Float_t v_molrad_vert_9s_2pe;
   Float_t v_mean_9s_2pe;
   Float_t v_rms_9s_2pe;
   Float_t v_skew_9s_2pe;
   Float_t v_kurt_9s_2pe;
                                                                                
   //uvtrans parameters
   Float_t uv_asym_peak_9s_2pe;
   Float_t uv_asym_vert_9s_2pe;
   Float_t uv_molrad_peak_9s_2pe;
   Float_t uv_molrad_vert_9s_2pe;
   Float_t uv_mean_9s_2pe;
   Float_t uv_rms_9s_2pe;
   Float_t uv_skew_9s_2pe;
   Float_t uv_kurt_9s_2pe;
   Float_t uv_ratio_9s_2pe;

   //9s cut 2pe dw
   //utrans parameters
   Float_t u_asym_peak_9s_2pe_dw;
   Float_t u_asym_vert_9s_2pe_dw;
   Float_t u_molrad_peak_9s_2pe_dw;
   Float_t u_molrad_vert_9s_2pe_dw;
   Float_t u_mean_9s_2pe_dw;
   Float_t u_rms_9s_2pe_dw;
   Float_t u_skew_9s_2pe_dw;
   Float_t u_kurt_9s_2pe_dw;
                                                                                
                                                                                
   //vtrans parameters
   Float_t v_asym_peak_9s_2pe_dw;
   Float_t v_asym_vert_9s_2pe_dw;
   Float_t v_molrad_peak_9s_2pe_dw;
   Float_t v_molrad_vert_9s_2pe_dw;
   Float_t v_mean_9s_2pe_dw;
   Float_t v_rms_9s_2pe_dw;
   Float_t v_skew_9s_2pe_dw;
   Float_t v_kurt_9s_2pe_dw;
                                                                                
   //uvtrans parameters
   Float_t uv_asym_peak_9s_2pe_dw;
   Float_t uv_asym_vert_9s_2pe_dw;
   Float_t uv_molrad_peak_9s_2pe_dw;
   Float_t uv_molrad_vert_9s_2pe_dw;
   Float_t uv_mean_9s_2pe_dw;
   Float_t uv_rms_9s_2pe_dw;
   Float_t uv_skew_9s_2pe_dw;
   Float_t uv_kurt_9s_2pe_dw;
   Float_t uv_ratio_9s_2pe_dw;

    TH1F *lenepl; //!  don't write this to the tree
    TH1F *ph_hist; //!
    TH1F *tenestu; //!  don't write this to the tree
    TH1F *tenestv; //!  don't write this to the tree

    TH1F *tenestu_9s_2pe; //!  don't write this to the tree
    TH1F *tenestv_9s_2pe; //!  don't write this to the tree
                                                                                
    TH1F *tenestu_9s_2pe_dw; //!  don't write this to the tree
    TH1F *tenestv_9s_2pe_dw; //!  don't write this to the tree
                                                                                


    TF1 *efit; //!  don't write this to the tree
    TF1 *hfit; //!  don't write this to the tree

    TF1 *efit_maxwell; //!  don't write this to the tree
    TF1 *efit_maxwell3; //!  don't write this to the tree
    TF1 *ufit; //!  don't write this to the tree
    TF1 *vfit; //!  don't write this to the tree

    TPaveText *info1; //! don't write this to the tree

private:

ClassDef(Shwfit,7)
};

#endif// SHWFIT_H
