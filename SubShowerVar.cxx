#include "NueAna/SubShowerVar.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(SubShowerVar)

SubShowerVar::SubShowerVar():
  ncluster(ANtpDefVal::kInt), //# of clusters in primary shower
  nclusterU(ANtpDefVal::kInt), //# of clusters in primary shower in U
  nclusterV(ANtpDefVal::kInt), //# of clusters in primary shower in V
  nPhysClusterU(ANtpDefVal::kInt), //# of "physics" clusters in U
  nPhysClusterV(ANtpDefVal::kInt), //# of "physics" clusters in V
  nstp0U(ANtpDefVal::kInt), //# of strips in primary SubShower in U
  nstp0V(ANtpDefVal::kInt), //# of strips in primary SubShower in V
  E2to1U(ANtpDefVal::kDouble), //PH ratio of 2ndary SubShower over primary SubShower in U
  E2to1V(ANtpDefVal::kDouble), //PH ratio of 2ndary SubShower over primary SubShower in V
  PHAvgIDU(ANtpDefVal::kDouble), //PH averaged clu ID in U
  PHAvgIDV(ANtpDefVal::kDouble), //in V
  PHAvgProbEMU(ANtpDefVal::kDouble), //PH averaged clu probem in U
  PHAvgProbEMV(ANtpDefVal::kDouble), //in V
  PHFracRMSU(ANtpDefVal::kDouble), //RMS of PH fraction of clusters in U
  PHFracRMSV(ANtpDefVal::kDouble), //in V
  PHAvgDevU(ANtpDefVal::kDouble), //PH weighted average of SubShower spread in U
  PHAvgDevV(ANtpDefVal::kDouble), //in V
  pid(ANtpDefVal::kDouble), //pid from NN trained on SubShower variables
  nclusterU_he(ANtpDefVal::kInt),
  nclusterV_he(ANtpDefVal::kInt),
  nEMU(ANtpDefVal::kInt),
  nEMV(ANtpDefVal::kInt),
  nHadU(ANtpDefVal::kInt),
  nHadV(ANtpDefVal::kInt),
  nTrkU(ANtpDefVal::kInt),
  nTrkV(ANtpDefVal::kInt),
  nXTkU(ANtpDefVal::kInt),
  nXTkV(ANtpDefVal::kInt),
  nHalU(ANtpDefVal::kInt),
  nHalV(ANtpDefVal::kInt),
  nEMU_he(ANtpDefVal::kInt),
  nEMV_he(ANtpDefVal::kInt),
  nHadU_he(ANtpDefVal::kInt),
  nHadV_he(ANtpDefVal::kInt),
  nTrkU_he(ANtpDefVal::kInt),
  nTrkV_he(ANtpDefVal::kInt),
  nXTkU_he(ANtpDefVal::kInt),
  nXTkV_he(ANtpDefVal::kInt),
  nHalU_he(ANtpDefVal::kInt),
  nHalV_he(ANtpDefVal::kInt),
  eEMU(ANtpDefVal::kDouble),
  eEMV(ANtpDefVal::kDouble),
  eHadU(ANtpDefVal::kDouble),
  eHadV(ANtpDefVal::kDouble),
  eTrkU(ANtpDefVal::kDouble),
  eTrkV(ANtpDefVal::kDouble),
  eXTkU(ANtpDefVal::kDouble),
  eXTkV(ANtpDefVal::kDouble),
  eHalU(ANtpDefVal::kDouble),
  eHalV(ANtpDefVal::kDouble),  
  wEMU(ANtpDefVal::kDouble),
  wEMV(ANtpDefVal::kDouble),
  wHadU(ANtpDefVal::kDouble),
  wHadV(ANtpDefVal::kDouble),
  wTrkU(ANtpDefVal::kDouble),
  wTrkV(ANtpDefVal::kDouble),
  wXTkU(ANtpDefVal::kDouble),
  wXTkV(ANtpDefVal::kDouble),
  wHalU(ANtpDefVal::kDouble),
  wHalV(ANtpDefVal::kDouble)
{}

SubShowerVar::~SubShowerVar(){}

void SubShowerVar::Zero()
{
  ncluster = 0;
  nclusterU = 0;
  nclusterV = 0;
  nPhysClusterU = 0;
  nPhysClusterV = 0;
  nstp0U = 0;
  nstp0V = 0;
  E2to1U = 0;
  E2to1V = 0;
  PHAvgIDU = 0;
  PHAvgIDV = 0;
  PHAvgProbEMU = 0;
  PHAvgProbEMV = 0;
  PHFracRMSU = 0;
  PHFracRMSV = 0;
  PHAvgDevU = 0;
  PHAvgDevV = 0;
  pid = 0;
  nclusterU_he = 0;
  nclusterV_he = 0;
  nEMU = 0;
  nEMV = 0;
  nHadU = 0;
  nHadV = 0;
  nTrkU = 0;
  nTrkV = 0;
  nXTkU = 0;
  nXTkV = 0;
  nHalU = 0;
  nHalV = 0;
  nEMU_he = 0;
  nEMV_he = 0;
  nHadU_he = 0;
  nHadV_he = 0;
  nTrkU_he = 0;
  nTrkV_he = 0;
  nXTkU_he = 0;
  nXTkV_he = 0;
  nHalU_he = 0;
  nHalV_he = 0;
  eEMU = 0;
  eEMV = 0;
  eHadU = 0;
  eHadV = 0;
  eTrkU = 0;
  eTrkV = 0;
  eXTkU = 0;
  eXTkV = 0;
  eHalU = 0;
  eHalV = 0;  
  wEMU = 0;
  wEMV = 0;
  wHadU = 0;
  wHadV = 0;
  wTrkU = 0;
  wTrkV = 0;
  wXTkU = 0;
  wXTkV = 0;
  wHalU = 0;
  wHalV = 0;
}

void SubShowerVar::Reset()
{
  ncluster = ANtpDefVal::kInt;
  nclusterU = ANtpDefVal::kInt;
  nclusterV = ANtpDefVal::kInt;
  nPhysClusterU = ANtpDefVal::kInt;
  nPhysClusterV = ANtpDefVal::kInt;
  nstp0U = ANtpDefVal::kInt;
  nstp0V = ANtpDefVal::kInt;
  E2to1U = ANtpDefVal::kDouble;
  E2to1V = ANtpDefVal::kDouble;
  PHAvgIDU = ANtpDefVal::kDouble;
  PHAvgIDV = ANtpDefVal::kDouble;
  PHAvgProbEMU = ANtpDefVal::kDouble;
  PHAvgProbEMV = ANtpDefVal::kDouble;
  PHFracRMSU = ANtpDefVal::kDouble;
  PHFracRMSV = ANtpDefVal::kDouble;
  PHAvgDevU = ANtpDefVal::kDouble;
  PHAvgDevV = ANtpDefVal::kDouble;
  pid = ANtpDefVal::kDouble;
  nclusterU_he = ANtpDefVal::kInt;
  nclusterV_he = ANtpDefVal::kInt;
  nEMU = ANtpDefVal::kInt;
  nEMV = ANtpDefVal::kInt;
  nHadU = ANtpDefVal::kInt;
  nHadV = ANtpDefVal::kInt;
  nTrkU = ANtpDefVal::kInt;
  nTrkV = ANtpDefVal::kInt;
  nXTkU = ANtpDefVal::kInt;
  nXTkV = ANtpDefVal::kInt;
  nHalU = ANtpDefVal::kInt;
  nHalV = ANtpDefVal::kInt;
  nEMU_he = ANtpDefVal::kInt;
  nEMV_he = ANtpDefVal::kInt;
  nHadU_he = ANtpDefVal::kInt;
  nHadV_he = ANtpDefVal::kInt;
  nTrkU_he = ANtpDefVal::kInt;
  nTrkV_he = ANtpDefVal::kInt;
  nXTkU_he = ANtpDefVal::kInt;
  nXTkV_he = ANtpDefVal::kInt;
  nHalU_he = ANtpDefVal::kInt;
  nHalV_he = ANtpDefVal::kInt;
  eEMU = ANtpDefVal::kDouble;
  eEMV = ANtpDefVal::kDouble;
  eHadU = ANtpDefVal::kDouble;
  eHadV = ANtpDefVal::kDouble;
  eTrkU = ANtpDefVal::kDouble;
  eTrkV = ANtpDefVal::kDouble;
  eXTkU = ANtpDefVal::kDouble;
  eXTkV = ANtpDefVal::kDouble;
  eHalU = ANtpDefVal::kDouble;
  eHalV = ANtpDefVal::kDouble;  
  wEMU = ANtpDefVal::kDouble;
  wEMV = ANtpDefVal::kDouble;
  wHadU = ANtpDefVal::kDouble;
  wHadV = ANtpDefVal::kDouble;
  wTrkU = ANtpDefVal::kDouble;
  wTrkV = ANtpDefVal::kDouble;
  wXTkU = ANtpDefVal::kDouble;
  wXTkV = ANtpDefVal::kDouble;
  wHalU = ANtpDefVal::kDouble;
  wHalV = ANtpDefVal::kDouble;
}
