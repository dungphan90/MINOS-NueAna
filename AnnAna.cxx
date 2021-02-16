#include "NueAna/AnnAna.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NueStandard.h"
#include "math.h"
#include "TFile.h"
#include <iostream>
#include <fstream>

using namespace std;
ClassImp(AnnAna)

TMultiLayerPerceptron* AnnAna::fneuralNet_6inp = 0;
TMultiLayerPerceptron* AnnAna::fneuralNet_30inp = 0;
TMultiLayerPerceptron* AnnAna::fneuralNet_11inp = 0;
TMultiLayerPerceptron* AnnAna::fneuralNet_11inp_daikon04 = 0;
TMultiLayerPerceptron* AnnAna::fneuralNet_14inp_daikon04 = 0;

AnnAna::AnnAna(NueRecord &nr):
  nuerec(nr)
{
  static bool first = true;
  char *srt_dir = getenv("SRT_PRIVATE_CONTEXT");
  char annfile[10000];
  char annfile2[10000];
  char annfile3[10000];
  char annfile4[10000];
  char annfile5[10000];

  if(first){
    sprintf(annfile,"%s/NueAna/data/ann040907_6inp.root",srt_dir);
    ifstream Test(annfile);
    if (!Test){
      srt_dir = getenv("SRT_PUBLIC_CONTEXT");
      sprintf(annfile,"%s/NueAna/data/ann040907_6inp.root",srt_dir);
      ifstream Test_again(annfile);
      if (!Test_again){
        cout<<"Couldn't find ANN object, blame Tingjun"<<endl;
        exit(0);
      }
    }
    static TFile *f1 = TFile::Open(annfile);
    fneuralNet_6inp = (TMultiLayerPerceptron*) f1->Get("mlp");
  }

  if (!fneuralNet_6inp) {
    cout<<"Couldn't find ANN object for FracVarAna, blame Tingjun"<<endl;
    exit(0);
  }

  if(first){
    srt_dir = getenv("SRT_PRIVATE_CONTEXT");
    sprintf(annfile2,"%s/NueAna/data/ann050107_30inp.root",srt_dir);
    ifstream Test2(annfile2);
    if (!Test2){
      srt_dir = getenv("SRT_PUBLIC_CONTEXT");
      sprintf(annfile2,"%s/NueAna/data/ann050107_30inp.root",srt_dir);
      ifstream Test_again2(annfile2);
      if (!Test_again2){
        cout<<"Couldn't find ANN object, blame Tingjun"<<endl;
        exit(0);
      }
    }
    static TFile *f2 = TFile::Open(annfile2);
    fneuralNet_30inp = (TMultiLayerPerceptron*) f2->Get("mlp");
  }

  if (!fneuralNet_30inp) {
    cout<<"Couldn't find ANN object, blame Tingjun"<<endl;
    exit(0);
  }

  if(first){
    srt_dir = getenv("SRT_PRIVATE_CONTEXT");
    sprintf(annfile3,"%s/NueAna/data/ann082608_11inp.root",srt_dir);
    ifstream Test2(annfile3);
    if (!Test2){
      srt_dir = getenv("SRT_PUBLIC_CONTEXT");
      sprintf(annfile3,"%s/NueAna/data/ann082608_11inp.root",srt_dir);
      ifstream Test_again2(annfile3);
      if (!Test_again2){
        cout<<"Couldn't find ANN object, blame Tingjun"<<endl;
        exit(0);
      }
    }
    static TFile *f3 = TFile::Open(annfile3);
    fneuralNet_11inp = (TMultiLayerPerceptron*) f3->Get("mlp");
  }

  if (!fneuralNet_11inp) {
    cout<<"Couldn't find ANN object, blame Tingjun"<<endl;
    exit(0);
  }

  if (first) {
    srt_dir = getenv("SRT_PRIVATE_CONTEXT");
    sprintf(annfile4,"%s/NueAna/data/ann11_400_14_9.root",srt_dir);
    ifstream Test2(annfile4); 
    if (!Test2){
      srt_dir = getenv("SRT_PUBLIC_CONTEXT");
      sprintf(annfile4,"%s/NueAna/data/ann11_400_14_9.root",srt_dir);
      ifstream Test_again2(annfile4);
      if (!Test_again2){
        cout<<"Couldn't find ANN object, blame Jiajie"<<endl;
        exit(0);
      }
    }
    static TFile *f4 = TFile::Open(annfile4);
    fneuralNet_11inp_daikon04 = (TMultiLayerPerceptron*) f4->Get("mlp");
  }

  if (!fneuralNet_11inp_daikon04) {
    cout<<"Couldn't find ANN object, blame Jiajie"<<endl;
    exit(0);
  }

  if (first) {
    srt_dir = getenv("SRT_PRIVATE_CONTEXT");
    sprintf(annfile5,"%s/NueAna/data/ann14_400_14_9.root",srt_dir);
    ifstream Test2(annfile5);
    if (!Test2){
      srt_dir = getenv("SRT_PUBLIC_CONTEXT");
      sprintf(annfile5,"%s/NueAna/data/ann14_400_14_9.root",srt_dir);
      ifstream Test_again2(annfile5);
      if (!Test_again2){
        cout<<"Couldn't find ANN object, blame Jiajie"<<endl;
        exit(0);
      }
    }
    static TFile *f5 = TFile::Open(annfile5);
    fneuralNet_14inp_daikon04 = (TMultiLayerPerceptron*) f5->Get("mlp");
  }

  if (!fneuralNet_14inp_daikon04) {
    cout<<"Couldn't find ANN object, blame Jiajie"<<endl;
    exit(0);
  }

  if (first){
    cout<<"Reading ann_6inp from : "<<annfile<<endl;
    cout<<"Reading ann_30inp from : "<<annfile2<<endl;
    cout<<"Reading ann_11inp from : "<<annfile3<<endl;
    cout<<"Reading ann_11inp_daikon04 from : "<<annfile4<<endl;
    cout<<"Reading ann_14inp_daikon04 from : "<<annfile5<<endl;
    first = false;
  }
}

AnnAna::~AnnAna()
{}

void AnnAna::Analyze(){


  int pass = 1;
  if (nuerec.shwfit.par_a<-1000) pass = 0;
  if (nuerec.shwfit.par_b<-1000) pass = 0;
  if (nuerec.shwfit.uv_molrad_peak<-1000) pass = 0;
  if (nuerec.shwfit.uv_rms<-1000) pass = 0;
  if (nuerec.angcluster.fACluRmsShwAxis<-1000) pass = 0;
  if (nuerec.angcluster.fACluRmsZAxis<-1000) pass = 0;
  if (nuerec.hitcalc.fHitLongEnergy<-1000) pass = 0;
  if (nuerec.hitcalc.fHitLongEnergy>1000) pass = 0;
  if (nuerec.hitcalc.fHitTransCMEnergy<-1000) pass = 0;
  if (nuerec.hitcalc.fHitTransCMEnergy>2000) pass = 0;
  if (nuerec.hitcalc.fHitTransEnergyRatio<-1000) pass = 0;
  if (nuerec.hitcalc.fHitLongEnergyRatio<-1000) pass = 0;
  if (nuerec.hitcalc.fHitFarAngle<-1000) pass = 0;
  if (nuerec.hitcalc.fHitPeakAngle<-1000) pass = 0;
  if (nuerec.mstvars.e4w<-1000) pass = 0;
  if (nuerec.mstvars.e4w>500) pass = 0;
  if (nuerec.mstvars.o4w<-1000) pass = 0;
  if (nuerec.mstvars.o4w>500) pass = 0;
  if (nuerec.mstvars.osmtot<-1000) pass = 0;
  if (nuerec.mstvars.osmtot>500) pass = 0;
  if (nuerec.mstvars.eeprob<-1000) pass = 0;
  if (nuerec.mstvars.oeprob<-1000) pass = 0;
  if (nuerec.fracvars.fract_2_planes<-1000) pass = 0;
  if (nuerec.fracvars.fract_4_planes<-1000) pass = 0;
  if (nuerec.fracvars.fract_6_planes<-1000) pass = 0;
  if (nuerec.fracvars.fract_4_counters<-1000) pass = 0;
  if (nuerec.fracvars.fract_8_counters<-1000) pass = 0;
  if (nuerec.fracvars.fract_12_counters<-1000) pass = 0;
  if (nuerec.fracvars.fract_road<-1000) pass = 0;
  if (nuerec.subshowervars.E2to1U<-1000) pass = 0;
  if (nuerec.subshowervars.E2to1V<-1000) pass = 0;
  if (nuerec.subshowervars.PHFracRMSU<-1000) pass = 0;
  if (nuerec.subshowervars.PHFracRMSV<-1000) pass = 0;
  if (nuerec.subshowervars.PHAvgDevU<-1000) pass = 0;
  if (nuerec.subshowervars.PHAvgDevV<-1000) pass = 0;

  if (!pass) {
    nuerec.ann.pid_30inp = ANtpDefaultValue::kDouble;
  }
  else {
    
    Double_t params[30];
    params[0] = nuerec.shwfit.par_a;
    params[1] = nuerec.shwfit.par_b;
    params[2] = nuerec.shwfit.uv_molrad_peak;
    params[3] = nuerec.shwfit.uv_rms;
    params[4] = nuerec.angcluster.fACluRmsShwAxis;
    params[5] = nuerec.angcluster.fACluRmsZAxis;
    params[6] = nuerec.hitcalc.fHitLongEnergy;
    params[7] = nuerec.hitcalc.fHitTransCMEnergy;
    params[8] = nuerec.hitcalc.fHitTransEnergyRatio;
    params[9] = nuerec.hitcalc.fHitLongEnergyRatio;
    params[10] = nuerec.hitcalc.fHitFarAngle;
    params[11] = nuerec.hitcalc.fHitPeakAngle;
    params[12] = nuerec.mstvars.e4w;
    params[13] = nuerec.mstvars.o4w;
    params[14] = nuerec.mstvars.osmtot;
    params[15] = nuerec.mstvars.eeprob;
    params[16] = nuerec.mstvars.oeprob;
    params[17] = nuerec.fracvars.fract_2_planes;
    params[18] = nuerec.fracvars.fract_4_planes;
    params[19] = nuerec.fracvars.fract_6_planes;
    params[20] = nuerec.fracvars.fract_4_counters;
    params[21] = nuerec.fracvars.fract_8_counters;
    params[22] = nuerec.fracvars.fract_12_counters;
    params[23] = nuerec.fracvars.fract_road;
    params[24] = nuerec.subshowervars.E2to1U;
    params[25] = nuerec.subshowervars.E2to1V;
    params[26] = nuerec.subshowervars.PHFracRMSU;
    params[27] = nuerec.subshowervars.PHFracRMSV;
    params[28] = nuerec.subshowervars.PHAvgDevU;
    params[29] = nuerec.subshowervars.PHAvgDevV;   

    nuerec.ann.pid_30inp = fneuralNet_30inp->Evaluate(0,params);
  }
  
  double par[6];
  par[0] = nuerec.fracvars.fract_road;
  par[1] = nuerec.shwfit.par_a;
  par[2] = nuerec.shwfit.par_b;
  par[3] = nuerec.hitcalc.fHitLongEnergy;
  par[4] = nuerec.mstvars.e4w+nuerec.mstvars.o4w;
  par[5] = nuerec.subshowervars.PHFracRMSU+nuerec.subshowervars.PHFracRMSV;
  bool inputok = true;
  for (int idx = 0; idx<6; idx++){
    inputok = inputok&&par[idx]>-100;
  }
  if (inputok) nuerec.ann.pid_6inp = fneuralNet_6inp->Evaluate(0,par);
  else nuerec.ann.pid_6inp = ANtpDefaultValue::kDouble;

  // retuned ANN11 PID - with dogwood reconstruction and daikon04 MC 

  double pars[11];

  pars[0] = nuerec.shwfit.par_a;
  pars[1] = nuerec.shwfit.par_b;
  pars[2] = nuerec.shwfit.uv_molrad_peak_9s_2pe_dw;
  pars[3] = nuerec.shwfit.uv_rms_9s_2pe_dw;
  pars[4] = nuerec.mstvars.e4w+nuerec.mstvars.o4w;
  pars[5] = nuerec.fracvars.fract_2_planes;
  pars[6] = nuerec.fracvars.fract_4_planes;
  pars[7] = nuerec.fracvars.fract_6_planes;
  pars[8] = nuerec.fracvars.fract_8_counters;
  pars[9] = nuerec.fracvars.fract_road;
  pars[10]=nuerec.shwfit.LongE;
  pass = 1;
  if (nuerec.shwfit.par_a<-1000) pass = 0;
  if (nuerec.shwfit.par_b<-1000) pass = 0;
  if (nuerec.shwfit.uv_molrad_peak_9s_2pe_dw<-1000) pass = 0;
  if (nuerec.shwfit.uv_rms_9s_2pe_dw<-1000) pass = 0;
  if (nuerec.mstvars.e4w<-1000) pass = 0;
  if (nuerec.mstvars.e4w>500) pass = 0;
  if (nuerec.mstvars.o4w<-1000) pass = 0;
  if (nuerec.mstvars.o4w>500) pass = 0;
  if (nuerec.fracvars.fract_2_planes<-1000) pass = 0;
  if (nuerec.fracvars.fract_4_planes<-1000) pass = 0;
  if (nuerec.fracvars.fract_6_planes<-1000) pass = 0;
  if (nuerec.fracvars.fract_8_counters<-1000) pass = 0;
  if (nuerec.fracvars.fract_road<-1000) pass = 0;
  if (nuerec.shwfit.LongE<-1000) pass = 0;
  if (nuerec.shwfit.LongE>1000) pass = 0;
  
  if (pass) nuerec.ann.pid_11inp = fneuralNet_11inp->Evaluate(0,pars);
  else nuerec.ann.pid_11inp = ANtpDefaultValue::kDouble;

  if (pass) nuerec.ann.pid_11inp_daikon04 = fneuralNet_11inp_daikon04->Evaluate(0,pars);
  else nuerec.ann.pid_11inp_daikon04 = ANtpDefaultValue::kDouble;

  // ANN14 PID - with dogwood reconstruction and daikon04 MC
  // Since it includes event energy, we should call the energy correction first 

  NueConvention::NueEnergyCorrection(&nuerec);
 
  double paras[14];
  paras[0] = nuerec.shwfit.par_a;
  paras[1] = nuerec.shwfit.par_b;
  paras[2] = nuerec.shwfit.uv_rms_9s_2pe_dw;
  paras[3] = nuerec.fracvars.fract_road;
  paras[4] = nuerec.shwfit.LongE;
  paras[5] = nuerec.fracvars.shw_max;
  paras[6] = nuerec.fracvars.shw_slp;
  paras[7] = nuerec.srevent.phNueGeV;
  paras[8] = nuerec.fracvars.dis2stp;
  paras[9] = nuerec.fracvars.fract_1_plane;
  paras[10] = nuerec.fracvars.fract_5_planes;
  paras[11] = nuerec.fracvars.fract_6_counters;
  paras[12] = nuerec.fracvars.fract_20_counters;
  paras[13] = TMath::Abs(nuerec.srevent.endPlane - nuerec.srevent.begPlane);
   
  pass = 1;

  if (nuerec.shwfit.par_a<-1000) pass = 0;
  if (nuerec.shwfit.par_b<-1000) pass = 0;
  if (nuerec.shwfit.uv_rms_9s_2pe_dw<-1000) pass = 0;
  if (nuerec.fracvars.fract_1_plane<-1000) pass = 0;
  if (nuerec.fracvars.fract_5_planes<-1000) pass = 0;
  if (nuerec.fracvars.fract_6_counters<-1000) pass = 0;
  if (nuerec.fracvars.fract_20_counters<-1000) pass = 0;
  if (nuerec.fracvars.fract_road<-1000) pass = 0;
  if (nuerec.shwfit.LongE<-1000) pass = 0;
  if (nuerec.shwfit.LongE>1000) pass = 0;
  if (nuerec.fracvars.shw_max<-1000) pass = 0;
  if (nuerec.fracvars.dis2stp<-1000) pass = 0;
  if (nuerec.fracvars.shw_slp<0) pass = 0;
  if (nuerec.srevent.endPlane<-1000) pass = 0;
  if (nuerec.srevent.begPlane<-1000) pass = 0;

  if (pass) {
    nuerec.ann.pid_14inp_daikon04 = fneuralNet_14inp_daikon04->Evaluate(0,paras);
    nuerec.ann.pid = fneuralNet_14inp_daikon04->Evaluate(0,paras);
  } else {
    nuerec.ann.pid_14inp_daikon04 = ANtpDefaultValue::kDouble;
    nuerec.ann.pid = ANtpDefaultValue::kDouble;
  }
}

double AnnAna::value(int index,double in0,double in1,double in2,double in3,double in4,double in5,double in6,double in7,double in8,double in9,double in10,double in11,double in12,double in13,double in14,double in15,double in16,double in17,double in18,double in19,double in20,double in21,double in22,double in23,double in24,double in25,double in26,double in27,double in28,double in29) {
   input0 = (in0 - 3.26992)/2.06389;
   input1 = (in1 - 0.642738)/0.531793;
   input2 = (in2 - 6.09479)/3.08303;
   input3 = (in3 - 3.25365)/1.08684;
   input4 = (in4 - 0.0736766)/0.0476116;
   input5 = (in5 - 0.0686119)/0.0525439;
   input6 = (in6 - 119.134)/78.7005;
   input7 = (in7 - 150.227)/151.711;
   input8 = (in8 - 0.140695)/0.135734;
   input9 = (in9 - 0.79371)/0.126991;
   input10 = (in10 - 0.570399)/0.325531;
   input11 = (in11 - 0.407685)/0.328878;
   input12 = (in12 - 36.1366)/31.5127;
   input13 = (in13 - 35.6709)/31.1802;
   input14 = (in14 - 32.4559)/15.8662;
   input15 = (in15 - 5.70368)/5.50635;
   input16 = (in16 - 5.59754)/5.46595;
   input17 = (in17 - 0.462635)/0.118817;
   input18 = (in18 - 0.694673)/0.133558;
   input19 = (in19 - 0.834602)/0.124673;
   input20 = (in20 - 0.482597)/0.107022;
   input21 = (in21 - 0.674527)/0.106117;
   input22 = (in22 - 0.781535)/0.0969513;
   input23 = (in23 - 0.721024)/0.162524;
   input24 = (in24 - 0.255528)/0.265611;
   input25 = (in25 - 0.251731)/0.265088;
   input26 = (in26 - 0.289288)/0.0953619;
   input27 = (in27 - 0.291994)/0.0972296;
   input28 = (in28 - 0.0249709)/0.0264772;
   input29 = (in29 - 0.0248471)/0.0278384;
   switch(index) {
     case 0:
         return neuron0xa0068f8();
     default:
         return 0.;
   }
}

double AnnAna::neuron0xa0253a8() {
   return input0;
}

double AnnAna::neuron0xa025510() {
   return input1;
}

double AnnAna::neuron0xa025678() {
   return input2;
}

double AnnAna::neuron0xa0257e0() {
   return input3;
}

double AnnAna::neuron0xa025948() {
   return input4;
}

double AnnAna::neuron0xa025ab8() {
   return input5;
}

double AnnAna::neuron0xa025c20() {
   return input6;
}

double AnnAna::neuron0xa025d88() {
   return input7;
}

double AnnAna::neuron0xa025ef0() {
   return input8;
}

double AnnAna::neuron0xa026058() {
   return input9;
}

double AnnAna::neuron0xa0261c0() {
   return input10;
}

double AnnAna::neuron0xa026328() {
   return input11;
}

double AnnAna::neuron0xa026490() {
   return input12;
}

double AnnAna::neuron0xa0265f8() {
   return input13;
}

double AnnAna::neuron0xa026760() {
   return input14;
}

double AnnAna::neuron0xa0268c8() {
   return input15;
}

double AnnAna::neuron0xa026a30() {
   return input16;
}

double AnnAna::neuron0xa026ca8() {
   return input17;
}

double AnnAna::neuron0xa026d80() {
   return input18;
}

double AnnAna::neuron0xa026ee8() {
   return input19;
}

double AnnAna::neuron0xa027050() {
   return input20;
}

double AnnAna::neuron0xa0271b8() {
   return input21;
}

double AnnAna::neuron0xa027320() {
   return input22;
}

double AnnAna::neuron0xa027488() {
   return input23;
}

double AnnAna::neuron0xa0275f0() {
   return input24;
}

double AnnAna::neuron0xa027758() {
   return input25;
}

double AnnAna::neuron0xa0278c0() {
   return input26;
}

double AnnAna::neuron0xa0064a0() {
   return input27;
}

double AnnAna::neuron0xa006608() {
   return input28;
}

double AnnAna::neuron0xa006770() {
   return input29;
}

double AnnAna::neuron0xa0069f0() {
   double input = -1.8388;
   input += synapse0xa006b40();
   input += synapse0xa006b68();
   input += synapse0xa006b90();
   input += synapse0xa006bb8();
   input += synapse0xa006be0();
   input += synapse0xa006c08();
   input += synapse0xa006c30();
   input += synapse0xa006c58();
   input += synapse0xa006c80();
   input += synapse0xa006ca8();
   input += synapse0xa006cd0();
   input += synapse0xa006cf8();
   input += synapse0xa006d20();
   input += synapse0xa006d48();
   input += synapse0xa006d70();
   input += synapse0xa006d98();
   input += synapse0xa006dc0();
   input += synapse0xa006ef8();
   input += synapse0xa006f20();
   input += synapse0xa006f48();
   input += synapse0xa006f70();
   input += synapse0xa006f98();
   input += synapse0xa006fc0();
   input += synapse0xa006fe8();
   input += synapse0xa007010();
   input += synapse0xa007038();
   input += synapse0xa007060();
   input += synapse0xa007088();
   input += synapse0xa0070b0();
   input += synapse0xa0070d8();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa007100() {
   double input = 0.592808;
   input += synapse0xa0071e0();
   input += synapse0xa007208();
   input += synapse0xa007230();
   input += synapse0xa0279e0();
   input += synapse0xa027a08();
   input += synapse0xa0315d0();
   input += synapse0xa0315f8();
   input += synapse0xa006e70();
   input += synapse0xa006e98();
   input += synapse0xa006ec0();
   input += synapse0xa007360();
   input += synapse0xa007388();
   input += synapse0xa0073b0();
   input += synapse0xa0073d8();
   input += synapse0xa007400();
   input += synapse0xa007428();
   input += synapse0xa007450();
   input += synapse0xa007500();
   input += synapse0xa007528();
   input += synapse0xa007550();
   input += synapse0xa007578();
   input += synapse0xa0075a0();
   input += synapse0xa0075c8();
   input += synapse0xa0075f0();
   input += synapse0xa007618();
   input += synapse0xa007640();
   input += synapse0xa007668();
   input += synapse0xa007690();
   input += synapse0xa0076b8();
   input += synapse0xa0076e0();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa007708() {
   double input = -0.303101;
   input += synapse0xa026c20();
   input += synapse0xa026c48();
   input += synapse0xa026c70();
   input += synapse0xa007938();
   input += synapse0xa007960();
   input += synapse0xa007258();
   input += synapse0xa007280();
   input += synapse0xa0072a8();
   input += synapse0xa0072d0();
   input += synapse0xa0072f8();
   input += synapse0xa007320();
   input += synapse0xa007b90();
   input += synapse0xa007bb8();
   input += synapse0xa007be0();
   input += synapse0xa007c08();
   input += synapse0xa007c30();
   input += synapse0xa007c58();
   input += synapse0xa007d08();
   input += synapse0xa007d30();
   input += synapse0xa007d58();
   input += synapse0xa007d80();
   input += synapse0xa007da8();
   input += synapse0xa007dd0();
   input += synapse0xa007df8();
   input += synapse0xa007e20();
   input += synapse0xa007e48();
   input += synapse0xa007e70();
   input += synapse0xa007e98();
   input += synapse0xa007ec0();
   input += synapse0xa007ee8();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa007f10() {
   double input = 0.160484;
   input += synapse0xa008038();
   input += synapse0xa008060();
   input += synapse0xa008088();
   input += synapse0xa0080b0();
   input += synapse0xa0080d8();
   input += synapse0xa008100();
   input += synapse0xa008128();
   input += synapse0xa008150();
   input += synapse0xa008178();
   input += synapse0xa0081a0();
   input += synapse0xa0081c8();
   input += synapse0xa0081f0();
   input += synapse0xa008218();
   input += synapse0xa008240();
   input += synapse0xa008268();
   input += synapse0xa008290();
   input += synapse0xa0082b8();
   input += synapse0xa008368();
   input += synapse0xa008390();
   input += synapse0xa0083b8();
   input += synapse0xa0083e0();
   input += synapse0xa008408();
   input += synapse0xa008430();
   input += synapse0xa008458();
   input += synapse0xa008480();
   input += synapse0xa0084a8();
   input += synapse0xa0084d0();
   input += synapse0xa0084f8();
   input += synapse0xa008520();
   input += synapse0xa008548();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa008570() {
   double input = 0.034072;
   input += synapse0xa008698();
   input += synapse0xa0086c0();
   input += synapse0xa0086e8();
   input += synapse0xa008710();
   input += synapse0xa008738();
   input += synapse0xa008760();
   input += synapse0xa008788();
   input += synapse0xa0087b0();
   input += synapse0xa0087d8();
   input += synapse0xa007988();
   input += synapse0xa0079b0();
   input += synapse0xa0079d8();
   input += synapse0xa007a00();
   input += synapse0xa007a28();
   input += synapse0xa007a50();
   input += synapse0xa007a78();
   input += synapse0xa007aa0();
   input += synapse0xa007b50();
   input += synapse0xa008c08();
   input += synapse0xa008c30();
   input += synapse0xa008c58();
   input += synapse0xa008c80();
   input += synapse0xa008ca8();
   input += synapse0xa008cd0();
   input += synapse0xa008cf8();
   input += synapse0xa008d20();
   input += synapse0xa008d48();
   input += synapse0xa008d70();
   input += synapse0xa008d98();
   input += synapse0xa008dc0();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa008de8() {
   double input = 0.0177348;
   input += synapse0xa008f10();
   input += synapse0xa008f38();
   input += synapse0xa008f60();
   input += synapse0xa008f88();
   input += synapse0xa008fb0();
   input += synapse0xa008fd8();
   input += synapse0xa009000();
   input += synapse0xa009028();
   input += synapse0xa009050();
   input += synapse0xa009078();
   input += synapse0xa0090a0();
   input += synapse0xa0090c8();
   input += synapse0xa0090f0();
   input += synapse0xa009118();
   input += synapse0xa009140();
   input += synapse0xa009168();
   input += synapse0xa009190();
   input += synapse0xa009240();
   input += synapse0xa009268();
   input += synapse0xa009290();
   input += synapse0xa0092b8();
   input += synapse0xa0092e0();
   input += synapse0xa009308();
   input += synapse0xa009330();
   input += synapse0xa009358();
   input += synapse0xa009380();
   input += synapse0xa0093a8();
   input += synapse0xa0093d0();
   input += synapse0xa0093f8();
   input += synapse0xa009420();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa009448() {
   double input = -0.490912;
   input += synapse0xa009570();
   input += synapse0xa009598();
   input += synapse0xa0095c0();
   input += synapse0xa0095e8();
   input += synapse0xa009610();
   input += synapse0xa009638();
   input += synapse0xa009660();
   input += synapse0xa009688();
   input += synapse0xa0096b0();
   input += synapse0xa0096d8();
   input += synapse0xa009700();
   input += synapse0xa009728();
   input += synapse0xa009750();
   input += synapse0xa009778();
   input += synapse0xa0097a0();
   input += synapse0xa0097c8();
   input += synapse0xa0097f0();
   input += synapse0xa0098a0();
   input += synapse0xa0098c8();
   input += synapse0xa0098f0();
   input += synapse0xa009918();
   input += synapse0xa009940();
   input += synapse0xa009968();
   input += synapse0xa009990();
   input += synapse0xa0099b8();
   input += synapse0xa0099e0();
   input += synapse0xa009a08();
   input += synapse0xa009a30();
   input += synapse0xa009a58();
   input += synapse0xa009a80();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa009aa8() {
   double input = -0.420931;
   input += synapse0xa009bd0();
   input += synapse0xa009bf8();
   input += synapse0xa009c20();
   input += synapse0xa009c48();
   input += synapse0xa009c70();
   input += synapse0xa009c98();
   input += synapse0xa009cc0();
   input += synapse0xa009ce8();
   input += synapse0xa009d10();
   input += synapse0xa009d38();
   input += synapse0xa009d60();
   input += synapse0xa009d88();
   input += synapse0xa009db0();
   input += synapse0xa009dd8();
   input += synapse0xa009e00();
   input += synapse0xa009e28();
   input += synapse0xa009e50();
   input += synapse0xa009f00();
   input += synapse0xa009f28();
   input += synapse0xa009f50();
   input += synapse0xa009f78();
   input += synapse0xa009fa0();
   input += synapse0xa009fc8();
   input += synapse0xa009ff0();
   input += synapse0xa00a018();
   input += synapse0xa00a040();
   input += synapse0xa00a068();
   input += synapse0xa00a090();
   input += synapse0xa00a0b8();
   input += synapse0xa00a0e0();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00a108() {
   double input = -0.189881;
   input += synapse0xa00a230();
   input += synapse0xa00a258();
   input += synapse0xa00a280();
   input += synapse0xa00a2a8();
   input += synapse0xa00a2d0();
   input += synapse0xa00a2f8();
   input += synapse0xa00a320();
   input += synapse0xa00a348();
   input += synapse0xa00a370();
   input += synapse0xa00a398();
   input += synapse0xa00a3c0();
   input += synapse0xa00a3e8();
   input += synapse0xa00a410();
   input += synapse0xa00a438();
   input += synapse0xa00a460();
   input += synapse0xa00a488();
   input += synapse0xa00a4b0();
   input += synapse0xa024cf8();
   input += synapse0xa008800();
   input += synapse0xa008828();
   input += synapse0xa008850();
   input += synapse0xa008878();
   input += synapse0xa0088a0();
   input += synapse0xa0088c8();
   input += synapse0xa0088f0();
   input += synapse0xa008918();
   input += synapse0xa008940();
   input += synapse0xa008968();
   input += synapse0xa008990();
   input += synapse0xa0089b8();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa0089e0() {
   double input = 0.447111;
   input += synapse0xa008b30();
   input += synapse0xa008b58();
   input += synapse0xa008b80();
   input += synapse0xa008ba8();
   input += synapse0xa008bd0();
   input += synapse0xa00ad68();
   input += synapse0xa00ad90();
   input += synapse0xa00adb8();
   input += synapse0xa00ade0();
   input += synapse0xa00ae08();
   input += synapse0xa00ae30();
   input += synapse0xa00ae58();
   input += synapse0xa00ae80();
   input += synapse0xa00aea8();
   input += synapse0xa00aed0();
   input += synapse0xa00aef8();
   input += synapse0xa00af20();
   input += synapse0xa00afd0();
   input += synapse0xa00aff8();
   input += synapse0xa00b020();
   input += synapse0xa00b048();
   input += synapse0xa00b070();
   input += synapse0xa00b098();
   input += synapse0xa00b0c0();
   input += synapse0xa00b0e8();
   input += synapse0xa00b110();
   input += synapse0xa00b138();
   input += synapse0xa00b160();
   input += synapse0xa00b188();
   input += synapse0xa00b1b0();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00b1d8() {
   double input = 0.717482;
   input += synapse0xa00b300();
   input += synapse0xa00b328();
   input += synapse0xa00b350();
   input += synapse0xa00b378();
   input += synapse0xa00b3a0();
   input += synapse0xa00b3c8();
   input += synapse0xa00b3f0();
   input += synapse0xa00b418();
   input += synapse0xa00b440();
   input += synapse0xa00b468();
   input += synapse0xa00b490();
   input += synapse0xa00b4b8();
   input += synapse0xa00b4e0();
   input += synapse0xa00b508();
   input += synapse0xa00b530();
   input += synapse0xa00b558();
   input += synapse0xa00b580();
   input += synapse0xa00b630();
   input += synapse0xa00b658();
   input += synapse0xa00b680();
   input += synapse0xa00b6a8();
   input += synapse0xa00b6d0();
   input += synapse0xa00b6f8();
   input += synapse0xa00b720();
   input += synapse0xa00b748();
   input += synapse0xa00b770();
   input += synapse0xa00b798();
   input += synapse0xa00b7c0();
   input += synapse0xa00b7e8();
   input += synapse0xa00b810();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00b838() {
   double input = -0.575772;
   input += synapse0xa00b960();
   input += synapse0xa00b988();
   input += synapse0xa00b9b0();
   input += synapse0xa00b9d8();
   input += synapse0xa00ba00();
   input += synapse0xa00ba28();
   input += synapse0xa00ba50();
   input += synapse0xa00ba78();
   input += synapse0xa00baa0();
   input += synapse0xa00bac8();
   input += synapse0xa00baf0();
   input += synapse0xa00bb18();
   input += synapse0xa00bb40();
   input += synapse0xa00bb68();
   input += synapse0xa00bb90();
   input += synapse0xa00bbb8();
   input += synapse0xa00bbe0();
   input += synapse0xa00bc90();
   input += synapse0xa00bcb8();
   input += synapse0xa00bce0();
   input += synapse0xa00bd08();
   input += synapse0xa00bd30();
   input += synapse0xa00bd58();
   input += synapse0xa00bd80();
   input += synapse0xa00bda8();
   input += synapse0xa00bdd0();
   input += synapse0xa00bdf8();
   input += synapse0xa00be20();
   input += synapse0xa00be48();
   input += synapse0xa00be70();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00be98() {
   double input = 0.176399;
   input += synapse0xa00bfc0();
   input += synapse0xa00bfe8();
   input += synapse0xa00c010();
   input += synapse0xa00c038();
   input += synapse0xa00c060();
   input += synapse0xa00c088();
   input += synapse0xa00c0b0();
   input += synapse0xa00c0d8();
   input += synapse0xa00c100();
   input += synapse0xa00c128();
   input += synapse0xa00c150();
   input += synapse0xa00c178();
   input += synapse0xa00c1a0();
   input += synapse0xa00c1c8();
   input += synapse0xa00c1f0();
   input += synapse0xa00c218();
   input += synapse0xa00c240();
   input += synapse0xa00c2f0();
   input += synapse0xa00c318();
   input += synapse0xa00c340();
   input += synapse0xa00c368();
   input += synapse0xa00c390();
   input += synapse0xa00c3b8();
   input += synapse0xa00c3e0();
   input += synapse0xa00c408();
   input += synapse0xa00c430();
   input += synapse0xa00c458();
   input += synapse0xa00c480();
   input += synapse0xa00c4a8();
   input += synapse0xa00c4d0();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00c4f8() {
   double input = -0.727007;
   input += synapse0xa00c620();
   input += synapse0xa00c648();
   input += synapse0xa00c670();
   input += synapse0xa00c698();
   input += synapse0xa00c6c0();
   input += synapse0xa00c6e8();
   input += synapse0xa00c710();
   input += synapse0xa00c738();
   input += synapse0xa00c760();
   input += synapse0xa00c788();
   input += synapse0xa00c7b0();
   input += synapse0xa00c7d8();
   input += synapse0xa00c800();
   input += synapse0xa00c828();
   input += synapse0xa00c850();
   input += synapse0xa00c878();
   input += synapse0xa00c8a0();
   input += synapse0xa00c950();
   input += synapse0xa00c978();
   input += synapse0xa00c9a0();
   input += synapse0xa00c9c8();
   input += synapse0xa00c9f0();
   input += synapse0xa00ca18();
   input += synapse0xa00ca40();
   input += synapse0xa00ca68();
   input += synapse0xa00ca90();
   input += synapse0xa00cab8();
   input += synapse0xa00cae0();
   input += synapse0xa00cb08();
   input += synapse0xa00cb30();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00cb58() {
   double input = 0.174627;
   input += synapse0xa00cc80();
   input += synapse0xa00cca8();
   input += synapse0xa00ccd0();
   input += synapse0xa00ccf8();
   input += synapse0xa00cd20();
   input += synapse0xa00cd48();
   input += synapse0xa00cd70();
   input += synapse0xa00cd98();
   input += synapse0xa00cdc0();
   input += synapse0xa00cde8();
   input += synapse0xa00ce10();
   input += synapse0xa00ce38();
   input += synapse0xa00ce60();
   input += synapse0xa00ce88();
   input += synapse0xa00ceb0();
   input += synapse0xa00ced8();
   input += synapse0xa00cf00();
   input += synapse0xa00cfb0();
   input += synapse0xa00cfd8();
   input += synapse0xa00d000();
   input += synapse0xa00d028();
   input += synapse0xa00d050();
   input += synapse0xa00d078();
   input += synapse0xa00d0a0();
   input += synapse0xa00d0c8();
   input += synapse0xa00d0f0();
   input += synapse0xa00d118();
   input += synapse0xa00d140();
   input += synapse0xa00d168();
   input += synapse0xa00d190();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00d1d8() {
   double input = -0.308594;
   input += synapse0xa00d2e0();
   input += synapse0xa00d308();
   input += synapse0xa00d330();
   input += synapse0xa00d358();
   input += synapse0xa00d380();
   input += synapse0xa00d3a8();
   input += synapse0xa00d3d0();
   input += synapse0xa00d3f8();
   input += synapse0xa00d420();
   input += synapse0xa00d448();
   input += synapse0xa00d470();
   input += synapse0xa00d498();
   input += synapse0xa00d4c0();
   input += synapse0xa00d4e8();
   input += synapse0xa00d510();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00d538() {
   double input = 0.447124;
   input += synapse0xa00d6a8();
   input += synapse0xa00d6d0();
   input += synapse0xa00d6f8();
   input += synapse0xa00d720();
   input += synapse0xa00d748();
   input += synapse0xa00d770();
   input += synapse0xa00d798();
   input += synapse0xa00d7c0();
   input += synapse0xa00d7e8();
   input += synapse0xa00d810();
   input += synapse0xa00d838();
   input += synapse0xa00d860();
   input += synapse0xa00d888();
   input += synapse0xa00d8b0();
   input += synapse0xa00d8d8();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00d900() {
   double input = -0.609321;
   input += synapse0xa00da70();
   input += synapse0xa00da98();
   input += synapse0xa00dac0();
   input += synapse0xa00dae8();
   input += synapse0xa00db10();
   input += synapse0xa00db38();
   input += synapse0xa00db60();
   input += synapse0xa00db88();
   input += synapse0xa00dbb0();
   input += synapse0xa00dbd8();
   input += synapse0xa00dc00();
   input += synapse0xa00dc28();
   input += synapse0xa00dc50();
   input += synapse0xa00dc78();
   input += synapse0xa00dca0();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00dcc8() {
   double input = -0.386222;
   input += synapse0xa00de38();
   input += synapse0xa00de60();
   input += synapse0xa00de88();
   input += synapse0xa00deb0();
   input += synapse0xa00ded8();
   input += synapse0xa00df00();
   input += synapse0xa00df28();
   input += synapse0xa00df50();
   input += synapse0xa00df78();
   input += synapse0xa00dfa0();
   input += synapse0xa00dfc8();
   input += synapse0xa00dff0();
   input += synapse0xa00e018();
   input += synapse0xa00e040();
   input += synapse0xa00e068();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00e090() {
   double input = -0.473866;
   input += synapse0xa00e200();
   input += synapse0xa00e228();
   input += synapse0xa00e250();
   input += synapse0xa00a560();
   input += synapse0xa00a588();
   input += synapse0xa00a5b0();
   input += synapse0xa00a5d8();
   input += synapse0xa00a600();
   input += synapse0xa00a628();
   input += synapse0xa00a650();
   input += synapse0xa00a678();
   input += synapse0xa00a6a0();
   input += synapse0xa00a6c8();
   input += synapse0xa00a6f0();
   input += synapse0xa00a718();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00a740() {
   double input = -0.103048;
   input += synapse0xa00a8b0();
   input += synapse0xa00a8d8();
   input += synapse0xa00a900();
   input += synapse0xa00a928();
   input += synapse0xa00a950();
   input += synapse0xa00a978();
   input += synapse0xa00a9a0();
   input += synapse0xa00a9c8();
   input += synapse0xa00a9f0();
   input += synapse0xa00aa18();
   input += synapse0xa00aa40();
   input += synapse0xa00aa68();
   input += synapse0xa00aa90();
   input += synapse0xa00aab8();
   input += synapse0xa00aae0();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00ab08() {
   double input = 0.472819;
   input += synapse0xa00ac78();
   input += synapse0xa00aca0();
   input += synapse0xa00acc8();
   input += synapse0xa00acf0();
   input += synapse0xa00ad18();
   input += synapse0xa00ad40();
   input += synapse0xa00f280();
   input += synapse0xa00f2a8();
   input += synapse0xa00f2d0();
   input += synapse0xa00f2f8();
   input += synapse0xa00f320();
   input += synapse0xa00f348();
   input += synapse0xa00f370();
   input += synapse0xa00f398();
   input += synapse0xa00f3c0();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00f3e8() {
   double input = 0.39081;
   input += synapse0xa00f558();
   input += synapse0xa00f580();
   input += synapse0xa00f5a8();
   input += synapse0xa00f5d0();
   input += synapse0xa00f5f8();
   input += synapse0xa00f620();
   input += synapse0xa00f648();
   input += synapse0xa00f670();
   input += synapse0xa00f698();
   input += synapse0xa00f6c0();
   input += synapse0xa00f6e8();
   input += synapse0xa00f710();
   input += synapse0xa00f738();
   input += synapse0xa00f760();
   input += synapse0xa00f788();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00f7b0() {
   double input = -0.244167;
   input += synapse0xa00f920();
   input += synapse0xa00f948();
   input += synapse0xa00f970();
   input += synapse0xa00f998();
   input += synapse0xa00f9c0();
   input += synapse0xa00f9e8();
   input += synapse0xa00fa10();
   input += synapse0xa00fa38();
   input += synapse0xa00fa60();
   input += synapse0xa00fa88();
   input += synapse0xa00fab0();
   input += synapse0xa00fad8();
   input += synapse0xa00fb00();
   input += synapse0xa00fb28();
   input += synapse0xa00fb50();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00fb78() {
   double input = 0.314979;
   input += synapse0xa00fce8();
   input += synapse0xa00fd10();
   input += synapse0xa00fd38();
   input += synapse0xa00fd60();
   input += synapse0xa00fd88();
   input += synapse0xa00fdb0();
   input += synapse0xa00fdd8();
   input += synapse0xa00fe00();
   input += synapse0xa00fe28();
   input += synapse0xa00fe50();
   input += synapse0xa00fe78();
   input += synapse0xa00fea0();
   input += synapse0xa00fec8();
   input += synapse0xa00fef0();
   input += synapse0xa00ff18();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa00ff40() {
   double input = 0.446355;
   input += synapse0xa0100b0();
   input += synapse0xa0100d8();
   input += synapse0xa010100();
   input += synapse0xa010128();
   input += synapse0xa010150();
   input += synapse0xa010178();
   input += synapse0xa0101a0();
   input += synapse0xa0101c8();
   input += synapse0xa0101f0();
   input += synapse0xa010218();
   input += synapse0xa010240();
   input += synapse0xa010268();
   input += synapse0xa010290();
   input += synapse0xa0102b8();
   input += synapse0xa0102e0();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa010308() {
   double input = 0.276952;
   input += synapse0xa010478();
   input += synapse0xa0104a0();
   input += synapse0xa0104c8();
   input += synapse0xa0104f0();
   input += synapse0xa010518();
   input += synapse0xa010540();
   input += synapse0xa010568();
   input += synapse0xa010590();
   input += synapse0xa0105b8();
   input += synapse0xa0105e0();
   input += synapse0xa010608();
   input += synapse0xa010630();
   input += synapse0xa010658();
   input += synapse0xa010680();
   input += synapse0xa0106a8();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa0106d0() {
   double input = 0.10218;
   input += synapse0xa010840();
   input += synapse0xa010868();
   input += synapse0xa010890();
   input += synapse0xa0108b8();
   input += synapse0xa0108e0();
   input += synapse0xa010908();
   input += synapse0xa010930();
   input += synapse0xa010958();
   input += synapse0xa010980();
   input += synapse0xa0109a8();
   input += synapse0xa0109d0();
   input += synapse0xa0109f8();
   input += synapse0xa010a20();
   input += synapse0xa010a48();
   input += synapse0xa010a70();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa010a98() {
   double input = -0.399362;
   input += synapse0xa010c08();
   input += synapse0xa010c30();
   input += synapse0xa010c58();
   input += synapse0xa010c80();
   input += synapse0xa010ca8();
   input += synapse0xa010cd0();
   input += synapse0xa010cf8();
   input += synapse0xa010d20();
   input += synapse0xa010d48();
   input += synapse0xa010d70();
   input += synapse0xa010d98();
   input += synapse0xa010dc0();
   input += synapse0xa010de8();
   input += synapse0xa010e10();
   input += synapse0xa010e38();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa010e60() {
   double input = -0.550638;
   input += synapse0xa010fd0();
   input += synapse0xa010ff8();
   input += synapse0xa011020();
   input += synapse0xa011048();
   input += synapse0xa011070();
   input += synapse0xa011098();
   input += synapse0xa0110c0();
   input += synapse0xa0110e8();
   input += synapse0xa011110();
   input += synapse0xa011138();
   input += synapse0xa011160();
   input += synapse0xa011188();
   input += synapse0xa0111b0();
   input += synapse0xa0111d8();
   input += synapse0xa011200();
   return ((1/(1+exp(-input)))*1)+0;
}

double AnnAna::neuron0xa0068f8() {
   double input = -0.0428478;
   input += synapse0xa0112b8();
   input += synapse0xa0112e0();
   input += synapse0xa011308();
   input += synapse0xa011330();
   input += synapse0xa011358();
   input += synapse0xa011380();
   input += synapse0xa0113a8();
   input += synapse0xa0113d0();
   input += synapse0xa0113f8();
   input += synapse0xa011420();
   input += synapse0xa011448();
   input += synapse0xa011470();
   input += synapse0xa011498();
   input += synapse0xa0114c0();
   input += synapse0xa0114e8();
   return input;
}

double AnnAna::synapse0xa006b40() {
   return (neuron0xa0253a8()*0.436848);
}

double AnnAna::synapse0xa006b68() {
   return (neuron0xa025510()*0.774105);
}

double AnnAna::synapse0xa006b90() {
   return (neuron0xa025678()*-0.401912);
}

double AnnAna::synapse0xa006bb8() {
   return (neuron0xa0257e0()*-0.950019);
}

double AnnAna::synapse0xa006be0() {
   return (neuron0xa025948()*-0.185474);
}

double AnnAna::synapse0xa006c08() {
   return (neuron0xa025ab8()*0.29392);
}

double AnnAna::synapse0xa006c30() {
   return (neuron0xa025c20()*0.227605);
}

double AnnAna::synapse0xa006c58() {
   return (neuron0xa025d88()*-0.105657);
}

double AnnAna::synapse0xa006c80() {
   return (neuron0xa025ef0()*-0.334641);
}

double AnnAna::synapse0xa006ca8() {
   return (neuron0xa026058()*-0.545888);
}

double AnnAna::synapse0xa006cd0() {
   return (neuron0xa0261c0()*0.00351596);
}

double AnnAna::synapse0xa006cf8() {
   return (neuron0xa026328()*0.602204);
}

double AnnAna::synapse0xa006d20() {
   return (neuron0xa026490()*-0.10411);
}

double AnnAna::synapse0xa006d48() {
   return (neuron0xa0265f8()*0.0760779);
}

double AnnAna::synapse0xa006d70() {
   return (neuron0xa026760()*0.131848);
}

double AnnAna::synapse0xa006d98() {
   return (neuron0xa0268c8()*0.302881);
}

double AnnAna::synapse0xa006dc0() {
   return (neuron0xa026a30()*0.48811);
}

double AnnAna::synapse0xa006ef8() {
   return (neuron0xa026ca8()*0.568274);
}

double AnnAna::synapse0xa006f20() {
   return (neuron0xa026d80()*0.563196);
}

double AnnAna::synapse0xa006f48() {
   return (neuron0xa026ee8()*0.228305);
}

double AnnAna::synapse0xa006f70() {
   return (neuron0xa027050()*-0.574875);
}

double AnnAna::synapse0xa006f98() {
   return (neuron0xa0271b8()*-0.728777);
}

double AnnAna::synapse0xa006fc0() {
   return (neuron0xa027320()*-1.00707);
}

double AnnAna::synapse0xa006fe8() {
   return (neuron0xa027488()*-1.58895);
}

double AnnAna::synapse0xa007010() {
   return (neuron0xa0275f0()*0.325033);
}

double AnnAna::synapse0xa007038() {
   return (neuron0xa027758()*-0.0570835);
}

double AnnAna::synapse0xa007060() {
   return (neuron0xa0278c0()*-0.414386);
}

double AnnAna::synapse0xa007088() {
   return (neuron0xa0064a0()*-0.220013);
}

double AnnAna::synapse0xa0070b0() {
   return (neuron0xa006608()*1.02402);
}

double AnnAna::synapse0xa0070d8() {
   return (neuron0xa006770()*0.634818);
}

double AnnAna::synapse0xa0071e0() {
   return (neuron0xa0253a8()*-0.735289);
}

double AnnAna::synapse0xa007208() {
   return (neuron0xa025510()*-0.256214);
}

double AnnAna::synapse0xa007230() {
   return (neuron0xa025678()*0.327874);
}

double AnnAna::synapse0xa0279e0() {
   return (neuron0xa0257e0()*0.113329);
}

double AnnAna::synapse0xa027a08() {
   return (neuron0xa025948()*-0.000581991);
}

double AnnAna::synapse0xa0315d0() {
   return (neuron0xa025ab8()*0.146604);
}

double AnnAna::synapse0xa0315f8() {
   return (neuron0xa025c20()*0.380137);
}

double AnnAna::synapse0xa006e70() {
   return (neuron0xa025d88()*-0.70031);
}

double AnnAna::synapse0xa006e98() {
   return (neuron0xa025ef0()*-0.210992);
}

double AnnAna::synapse0xa006ec0() {
   return (neuron0xa026058()*-0.822514);
}

double AnnAna::synapse0xa007360() {
   return (neuron0xa0261c0()*0.393659);
}

double AnnAna::synapse0xa007388() {
   return (neuron0xa026328()*0.0910783);
}

double AnnAna::synapse0xa0073b0() {
   return (neuron0xa026490()*0.256777);
}

double AnnAna::synapse0xa0073d8() {
   return (neuron0xa0265f8()*0.521723);
}

double AnnAna::synapse0xa007400() {
   return (neuron0xa026760()*-0.00513969);
}

double AnnAna::synapse0xa007428() {
   return (neuron0xa0268c8()*0.222863);
}

double AnnAna::synapse0xa007450() {
   return (neuron0xa026a30()*0.1745);
}

double AnnAna::synapse0xa007500() {
   return (neuron0xa026ca8()*0.575282);
}

double AnnAna::synapse0xa007528() {
   return (neuron0xa026d80()*-0.653654);
}

double AnnAna::synapse0xa007550() {
   return (neuron0xa026ee8()*-0.391515);
}

double AnnAna::synapse0xa007578() {
   return (neuron0xa027050()*0.325741);
}

double AnnAna::synapse0xa0075a0() {
   return (neuron0xa0271b8()*-0.439108);
}

double AnnAna::synapse0xa0075c8() {
   return (neuron0xa027320()*-1.48929);
}

double AnnAna::synapse0xa0075f0() {
   return (neuron0xa027488()*-1.94498);
}

double AnnAna::synapse0xa007618() {
   return (neuron0xa0275f0()*0.160284);
}

double AnnAna::synapse0xa007640() {
   return (neuron0xa027758()*0.368076);
}

double AnnAna::synapse0xa007668() {
   return (neuron0xa0278c0()*-1.20931);
}

double AnnAna::synapse0xa007690() {
   return (neuron0xa0064a0()*-0.965336);
}

double AnnAna::synapse0xa0076b8() {
   return (neuron0xa006608()*1.00126);
}

double AnnAna::synapse0xa0076e0() {
   return (neuron0xa006770()*0.910948);
}

double AnnAna::synapse0xa026c20() {
   return (neuron0xa0253a8()*0.260956);
}

double AnnAna::synapse0xa026c48() {
   return (neuron0xa025510()*-0.546665);
}

double AnnAna::synapse0xa026c70() {
   return (neuron0xa025678()*0.526897);
}

double AnnAna::synapse0xa007938() {
   return (neuron0xa0257e0()*-0.403909);
}

double AnnAna::synapse0xa007960() {
   return (neuron0xa025948()*0.0840154);
}

double AnnAna::synapse0xa007258() {
   return (neuron0xa025ab8()*-0.181505);
}

double AnnAna::synapse0xa007280() {
   return (neuron0xa025c20()*-0.616972);
}

double AnnAna::synapse0xa0072a8() {
   return (neuron0xa025d88()*-1.15893);
}

double AnnAna::synapse0xa0072d0() {
   return (neuron0xa025ef0()*-0.169514);
}

double AnnAna::synapse0xa0072f8() {
   return (neuron0xa026058()*-0.0655928);
}

double AnnAna::synapse0xa007320() {
   return (neuron0xa0261c0()*-0.0486667);
}

double AnnAna::synapse0xa007b90() {
   return (neuron0xa026328()*0.0543502);
}

double AnnAna::synapse0xa007bb8() {
   return (neuron0xa026490()*0.0244697);
}

double AnnAna::synapse0xa007be0() {
   return (neuron0xa0265f8()*0.249698);
}

double AnnAna::synapse0xa007c08() {
   return (neuron0xa026760()*-0.179008);
}

double AnnAna::synapse0xa007c30() {
   return (neuron0xa0268c8()*0.445441);
}

double AnnAna::synapse0xa007c58() {
   return (neuron0xa026a30()*0.214446);
}

double AnnAna::synapse0xa007d08() {
   return (neuron0xa026ca8()*-0.253676);
}

double AnnAna::synapse0xa007d30() {
   return (neuron0xa026d80()*0.0694983);
}

double AnnAna::synapse0xa007d58() {
   return (neuron0xa026ee8()*0.517132);
}

double AnnAna::synapse0xa007d80() {
   return (neuron0xa027050()*0.244857);
}

double AnnAna::synapse0xa007da8() {
   return (neuron0xa0271b8()*0.581814);
}

double AnnAna::synapse0xa007dd0() {
   return (neuron0xa027320()*0.473134);
}

double AnnAna::synapse0xa007df8() {
   return (neuron0xa027488()*-0.401601);
}

double AnnAna::synapse0xa007e20() {
   return (neuron0xa0275f0()*0.259251);
}

double AnnAna::synapse0xa007e48() {
   return (neuron0xa027758()*-0.291009);
}

double AnnAna::synapse0xa007e70() {
   return (neuron0xa0278c0()*0.160222);
}

double AnnAna::synapse0xa007e98() {
   return (neuron0xa0064a0()*0.386337);
}

double AnnAna::synapse0xa007ec0() {
   return (neuron0xa006608()*0.192261);
}

double AnnAna::synapse0xa007ee8() {
   return (neuron0xa006770()*0.0680474);
}

double AnnAna::synapse0xa008038() {
   return (neuron0xa0253a8()*-0.134303);
}

double AnnAna::synapse0xa008060() {
   return (neuron0xa025510()*0.107171);
}

double AnnAna::synapse0xa008088() {
   return (neuron0xa025678()*-0.398607);
}

double AnnAna::synapse0xa0080b0() {
   return (neuron0xa0257e0()*0.175908);
}

double AnnAna::synapse0xa0080d8() {
   return (neuron0xa025948()*-0.00423978);
}

double AnnAna::synapse0xa008100() {
   return (neuron0xa025ab8()*0.0392909);
}

double AnnAna::synapse0xa008128() {
   return (neuron0xa025c20()*0.303689);
}

double AnnAna::synapse0xa008150() {
   return (neuron0xa025d88()*0.0661471);
}

double AnnAna::synapse0xa008178() {
   return (neuron0xa025ef0()*-0.261519);
}

double AnnAna::synapse0xa0081a0() {
   return (neuron0xa026058()*-0.349336);
}

double AnnAna::synapse0xa0081c8() {
   return (neuron0xa0261c0()*0.66851);
}

double AnnAna::synapse0xa0081f0() {
   return (neuron0xa026328()*0.794694);
}

double AnnAna::synapse0xa008218() {
   return (neuron0xa026490()*0.388376);
}

double AnnAna::synapse0xa008240() {
   return (neuron0xa0265f8()*0.0197791);
}

double AnnAna::synapse0xa008268() {
   return (neuron0xa026760()*-0.133443);
}

double AnnAna::synapse0xa008290() {
   return (neuron0xa0268c8()*0.240394);
}

double AnnAna::synapse0xa0082b8() {
   return (neuron0xa026a30()*0.447989);
}

double AnnAna::synapse0xa008368() {
   return (neuron0xa026ca8()*0.707603);
}

double AnnAna::synapse0xa008390() {
   return (neuron0xa026d80()*0.348355);
}

double AnnAna::synapse0xa0083b8() {
   return (neuron0xa026ee8()*0.0764537);
}

double AnnAna::synapse0xa0083e0() {
   return (neuron0xa027050()*0.17472);
}

double AnnAna::synapse0xa008408() {
   return (neuron0xa0271b8()*0.655574);
}

double AnnAna::synapse0xa008430() {
   return (neuron0xa027320()*-0.348972);
}

double AnnAna::synapse0xa008458() {
   return (neuron0xa027488()*-0.897831);
}

double AnnAna::synapse0xa008480() {
   return (neuron0xa0275f0()*-0.0809154);
}

double AnnAna::synapse0xa0084a8() {
   return (neuron0xa027758()*-0.121484);
}

double AnnAna::synapse0xa0084d0() {
   return (neuron0xa0278c0()*0.291519);
}

double AnnAna::synapse0xa0084f8() {
   return (neuron0xa0064a0()*-0.51407);
}

double AnnAna::synapse0xa008520() {
   return (neuron0xa006608()*0.363543);
}

double AnnAna::synapse0xa008548() {
   return (neuron0xa006770()*0.537082);
}

double AnnAna::synapse0xa008698() {
   return (neuron0xa0253a8()*-0.671042);
}

double AnnAna::synapse0xa0086c0() {
   return (neuron0xa025510()*-0.590131);
}

double AnnAna::synapse0xa0086e8() {
   return (neuron0xa025678()*-0.496622);
}

double AnnAna::synapse0xa008710() {
   return (neuron0xa0257e0()*-0.355808);
}

double AnnAna::synapse0xa008738() {
   return (neuron0xa025948()*0.0679145);
}

double AnnAna::synapse0xa008760() {
   return (neuron0xa025ab8()*-0.296694);
}

double AnnAna::synapse0xa008788() {
   return (neuron0xa025c20()*0.216053);
}

double AnnAna::synapse0xa0087b0() {
   return (neuron0xa025d88()*0.273578);
}

double AnnAna::synapse0xa0087d8() {
   return (neuron0xa025ef0()*-0.130873);
}

double AnnAna::synapse0xa007988() {
   return (neuron0xa026058()*0.224404);
}

double AnnAna::synapse0xa0079b0() {
   return (neuron0xa0261c0()*-0.110817);
}

double AnnAna::synapse0xa0079d8() {
   return (neuron0xa026328()*-0.120695);
}

double AnnAna::synapse0xa007a00() {
   return (neuron0xa026490()*0.126341);
}

double AnnAna::synapse0xa007a28() {
   return (neuron0xa0265f8()*0.0950813);
}

double AnnAna::synapse0xa007a50() {
   return (neuron0xa026760()*0.313388);
}

double AnnAna::synapse0xa007a78() {
   return (neuron0xa0268c8()*-0.235995);
}

double AnnAna::synapse0xa007aa0() {
   return (neuron0xa026a30()*-0.298999);
}

double AnnAna::synapse0xa007b50() {
   return (neuron0xa026ca8()*-0.512071);
}

double AnnAna::synapse0xa008c08() {
   return (neuron0xa026d80()*0.0958581);
}

double AnnAna::synapse0xa008c30() {
   return (neuron0xa026ee8()*0.0440837);
}

double AnnAna::synapse0xa008c58() {
   return (neuron0xa027050()*0.258726);
}

double AnnAna::synapse0xa008c80() {
   return (neuron0xa0271b8()*0.0441915);
}

double AnnAna::synapse0xa008ca8() {
   return (neuron0xa027320()*-0.220244);
}

double AnnAna::synapse0xa008cd0() {
   return (neuron0xa027488()*-0.111962);
}

double AnnAna::synapse0xa008cf8() {
   return (neuron0xa0275f0()*0.069718);
}

double AnnAna::synapse0xa008d20() {
   return (neuron0xa027758()*-0.621424);
}

double AnnAna::synapse0xa008d48() {
   return (neuron0xa0278c0()*-0.256212);
}

double AnnAna::synapse0xa008d70() {
   return (neuron0xa0064a0()*-0.198411);
}

double AnnAna::synapse0xa008d98() {
   return (neuron0xa006608()*-0.183372);
}

double AnnAna::synapse0xa008dc0() {
   return (neuron0xa006770()*0.483285);
}

double AnnAna::synapse0xa008f10() {
   return (neuron0xa0253a8()*-0.32979);
}

double AnnAna::synapse0xa008f38() {
   return (neuron0xa025510()*-0.243305);
}

double AnnAna::synapse0xa008f60() {
   return (neuron0xa025678()*-0.19613);
}

double AnnAna::synapse0xa008f88() {
   return (neuron0xa0257e0()*0.699315);
}

double AnnAna::synapse0xa008fb0() {
   return (neuron0xa025948()*-0.238795);
}

double AnnAna::synapse0xa008fd8() {
   return (neuron0xa025ab8()*-0.495275);
}

double AnnAna::synapse0xa009000() {
   return (neuron0xa025c20()*0.343178);
}

double AnnAna::synapse0xa009028() {
   return (neuron0xa025d88()*1.969);
}

double AnnAna::synapse0xa009050() {
   return (neuron0xa025ef0()*0.300509);
}

double AnnAna::synapse0xa009078() {
   return (neuron0xa026058()*0.350022);
}

double AnnAna::synapse0xa0090a0() {
   return (neuron0xa0261c0()*-0.0383063);
}

double AnnAna::synapse0xa0090c8() {
   return (neuron0xa026328()*-0.353061);
}

double AnnAna::synapse0xa0090f0() {
   return (neuron0xa026490()*0.00480291);
}

double AnnAna::synapse0xa009118() {
   return (neuron0xa0265f8()*-0.176569);
}

double AnnAna::synapse0xa009140() {
   return (neuron0xa026760()*0.781756);
}

double AnnAna::synapse0xa009168() {
   return (neuron0xa0268c8()*-0.258833);
}

double AnnAna::synapse0xa009190() {
   return (neuron0xa026a30()*-0.624973);
}

double AnnAna::synapse0xa009240() {
   return (neuron0xa026ca8()*-1.14503);
}

double AnnAna::synapse0xa009268() {
   return (neuron0xa026d80()*0.766695);
}

double AnnAna::synapse0xa009290() {
   return (neuron0xa026ee8()*1.62363);
}

double AnnAna::synapse0xa0092b8() {
   return (neuron0xa027050()*-1.07598);
}

double AnnAna::synapse0xa0092e0() {
   return (neuron0xa0271b8()*0.366591);
}

double AnnAna::synapse0xa009308() {
   return (neuron0xa027320()*0.538998);
}

double AnnAna::synapse0xa009330() {
   return (neuron0xa027488()*1.65703);
}

double AnnAna::synapse0xa009358() {
   return (neuron0xa0275f0()*-0.12476);
}

double AnnAna::synapse0xa009380() {
   return (neuron0xa027758()*-0.140061);
}

double AnnAna::synapse0xa0093a8() {
   return (neuron0xa0278c0()*0.903589);
}

double AnnAna::synapse0xa0093d0() {
   return (neuron0xa0064a0()*0.687103);
}

double AnnAna::synapse0xa0093f8() {
   return (neuron0xa006608()*-0.236736);
}

double AnnAna::synapse0xa009420() {
   return (neuron0xa006770()*-0.887707);
}

double AnnAna::synapse0xa009570() {
   return (neuron0xa0253a8()*0.342888);
}

double AnnAna::synapse0xa009598() {
   return (neuron0xa025510()*-0.15214);
}

double AnnAna::synapse0xa0095c0() {
   return (neuron0xa025678()*-0.140502);
}

double AnnAna::synapse0xa0095e8() {
   return (neuron0xa0257e0()*0.144491);
}

double AnnAna::synapse0xa009610() {
   return (neuron0xa025948()*-0.250981);
}

double AnnAna::synapse0xa009638() {
   return (neuron0xa025ab8()*-0.00509901);
}

double AnnAna::synapse0xa009660() {
   return (neuron0xa025c20()*0.0395716);
}

double AnnAna::synapse0xa009688() {
   return (neuron0xa025d88()*0.508624);
}

double AnnAna::synapse0xa0096b0() {
   return (neuron0xa025ef0()*0.163077);
}

double AnnAna::synapse0xa0096d8() {
   return (neuron0xa026058()*-0.648062);
}

double AnnAna::synapse0xa009700() {
   return (neuron0xa0261c0()*0.275792);
}

double AnnAna::synapse0xa009728() {
   return (neuron0xa026328()*0.172029);
}

double AnnAna::synapse0xa009750() {
   return (neuron0xa026490()*0.0777463);
}

double AnnAna::synapse0xa009778() {
   return (neuron0xa0265f8()*-0.0863003);
}

double AnnAna::synapse0xa0097a0() {
   return (neuron0xa026760()*0.180699);
}

double AnnAna::synapse0xa0097c8() {
   return (neuron0xa0268c8()*-0.165746);
}

double AnnAna::synapse0xa0097f0() {
   return (neuron0xa026a30()*0.257052);
}

double AnnAna::synapse0xa0098a0() {
   return (neuron0xa026ca8()*-0.439523);
}

double AnnAna::synapse0xa0098c8() {
   return (neuron0xa026d80()*0.111393);
}

double AnnAna::synapse0xa0098f0() {
   return (neuron0xa026ee8()*0.337367);
}

double AnnAna::synapse0xa009918() {
   return (neuron0xa027050()*0.326849);
}

double AnnAna::synapse0xa009940() {
   return (neuron0xa0271b8()*0.115823);
}

double AnnAna::synapse0xa009968() {
   return (neuron0xa027320()*-0.100338);
}

double AnnAna::synapse0xa009990() {
   return (neuron0xa027488()*0.333152);
}

double AnnAna::synapse0xa0099b8() {
   return (neuron0xa0275f0()*0.302921);
}

double AnnAna::synapse0xa0099e0() {
   return (neuron0xa027758()*-0.182395);
}

double AnnAna::synapse0xa009a08() {
   return (neuron0xa0278c0()*0.247547);
}

double AnnAna::synapse0xa009a30() {
   return (neuron0xa0064a0()*-0.101442);
}

double AnnAna::synapse0xa009a58() {
   return (neuron0xa006608()*0.172986);
}

double AnnAna::synapse0xa009a80() {
   return (neuron0xa006770()*-0.0835321);
}

double AnnAna::synapse0xa009bd0() {
   return (neuron0xa0253a8()*0.271695);
}

double AnnAna::synapse0xa009bf8() {
   return (neuron0xa025510()*-0.15443);
}

double AnnAna::synapse0xa009c20() {
   return (neuron0xa025678()*0.147905);
}

double AnnAna::synapse0xa009c48() {
   return (neuron0xa0257e0()*-0.106713);
}

double AnnAna::synapse0xa009c70() {
   return (neuron0xa025948()*-0.230814);
}

double AnnAna::synapse0xa009c98() {
   return (neuron0xa025ab8()*0.0309136);
}

double AnnAna::synapse0xa009cc0() {
   return (neuron0xa025c20()*-0.8201);
}

double AnnAna::synapse0xa009ce8() {
   return (neuron0xa025d88()*-0.390933);
}

double AnnAna::synapse0xa009d10() {
   return (neuron0xa025ef0()*0.0924757);
}

double AnnAna::synapse0xa009d38() {
   return (neuron0xa026058()*-0.3874);
}

double AnnAna::synapse0xa009d60() {
   return (neuron0xa0261c0()*0.384426);
}

double AnnAna::synapse0xa009d88() {
   return (neuron0xa026328()*0.169974);
}

double AnnAna::synapse0xa009db0() {
   return (neuron0xa026490()*-0.199647);
}

double AnnAna::synapse0xa009dd8() {
   return (neuron0xa0265f8()*0.23354);
}

double AnnAna::synapse0xa009e00() {
   return (neuron0xa026760()*-0.284328);
}

double AnnAna::synapse0xa009e28() {
   return (neuron0xa0268c8()*0.211814);
}

double AnnAna::synapse0xa009e50() {
   return (neuron0xa026a30()*0.4718);
}

double AnnAna::synapse0xa009f00() {
   return (neuron0xa026ca8()*0.0913299);
}

double AnnAna::synapse0xa009f28() {
   return (neuron0xa026d80()*-0.360952);
}

double AnnAna::synapse0xa009f50() {
   return (neuron0xa026ee8()*-0.169745);
}

double AnnAna::synapse0xa009f78() {
   return (neuron0xa027050()*0.897452);
}

double AnnAna::synapse0xa009fa0() {
   return (neuron0xa0271b8()*0.00719067);
}

double AnnAna::synapse0xa009fc8() {
   return (neuron0xa027320()*0.0390898);
}

double AnnAna::synapse0xa009ff0() {
   return (neuron0xa027488()*-0.239637);
}

double AnnAna::synapse0xa00a018() {
   return (neuron0xa0275f0()*-0.00529673);
}

double AnnAna::synapse0xa00a040() {
   return (neuron0xa027758()*0.581812);
}

double AnnAna::synapse0xa00a068() {
   return (neuron0xa0278c0()*0.171289);
}

double AnnAna::synapse0xa00a090() {
   return (neuron0xa0064a0()*0.252467);
}

double AnnAna::synapse0xa00a0b8() {
   return (neuron0xa006608()*-0.1713);
}

double AnnAna::synapse0xa00a0e0() {
   return (neuron0xa006770()*-0.112319);
}

double AnnAna::synapse0xa00a230() {
   return (neuron0xa0253a8()*-0.130216);
}

double AnnAna::synapse0xa00a258() {
   return (neuron0xa025510()*-0.432291);
}

double AnnAna::synapse0xa00a280() {
   return (neuron0xa025678()*0.650963);
}

double AnnAna::synapse0xa00a2a8() {
   return (neuron0xa0257e0()*0.0910798);
}

double AnnAna::synapse0xa00a2d0() {
   return (neuron0xa025948()*-0.449117);
}

double AnnAna::synapse0xa00a2f8() {
   return (neuron0xa025ab8()*0.110845);
}

double AnnAna::synapse0xa00a320() {
   return (neuron0xa025c20()*-0.0629093);
}

double AnnAna::synapse0xa00a348() {
   return (neuron0xa025d88()*-0.178477);
}

double AnnAna::synapse0xa00a370() {
   return (neuron0xa025ef0()*-0.0833317);
}

double AnnAna::synapse0xa00a398() {
   return (neuron0xa026058()*-0.278471);
}

double AnnAna::synapse0xa00a3c0() {
   return (neuron0xa0261c0()*0.615223);
}

double AnnAna::synapse0xa00a3e8() {
   return (neuron0xa026328()*-0.0790409);
}

double AnnAna::synapse0xa00a410() {
   return (neuron0xa026490()*0.432635);
}

double AnnAna::synapse0xa00a438() {
   return (neuron0xa0265f8()*-0.231539);
}

double AnnAna::synapse0xa00a460() {
   return (neuron0xa026760()*-0.0148876);
}

double AnnAna::synapse0xa00a488() {
   return (neuron0xa0268c8()*0.263407);
}

double AnnAna::synapse0xa00a4b0() {
   return (neuron0xa026a30()*0.134533);
}

double AnnAna::synapse0xa024cf8() {
   return (neuron0xa026ca8()*0.458897);
}

double AnnAna::synapse0xa008800() {
   return (neuron0xa026d80()*-0.282452);
}

double AnnAna::synapse0xa008828() {
   return (neuron0xa026ee8()*0.147223);
}

double AnnAna::synapse0xa008850() {
   return (neuron0xa027050()*0.751823);
}

double AnnAna::synapse0xa008878() {
   return (neuron0xa0271b8()*0.682027);
}

double AnnAna::synapse0xa0088a0() {
   return (neuron0xa027320()*-0.224508);
}

double AnnAna::synapse0xa0088c8() {
   return (neuron0xa027488()*-0.339937);
}

double AnnAna::synapse0xa0088f0() {
   return (neuron0xa0275f0()*0.522873);
}

double AnnAna::synapse0xa008918() {
   return (neuron0xa027758()*0.612153);
}

double AnnAna::synapse0xa008940() {
   return (neuron0xa0278c0()*-0.114982);
}

double AnnAna::synapse0xa008968() {
   return (neuron0xa0064a0()*-0.318636);
}

double AnnAna::synapse0xa008990() {
   return (neuron0xa006608()*0.18613);
}

double AnnAna::synapse0xa0089b8() {
   return (neuron0xa006770()*0.401373);
}

double AnnAna::synapse0xa008b30() {
   return (neuron0xa0253a8()*-0.0293046);
}

double AnnAna::synapse0xa008b58() {
   return (neuron0xa025510()*0.987525);
}

double AnnAna::synapse0xa008b80() {
   return (neuron0xa025678()*0.167438);
}

double AnnAna::synapse0xa008ba8() {
   return (neuron0xa0257e0()*0.727673);
}

double AnnAna::synapse0xa008bd0() {
   return (neuron0xa025948()*0.477319);
}

double AnnAna::synapse0xa00ad68() {
   return (neuron0xa025ab8()*0.0127759);
}

double AnnAna::synapse0xa00ad90() {
   return (neuron0xa025c20()*0.916618);
}

double AnnAna::synapse0xa00adb8() {
   return (neuron0xa025d88()*0.720388);
}

double AnnAna::synapse0xa00ade0() {
   return (neuron0xa025ef0()*0.349514);
}

double AnnAna::synapse0xa00ae08() {
   return (neuron0xa026058()*-0.193762);
}

double AnnAna::synapse0xa00ae30() {
   return (neuron0xa0261c0()*-0.220174);
}

double AnnAna::synapse0xa00ae58() {
   return (neuron0xa026328()*-0.652297);
}

double AnnAna::synapse0xa00ae80() {
   return (neuron0xa026490()*-0.226863);
}

double AnnAna::synapse0xa00aea8() {
   return (neuron0xa0265f8()*0.0820539);
}

double AnnAna::synapse0xa00aed0() {
   return (neuron0xa026760()*0.603436);
}

double AnnAna::synapse0xa00aef8() {
   return (neuron0xa0268c8()*-0.36843);
}

double AnnAna::synapse0xa00af20() {
   return (neuron0xa026a30()*-0.299138);
}

double AnnAna::synapse0xa00afd0() {
   return (neuron0xa026ca8()*-0.156747);
}

double AnnAna::synapse0xa00aff8() {
   return (neuron0xa026d80()*1.1197);
}

double AnnAna::synapse0xa00b020() {
   return (neuron0xa026ee8()*0.35503);
}

double AnnAna::synapse0xa00b048() {
   return (neuron0xa027050()*-0.395081);
}

double AnnAna::synapse0xa00b070() {
   return (neuron0xa0271b8()*0.0385539);
}

double AnnAna::synapse0xa00b098() {
   return (neuron0xa027320()*1.00859);
}

double AnnAna::synapse0xa00b0c0() {
   return (neuron0xa027488()*1.40315);
}

double AnnAna::synapse0xa00b0e8() {
   return (neuron0xa0275f0()*-0.150149);
}

double AnnAna::synapse0xa00b110() {
   return (neuron0xa027758()*-0.347005);
}

double AnnAna::synapse0xa00b138() {
   return (neuron0xa0278c0()*0.527464);
}

double AnnAna::synapse0xa00b160() {
   return (neuron0xa0064a0()*0.597);
}

double AnnAna::synapse0xa00b188() {
   return (neuron0xa006608()*-0.901939);
}

double AnnAna::synapse0xa00b1b0() {
   return (neuron0xa006770()*-0.041256);
}

double AnnAna::synapse0xa00b300() {
   return (neuron0xa0253a8()*0.444164);
}

double AnnAna::synapse0xa00b328() {
   return (neuron0xa025510()*-0.0558126);
}

double AnnAna::synapse0xa00b350() {
   return (neuron0xa025678()*-0.171957);
}

double AnnAna::synapse0xa00b378() {
   return (neuron0xa0257e0()*0.0139001);
}

double AnnAna::synapse0xa00b3a0() {
   return (neuron0xa025948()*0.113559);
}

double AnnAna::synapse0xa00b3c8() {
   return (neuron0xa025ab8()*-0.226963);
}

double AnnAna::synapse0xa00b3f0() {
   return (neuron0xa025c20()*-0.35174);
}

double AnnAna::synapse0xa00b418() {
   return (neuron0xa025d88()*0.285569);
}

double AnnAna::synapse0xa00b440() {
   return (neuron0xa025ef0()*-0.332591);
}

double AnnAna::synapse0xa00b468() {
   return (neuron0xa026058()*0.586186);
}

double AnnAna::synapse0xa00b490() {
   return (neuron0xa0261c0()*-0.172417);
}

double AnnAna::synapse0xa00b4b8() {
   return (neuron0xa026328()*-0.187434);
}

double AnnAna::synapse0xa00b4e0() {
   return (neuron0xa026490()*0.0681323);
}

double AnnAna::synapse0xa00b508() {
   return (neuron0xa0265f8()*-0.00406189);
}

double AnnAna::synapse0xa00b530() {
   return (neuron0xa026760()*-0.461525);
}

double AnnAna::synapse0xa00b558() {
   return (neuron0xa0268c8()*-0.34174);
}

double AnnAna::synapse0xa00b580() {
   return (neuron0xa026a30()*0.122852);
}

double AnnAna::synapse0xa00b630() {
   return (neuron0xa026ca8()*-0.700054);
}

double AnnAna::synapse0xa00b658() {
   return (neuron0xa026d80()*0.221246);
}

double AnnAna::synapse0xa00b680() {
   return (neuron0xa026ee8()*0.270049);
}

double AnnAna::synapse0xa00b6a8() {
   return (neuron0xa027050()*-0.655255);
}

double AnnAna::synapse0xa00b6d0() {
   return (neuron0xa0271b8()*-0.225891);
}

double AnnAna::synapse0xa00b6f8() {
   return (neuron0xa027320()*0.674855);
}

double AnnAna::synapse0xa00b720() {
   return (neuron0xa027488()*0.377815);
}

double AnnAna::synapse0xa00b748() {
   return (neuron0xa0275f0()*0.364611);
}

double AnnAna::synapse0xa00b770() {
   return (neuron0xa027758()*-0.201067);
}

double AnnAna::synapse0xa00b798() {
   return (neuron0xa0278c0()*0.266676);
}

double AnnAna::synapse0xa00b7c0() {
   return (neuron0xa0064a0()*0.441911);
}

double AnnAna::synapse0xa00b7e8() {
   return (neuron0xa006608()*-0.0838508);
}

double AnnAna::synapse0xa00b810() {
   return (neuron0xa006770()*-0.638238);
}

double AnnAna::synapse0xa00b960() {
   return (neuron0xa0253a8()*-0.0695893);
}

double AnnAna::synapse0xa00b988() {
   return (neuron0xa025510()*0.0217867);
}

double AnnAna::synapse0xa00b9b0() {
   return (neuron0xa025678()*0.780309);
}

double AnnAna::synapse0xa00b9d8() {
   return (neuron0xa0257e0()*-0.749993);
}

double AnnAna::synapse0xa00ba00() {
   return (neuron0xa025948()*0.326615);
}

double AnnAna::synapse0xa00ba28() {
   return (neuron0xa025ab8()*-0.310077);
}

double AnnAna::synapse0xa00ba50() {
   return (neuron0xa025c20()*-1.34546);
}

double AnnAna::synapse0xa00ba78() {
   return (neuron0xa025d88()*-1.29269);
}

double AnnAna::synapse0xa00baa0() {
   return (neuron0xa025ef0()*0.430433);
}

double AnnAna::synapse0xa00bac8() {
   return (neuron0xa026058()*-0.475169);
}

double AnnAna::synapse0xa00baf0() {
   return (neuron0xa0261c0()*0.205748);
}

double AnnAna::synapse0xa00bb18() {
   return (neuron0xa026328()*0.250818);
}

double AnnAna::synapse0xa00bb40() {
   return (neuron0xa026490()*0.35009);
}

double AnnAna::synapse0xa00bb68() {
   return (neuron0xa0265f8()*0.423148);
}

double AnnAna::synapse0xa00bb90() {
   return (neuron0xa026760()*-0.229928);
}

double AnnAna::synapse0xa00bbb8() {
   return (neuron0xa0268c8()*0.54458);
}

double AnnAna::synapse0xa00bbe0() {
   return (neuron0xa026a30()*0.487381);
}

double AnnAna::synapse0xa00bc90() {
   return (neuron0xa026ca8()*0.915117);
}

double AnnAna::synapse0xa00bcb8() {
   return (neuron0xa026d80()*-0.925384);
}

double AnnAna::synapse0xa00bce0() {
   return (neuron0xa026ee8()*-0.590195);
}

double AnnAna::synapse0xa00bd08() {
   return (neuron0xa027050()*1.05254);
}

double AnnAna::synapse0xa00bd30() {
   return (neuron0xa0271b8()*-0.112009);
}

double AnnAna::synapse0xa00bd58() {
   return (neuron0xa027320()*-0.813265);
}

double AnnAna::synapse0xa00bd80() {
   return (neuron0xa027488()*-0.800471);
}

double AnnAna::synapse0xa00bda8() {
   return (neuron0xa0275f0()*-0.133386);
}

double AnnAna::synapse0xa00bdd0() {
   return (neuron0xa027758()*0.0249995);
}

double AnnAna::synapse0xa00bdf8() {
   return (neuron0xa0278c0()*-0.241862);
}

double AnnAna::synapse0xa00be20() {
   return (neuron0xa0064a0()*0.214946);
}

double AnnAna::synapse0xa00be48() {
   return (neuron0xa006608()*0.21264);
}

double AnnAna::synapse0xa00be70() {
   return (neuron0xa006770()*0.807581);
}

double AnnAna::synapse0xa00bfc0() {
   return (neuron0xa0253a8()*0.286766);
}

double AnnAna::synapse0xa00bfe8() {
   return (neuron0xa025510()*-0.544573);
}

double AnnAna::synapse0xa00c010() {
   return (neuron0xa025678()*0.835343);
}

double AnnAna::synapse0xa00c038() {
   return (neuron0xa0257e0()*0.035209);
}

double AnnAna::synapse0xa00c060() {
   return (neuron0xa025948()*0.35079);
}

double AnnAna::synapse0xa00c088() {
   return (neuron0xa025ab8()*0.748622);
}

double AnnAna::synapse0xa00c0b0() {
   return (neuron0xa025c20()*-1.10798);
}

double AnnAna::synapse0xa00c0d8() {
   return (neuron0xa025d88()*-0.933313);
}

double AnnAna::synapse0xa00c100() {
   return (neuron0xa025ef0()*-0.474041);
}

double AnnAna::synapse0xa00c128() {
   return (neuron0xa026058()*-0.228029);
}

double AnnAna::synapse0xa00c150() {
   return (neuron0xa0261c0()*0.368422);
}

double AnnAna::synapse0xa00c178() {
   return (neuron0xa026328()*0.58559);
}

double AnnAna::synapse0xa00c1a0() {
   return (neuron0xa026490()*-0.45313);
}

double AnnAna::synapse0xa00c1c8() {
   return (neuron0xa0265f8()*0.12924);
}

double AnnAna::synapse0xa00c1f0() {
   return (neuron0xa026760()*-0.44446);
}

double AnnAna::synapse0xa00c218() {
   return (neuron0xa0268c8()*-0.130893);
}

double AnnAna::synapse0xa00c240() {
   return (neuron0xa026a30()*0.000260408);
}

double AnnAna::synapse0xa00c2f0() {
   return (neuron0xa026ca8()*0.475652);
}

double AnnAna::synapse0xa00c318() {
   return (neuron0xa026d80()*-0.35965);
}

double AnnAna::synapse0xa00c340() {
   return (neuron0xa026ee8()*-0.240043);
}

double AnnAna::synapse0xa00c368() {
   return (neuron0xa027050()*0.646617);
}

double AnnAna::synapse0xa00c390() {
   return (neuron0xa0271b8()*0.139155);
}

double AnnAna::synapse0xa00c3b8() {
   return (neuron0xa027320()*0.0907881);
}

double AnnAna::synapse0xa00c3e0() {
   return (neuron0xa027488()*-0.566367);
}

double AnnAna::synapse0xa00c408() {
   return (neuron0xa0275f0()*0.284069);
}

double AnnAna::synapse0xa00c430() {
   return (neuron0xa027758()*0.790254);
}

double AnnAna::synapse0xa00c458() {
   return (neuron0xa0278c0()*-0.79975);
}

double AnnAna::synapse0xa00c480() {
   return (neuron0xa0064a0()*-0.261207);
}

double AnnAna::synapse0xa00c4a8() {
   return (neuron0xa006608()*0.889924);
}

double AnnAna::synapse0xa00c4d0() {
   return (neuron0xa006770()*-0.0449055);
}

double AnnAna::synapse0xa00c620() {
   return (neuron0xa0253a8()*0.106782);
}

double AnnAna::synapse0xa00c648() {
   return (neuron0xa025510()*-0.172299);
}

double AnnAna::synapse0xa00c670() {
   return (neuron0xa025678()*0.535113);
}

double AnnAna::synapse0xa00c698() {
   return (neuron0xa0257e0()*-0.235005);
}

double AnnAna::synapse0xa00c6c0() {
   return (neuron0xa025948()*0.0945094);
}

double AnnAna::synapse0xa00c6e8() {
   return (neuron0xa025ab8()*-0.430856);
}

double AnnAna::synapse0xa00c710() {
   return (neuron0xa025c20()*-0.261229);
}

double AnnAna::synapse0xa00c738() {
   return (neuron0xa025d88()*-0.000573766);
}

double AnnAna::synapse0xa00c760() {
   return (neuron0xa025ef0()*-0.144266);
}

double AnnAna::synapse0xa00c788() {
   return (neuron0xa026058()*0.318798);
}

double AnnAna::synapse0xa00c7b0() {
   return (neuron0xa0261c0()*0.0272416);
}

double AnnAna::synapse0xa00c7d8() {
   return (neuron0xa026328()*-0.10133);
}

double AnnAna::synapse0xa00c800() {
   return (neuron0xa026490()*0.507506);
}

double AnnAna::synapse0xa00c828() {
   return (neuron0xa0265f8()*-0.227377);
}

double AnnAna::synapse0xa00c850() {
   return (neuron0xa026760()*-0.804315);
}

double AnnAna::synapse0xa00c878() {
   return (neuron0xa0268c8()*0.37846);
}

double AnnAna::synapse0xa00c8a0() {
   return (neuron0xa026a30()*0.323508);
}

double AnnAna::synapse0xa00c950() {
   return (neuron0xa026ca8()*-0.0343451);
}

double AnnAna::synapse0xa00c978() {
   return (neuron0xa026d80()*-0.164686);
}

double AnnAna::synapse0xa00c9a0() {
   return (neuron0xa026ee8()*-0.229259);
}

double AnnAna::synapse0xa00c9c8() {
   return (neuron0xa027050()*-0.735297);
}

double AnnAna::synapse0xa00c9f0() {
   return (neuron0xa0271b8()*-0.49187);
}

double AnnAna::synapse0xa00ca18() {
   return (neuron0xa027320()*-0.847951);
}

double AnnAna::synapse0xa00ca40() {
   return (neuron0xa027488()*-0.662339);
}

double AnnAna::synapse0xa00ca68() {
   return (neuron0xa0275f0()*0.470416);
}

double AnnAna::synapse0xa00ca90() {
   return (neuron0xa027758()*0.0106623);
}

double AnnAna::synapse0xa00cab8() {
   return (neuron0xa0278c0()*-0.331726);
}

double AnnAna::synapse0xa00cae0() {
   return (neuron0xa0064a0()*-1.05366);
}

double AnnAna::synapse0xa00cb08() {
   return (neuron0xa006608()*-0.209267);
}

double AnnAna::synapse0xa00cb30() {
   return (neuron0xa006770()*0.397805);
}

double AnnAna::synapse0xa00cc80() {
   return (neuron0xa0253a8()*-0.194431);
}

double AnnAna::synapse0xa00cca8() {
   return (neuron0xa025510()*0.50544);
}

double AnnAna::synapse0xa00ccd0() {
   return (neuron0xa025678()*0.529926);
}

double AnnAna::synapse0xa00ccf8() {
   return (neuron0xa0257e0()*0.597738);
}

double AnnAna::synapse0xa00cd20() {
   return (neuron0xa025948()*0.284496);
}

double AnnAna::synapse0xa00cd48() {
   return (neuron0xa025ab8()*0.182485);
}

double AnnAna::synapse0xa00cd70() {
   return (neuron0xa025c20()*-0.549206);
}

double AnnAna::synapse0xa00cd98() {
   return (neuron0xa025d88()*0.529639);
}

double AnnAna::synapse0xa00cdc0() {
   return (neuron0xa025ef0()*-0.0321054);
}

double AnnAna::synapse0xa00cde8() {
   return (neuron0xa026058()*-0.366778);
}

double AnnAna::synapse0xa00ce10() {
   return (neuron0xa0261c0()*-0.0941685);
}

double AnnAna::synapse0xa00ce38() {
   return (neuron0xa026328()*0.555763);
}

double AnnAna::synapse0xa00ce60() {
   return (neuron0xa026490()*-0.151239);
}

double AnnAna::synapse0xa00ce88() {
   return (neuron0xa0265f8()*0.300054);
}

double AnnAna::synapse0xa00ceb0() {
   return (neuron0xa026760()*0.0890819);
}

double AnnAna::synapse0xa00ced8() {
   return (neuron0xa0268c8()*0.506204);
}

double AnnAna::synapse0xa00cf00() {
   return (neuron0xa026a30()*-0.133234);
}

double AnnAna::synapse0xa00cfb0() {
   return (neuron0xa026ca8()*0.287972);
}

double AnnAna::synapse0xa00cfd8() {
   return (neuron0xa026d80()*0.242932);
}

double AnnAna::synapse0xa00d000() {
   return (neuron0xa026ee8()*-0.0751116);
}

double AnnAna::synapse0xa00d028() {
   return (neuron0xa027050()*-0.258477);
}

double AnnAna::synapse0xa00d050() {
   return (neuron0xa0271b8()*-0.370537);
}

double AnnAna::synapse0xa00d078() {
   return (neuron0xa027320()*-0.173043);
}

double AnnAna::synapse0xa00d0a0() {
   return (neuron0xa027488()*-0.430735);
}

double AnnAna::synapse0xa00d0c8() {
   return (neuron0xa0275f0()*0.659618);
}

double AnnAna::synapse0xa00d0f0() {
   return (neuron0xa027758()*0.42113);
}

double AnnAna::synapse0xa00d118() {
   return (neuron0xa0278c0()*0.181532);
}

double AnnAna::synapse0xa00d140() {
   return (neuron0xa0064a0()*0.0047959);
}

double AnnAna::synapse0xa00d168() {
   return (neuron0xa006608()*-0.177799);
}

double AnnAna::synapse0xa00d190() {
   return (neuron0xa006770()*0.456877);
}

double AnnAna::synapse0xa00d2e0() {
   return (neuron0xa0069f0()*-0.384874);
}

double AnnAna::synapse0xa00d308() {
   return (neuron0xa007100()*0.0940238);
}

double AnnAna::synapse0xa00d330() {
   return (neuron0xa007708()*-0.155807);
}

double AnnAna::synapse0xa00d358() {
   return (neuron0xa007f10()*0.473196);
}

double AnnAna::synapse0xa00d380() {
   return (neuron0xa008570()*-0.216166);
}

double AnnAna::synapse0xa00d3a8() {
   return (neuron0xa008de8()*0.33456);
}

double AnnAna::synapse0xa00d3d0() {
   return (neuron0xa009448()*-0.306346);
}

double AnnAna::synapse0xa00d3f8() {
   return (neuron0xa009aa8()*0.186838);
}

double AnnAna::synapse0xa00d420() {
   return (neuron0xa00a108()*0.337553);
}

double AnnAna::synapse0xa00d448() {
   return (neuron0xa0089e0()*0.0421072);
}

double AnnAna::synapse0xa00d470() {
   return (neuron0xa00b1d8()*0.447237);
}

double AnnAna::synapse0xa00d498() {
   return (neuron0xa00b838()*-0.120969);
}

double AnnAna::synapse0xa00d4c0() {
   return (neuron0xa00be98()*0.0539114);
}

double AnnAna::synapse0xa00d4e8() {
   return (neuron0xa00c4f8()*0.383094);
}

double AnnAna::synapse0xa00d510() {
   return (neuron0xa00cb58()*-0.199684);
}

double AnnAna::synapse0xa00d6a8() {
   return (neuron0xa0069f0()*0.484647);
}

double AnnAna::synapse0xa00d6d0() {
   return (neuron0xa007100()*0.105443);
}

double AnnAna::synapse0xa00d6f8() {
   return (neuron0xa007708()*0.0275635);
}

double AnnAna::synapse0xa00d720() {
   return (neuron0xa007f10()*-0.284578);
}

double AnnAna::synapse0xa00d748() {
   return (neuron0xa008570()*0.251811);
}

double AnnAna::synapse0xa00d770() {
   return (neuron0xa008de8()*0.251579);
}

double AnnAna::synapse0xa00d798() {
   return (neuron0xa009448()*-0.380857);
}

double AnnAna::synapse0xa00d7c0() {
   return (neuron0xa009aa8()*0.259627);
}

double AnnAna::synapse0xa00d7e8() {
   return (neuron0xa00a108()*0.462864);
}

double AnnAna::synapse0xa00d810() {
   return (neuron0xa0089e0()*0.413905);
}

double AnnAna::synapse0xa00d838() {
   return (neuron0xa00b1d8()*0.0258886);
}

double AnnAna::synapse0xa00d860() {
   return (neuron0xa00b838()*0.303609);
}

double AnnAna::synapse0xa00d888() {
   return (neuron0xa00be98()*0.387216);
}

double AnnAna::synapse0xa00d8b0() {
   return (neuron0xa00c4f8()*0.0782564);
}

double AnnAna::synapse0xa00d8d8() {
   return (neuron0xa00cb58()*0.402414);
}

double AnnAna::synapse0xa00da70() {
   return (neuron0xa0069f0()*-0.348239);
}

double AnnAna::synapse0xa00da98() {
   return (neuron0xa007100()*-0.604211);
}

double AnnAna::synapse0xa00dac0() {
   return (neuron0xa007708()*-0.439808);
}

double AnnAna::synapse0xa00dae8() {
   return (neuron0xa007f10()*-0.666902);
}

double AnnAna::synapse0xa00db10() {
   return (neuron0xa008570()*-0.129008);
}

double AnnAna::synapse0xa00db38() {
   return (neuron0xa008de8()*0.351192);
}

double AnnAna::synapse0xa00db60() {
   return (neuron0xa009448()*-0.0624082);
}

double AnnAna::synapse0xa00db88() {
   return (neuron0xa009aa8()*-0.712092);
}

double AnnAna::synapse0xa00dbb0() {
   return (neuron0xa00a108()*-0.010182);
}

double AnnAna::synapse0xa00dbd8() {
   return (neuron0xa0089e0()*-0.0257799);
}

double AnnAna::synapse0xa00dc00() {
   return (neuron0xa00b1d8()*0.158289);
}

double AnnAna::synapse0xa00dc28() {
   return (neuron0xa00b838()*-0.330146);
}

double AnnAna::synapse0xa00dc50() {
   return (neuron0xa00be98()*-0.153047);
}

double AnnAna::synapse0xa00dc78() {
   return (neuron0xa00c4f8()*-0.0921071);
}

double AnnAna::synapse0xa00dca0() {
   return (neuron0xa00cb58()*-0.0358956);
}

double AnnAna::synapse0xa00de38() {
   return (neuron0xa0069f0()*-0.488304);
}

double AnnAna::synapse0xa00de60() {
   return (neuron0xa007100()*0.144286);
}

double AnnAna::synapse0xa00de88() {
   return (neuron0xa007708()*-0.263408);
}

double AnnAna::synapse0xa00deb0() {
   return (neuron0xa007f10()*0.28948);
}

double AnnAna::synapse0xa00ded8() {
   return (neuron0xa008570()*-0.178015);
}

double AnnAna::synapse0xa00df00() {
   return (neuron0xa008de8()*-0.308244);
}

double AnnAna::synapse0xa00df28() {
   return (neuron0xa009448()*-0.313449);
}

double AnnAna::synapse0xa00df50() {
   return (neuron0xa009aa8()*-0.20763);
}

double AnnAna::synapse0xa00df78() {
   return (neuron0xa00a108()*-0.590426);
}

double AnnAna::synapse0xa00dfa0() {
   return (neuron0xa0089e0()*-0.0908632);
}

double AnnAna::synapse0xa00dfc8() {
   return (neuron0xa00b1d8()*0.153761);
}

double AnnAna::synapse0xa00dff0() {
   return (neuron0xa00b838()*0.107364);
}

double AnnAna::synapse0xa00e018() {
   return (neuron0xa00be98()*-0.395926);
}

double AnnAna::synapse0xa00e040() {
   return (neuron0xa00c4f8()*0.0233321);
}

double AnnAna::synapse0xa00e068() {
   return (neuron0xa00cb58()*-0.0865059);
}

double AnnAna::synapse0xa00e200() {
   return (neuron0xa0069f0()*0.204983);
}

double AnnAna::synapse0xa00e228() {
   return (neuron0xa007100()*-0.539987);
}

double AnnAna::synapse0xa00e250() {
   return (neuron0xa007708()*0.179905);
}

double AnnAna::synapse0xa00a560() {
   return (neuron0xa007f10()*0.202662);
}

double AnnAna::synapse0xa00a588() {
   return (neuron0xa008570()*0.0387127);
}

double AnnAna::synapse0xa00a5b0() {
   return (neuron0xa008de8()*-0.429266);
}

double AnnAna::synapse0xa00a5d8() {
   return (neuron0xa009448()*-0.0717374);
}

double AnnAna::synapse0xa00a600() {
   return (neuron0xa009aa8()*-0.324423);
}

double AnnAna::synapse0xa00a628() {
   return (neuron0xa00a108()*-0.129658);
}

double AnnAna::synapse0xa00a650() {
   return (neuron0xa0089e0()*-0.171942);
}

double AnnAna::synapse0xa00a678() {
   return (neuron0xa00b1d8()*-0.128335);
}

double AnnAna::synapse0xa00a6a0() {
   return (neuron0xa00b838()*-0.585285);
}

double AnnAna::synapse0xa00a6c8() {
   return (neuron0xa00be98()*0.14541);
}

double AnnAna::synapse0xa00a6f0() {
   return (neuron0xa00c4f8()*0.116967);
}

double AnnAna::synapse0xa00a718() {
   return (neuron0xa00cb58()*-0.577416);
}

double AnnAna::synapse0xa00a8b0() {
   return (neuron0xa0069f0()*0.133004);
}

double AnnAna::synapse0xa00a8d8() {
   return (neuron0xa007100()*0.246244);
}

double AnnAna::synapse0xa00a900() {
   return (neuron0xa007708()*-0.115466);
}

double AnnAna::synapse0xa00a928() {
   return (neuron0xa007f10()*-0.0188621);
}

double AnnAna::synapse0xa00a950() {
   return (neuron0xa008570()*-0.0144285);
}

double AnnAna::synapse0xa00a978() {
   return (neuron0xa008de8()*-0.053055);
}

double AnnAna::synapse0xa00a9a0() {
   return (neuron0xa009448()*-0.17596);
}

double AnnAna::synapse0xa00a9c8() {
   return (neuron0xa009aa8()*-0.0491296);
}

double AnnAna::synapse0xa00a9f0() {
   return (neuron0xa00a108()*0.322741);
}

double AnnAna::synapse0xa00aa18() {
   return (neuron0xa0089e0()*-0.243188);
}

double AnnAna::synapse0xa00aa40() {
   return (neuron0xa00b1d8()*0.207166);
}

double AnnAna::synapse0xa00aa68() {
   return (neuron0xa00b838()*-0.628805);
}

double AnnAna::synapse0xa00aa90() {
   return (neuron0xa00be98()*-0.0263036);
}

double AnnAna::synapse0xa00aab8() {
   return (neuron0xa00c4f8()*-0.40081);
}

double AnnAna::synapse0xa00aae0() {
   return (neuron0xa00cb58()*0.418532);
}

double AnnAna::synapse0xa00ac78() {
   return (neuron0xa0069f0()*-2.40624);
}

double AnnAna::synapse0xa00aca0() {
   return (neuron0xa007100()*-1.33207);
}

double AnnAna::synapse0xa00acc8() {
   return (neuron0xa007708()*-0.934918);
}

double AnnAna::synapse0xa00acf0() {
   return (neuron0xa007f10()*-0.605201);
}

double AnnAna::synapse0xa00ad18() {
   return (neuron0xa008570()*-0.00896515);
}

double AnnAna::synapse0xa00ad40() {
   return (neuron0xa008de8()*0.856137);
}

double AnnAna::synapse0xa00f280() {
   return (neuron0xa009448()*-0.436424);
}

double AnnAna::synapse0xa00f2a8() {
   return (neuron0xa009aa8()*0.0240858);
}

double AnnAna::synapse0xa00f2d0() {
   return (neuron0xa00a108()*-0.589296);
}

double AnnAna::synapse0xa00f2f8() {
   return (neuron0xa0089e0()*0.696249);
}

double AnnAna::synapse0xa00f320() {
   return (neuron0xa00b1d8()*1.19113);
}

double AnnAna::synapse0xa00f348() {
   return (neuron0xa00b838()*-0.549915);
}

double AnnAna::synapse0xa00f370() {
   return (neuron0xa00be98()*-1.24023);
}

double AnnAna::synapse0xa00f398() {
   return (neuron0xa00c4f8()*-0.169482);
}

double AnnAna::synapse0xa00f3c0() {
   return (neuron0xa00cb58()*0.405023);
}

double AnnAna::synapse0xa00f558() {
   return (neuron0xa0069f0()*0.396361);
}

double AnnAna::synapse0xa00f580() {
   return (neuron0xa007100()*0.464048);
}

double AnnAna::synapse0xa00f5a8() {
   return (neuron0xa007708()*0.530724);
}

double AnnAna::synapse0xa00f5d0() {
   return (neuron0xa007f10()*0.295795);
}

double AnnAna::synapse0xa00f5f8() {
   return (neuron0xa008570()*0.334223);
}

double AnnAna::synapse0xa00f620() {
   return (neuron0xa008de8()*-0.437479);
}

double AnnAna::synapse0xa00f648() {
   return (neuron0xa009448()*0.148524);
}

double AnnAna::synapse0xa00f670() {
   return (neuron0xa009aa8()*-0.0329386);
}

double AnnAna::synapse0xa00f698() {
   return (neuron0xa00a108()*0.0252872);
}

double AnnAna::synapse0xa00f6c0() {
   return (neuron0xa0089e0()*-0.285163);
}

double AnnAna::synapse0xa00f6e8() {
   return (neuron0xa00b1d8()*-0.0126029);
}

double AnnAna::synapse0xa00f710() {
   return (neuron0xa00b838()*0.456445);
}

double AnnAna::synapse0xa00f738() {
   return (neuron0xa00be98()*0.114465);
}

double AnnAna::synapse0xa00f760() {
   return (neuron0xa00c4f8()*0.0246948);
}

double AnnAna::synapse0xa00f788() {
   return (neuron0xa00cb58()*-0.835319);
}

double AnnAna::synapse0xa00f920() {
   return (neuron0xa0069f0()*0.0626803);
}

double AnnAna::synapse0xa00f948() {
   return (neuron0xa007100()*-0.527751);
}

double AnnAna::synapse0xa00f970() {
   return (neuron0xa007708()*-0.140207);
}

double AnnAna::synapse0xa00f998() {
   return (neuron0xa007f10()*-0.423579);
}

double AnnAna::synapse0xa00f9c0() {
   return (neuron0xa008570()*0.0312344);
}

double AnnAna::synapse0xa00f9e8() {
   return (neuron0xa008de8()*-0.0753875);
}

double AnnAna::synapse0xa00fa10() {
   return (neuron0xa009448()*0.420892);
}

double AnnAna::synapse0xa00fa38() {
   return (neuron0xa009aa8()*0.155242);
}

double AnnAna::synapse0xa00fa60() {
   return (neuron0xa00a108()*0.134481);
}

double AnnAna::synapse0xa00fa88() {
   return (neuron0xa0089e0()*0.399717);
}

double AnnAna::synapse0xa00fab0() {
   return (neuron0xa00b1d8()*-0.123281);
}

double AnnAna::synapse0xa00fad8() {
   return (neuron0xa00b838()*0.288251);
}

double AnnAna::synapse0xa00fb00() {
   return (neuron0xa00be98()*-0.341913);
}

double AnnAna::synapse0xa00fb28() {
   return (neuron0xa00c4f8()*-0.379045);
}

double AnnAna::synapse0xa00fb50() {
   return (neuron0xa00cb58()*0.232816);
}

double AnnAna::synapse0xa00fce8() {
   return (neuron0xa0069f0()*0.0699192);
}

double AnnAna::synapse0xa00fd10() {
   return (neuron0xa007100()*-0.0457207);
}

double AnnAna::synapse0xa00fd38() {
   return (neuron0xa007708()*0.56123);
}

double AnnAna::synapse0xa00fd60() {
   return (neuron0xa007f10()*0.214979);
}

double AnnAna::synapse0xa00fd88() {
   return (neuron0xa008570()*-0.321268);
}

double AnnAna::synapse0xa00fdb0() {
   return (neuron0xa008de8()*0.165894);
}

double AnnAna::synapse0xa00fdd8() {
   return (neuron0xa009448()*-0.357282);
}

double AnnAna::synapse0xa00fe00() {
   return (neuron0xa009aa8()*-0.21962);
}

double AnnAna::synapse0xa00fe28() {
   return (neuron0xa00a108()*-0.311705);
}

double AnnAna::synapse0xa00fe50() {
   return (neuron0xa0089e0()*-0.047297);
}

double AnnAna::synapse0xa00fe78() {
   return (neuron0xa00b1d8()*0.0866319);
}

double AnnAna::synapse0xa00fea0() {
   return (neuron0xa00b838()*0.13174);
}

double AnnAna::synapse0xa00fec8() {
   return (neuron0xa00be98()*-0.145428);
}

double AnnAna::synapse0xa00fef0() {
   return (neuron0xa00c4f8()*0.389642);
}

double AnnAna::synapse0xa00ff18() {
   return (neuron0xa00cb58()*-0.170325);
}

double AnnAna::synapse0xa0100b0() {
   return (neuron0xa0069f0()*-0.655344);
}

double AnnAna::synapse0xa0100d8() {
   return (neuron0xa007100()*-0.0593815);
}

double AnnAna::synapse0xa010100() {
   return (neuron0xa007708()*0.0528425);
}

double AnnAna::synapse0xa010128() {
   return (neuron0xa007f10()*0.349052);
}

double AnnAna::synapse0xa010150() {
   return (neuron0xa008570()*-0.0697263);
}

double AnnAna::synapse0xa010178() {
   return (neuron0xa008de8()*0.229337);
}

double AnnAna::synapse0xa0101a0() {
   return (neuron0xa009448()*0.109145);
}

double AnnAna::synapse0xa0101c8() {
   return (neuron0xa009aa8()*-0.0823547);
}

double AnnAna::synapse0xa0101f0() {
   return (neuron0xa00a108()*0.698658);
}

double AnnAna::synapse0xa010218() {
   return (neuron0xa0089e0()*0.363601);
}

double AnnAna::synapse0xa010240() {
   return (neuron0xa00b1d8()*0.110295);
}

double AnnAna::synapse0xa010268() {
   return (neuron0xa00b838()*0.178005);
}

double AnnAna::synapse0xa010290() {
   return (neuron0xa00be98()*0.640629);
}

double AnnAna::synapse0xa0102b8() {
   return (neuron0xa00c4f8()*-0.776534);
}

double AnnAna::synapse0xa0102e0() {
   return (neuron0xa00cb58()*0.0330286);
}

double AnnAna::synapse0xa010478() {
   return (neuron0xa0069f0()*0.00257742);
}

double AnnAna::synapse0xa0104a0() {
   return (neuron0xa007100()*0.599457);
}

double AnnAna::synapse0xa0104c8() {
   return (neuron0xa007708()*0.285182);
}

double AnnAna::synapse0xa0104f0() {
   return (neuron0xa007f10()*-0.177449);
}

double AnnAna::synapse0xa010518() {
   return (neuron0xa008570()*-0.209605);
}

double AnnAna::synapse0xa010540() {
   return (neuron0xa008de8()*-0.0524468);
}

double AnnAna::synapse0xa010568() {
   return (neuron0xa009448()*-0.269154);
}

double AnnAna::synapse0xa010590() {
   return (neuron0xa009aa8()*0.261944);
}

double AnnAna::synapse0xa0105b8() {
   return (neuron0xa00a108()*0.366011);
}

double AnnAna::synapse0xa0105e0() {
   return (neuron0xa0089e0()*0.0822091);
}

double AnnAna::synapse0xa010608() {
   return (neuron0xa00b1d8()*-0.206204);
}

double AnnAna::synapse0xa010630() {
   return (neuron0xa00b838()*-0.401717);
}

double AnnAna::synapse0xa010658() {
   return (neuron0xa00be98()*0.543885);
}

double AnnAna::synapse0xa010680() {
   return (neuron0xa00c4f8()*0.441203);
}

double AnnAna::synapse0xa0106a8() {
   return (neuron0xa00cb58()*0.238429);
}

double AnnAna::synapse0xa010840() {
   return (neuron0xa0069f0()*-0.228588);
}

double AnnAna::synapse0xa010868() {
   return (neuron0xa007100()*-0.154059);
}

double AnnAna::synapse0xa010890() {
   return (neuron0xa007708()*-0.354672);
}

double AnnAna::synapse0xa0108b8() {
   return (neuron0xa007f10()*0.00283532);
}

double AnnAna::synapse0xa0108e0() {
   return (neuron0xa008570()*0.346283);
}

double AnnAna::synapse0xa010908() {
   return (neuron0xa008de8()*-0.386103);
}

double AnnAna::synapse0xa010930() {
   return (neuron0xa009448()*0.192452);
}

double AnnAna::synapse0xa010958() {
   return (neuron0xa009aa8()*0.0429158);
}

double AnnAna::synapse0xa010980() {
   return (neuron0xa00a108()*0.224287);
}

double AnnAna::synapse0xa0109a8() {
   return (neuron0xa0089e0()*-0.248358);
}

double AnnAna::synapse0xa0109d0() {
   return (neuron0xa00b1d8()*-0.546647);
}

double AnnAna::synapse0xa0109f8() {
   return (neuron0xa00b838()*0.281783);
}

double AnnAna::synapse0xa010a20() {
   return (neuron0xa00be98()*-0.263914);
}

double AnnAna::synapse0xa010a48() {
   return (neuron0xa00c4f8()*0.36657);
}

double AnnAna::synapse0xa010a70() {
   return (neuron0xa00cb58()*-0.462262);
}

double AnnAna::synapse0xa010c08() {
   return (neuron0xa0069f0()*0.413605);
}

double AnnAna::synapse0xa010c30() {
   return (neuron0xa007100()*-0.46025);
}

double AnnAna::synapse0xa010c58() {
   return (neuron0xa007708()*-0.349511);
}

double AnnAna::synapse0xa010c80() {
   return (neuron0xa007f10()*0.227848);
}

double AnnAna::synapse0xa010ca8() {
   return (neuron0xa008570()*0.324642);
}

double AnnAna::synapse0xa010cd0() {
   return (neuron0xa008de8()*-0.197649);
}

double AnnAna::synapse0xa010cf8() {
   return (neuron0xa009448()*-0.255278);
}

double AnnAna::synapse0xa010d20() {
   return (neuron0xa009aa8()*0.180081);
}

double AnnAna::synapse0xa010d48() {
   return (neuron0xa00a108()*-0.236719);
}

double AnnAna::synapse0xa010d70() {
   return (neuron0xa0089e0()*-0.257446);
}

double AnnAna::synapse0xa010d98() {
   return (neuron0xa00b1d8()*0.369206);
}

double AnnAna::synapse0xa010dc0() {
   return (neuron0xa00b838()*0.249322);
}

double AnnAna::synapse0xa010de8() {
   return (neuron0xa00be98()*-0.459041);
}

double AnnAna::synapse0xa010e10() {
   return (neuron0xa00c4f8()*-0.107539);
}

double AnnAna::synapse0xa010e38() {
   return (neuron0xa00cb58()*-0.337229);
}

double AnnAna::synapse0xa010fd0() {
   return (neuron0xa0069f0()*-0.31218);
}

double AnnAna::synapse0xa010ff8() {
   return (neuron0xa007100()*-0.41505);
}

double AnnAna::synapse0xa011020() {
   return (neuron0xa007708()*0.288695);
}

double AnnAna::synapse0xa011048() {
   return (neuron0xa007f10()*0.192245);
}

double AnnAna::synapse0xa011070() {
   return (neuron0xa008570()*-0.468671);
}

double AnnAna::synapse0xa011098() {
   return (neuron0xa008de8()*0.298661);
}

double AnnAna::synapse0xa0110c0() {
   return (neuron0xa009448()*-0.176758);
}

double AnnAna::synapse0xa0110e8() {
   return (neuron0xa009aa8()*0.160656);
}

double AnnAna::synapse0xa011110() {
   return (neuron0xa00a108()*-0.338612);
}

double AnnAna::synapse0xa011138() {
   return (neuron0xa0089e0()*-0.0539761);
}

double AnnAna::synapse0xa011160() {
   return (neuron0xa00b1d8()*0.0142956);
}

double AnnAna::synapse0xa011188() {
   return (neuron0xa00b838()*-0.0898725);
}

double AnnAna::synapse0xa0111b0() {
   return (neuron0xa00be98()*-0.466561);
}

double AnnAna::synapse0xa0111d8() {
   return (neuron0xa00c4f8()*-0.196646);
}

double AnnAna::synapse0xa011200() {
   return (neuron0xa00cb58()*-0.239142);
}

double AnnAna::synapse0xa0112b8() {
   return (neuron0xa00d1d8()*-0.177715);
}

double AnnAna::synapse0xa0112e0() {
   return (neuron0xa00d538()*-0.102382);
}

double AnnAna::synapse0xa011308() {
   return (neuron0xa00d900()*-0.340363);
}

double AnnAna::synapse0xa011330() {
   return (neuron0xa00dcc8()*0.198273);
}

double AnnAna::synapse0xa011358() {
   return (neuron0xa00e090()*0.0571124);
}

double AnnAna::synapse0xa011380() {
   return (neuron0xa00a740()*0.608415);
}

double AnnAna::synapse0xa0113a8() {
   return (neuron0xa00ab08()*1.23687);
}

double AnnAna::synapse0xa0113d0() {
   return (neuron0xa00f3e8()*0.0815651);
}

double AnnAna::synapse0xa0113f8() {
   return (neuron0xa00f7b0()*-0.338155);
}

double AnnAna::synapse0xa011420() {
   return (neuron0xa00fb78()*-0.343967);
}

double AnnAna::synapse0xa011448() {
   return (neuron0xa00ff40()*0.434548);
}

double AnnAna::synapse0xa011470() {
   return (neuron0xa010308()*0.0223337);
}

double AnnAna::synapse0xa011498() {
   return (neuron0xa0106d0()*-0.277661);
}

double AnnAna::synapse0xa0114c0() {
   return (neuron0xa010a98()*0.0316252);
}

double AnnAna::synapse0xa0114e8() {
   return (neuron0xa010e60()*-0.138036);
}

