#include "NueAna/mlpANN.h"
#include "math.h"

double mlpANN::value(int index,double in0,double in1,double in2,double in3,double in4,double in5,double in6,double in7,double in8,double in9,double in10,double in11,double in12,double in13) {
   input0 = (in0 - 0.310947)/0.116781;
   input1 = (in1 - 0.46791)/0.129756;
   input2 = (in2 - 0.59676)/0.140262;
   input3 = (in3 - 0.69552)/0.143586;
   input4 = (in4 - 0.772466)/0.140853;
   input5 = (in5 - 0.829331)/0.134609;
   input6 = (in6 - 0.305036)/0.107442;
   input7 = (in7 - 0.465258)/0.120199;
   input8 = (in8 - 0.574227)/0.12166;
   input9 = (in9 - 0.654448)/0.119232;
   input10 = (in10 - 0.715649)/0.115052;
   input11 = (in11 - 0.763165)/0.110242;
   input12 = (in12 - 0.691771)/0.15586;
   input13 = (in13 - 3.88725)/2.76409;
   switch(index) {
     case 0:
         return neuron0x9e28dd0();
     default:
         return 0.;
   }
}

double mlpANN::neuron0x95a0f78() {
   return input0;
}

double mlpANN::neuron0x931f328() {
   return input1;
}

double mlpANN::neuron0x9b48b50() {
   return input2;
}

double mlpANN::neuron0x9b48cc8() {
   return input3;
}

double mlpANN::neuron0x9b48e40() {
   return input4;
}

double mlpANN::neuron0x9b48fb8() {
   return input5;
}

double mlpANN::neuron0x9e282a0() {
   return input6;
}

double mlpANN::neuron0x9e283d0() {
   return input7;
}

double mlpANN::neuron0x9e28548() {
   return input8;
}

double mlpANN::neuron0x9e286c0() {
   return input9;
}

double mlpANN::neuron0x9e28838() {
   return input10;
}

double mlpANN::neuron0x9e28988() {
   return input11;
}

double mlpANN::neuron0x9e28b08() {
   return input12;
}

double mlpANN::neuron0x9e28c80() {
   return input13;
}

double mlpANN::neuron0x9e28ee8() {
   double input = 0.478999;
   input += synapse0x9e29038();
   input += synapse0x9e29060();
   input += synapse0x9e29088();
   input += synapse0x9e290b0();
   input += synapse0x9e290d8();
   input += synapse0x9e29100();
   input += synapse0x9e29128();
   input += synapse0x9e29150();
   input += synapse0x9e29178();
   input += synapse0x9e291a0();
   input += synapse0x9e291c8();
   input += synapse0x9e291f0();
   input += synapse0x9e29218();
   input += synapse0x9e29240();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e29268() {
   double input = -1.47329;
   input += synapse0x9e293d8();
   input += synapse0x9e29400();
   input += synapse0x9e29428();
   input += synapse0x9e294d8();
   input += synapse0x9e29500();
   input += synapse0x9e29528();
   input += synapse0x9e29550();
   input += synapse0x9e29578();
   input += synapse0x9e295a0();
   input += synapse0x9e295c8();
   input += synapse0x9e295f0();
   input += synapse0x9e29618();
   input += synapse0x9e29640();
   input += synapse0x9e29668();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e29690() {
   double input = 0.930365;
   input += synapse0x90ac248();
   input += synapse0x9645e38();
   input += synapse0x9560d78();
   input += synapse0x9560da0();
   input += synapse0x9560dc8();
   input += synapse0x9560df0();
   input += synapse0x9e29450();
   input += synapse0x9e29478();
   input += synapse0x9e294a0();
   input += synapse0x9e87e70();
   input += synapse0x9e87e98();
   input += synapse0x9e87ec0();
   input += synapse0x9e87ee8();
   input += synapse0x9e87f10();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e87f38() {
   double input = -2.52937;
   input += synapse0x9e88040();
   input += synapse0x9e88068();
   input += synapse0x9e88090();
   input += synapse0x9e880b8();
   input += synapse0x9e880e0();
   input += synapse0x9e88108();
   input += synapse0x9e88130();
   input += synapse0x9e88158();
   input += synapse0x9e88180();
   input += synapse0x9e881a8();
   input += synapse0x9e881d0();
   input += synapse0x9e881f8();
   input += synapse0x9e88220();
   input += synapse0x9e88248();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e88270() {
   double input = 2.0951;
   input += synapse0x9e883c0();
   input += synapse0x9e883e8();
   input += synapse0x9e88410();
   input += synapse0x9e88438();
   input += synapse0x9e88460();
   input += synapse0x9e88488();
   input += synapse0x9e884b0();
   input += synapse0x9e884d8();
   input += synapse0x9e88500();
   input += synapse0x90ac218();
   input += synapse0x9e87d68();
   input += synapse0x9e87d90();
   input += synapse0x9e87db8();
   input += synapse0x9e87de0();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e88730() {
   double input = 3.26255;
   input += synapse0x9e888a0();
   input += synapse0x9e888c8();
   input += synapse0x9e888f0();
   input += synapse0x9e88918();
   input += synapse0x9e88940();
   input += synapse0x9e88968();
   input += synapse0x9e88990();
   input += synapse0x9e889b8();
   input += synapse0x9e889e0();
   input += synapse0x9e88a08();
   input += synapse0x9e88a30();
   input += synapse0x9e88a58();
   input += synapse0x9e88a80();
   input += synapse0x9e88aa8();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e88ad0() {
   double input = -0.455949;
   input += synapse0x9e88c40();
   input += synapse0x9e88c68();
   input += synapse0x9e88c90();
   input += synapse0x9e88cb8();
   input += synapse0x9e88ce0();
   input += synapse0x9e88d08();
   input += synapse0x9e88d30();
   input += synapse0x9e88d58();
   input += synapse0x9e88d80();
   input += synapse0x9e88da8();
   input += synapse0x9e88dd0();
   input += synapse0x9e88df8();
   input += synapse0x9e88e20();
   input += synapse0x9e88e48();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e88e70() {
   double input = 1.41615;
   input += synapse0x9e88fe0();
   input += synapse0x9e89008();
   input += synapse0x9e89030();
   input += synapse0x9e89058();
   input += synapse0x9e89080();
   input += synapse0x9e890a8();
   input += synapse0x9e890d0();
   input += synapse0x9e890f8();
   input += synapse0x9e89120();
   input += synapse0x9e89148();
   input += synapse0x9e89170();
   input += synapse0x9e89198();
   input += synapse0x9e891c0();
   input += synapse0x9e891e8();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e89210() {
   double input = -2.01413;
   input += synapse0x9e89380();
   input += synapse0x9e893a8();
   input += synapse0x9e893d0();
   input += synapse0x9e893f8();
   input += synapse0x9e89420();
   input += synapse0x9e89448();
   input += synapse0x9e89470();
   input += synapse0x9e89498();
   input += synapse0x9e894c0();
   input += synapse0x9e894e8();
   input += synapse0x9e89510();
   input += synapse0x9e89538();
   input += synapse0x9e89560();
   input += synapse0x9e89588();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e895b0() {
   double input = -0.114291;
   input += synapse0x9e87e08();
   input += synapse0x9e87e30();
   input += synapse0x9aebc28();
   input += synapse0x9e88528();
   input += synapse0x9e88550();
   input += synapse0x9e88578();
   input += synapse0x9e885a0();
   input += synapse0x9e885c8();
   input += synapse0x9e885f0();
   input += synapse0x9e88618();
   input += synapse0x9e88640();
   input += synapse0x9e88668();
   input += synapse0x9e88690();
   input += synapse0x9e886b8();
   return ((1/(1+exp(-input)))*1)+0;
}

double mlpANN::neuron0x9e28dd0() {
   double input = -0.622413;
   input += synapse0x9aec0a0();
   input += synapse0x9aec0c8();
   input += synapse0x9aec0f0();
   input += synapse0x9aec118();
   input += synapse0x9aec140();
   input += synapse0x9aec168();
   input += synapse0x9aec190();
   input += synapse0x9aec1b8();
   input += synapse0x9aec1e0();
   input += synapse0x9aec208();
   return input;
}

double mlpANN::synapse0x9e29038() {
   return (neuron0x95a0f78()*-0.985666);
}

double mlpANN::synapse0x9e29060() {
   return (neuron0x931f328()*0.122559);
}

double mlpANN::synapse0x9e29088() {
   return (neuron0x9b48b50()*-0.0452834);
}

double mlpANN::synapse0x9e290b0() {
   return (neuron0x9b48cc8()*0.750915);
}

double mlpANN::synapse0x9e290d8() {
   return (neuron0x9b48e40()*0.103577);
}

double mlpANN::synapse0x9e29100() {
   return (neuron0x9b48fb8()*0.170561);
}

double mlpANN::synapse0x9e29128() {
   return (neuron0x9e282a0()*-0.59644);
}

double mlpANN::synapse0x9e29150() {
   return (neuron0x9e283d0()*0.613393);
}

double mlpANN::synapse0x9e29178() {
   return (neuron0x9e28548()*0.53856);
}

double mlpANN::synapse0x9e291a0() {
   return (neuron0x9e286c0()*0.0239845);
}

double mlpANN::synapse0x9e291c8() {
   return (neuron0x9e28838()*-1.58596);
}

double mlpANN::synapse0x9e291f0() {
   return (neuron0x9e28988()*-0.224198);
}

double mlpANN::synapse0x9e29218() {
   return (neuron0x9e28b08()*1.47166);
}

double mlpANN::synapse0x9e29240() {
   return (neuron0x9e28c80()*-0.735476);
}

double mlpANN::synapse0x9e293d8() {
   return (neuron0x95a0f78()*-0.850271);
}

double mlpANN::synapse0x9e29400() {
   return (neuron0x931f328()*1.7883);
}

double mlpANN::synapse0x9e29428() {
   return (neuron0x9b48b50()*1.21318);
}

double mlpANN::synapse0x9e294d8() {
   return (neuron0x9b48cc8()*-0.474948);
}

double mlpANN::synapse0x9e29500() {
   return (neuron0x9b48e40()*0.61067);
}

double mlpANN::synapse0x9e29528() {
   return (neuron0x9b48fb8()*-1.6303);
}

double mlpANN::synapse0x9e29550() {
   return (neuron0x9e282a0()*0.481609);
}

double mlpANN::synapse0x9e29578() {
   return (neuron0x9e283d0()*-1.34116);
}

double mlpANN::synapse0x9e295a0() {
   return (neuron0x9e28548()*-1.4707);
}

double mlpANN::synapse0x9e295c8() {
   return (neuron0x9e286c0()*-0.828891);
}

double mlpANN::synapse0x9e295f0() {
   return (neuron0x9e28838()*-0.237887);
}

double mlpANN::synapse0x9e29618() {
   return (neuron0x9e28988()*0.459936);
}

double mlpANN::synapse0x9e29640() {
   return (neuron0x9e28b08()*-0.431486);
}

double mlpANN::synapse0x9e29668() {
   return (neuron0x9e28c80()*0.11127);
}

double mlpANN::synapse0x90ac248() {
   return (neuron0x95a0f78()*-1.65172);
}

double mlpANN::synapse0x9645e38() {
   return (neuron0x931f328()*1.55661);
}

double mlpANN::synapse0x9560d78() {
   return (neuron0x9b48b50()*1.10513);
}

double mlpANN::synapse0x9560da0() {
   return (neuron0x9b48cc8()*-1.72435);
}

double mlpANN::synapse0x9560dc8() {
   return (neuron0x9b48e40()*3.07762);
}

double mlpANN::synapse0x9560df0() {
   return (neuron0x9b48fb8()*-2.33158);
}

double mlpANN::synapse0x9e29450() {
   return (neuron0x9e282a0()*-1.43899);
}

double mlpANN::synapse0x9e29478() {
   return (neuron0x9e283d0()*0.177617);
}

double mlpANN::synapse0x9e294a0() {
   return (neuron0x9e28548()*0.218684);
}

double mlpANN::synapse0x9e87e70() {
   return (neuron0x9e286c0()*0.415998);
}

double mlpANN::synapse0x9e87e98() {
   return (neuron0x9e28838()*0.506269);
}

double mlpANN::synapse0x9e87ec0() {
   return (neuron0x9e28988()*-0.463399);
}

double mlpANN::synapse0x9e87ee8() {
   return (neuron0x9e28b08()*0.911479);
}

double mlpANN::synapse0x9e87f10() {
   return (neuron0x9e28c80()*-1.86222);
}

double mlpANN::synapse0x9e88040() {
   return (neuron0x95a0f78()*-1.12071);
}

double mlpANN::synapse0x9e88068() {
   return (neuron0x931f328()*-0.384408);
}

double mlpANN::synapse0x9e88090() {
   return (neuron0x9b48b50()*0.016917);
}

double mlpANN::synapse0x9e880b8() {
   return (neuron0x9b48cc8()*-0.782961);
}

double mlpANN::synapse0x9e880e0() {
   return (neuron0x9b48e40()*1.00143);
}

double mlpANN::synapse0x9e88108() {
   return (neuron0x9b48fb8()*0.452493);
}

double mlpANN::synapse0x9e88130() {
   return (neuron0x9e282a0()*-0.0605043);
}

double mlpANN::synapse0x9e88158() {
   return (neuron0x9e283d0()*-1.17155);
}

double mlpANN::synapse0x9e88180() {
   return (neuron0x9e28548()*0.931899);
}

double mlpANN::synapse0x9e881a8() {
   return (neuron0x9e286c0()*1.48143);
}

double mlpANN::synapse0x9e881d0() {
   return (neuron0x9e28838()*0.289564);
}

double mlpANN::synapse0x9e881f8() {
   return (neuron0x9e28988()*-2.12976);
}

double mlpANN::synapse0x9e88220() {
   return (neuron0x9e28b08()*2.59793);
}

double mlpANN::synapse0x9e88248() {
   return (neuron0x9e28c80()*0.204235);
}

double mlpANN::synapse0x9e883c0() {
   return (neuron0x95a0f78()*0.450529);
}

double mlpANN::synapse0x9e883e8() {
   return (neuron0x931f328()*-1.99407);
}

double mlpANN::synapse0x9e88410() {
   return (neuron0x9b48b50()*1.67663);
}

double mlpANN::synapse0x9e88438() {
   return (neuron0x9b48cc8()*0.734318);
}

double mlpANN::synapse0x9e88460() {
   return (neuron0x9b48e40()*-0.0403526);
}

double mlpANN::synapse0x9e88488() {
   return (neuron0x9b48fb8()*-0.958477);
}

double mlpANN::synapse0x9e884b0() {
   return (neuron0x9e282a0()*0.22483);
}

double mlpANN::synapse0x9e884d8() {
   return (neuron0x9e283d0()*-0.0411579);
}

double mlpANN::synapse0x9e88500() {
   return (neuron0x9e28548()*0.254688);
}

double mlpANN::synapse0x90ac218() {
   return (neuron0x9e286c0()*-0.946474);
}

double mlpANN::synapse0x9e87d68() {
   return (neuron0x9e28838()*0.456723);
}

double mlpANN::synapse0x9e87d90() {
   return (neuron0x9e28988()*0.247478);
}

double mlpANN::synapse0x9e87db8() {
   return (neuron0x9e28b08()*-0.104833);
}

double mlpANN::synapse0x9e87de0() {
   return (neuron0x9e28c80()*0.39331);
}

double mlpANN::synapse0x9e888a0() {
   return (neuron0x95a0f78()*0.599938);
}

double mlpANN::synapse0x9e888c8() {
   return (neuron0x931f328()*-0.26023);
}

double mlpANN::synapse0x9e888f0() {
   return (neuron0x9b48b50()*-0.674781);
}

double mlpANN::synapse0x9e88918() {
   return (neuron0x9b48cc8()*0.47904);
}

double mlpANN::synapse0x9e88940() {
   return (neuron0x9b48e40()*-0.0475278);
}

double mlpANN::synapse0x9e88968() {
   return (neuron0x9b48fb8()*-1.85906);
}

double mlpANN::synapse0x9e88990() {
   return (neuron0x9e282a0()*0.884859);
}

double mlpANN::synapse0x9e889b8() {
   return (neuron0x9e283d0()*1.19617);
}

double mlpANN::synapse0x9e889e0() {
   return (neuron0x9e28548()*-1.18885);
}

double mlpANN::synapse0x9e88a08() {
   return (neuron0x9e286c0()*-0.596585);
}

double mlpANN::synapse0x9e88a30() {
   return (neuron0x9e28838()*0.625185);
}

double mlpANN::synapse0x9e88a58() {
   return (neuron0x9e28988()*-0.927438);
}

double mlpANN::synapse0x9e88a80() {
   return (neuron0x9e28b08()*-0.404246);
}

double mlpANN::synapse0x9e88aa8() {
   return (neuron0x9e28c80()*-0.943553);
}

double mlpANN::synapse0x9e88c40() {
   return (neuron0x95a0f78()*-1.85163);
}

double mlpANN::synapse0x9e88c68() {
   return (neuron0x931f328()*1.17545);
}

double mlpANN::synapse0x9e88c90() {
   return (neuron0x9b48b50()*-1.64719);
}

double mlpANN::synapse0x9e88cb8() {
   return (neuron0x9b48cc8()*2.73542);
}

double mlpANN::synapse0x9e88ce0() {
   return (neuron0x9b48e40()*-0.00517145);
}

double mlpANN::synapse0x9e88d08() {
   return (neuron0x9b48fb8()*0.517336);
}

double mlpANN::synapse0x9e88d30() {
   return (neuron0x9e282a0()*-0.285551);
}

double mlpANN::synapse0x9e88d58() {
   return (neuron0x9e283d0()*-0.00915569);
}

double mlpANN::synapse0x9e88d80() {
   return (neuron0x9e28548()*0.672217);
}

double mlpANN::synapse0x9e88da8() {
   return (neuron0x9e286c0()*0.0802112);
}

double mlpANN::synapse0x9e88dd0() {
   return (neuron0x9e28838()*0.102955);
}

double mlpANN::synapse0x9e88df8() {
   return (neuron0x9e28988()*-1.6829);
}

double mlpANN::synapse0x9e88e20() {
   return (neuron0x9e28b08()*0.801354);
}

double mlpANN::synapse0x9e88e48() {
   return (neuron0x9e28c80()*-0.358761);
}

double mlpANN::synapse0x9e88fe0() {
   return (neuron0x95a0f78()*-0.0803179);
}

double mlpANN::synapse0x9e89008() {
   return (neuron0x931f328()*0.588231);
}

double mlpANN::synapse0x9e89030() {
   return (neuron0x9b48b50()*-0.952621);
}

double mlpANN::synapse0x9e89058() {
   return (neuron0x9b48cc8()*-0.323569);
}

double mlpANN::synapse0x9e89080() {
   return (neuron0x9b48e40()*-1.28158);
}

double mlpANN::synapse0x9e890a8() {
   return (neuron0x9b48fb8()*3.34688);
}

double mlpANN::synapse0x9e890d0() {
   return (neuron0x9e282a0()*-1.42211);
}

double mlpANN::synapse0x9e890f8() {
   return (neuron0x9e283d0()*1.23876);
}

double mlpANN::synapse0x9e89120() {
   return (neuron0x9e28548()*0.123066);
}

double mlpANN::synapse0x9e89148() {
   return (neuron0x9e286c0()*-0.378071);
}

double mlpANN::synapse0x9e89170() {
   return (neuron0x9e28838()*-0.182883);
}

double mlpANN::synapse0x9e89198() {
   return (neuron0x9e28988()*-2.24507);
}

double mlpANN::synapse0x9e891c0() {
   return (neuron0x9e28b08()*0.799389);
}

double mlpANN::synapse0x9e891e8() {
   return (neuron0x9e28c80()*0.303816);
}

double mlpANN::synapse0x9e89380() {
   return (neuron0x95a0f78()*-1.02521);
}

double mlpANN::synapse0x9e893a8() {
   return (neuron0x931f328()*-1.02596);
}

double mlpANN::synapse0x9e893d0() {
   return (neuron0x9b48b50()*-0.19172);
}

double mlpANN::synapse0x9e893f8() {
   return (neuron0x9b48cc8()*1.55273);
}

double mlpANN::synapse0x9e89420() {
   return (neuron0x9b48e40()*2.44123);
}

double mlpANN::synapse0x9e89448() {
   return (neuron0x9b48fb8()*-1.61579);
}

double mlpANN::synapse0x9e89470() {
   return (neuron0x9e282a0()*-0.682173);
}

double mlpANN::synapse0x9e89498() {
   return (neuron0x9e283d0()*0.711202);
}

double mlpANN::synapse0x9e894c0() {
   return (neuron0x9e28548()*-0.427465);
}

double mlpANN::synapse0x9e894e8() {
   return (neuron0x9e286c0()*0.769792);
}

double mlpANN::synapse0x9e89510() {
   return (neuron0x9e28838()*-1.84171);
}

double mlpANN::synapse0x9e89538() {
   return (neuron0x9e28988()*-0.305784);
}

double mlpANN::synapse0x9e89560() {
   return (neuron0x9e28b08()*0.442706);
}

double mlpANN::synapse0x9e89588() {
   return (neuron0x9e28c80()*-0.91867);
}

double mlpANN::synapse0x9e87e08() {
   return (neuron0x95a0f78()*1.42941);
}

double mlpANN::synapse0x9e87e30() {
   return (neuron0x931f328()*0.29531);
}

double mlpANN::synapse0x9aebc28() {
   return (neuron0x9b48b50()*-1.46342);
}

double mlpANN::synapse0x9e88528() {
   return (neuron0x9b48cc8()*0.173339);
}

double mlpANN::synapse0x9e88550() {
   return (neuron0x9b48e40()*0.0165939);
}

double mlpANN::synapse0x9e88578() {
   return (neuron0x9b48fb8()*-0.761793);
}

double mlpANN::synapse0x9e885a0() {
   return (neuron0x9e282a0()*0.758203);
}

double mlpANN::synapse0x9e885c8() {
   return (neuron0x9e283d0()*0.185506);
}

double mlpANN::synapse0x9e885f0() {
   return (neuron0x9e28548()*-0.498289);
}

double mlpANN::synapse0x9e88618() {
   return (neuron0x9e286c0()*0.231269);
}

double mlpANN::synapse0x9e88640() {
   return (neuron0x9e28838()*0.755853);
}

double mlpANN::synapse0x9e88668() {
   return (neuron0x9e28988()*1.54785);
}

double mlpANN::synapse0x9e88690() {
   return (neuron0x9e28b08()*-2.53174);
}

double mlpANN::synapse0x9e886b8() {
   return (neuron0x9e28c80()*1.55386);
}

double mlpANN::synapse0x9aec0a0() {
   return (neuron0x9e28ee8()*-1.13471);
}

double mlpANN::synapse0x9aec0c8() {
   return (neuron0x9e29268()*-0.241451);
}

double mlpANN::synapse0x9aec0f0() {
   return (neuron0x9e29690()*0.395077);
}

double mlpANN::synapse0x9aec118() {
   return (neuron0x9e87f38()*0.447769);
}

double mlpANN::synapse0x9aec140() {
   return (neuron0x9e88270()*0.850038);
}

double mlpANN::synapse0x9aec168() {
   return (neuron0x9e88730()*0.59626);
}

double mlpANN::synapse0x9aec190() {
   return (neuron0x9e88ad0()*0.693084);
}

double mlpANN::synapse0x9aec1b8() {
   return (neuron0x9e88e70()*0.417667);
}

double mlpANN::synapse0x9aec1e0() {
   return (neuron0x9e89210()*-0.480076);
}

double mlpANN::synapse0x9aec208() {
   return (neuron0x9e895b0()*-0.669786);
}

