#include "NueAna/SubShowerNN.h"
#include <cmath>

double SubShowerNN::Value(int index,double in0,double in1,double in2,double in3) {
   input0 = (in0 - 0)/1;
   input1 = (in1 - 0)/1;
   input2 = (in2 - 0)/1;
   input3 = (in3 - 0)/1;
   switch(index) {
     case 0:
         return neuron0x87e6f28();
     default:
         return 0.;
   }
}

double SubShowerNN::Value(int index, double* input) {
   input0 = (input[0] - 0)/1;
   input1 = (input[1] - 0)/1;
   input2 = (input[2] - 0)/1;
   input3 = (input[3] - 0)/1;
   switch(index) {
     case 0:
         return neuron0x87e6f28();
     default:
         return 0.;
   }
}

double SubShowerNN::neuron0x87ccb00() {
   return input0;
}

double SubShowerNN::neuron0x87e3780() {
   return input1;
}

double SubShowerNN::neuron0x87e3c30() {
   return input2;
}

double SubShowerNN::neuron0x87e40f8() {
   return input3;
}

double SubShowerNN::input0x87e3558() {
   double input = -2.51058;
   input += synapse0x87d8b98();
   input += synapse0x87e3758();
   input += synapse0x87e3e38();
   input += synapse0x87e42f8();
   return input;
}

double SubShowerNN::neuron0x87e3558() {
   double input = input0x87e3558();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e39c8() {
   double input = -0.276849;
   input += synapse0x87e3730();
   input += synapse0x87e39a0();
   input += synapse0x87e3c08();
   input += synapse0x87e4320();
   return input;
}

double SubShowerNN::neuron0x87e39c8() {
   double input = input0x87e39c8();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e3e88() {
   double input = -0.400797;
   input += synapse0x87e4080();
   input += synapse0x87e40a8();
   input += synapse0x87e3e60();
   input += synapse0x87e40d0();
   return input;
}

double SubShowerNN::neuron0x87e3e88() {
   double input = input0x87e3e88();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e4370() {
   double input = 2.16374;
   input += synapse0x87e45f0();
   input += synapse0x87e4618();
   input += synapse0x87e4640();
   input += synapse0x87e4348();
   return input;
}

double SubShowerNN::neuron0x87e4370() {
   double input = input0x87e4370();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e48e0() {
   double input = -0.157689;
   input += synapse0x87e4ad8();
   input += synapse0x87e4b00();
   input += synapse0x87e4b28();
   input += synapse0x87e4b50();
   return input;
}

double SubShowerNN::neuron0x87e48e0() {
   double input = input0x87e48e0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e4e60() {
   double input = -0.4857;
   input += synapse0x87e5058();
   input += synapse0x87e5080();
   input += synapse0x87e50a8();
   input += synapse0x87e50d0();
   return input;
}

double SubShowerNN::neuron0x87e4e60() {
   double input = input0x87e4e60();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e5430() {
   double input = 0.316371;
   input += synapse0x87e5628();
   input += synapse0x87e5650();
   input += synapse0x87e5678();
   input += synapse0x87e56a0();
   return input;
}

double SubShowerNN::neuron0x87e5430() {
   double input = input0x87e5430();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e5a50() {
   double input = -0.7789;
   input += synapse0x87e4568();
   input += synapse0x87e4590();
   input += synapse0x87e5d30();
   input += synapse0x87e5d58();
   return input;
}

double SubShowerNN::neuron0x87e5a50() {
   double input = input0x87e5a50();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e6158() {
   double input = 0.0102468;
   input += synapse0x87e6350();
   input += synapse0x87e6378();
   input += synapse0x87e63a0();
   input += synapse0x87e63c8();
   return input;
}

double SubShowerNN::neuron0x87e6158() {
   double input = input0x87e6158();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e6818() {
   double input = 0.888711;
   input += synapse0x87e6a10();
   input += synapse0x87e6a38();
   input += synapse0x87e6a60();
   input += synapse0x87e6a88();
   return input;
}

double SubShowerNN::neuron0x87e6818() {
   double input = input0x87e6818();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e4690() {
   double input = -1.42309;
   input += synapse0x87e4840();
   input += synapse0x87e4868();
   input += synapse0x87e4890();
   input += synapse0x87e4668();
   input += synapse0x87e48b8();
   input += synapse0x87e50f8();
   input += synapse0x87e56c8();
   input += synapse0x87e5d80();
   input += synapse0x87e63f0();
   input += synapse0x87e6ab0();
   return input;
}

double SubShowerNN::neuron0x87e4690() {
   double input = input0x87e4690();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e4ba0() {
   double input = 0.0430053;
   input += synapse0x87e4d98();
   input += synapse0x87e4dc0();
   input += synapse0x87e4de8();
   input += synapse0x87e4e10();
   input += synapse0x87e4b78();
   input += synapse0x87e4e38();
   input += synapse0x87e56f0();
   input += synapse0x87e5da8();
   input += synapse0x87e6418();
   input += synapse0x87e6ad8();
   return input;
}

double SubShowerNN::neuron0x87e4ba0() {
   double input = input0x87e4ba0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e5148() {
   double input = -0.633439;
   input += synapse0x87e5340();
   input += synapse0x87e5368();
   input += synapse0x87e5390();
   input += synapse0x87e53b8();
   input += synapse0x87e53e0();
   input += synapse0x87e5120();
   input += synapse0x87e5408();
   input += synapse0x87e5dd0();
   input += synapse0x87e6440();
   input += synapse0x87e6b00();
   return input;
}

double SubShowerNN::neuron0x87e5148() {
   double input = input0x87e5148();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e5740() {
   double input = -1.81575;
   input += synapse0x87e5938();
   input += synapse0x87e5960();
   input += synapse0x87e5988();
   input += synapse0x87e59b0();
   input += synapse0x87e59d8();
   input += synapse0x87e5a00();
   input += synapse0x87e5718();
   input += synapse0x87e5a28();
   input += synapse0x87e6468();
   input += synapse0x87e6b28();
   return input;
}

double SubShowerNN::neuron0x87e5740() {
   double input = input0x87e5740();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e5e20() {
   double input = -0.380953;
   input += synapse0x87e6018();
   input += synapse0x87e6040();
   input += synapse0x87e6068();
   input += synapse0x87e6090();
   input += synapse0x87e60b8();
   input += synapse0x87e60e0();
   input += synapse0x87e6108();
   input += synapse0x87e5df8();
   input += synapse0x87e6130();
   input += synapse0x87e6b50();
   return input;
}

double SubShowerNN::neuron0x87e5e20() {
   double input = input0x87e5e20();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e64b8() {
   double input = -0.054045;
   input += synapse0x87e66b0();
   input += synapse0x87e66d8();
   input += synapse0x87e6700();
   input += synapse0x87e6728();
   input += synapse0x87e6750();
   input += synapse0x87e6778();
   input += synapse0x87e67a0();
   input += synapse0x87e67c8();
   input += synapse0x87e6490();
   input += synapse0x87e67f0();
   return input;
}

double SubShowerNN::neuron0x87e64b8() {
   double input = input0x87e64b8();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e6ba0() {
   double input = 0.466166;
   input += synapse0x87e6d98();
   input += synapse0x87e6dc0();
   input += synapse0x87e6de8();
   input += synapse0x87e6e10();
   input += synapse0x87e6e38();
   input += synapse0x87e6e60();
   input += synapse0x87e6e88();
   input += synapse0x87e6eb0();
   input += synapse0x87e6ed8();
   input += synapse0x87e6b78();
   return input;
}

double SubShowerNN::neuron0x87e6ba0() {
   double input = input0x87e6ba0();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e7248() {
   double input = 0.247909;
   input += synapse0x87e7440();
   input += synapse0x87e7468();
   input += synapse0x87e7490();
   input += synapse0x87e74b8();
   input += synapse0x87e74e0();
   input += synapse0x87e7508();
   input += synapse0x87e7530();
   input += synapse0x87e7558();
   input += synapse0x87e7580();
   input += synapse0x87e75a8();
   return input;
}

double SubShowerNN::neuron0x87e7248() {
   double input = input0x87e7248();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e75f8() {
   double input = -0.353472;
   input += synapse0x87e77f0();
   input += synapse0x87e7818();
   input += synapse0x87e7840();
   input += synapse0x87e7868();
   input += synapse0x87e7890();
   input += synapse0x87e78b8();
   input += synapse0x87e78e0();
   input += synapse0x87e7908();
   input += synapse0x87e7930();
   input += synapse0x87e7958();
   return input;
}

double SubShowerNN::neuron0x87e75f8() {
   double input = input0x87e75f8();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e79a8() {
   double input = -0.364511;
   input += synapse0x87e7ba0();
   input += synapse0x87e7bc8();
   input += synapse0x87e7bf0();
   input += synapse0x87e7c18();
   input += synapse0x87e7c40();
   input += synapse0x87e7c68();
   input += synapse0x87e7c90();
   input += synapse0x87e7cb8();
   input += synapse0x87e7ce0();
   input += synapse0x87e7d08();
   return input;
}

double SubShowerNN::neuron0x87e79a8() {
   double input = input0x87e79a8();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double SubShowerNN::input0x87e6f28() {
   double input = -1.27279;
   input += synapse0x87e7130();
   input += synapse0x87e7158();
   input += synapse0x87e7180();
   input += synapse0x87e71a8();
   input += synapse0x87e71d0();
   input += synapse0x87e71f8();
   input += synapse0x87e6f00();
   input += synapse0x87e7220();
   input += synapse0x87e75d0();
   input += synapse0x87e7980();
   return input;
}

double SubShowerNN::neuron0x87e6f28() {
   double input = input0x87e6f28();
   return (input * 1)+0;
}

double SubShowerNN::synapse0x87d8b98() {
   return (neuron0x87ccb00()*-1.51964);
}

double SubShowerNN::synapse0x87e3758() {
   return (neuron0x87e3780()*-0.764593);
}

double SubShowerNN::synapse0x87e3e38() {
   return (neuron0x87e3c30()*6.05807);
}

double SubShowerNN::synapse0x87e42f8() {
   return (neuron0x87e40f8()*2.98751);
}

double SubShowerNN::synapse0x87e3730() {
   return (neuron0x87ccb00()*0.580341);
}

double SubShowerNN::synapse0x87e39a0() {
   return (neuron0x87e3780()*0.253629);
}

double SubShowerNN::synapse0x87e3c08() {
   return (neuron0x87e3c30()*0.962705);
}

double SubShowerNN::synapse0x87e4320() {
   return (neuron0x87e40f8()*-1.36657);
}

double SubShowerNN::synapse0x87e4080() {
   return (neuron0x87ccb00()*-2.02742);
}

double SubShowerNN::synapse0x87e40a8() {
   return (neuron0x87e3780()*1.82533);
}

double SubShowerNN::synapse0x87e3e60() {
   return (neuron0x87e3c30()*0.480191);
}

double SubShowerNN::synapse0x87e40d0() {
   return (neuron0x87e40f8()*-5.90924);
}

double SubShowerNN::synapse0x87e45f0() {
   return (neuron0x87ccb00()*-5.11713);
}

double SubShowerNN::synapse0x87e4618() {
   return (neuron0x87e3780()*1.76387);
}

double SubShowerNN::synapse0x87e4640() {
   return (neuron0x87e3c30()*0.699597);
}

double SubShowerNN::synapse0x87e4348() {
   return (neuron0x87e40f8()*5.47703);
}

double SubShowerNN::synapse0x87e4ad8() {
   return (neuron0x87ccb00()*-1.39427);
}

double SubShowerNN::synapse0x87e4b00() {
   return (neuron0x87e3780()*-1.33337);
}

double SubShowerNN::synapse0x87e4b28() {
   return (neuron0x87e3c30()*2.65988);
}

double SubShowerNN::synapse0x87e4b50() {
   return (neuron0x87e40f8()*-0.864033);
}

double SubShowerNN::synapse0x87e5058() {
   return (neuron0x87ccb00()*2.42384);
}

double SubShowerNN::synapse0x87e5080() {
   return (neuron0x87e3780()*1.59382);
}

double SubShowerNN::synapse0x87e50a8() {
   return (neuron0x87e3c30()*-7.0827);
}

double SubShowerNN::synapse0x87e50d0() {
   return (neuron0x87e40f8()*-1.62011);
}

double SubShowerNN::synapse0x87e5628() {
   return (neuron0x87ccb00()*-3.64683);
}

double SubShowerNN::synapse0x87e5650() {
   return (neuron0x87e3780()*-0.111531);
}

double SubShowerNN::synapse0x87e5678() {
   return (neuron0x87e3c30()*-5.87413);
}

double SubShowerNN::synapse0x87e56a0() {
   return (neuron0x87e40f8()*-0.888458);
}

double SubShowerNN::synapse0x87e4568() {
   return (neuron0x87ccb00()*3.70138);
}

double SubShowerNN::synapse0x87e4590() {
   return (neuron0x87e3780()*1.56974);
}

double SubShowerNN::synapse0x87e5d30() {
   return (neuron0x87e3c30()*-5.44227);
}

double SubShowerNN::synapse0x87e5d58() {
   return (neuron0x87e40f8()*1.24221);
}

double SubShowerNN::synapse0x87e6350() {
   return (neuron0x87ccb00()*-1.12404);
}

double SubShowerNN::synapse0x87e6378() {
   return (neuron0x87e3780()*0.581114);
}

double SubShowerNN::synapse0x87e63a0() {
   return (neuron0x87e3c30()*-0.0471626);
}

double SubShowerNN::synapse0x87e63c8() {
   return (neuron0x87e40f8()*-4.22058);
}

double SubShowerNN::synapse0x87e6a10() {
   return (neuron0x87ccb00()*0.329883);
}

double SubShowerNN::synapse0x87e6a38() {
   return (neuron0x87e3780()*2.01344);
}

double SubShowerNN::synapse0x87e6a60() {
   return (neuron0x87e3c30()*-1.29032);
}

double SubShowerNN::synapse0x87e6a88() {
   return (neuron0x87e40f8()*-1.46909);
}

double SubShowerNN::synapse0x87e4840() {
   return (neuron0x87e3558()*-1.53046);
}

double SubShowerNN::synapse0x87e4868() {
   return (neuron0x87e39c8()*-0.636854);
}

double SubShowerNN::synapse0x87e4890() {
   return (neuron0x87e3e88()*0.867331);
}

double SubShowerNN::synapse0x87e4668() {
   return (neuron0x87e4370()*-1.27728);
}

double SubShowerNN::synapse0x87e48b8() {
   return (neuron0x87e48e0()*0.0128366);
}

double SubShowerNN::synapse0x87e50f8() {
   return (neuron0x87e4e60()*-0.761449);
}

double SubShowerNN::synapse0x87e56c8() {
   return (neuron0x87e5430()*0.510321);
}

double SubShowerNN::synapse0x87e5d80() {
   return (neuron0x87e5a50()*-1.98975);
}

double SubShowerNN::synapse0x87e63f0() {
   return (neuron0x87e6158()*-0.716214);
}

double SubShowerNN::synapse0x87e6ab0() {
   return (neuron0x87e6818()*-0.589881);
}

double SubShowerNN::synapse0x87e4d98() {
   return (neuron0x87e3558()*0.174805);
}

double SubShowerNN::synapse0x87e4dc0() {
   return (neuron0x87e39c8()*0.609518);
}

double SubShowerNN::synapse0x87e4de8() {
   return (neuron0x87e3e88()*-0.227943);
}

double SubShowerNN::synapse0x87e4e10() {
   return (neuron0x87e4370()*0.259165);
}

double SubShowerNN::synapse0x87e4b78() {
   return (neuron0x87e48e0()*-0.896176);
}

double SubShowerNN::synapse0x87e4e38() {
   return (neuron0x87e4e60()*-0.0356925);
}

double SubShowerNN::synapse0x87e56f0() {
   return (neuron0x87e5430()*-0.0708626);
}

double SubShowerNN::synapse0x87e5da8() {
   return (neuron0x87e5a50()*0.816397);
}

double SubShowerNN::synapse0x87e6418() {
   return (neuron0x87e6158()*-0.127855);
}

double SubShowerNN::synapse0x87e6ad8() {
   return (neuron0x87e6818()*0.703495);
}

double SubShowerNN::synapse0x87e5340() {
   return (neuron0x87e3558()*-1.33098);
}

double SubShowerNN::synapse0x87e5368() {
   return (neuron0x87e39c8()*-0.362864);
}

double SubShowerNN::synapse0x87e5390() {
   return (neuron0x87e3e88()*2.40892);
}

double SubShowerNN::synapse0x87e53b8() {
   return (neuron0x87e4370()*-2.23644);
}

double SubShowerNN::synapse0x87e53e0() {
   return (neuron0x87e48e0()*-0.295975);
}

double SubShowerNN::synapse0x87e5120() {
   return (neuron0x87e4e60()*1.26264);
}

double SubShowerNN::synapse0x87e5408() {
   return (neuron0x87e5430()*1.15443);
}

double SubShowerNN::synapse0x87e5dd0() {
   return (neuron0x87e5a50()*-0.907784);
}

double SubShowerNN::synapse0x87e6440() {
   return (neuron0x87e6158()*1.70658);
}

double SubShowerNN::synapse0x87e6b00() {
   return (neuron0x87e6818()*1.10906);
}

double SubShowerNN::synapse0x87e5938() {
   return (neuron0x87e3558()*-2.55824);
}

double SubShowerNN::synapse0x87e5960() {
   return (neuron0x87e39c8()*0.6301);
}

double SubShowerNN::synapse0x87e5988() {
   return (neuron0x87e3e88()*5.41036);
}

double SubShowerNN::synapse0x87e59b0() {
   return (neuron0x87e4370()*-3.86037);
}

double SubShowerNN::synapse0x87e59d8() {
   return (neuron0x87e48e0()*1.94823);
}

double SubShowerNN::synapse0x87e5a00() {
   return (neuron0x87e4e60()*1.74082);
}

double SubShowerNN::synapse0x87e5718() {
   return (neuron0x87e5430()*2.92241);
}

double SubShowerNN::synapse0x87e5a28() {
   return (neuron0x87e5a50()*-1.84613);
}

double SubShowerNN::synapse0x87e6468() {
   return (neuron0x87e6158()*3.07201);
}

double SubShowerNN::synapse0x87e6b28() {
   return (neuron0x87e6818()*1.64651);
}

double SubShowerNN::synapse0x87e6018() {
   return (neuron0x87e3558()*-3.89571);
}

double SubShowerNN::synapse0x87e6040() {
   return (neuron0x87e39c8()*-0.383293);
}

double SubShowerNN::synapse0x87e6068() {
   return (neuron0x87e3e88()*3.41103);
}

double SubShowerNN::synapse0x87e6090() {
   return (neuron0x87e4370()*-2.20505);
}

double SubShowerNN::synapse0x87e60b8() {
   return (neuron0x87e48e0()*-0.274721);
}

double SubShowerNN::synapse0x87e60e0() {
   return (neuron0x87e4e60()*0.407058);
}

double SubShowerNN::synapse0x87e6108() {
   return (neuron0x87e5430()*0.820336);
}

double SubShowerNN::synapse0x87e5df8() {
   return (neuron0x87e5a50()*-0.217564);
}

double SubShowerNN::synapse0x87e6130() {
   return (neuron0x87e6158()*0.451884);
}

double SubShowerNN::synapse0x87e6b50() {
   return (neuron0x87e6818()*0.797602);
}

double SubShowerNN::synapse0x87e66b0() {
   return (neuron0x87e3558()*-1.84465);
}

double SubShowerNN::synapse0x87e66d8() {
   return (neuron0x87e39c8()*0.0325271);
}

double SubShowerNN::synapse0x87e6700() {
   return (neuron0x87e3e88()*2.50346);
}

double SubShowerNN::synapse0x87e6728() {
   return (neuron0x87e4370()*-1.6172);
}

double SubShowerNN::synapse0x87e6750() {
   return (neuron0x87e48e0()*0.0241894);
}

double SubShowerNN::synapse0x87e6778() {
   return (neuron0x87e4e60()*0.19091);
}

double SubShowerNN::synapse0x87e67a0() {
   return (neuron0x87e5430()*0.624186);
}

double SubShowerNN::synapse0x87e67c8() {
   return (neuron0x87e5a50()*-1.21348);
}

double SubShowerNN::synapse0x87e6490() {
   return (neuron0x87e6158()*1.38983);
}

double SubShowerNN::synapse0x87e67f0() {
   return (neuron0x87e6818()*1.35508);
}

double SubShowerNN::synapse0x87e6d98() {
   return (neuron0x87e3558()*-0.920274);
}

double SubShowerNN::synapse0x87e6dc0() {
   return (neuron0x87e39c8()*-0.58241);
}

double SubShowerNN::synapse0x87e6de8() {
   return (neuron0x87e3e88()*1.53953);
}

double SubShowerNN::synapse0x87e6e10() {
   return (neuron0x87e4370()*-0.457612);
}

double SubShowerNN::synapse0x87e6e38() {
   return (neuron0x87e48e0()*-0.268216);
}

double SubShowerNN::synapse0x87e6e60() {
   return (neuron0x87e4e60()*0.536563);
}

double SubShowerNN::synapse0x87e6e88() {
   return (neuron0x87e5430()*0.0858709);
}

double SubShowerNN::synapse0x87e6eb0() {
   return (neuron0x87e5a50()*-1.03249);
}

double SubShowerNN::synapse0x87e6ed8() {
   return (neuron0x87e6158()*1.14134);
}

double SubShowerNN::synapse0x87e6b78() {
   return (neuron0x87e6818()*0.162443);
}

double SubShowerNN::synapse0x87e7440() {
   return (neuron0x87e3558()*-1.44036);
}

double SubShowerNN::synapse0x87e7468() {
   return (neuron0x87e39c8()*-0.233911);
}

double SubShowerNN::synapse0x87e7490() {
   return (neuron0x87e3e88()*2.10351);
}

double SubShowerNN::synapse0x87e74b8() {
   return (neuron0x87e4370()*-1.49264);
}

double SubShowerNN::synapse0x87e74e0() {
   return (neuron0x87e48e0()*-0.0421463);
}

double SubShowerNN::synapse0x87e7508() {
   return (neuron0x87e4e60()*0.485113);
}

double SubShowerNN::synapse0x87e7530() {
   return (neuron0x87e5430()*0.0572048);
}

double SubShowerNN::synapse0x87e7558() {
   return (neuron0x87e5a50()*-1.32832);
}

double SubShowerNN::synapse0x87e7580() {
   return (neuron0x87e6158()*2.27081);
}

double SubShowerNN::synapse0x87e75a8() {
   return (neuron0x87e6818()*0.618502);
}

double SubShowerNN::synapse0x87e77f0() {
   return (neuron0x87e3558()*1.94999);
}

double SubShowerNN::synapse0x87e7818() {
   return (neuron0x87e39c8()*0.532869);
}

double SubShowerNN::synapse0x87e7840() {
   return (neuron0x87e3e88()*-2.03745);
}

double SubShowerNN::synapse0x87e7868() {
   return (neuron0x87e4370()*1.47153);
}

double SubShowerNN::synapse0x87e7890() {
   return (neuron0x87e48e0()*-0.0574973);
}

double SubShowerNN::synapse0x87e78b8() {
   return (neuron0x87e4e60()*-1.06016);
}

double SubShowerNN::synapse0x87e78e0() {
   return (neuron0x87e5430()*-1.36233);
}

double SubShowerNN::synapse0x87e7908() {
   return (neuron0x87e5a50()*1.34734);
}

double SubShowerNN::synapse0x87e7930() {
   return (neuron0x87e6158()*-1.07785);
}

double SubShowerNN::synapse0x87e7958() {
   return (neuron0x87e6818()*-1.39353);
}

double SubShowerNN::synapse0x87e7ba0() {
   return (neuron0x87e3558()*0.342644);
}

double SubShowerNN::synapse0x87e7bc8() {
   return (neuron0x87e39c8()*-0.0867636);
}

double SubShowerNN::synapse0x87e7bf0() {
   return (neuron0x87e3e88()*0.0258721);
}

double SubShowerNN::synapse0x87e7c18() {
   return (neuron0x87e4370()*0.510844);
}

double SubShowerNN::synapse0x87e7c40() {
   return (neuron0x87e48e0()*-0.472341);
}

double SubShowerNN::synapse0x87e7c68() {
   return (neuron0x87e4e60()*-0.369009);
}

double SubShowerNN::synapse0x87e7c90() {
   return (neuron0x87e5430()*-0.345007);
}

double SubShowerNN::synapse0x87e7cb8() {
   return (neuron0x87e5a50()*1.53966);
}

double SubShowerNN::synapse0x87e7ce0() {
   return (neuron0x87e6158()*-1.38018);
}

double SubShowerNN::synapse0x87e7d08() {
   return (neuron0x87e6818()*0.729445);
}

double SubShowerNN::synapse0x87e7130() {
   return (neuron0x87e4690()*2.37685);
}

double SubShowerNN::synapse0x87e7158() {
   return (neuron0x87e4ba0()*1.05212);
}

double SubShowerNN::synapse0x87e7180() {
   return (neuron0x87e5148()*0.244162);
}

double SubShowerNN::synapse0x87e71a8() {
   return (neuron0x87e5740()*3.53862);
}

double SubShowerNN::synapse0x87e71d0() {
   return (neuron0x87e5e20()*-2.01023);
}

double SubShowerNN::synapse0x87e71f8() {
   return (neuron0x87e64b8()*-0.708315);
}

double SubShowerNN::synapse0x87e6f00() {
   return (neuron0x87e6ba0()*-0.484825);
}

double SubShowerNN::synapse0x87e7220() {
   return (neuron0x87e7248()*-1.42554);
}

double SubShowerNN::synapse0x87e75d0() {
   return (neuron0x87e75f8()*0.227679);
}

double SubShowerNN::synapse0x87e7980() {
   return (neuron0x87e79a8()*1.75486);
}

