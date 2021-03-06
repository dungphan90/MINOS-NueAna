#ifndef ANNANA
#define ANNANA

#include "TObject.h"

#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "NueAna/NueAnaBase.h"
#include "TMultiLayerPerceptron.h"

class NueRecord;

class AnnAna : public TObject{ 

public:
   AnnAna(NueRecord &nr);
   ~AnnAna();

   void Analyze();

   double value(int index,double in0,double in1,double in2,double in3,double in4,double in5,double in6,double in7,double in8,double in9,double in10,double in11,double in12,double in13,double in14,double in15,double in16,double in17,double in18,double in19,double in20,double in21,double in22,double in23,double in24,double in25,double in26,double in27,double in28,double in29);
private:
   double input0;
   double input1;
   double input2;
   double input3;
   double input4;
   double input5;
   double input6;
   double input7;
   double input8;
   double input9;
   double input10;
   double input11;
   double input12;
   double input13;
   double input14;
   double input15;
   double input16;
   double input17;
   double input18;
   double input19;
   double input20;
   double input21;
   double input22;
   double input23;
   double input24;
   double input25;
   double input26;
   double input27;
   double input28;
   double input29;
   double neuron0xa0253a8();
   double neuron0xa025510();
   double neuron0xa025678();
   double neuron0xa0257e0();
   double neuron0xa025948();
   double neuron0xa025ab8();
   double neuron0xa025c20();
   double neuron0xa025d88();
   double neuron0xa025ef0();
   double neuron0xa026058();
   double neuron0xa0261c0();
   double neuron0xa026328();
   double neuron0xa026490();
   double neuron0xa0265f8();
   double neuron0xa026760();
   double neuron0xa0268c8();
   double neuron0xa026a30();
   double neuron0xa026ca8();
   double neuron0xa026d80();
   double neuron0xa026ee8();
   double neuron0xa027050();
   double neuron0xa0271b8();
   double neuron0xa027320();
   double neuron0xa027488();
   double neuron0xa0275f0();
   double neuron0xa027758();
   double neuron0xa0278c0();
   double neuron0xa0064a0();
   double neuron0xa006608();
   double neuron0xa006770();
   double neuron0xa0069f0();
   double neuron0xa007100();
   double neuron0xa007708();
   double neuron0xa007f10();
   double neuron0xa008570();
   double neuron0xa008de8();
   double neuron0xa009448();
   double neuron0xa009aa8();
   double neuron0xa00a108();
   double neuron0xa0089e0();
   double neuron0xa00b1d8();
   double neuron0xa00b838();
   double neuron0xa00be98();
   double neuron0xa00c4f8();
   double neuron0xa00cb58();
   double neuron0xa00d1d8();
   double neuron0xa00d538();
   double neuron0xa00d900();
   double neuron0xa00dcc8();
   double neuron0xa00e090();
   double neuron0xa00a740();
   double neuron0xa00ab08();
   double neuron0xa00f3e8();
   double neuron0xa00f7b0();
   double neuron0xa00fb78();
   double neuron0xa00ff40();
   double neuron0xa010308();
   double neuron0xa0106d0();
   double neuron0xa010a98();
   double neuron0xa010e60();
   double neuron0xa0068f8();
   double synapse0xa006b40();
   double synapse0xa006b68();
   double synapse0xa006b90();
   double synapse0xa006bb8();
   double synapse0xa006be0();
   double synapse0xa006c08();
   double synapse0xa006c30();
   double synapse0xa006c58();
   double synapse0xa006c80();
   double synapse0xa006ca8();
   double synapse0xa006cd0();
   double synapse0xa006cf8();
   double synapse0xa006d20();
   double synapse0xa006d48();
   double synapse0xa006d70();
   double synapse0xa006d98();
   double synapse0xa006dc0();
   double synapse0xa006ef8();
   double synapse0xa006f20();
   double synapse0xa006f48();
   double synapse0xa006f70();
   double synapse0xa006f98();
   double synapse0xa006fc0();
   double synapse0xa006fe8();
   double synapse0xa007010();
   double synapse0xa007038();
   double synapse0xa007060();
   double synapse0xa007088();
   double synapse0xa0070b0();
   double synapse0xa0070d8();
   double synapse0xa0071e0();
   double synapse0xa007208();
   double synapse0xa007230();
   double synapse0xa0279e0();
   double synapse0xa027a08();
   double synapse0xa0315d0();
   double synapse0xa0315f8();
   double synapse0xa006e70();
   double synapse0xa006e98();
   double synapse0xa006ec0();
   double synapse0xa007360();
   double synapse0xa007388();
   double synapse0xa0073b0();
   double synapse0xa0073d8();
   double synapse0xa007400();
   double synapse0xa007428();
   double synapse0xa007450();
   double synapse0xa007500();
   double synapse0xa007528();
   double synapse0xa007550();
   double synapse0xa007578();
   double synapse0xa0075a0();
   double synapse0xa0075c8();
   double synapse0xa0075f0();
   double synapse0xa007618();
   double synapse0xa007640();
   double synapse0xa007668();
   double synapse0xa007690();
   double synapse0xa0076b8();
   double synapse0xa0076e0();
   double synapse0xa026c20();
   double synapse0xa026c48();
   double synapse0xa026c70();
   double synapse0xa007938();
   double synapse0xa007960();
   double synapse0xa007258();
   double synapse0xa007280();
   double synapse0xa0072a8();
   double synapse0xa0072d0();
   double synapse0xa0072f8();
   double synapse0xa007320();
   double synapse0xa007b90();
   double synapse0xa007bb8();
   double synapse0xa007be0();
   double synapse0xa007c08();
   double synapse0xa007c30();
   double synapse0xa007c58();
   double synapse0xa007d08();
   double synapse0xa007d30();
   double synapse0xa007d58();
   double synapse0xa007d80();
   double synapse0xa007da8();
   double synapse0xa007dd0();
   double synapse0xa007df8();
   double synapse0xa007e20();
   double synapse0xa007e48();
   double synapse0xa007e70();
   double synapse0xa007e98();
   double synapse0xa007ec0();
   double synapse0xa007ee8();
   double synapse0xa008038();
   double synapse0xa008060();
   double synapse0xa008088();
   double synapse0xa0080b0();
   double synapse0xa0080d8();
   double synapse0xa008100();
   double synapse0xa008128();
   double synapse0xa008150();
   double synapse0xa008178();
   double synapse0xa0081a0();
   double synapse0xa0081c8();
   double synapse0xa0081f0();
   double synapse0xa008218();
   double synapse0xa008240();
   double synapse0xa008268();
   double synapse0xa008290();
   double synapse0xa0082b8();
   double synapse0xa008368();
   double synapse0xa008390();
   double synapse0xa0083b8();
   double synapse0xa0083e0();
   double synapse0xa008408();
   double synapse0xa008430();
   double synapse0xa008458();
   double synapse0xa008480();
   double synapse0xa0084a8();
   double synapse0xa0084d0();
   double synapse0xa0084f8();
   double synapse0xa008520();
   double synapse0xa008548();
   double synapse0xa008698();
   double synapse0xa0086c0();
   double synapse0xa0086e8();
   double synapse0xa008710();
   double synapse0xa008738();
   double synapse0xa008760();
   double synapse0xa008788();
   double synapse0xa0087b0();
   double synapse0xa0087d8();
   double synapse0xa007988();
   double synapse0xa0079b0();
   double synapse0xa0079d8();
   double synapse0xa007a00();
   double synapse0xa007a28();
   double synapse0xa007a50();
   double synapse0xa007a78();
   double synapse0xa007aa0();
   double synapse0xa007b50();
   double synapse0xa008c08();
   double synapse0xa008c30();
   double synapse0xa008c58();
   double synapse0xa008c80();
   double synapse0xa008ca8();
   double synapse0xa008cd0();
   double synapse0xa008cf8();
   double synapse0xa008d20();
   double synapse0xa008d48();
   double synapse0xa008d70();
   double synapse0xa008d98();
   double synapse0xa008dc0();
   double synapse0xa008f10();
   double synapse0xa008f38();
   double synapse0xa008f60();
   double synapse0xa008f88();
   double synapse0xa008fb0();
   double synapse0xa008fd8();
   double synapse0xa009000();
   double synapse0xa009028();
   double synapse0xa009050();
   double synapse0xa009078();
   double synapse0xa0090a0();
   double synapse0xa0090c8();
   double synapse0xa0090f0();
   double synapse0xa009118();
   double synapse0xa009140();
   double synapse0xa009168();
   double synapse0xa009190();
   double synapse0xa009240();
   double synapse0xa009268();
   double synapse0xa009290();
   double synapse0xa0092b8();
   double synapse0xa0092e0();
   double synapse0xa009308();
   double synapse0xa009330();
   double synapse0xa009358();
   double synapse0xa009380();
   double synapse0xa0093a8();
   double synapse0xa0093d0();
   double synapse0xa0093f8();
   double synapse0xa009420();
   double synapse0xa009570();
   double synapse0xa009598();
   double synapse0xa0095c0();
   double synapse0xa0095e8();
   double synapse0xa009610();
   double synapse0xa009638();
   double synapse0xa009660();
   double synapse0xa009688();
   double synapse0xa0096b0();
   double synapse0xa0096d8();
   double synapse0xa009700();
   double synapse0xa009728();
   double synapse0xa009750();
   double synapse0xa009778();
   double synapse0xa0097a0();
   double synapse0xa0097c8();
   double synapse0xa0097f0();
   double synapse0xa0098a0();
   double synapse0xa0098c8();
   double synapse0xa0098f0();
   double synapse0xa009918();
   double synapse0xa009940();
   double synapse0xa009968();
   double synapse0xa009990();
   double synapse0xa0099b8();
   double synapse0xa0099e0();
   double synapse0xa009a08();
   double synapse0xa009a30();
   double synapse0xa009a58();
   double synapse0xa009a80();
   double synapse0xa009bd0();
   double synapse0xa009bf8();
   double synapse0xa009c20();
   double synapse0xa009c48();
   double synapse0xa009c70();
   double synapse0xa009c98();
   double synapse0xa009cc0();
   double synapse0xa009ce8();
   double synapse0xa009d10();
   double synapse0xa009d38();
   double synapse0xa009d60();
   double synapse0xa009d88();
   double synapse0xa009db0();
   double synapse0xa009dd8();
   double synapse0xa009e00();
   double synapse0xa009e28();
   double synapse0xa009e50();
   double synapse0xa009f00();
   double synapse0xa009f28();
   double synapse0xa009f50();
   double synapse0xa009f78();
   double synapse0xa009fa0();
   double synapse0xa009fc8();
   double synapse0xa009ff0();
   double synapse0xa00a018();
   double synapse0xa00a040();
   double synapse0xa00a068();
   double synapse0xa00a090();
   double synapse0xa00a0b8();
   double synapse0xa00a0e0();
   double synapse0xa00a230();
   double synapse0xa00a258();
   double synapse0xa00a280();
   double synapse0xa00a2a8();
   double synapse0xa00a2d0();
   double synapse0xa00a2f8();
   double synapse0xa00a320();
   double synapse0xa00a348();
   double synapse0xa00a370();
   double synapse0xa00a398();
   double synapse0xa00a3c0();
   double synapse0xa00a3e8();
   double synapse0xa00a410();
   double synapse0xa00a438();
   double synapse0xa00a460();
   double synapse0xa00a488();
   double synapse0xa00a4b0();
   double synapse0xa024cf8();
   double synapse0xa008800();
   double synapse0xa008828();
   double synapse0xa008850();
   double synapse0xa008878();
   double synapse0xa0088a0();
   double synapse0xa0088c8();
   double synapse0xa0088f0();
   double synapse0xa008918();
   double synapse0xa008940();
   double synapse0xa008968();
   double synapse0xa008990();
   double synapse0xa0089b8();
   double synapse0xa008b30();
   double synapse0xa008b58();
   double synapse0xa008b80();
   double synapse0xa008ba8();
   double synapse0xa008bd0();
   double synapse0xa00ad68();
   double synapse0xa00ad90();
   double synapse0xa00adb8();
   double synapse0xa00ade0();
   double synapse0xa00ae08();
   double synapse0xa00ae30();
   double synapse0xa00ae58();
   double synapse0xa00ae80();
   double synapse0xa00aea8();
   double synapse0xa00aed0();
   double synapse0xa00aef8();
   double synapse0xa00af20();
   double synapse0xa00afd0();
   double synapse0xa00aff8();
   double synapse0xa00b020();
   double synapse0xa00b048();
   double synapse0xa00b070();
   double synapse0xa00b098();
   double synapse0xa00b0c0();
   double synapse0xa00b0e8();
   double synapse0xa00b110();
   double synapse0xa00b138();
   double synapse0xa00b160();
   double synapse0xa00b188();
   double synapse0xa00b1b0();
   double synapse0xa00b300();
   double synapse0xa00b328();
   double synapse0xa00b350();
   double synapse0xa00b378();
   double synapse0xa00b3a0();
   double synapse0xa00b3c8();
   double synapse0xa00b3f0();
   double synapse0xa00b418();
   double synapse0xa00b440();
   double synapse0xa00b468();
   double synapse0xa00b490();
   double synapse0xa00b4b8();
   double synapse0xa00b4e0();
   double synapse0xa00b508();
   double synapse0xa00b530();
   double synapse0xa00b558();
   double synapse0xa00b580();
   double synapse0xa00b630();
   double synapse0xa00b658();
   double synapse0xa00b680();
   double synapse0xa00b6a8();
   double synapse0xa00b6d0();
   double synapse0xa00b6f8();
   double synapse0xa00b720();
   double synapse0xa00b748();
   double synapse0xa00b770();
   double synapse0xa00b798();
   double synapse0xa00b7c0();
   double synapse0xa00b7e8();
   double synapse0xa00b810();
   double synapse0xa00b960();
   double synapse0xa00b988();
   double synapse0xa00b9b0();
   double synapse0xa00b9d8();
   double synapse0xa00ba00();
   double synapse0xa00ba28();
   double synapse0xa00ba50();
   double synapse0xa00ba78();
   double synapse0xa00baa0();
   double synapse0xa00bac8();
   double synapse0xa00baf0();
   double synapse0xa00bb18();
   double synapse0xa00bb40();
   double synapse0xa00bb68();
   double synapse0xa00bb90();
   double synapse0xa00bbb8();
   double synapse0xa00bbe0();
   double synapse0xa00bc90();
   double synapse0xa00bcb8();
   double synapse0xa00bce0();
   double synapse0xa00bd08();
   double synapse0xa00bd30();
   double synapse0xa00bd58();
   double synapse0xa00bd80();
   double synapse0xa00bda8();
   double synapse0xa00bdd0();
   double synapse0xa00bdf8();
   double synapse0xa00be20();
   double synapse0xa00be48();
   double synapse0xa00be70();
   double synapse0xa00bfc0();
   double synapse0xa00bfe8();
   double synapse0xa00c010();
   double synapse0xa00c038();
   double synapse0xa00c060();
   double synapse0xa00c088();
   double synapse0xa00c0b0();
   double synapse0xa00c0d8();
   double synapse0xa00c100();
   double synapse0xa00c128();
   double synapse0xa00c150();
   double synapse0xa00c178();
   double synapse0xa00c1a0();
   double synapse0xa00c1c8();
   double synapse0xa00c1f0();
   double synapse0xa00c218();
   double synapse0xa00c240();
   double synapse0xa00c2f0();
   double synapse0xa00c318();
   double synapse0xa00c340();
   double synapse0xa00c368();
   double synapse0xa00c390();
   double synapse0xa00c3b8();
   double synapse0xa00c3e0();
   double synapse0xa00c408();
   double synapse0xa00c430();
   double synapse0xa00c458();
   double synapse0xa00c480();
   double synapse0xa00c4a8();
   double synapse0xa00c4d0();
   double synapse0xa00c620();
   double synapse0xa00c648();
   double synapse0xa00c670();
   double synapse0xa00c698();
   double synapse0xa00c6c0();
   double synapse0xa00c6e8();
   double synapse0xa00c710();
   double synapse0xa00c738();
   double synapse0xa00c760();
   double synapse0xa00c788();
   double synapse0xa00c7b0();
   double synapse0xa00c7d8();
   double synapse0xa00c800();
   double synapse0xa00c828();
   double synapse0xa00c850();
   double synapse0xa00c878();
   double synapse0xa00c8a0();
   double synapse0xa00c950();
   double synapse0xa00c978();
   double synapse0xa00c9a0();
   double synapse0xa00c9c8();
   double synapse0xa00c9f0();
   double synapse0xa00ca18();
   double synapse0xa00ca40();
   double synapse0xa00ca68();
   double synapse0xa00ca90();
   double synapse0xa00cab8();
   double synapse0xa00cae0();
   double synapse0xa00cb08();
   double synapse0xa00cb30();
   double synapse0xa00cc80();
   double synapse0xa00cca8();
   double synapse0xa00ccd0();
   double synapse0xa00ccf8();
   double synapse0xa00cd20();
   double synapse0xa00cd48();
   double synapse0xa00cd70();
   double synapse0xa00cd98();
   double synapse0xa00cdc0();
   double synapse0xa00cde8();
   double synapse0xa00ce10();
   double synapse0xa00ce38();
   double synapse0xa00ce60();
   double synapse0xa00ce88();
   double synapse0xa00ceb0();
   double synapse0xa00ced8();
   double synapse0xa00cf00();
   double synapse0xa00cfb0();
   double synapse0xa00cfd8();
   double synapse0xa00d000();
   double synapse0xa00d028();
   double synapse0xa00d050();
   double synapse0xa00d078();
   double synapse0xa00d0a0();
   double synapse0xa00d0c8();
   double synapse0xa00d0f0();
   double synapse0xa00d118();
   double synapse0xa00d140();
   double synapse0xa00d168();
   double synapse0xa00d190();
   double synapse0xa00d2e0();
   double synapse0xa00d308();
   double synapse0xa00d330();
   double synapse0xa00d358();
   double synapse0xa00d380();
   double synapse0xa00d3a8();
   double synapse0xa00d3d0();
   double synapse0xa00d3f8();
   double synapse0xa00d420();
   double synapse0xa00d448();
   double synapse0xa00d470();
   double synapse0xa00d498();
   double synapse0xa00d4c0();
   double synapse0xa00d4e8();
   double synapse0xa00d510();
   double synapse0xa00d6a8();
   double synapse0xa00d6d0();
   double synapse0xa00d6f8();
   double synapse0xa00d720();
   double synapse0xa00d748();
   double synapse0xa00d770();
   double synapse0xa00d798();
   double synapse0xa00d7c0();
   double synapse0xa00d7e8();
   double synapse0xa00d810();
   double synapse0xa00d838();
   double synapse0xa00d860();
   double synapse0xa00d888();
   double synapse0xa00d8b0();
   double synapse0xa00d8d8();
   double synapse0xa00da70();
   double synapse0xa00da98();
   double synapse0xa00dac0();
   double synapse0xa00dae8();
   double synapse0xa00db10();
   double synapse0xa00db38();
   double synapse0xa00db60();
   double synapse0xa00db88();
   double synapse0xa00dbb0();
   double synapse0xa00dbd8();
   double synapse0xa00dc00();
   double synapse0xa00dc28();
   double synapse0xa00dc50();
   double synapse0xa00dc78();
   double synapse0xa00dca0();
   double synapse0xa00de38();
   double synapse0xa00de60();
   double synapse0xa00de88();
   double synapse0xa00deb0();
   double synapse0xa00ded8();
   double synapse0xa00df00();
   double synapse0xa00df28();
   double synapse0xa00df50();
   double synapse0xa00df78();
   double synapse0xa00dfa0();
   double synapse0xa00dfc8();
   double synapse0xa00dff0();
   double synapse0xa00e018();
   double synapse0xa00e040();
   double synapse0xa00e068();
   double synapse0xa00e200();
   double synapse0xa00e228();
   double synapse0xa00e250();
   double synapse0xa00a560();
   double synapse0xa00a588();
   double synapse0xa00a5b0();
   double synapse0xa00a5d8();
   double synapse0xa00a600();
   double synapse0xa00a628();
   double synapse0xa00a650();
   double synapse0xa00a678();
   double synapse0xa00a6a0();
   double synapse0xa00a6c8();
   double synapse0xa00a6f0();
   double synapse0xa00a718();
   double synapse0xa00a8b0();
   double synapse0xa00a8d8();
   double synapse0xa00a900();
   double synapse0xa00a928();
   double synapse0xa00a950();
   double synapse0xa00a978();
   double synapse0xa00a9a0();
   double synapse0xa00a9c8();
   double synapse0xa00a9f0();
   double synapse0xa00aa18();
   double synapse0xa00aa40();
   double synapse0xa00aa68();
   double synapse0xa00aa90();
   double synapse0xa00aab8();
   double synapse0xa00aae0();
   double synapse0xa00ac78();
   double synapse0xa00aca0();
   double synapse0xa00acc8();
   double synapse0xa00acf0();
   double synapse0xa00ad18();
   double synapse0xa00ad40();
   double synapse0xa00f280();
   double synapse0xa00f2a8();
   double synapse0xa00f2d0();
   double synapse0xa00f2f8();
   double synapse0xa00f320();
   double synapse0xa00f348();
   double synapse0xa00f370();
   double synapse0xa00f398();
   double synapse0xa00f3c0();
   double synapse0xa00f558();
   double synapse0xa00f580();
   double synapse0xa00f5a8();
   double synapse0xa00f5d0();
   double synapse0xa00f5f8();
   double synapse0xa00f620();
   double synapse0xa00f648();
   double synapse0xa00f670();
   double synapse0xa00f698();
   double synapse0xa00f6c0();
   double synapse0xa00f6e8();
   double synapse0xa00f710();
   double synapse0xa00f738();
   double synapse0xa00f760();
   double synapse0xa00f788();
   double synapse0xa00f920();
   double synapse0xa00f948();
   double synapse0xa00f970();
   double synapse0xa00f998();
   double synapse0xa00f9c0();
   double synapse0xa00f9e8();
   double synapse0xa00fa10();
   double synapse0xa00fa38();
   double synapse0xa00fa60();
   double synapse0xa00fa88();
   double synapse0xa00fab0();
   double synapse0xa00fad8();
   double synapse0xa00fb00();
   double synapse0xa00fb28();
   double synapse0xa00fb50();
   double synapse0xa00fce8();
   double synapse0xa00fd10();
   double synapse0xa00fd38();
   double synapse0xa00fd60();
   double synapse0xa00fd88();
   double synapse0xa00fdb0();
   double synapse0xa00fdd8();
   double synapse0xa00fe00();
   double synapse0xa00fe28();
   double synapse0xa00fe50();
   double synapse0xa00fe78();
   double synapse0xa00fea0();
   double synapse0xa00fec8();
   double synapse0xa00fef0();
   double synapse0xa00ff18();
   double synapse0xa0100b0();
   double synapse0xa0100d8();
   double synapse0xa010100();
   double synapse0xa010128();
   double synapse0xa010150();
   double synapse0xa010178();
   double synapse0xa0101a0();
   double synapse0xa0101c8();
   double synapse0xa0101f0();
   double synapse0xa010218();
   double synapse0xa010240();
   double synapse0xa010268();
   double synapse0xa010290();
   double synapse0xa0102b8();
   double synapse0xa0102e0();
   double synapse0xa010478();
   double synapse0xa0104a0();
   double synapse0xa0104c8();
   double synapse0xa0104f0();
   double synapse0xa010518();
   double synapse0xa010540();
   double synapse0xa010568();
   double synapse0xa010590();
   double synapse0xa0105b8();
   double synapse0xa0105e0();
   double synapse0xa010608();
   double synapse0xa010630();
   double synapse0xa010658();
   double synapse0xa010680();
   double synapse0xa0106a8();
   double synapse0xa010840();
   double synapse0xa010868();
   double synapse0xa010890();
   double synapse0xa0108b8();
   double synapse0xa0108e0();
   double synapse0xa010908();
   double synapse0xa010930();
   double synapse0xa010958();
   double synapse0xa010980();
   double synapse0xa0109a8();
   double synapse0xa0109d0();
   double synapse0xa0109f8();
   double synapse0xa010a20();
   double synapse0xa010a48();
   double synapse0xa010a70();
   double synapse0xa010c08();
   double synapse0xa010c30();
   double synapse0xa010c58();
   double synapse0xa010c80();
   double synapse0xa010ca8();
   double synapse0xa010cd0();
   double synapse0xa010cf8();
   double synapse0xa010d20();
   double synapse0xa010d48();
   double synapse0xa010d70();
   double synapse0xa010d98();
   double synapse0xa010dc0();
   double synapse0xa010de8();
   double synapse0xa010e10();
   double synapse0xa010e38();
   double synapse0xa010fd0();
   double synapse0xa010ff8();
   double synapse0xa011020();
   double synapse0xa011048();
   double synapse0xa011070();
   double synapse0xa011098();
   double synapse0xa0110c0();
   double synapse0xa0110e8();
   double synapse0xa011110();
   double synapse0xa011138();
   double synapse0xa011160();
   double synapse0xa011188();
   double synapse0xa0111b0();
   double synapse0xa0111d8();
   double synapse0xa011200();
   double synapse0xa0112b8();
   double synapse0xa0112e0();
   double synapse0xa011308();
   double synapse0xa011330();
   double synapse0xa011358();
   double synapse0xa011380();
   double synapse0xa0113a8();
   double synapse0xa0113d0();
   double synapse0xa0113f8();
   double synapse0xa011420();
   double synapse0xa011448();
   double synapse0xa011470();
   double synapse0xa011498();
   double synapse0xa0114c0();
   double synapse0xa0114e8();

   NueRecord &nuerec;
   
   static TMultiLayerPerceptron *fneuralNet_6inp;
   static TMultiLayerPerceptron *fneuralNet_30inp;
   static TMultiLayerPerceptron *fneuralNet_11inp;
   static TMultiLayerPerceptron *fneuralNet_11inp_daikon04;
   static TMultiLayerPerceptron *fneuralNet_14inp_daikon04;

   ClassDef(AnnAna,5)
   
};

#endif

