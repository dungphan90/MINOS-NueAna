#include <iostream>
#include "NueAna/MSTCalc.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"


MSTCalc::MSTCalc():
   enmsg(ANtpDefVal::kInt),
   ew1(ANtpDefVal::kFloat),
   enn1(ANtpDefVal::kInt),
   esm1(ANtpDefVal::kFloat),
   ewtot(ANtpDefVal::kFloat),
   enntot(ANtpDefVal::kInt),
   esmtot(ANtpDefVal::kFloat),
   e4w(ANtpDefVal::kFloat),
   e4sm(ANtpDefVal::kFloat),
   e4nn(ANtpDefVal::kInt),
   eb1(ANtpDefVal::kFloat),
   eb25(ANtpDefVal::kFloat),
   enbranch(ANtpDefVal::kInt),

   onmsg(ANtpDefVal::kInt),
   ow1(ANtpDefVal::kFloat),
   onn1(ANtpDefVal::kInt),
   osm1(ANtpDefVal::kFloat),
   owtot(ANtpDefVal::kFloat),
   onntot(ANtpDefVal::kInt),
   osmtot(ANtpDefVal::kFloat),
   o4w(ANtpDefVal::kFloat),
   o4sm(ANtpDefVal::kFloat),
   o4nn(ANtpDefVal::kInt),
   ob1(ANtpDefVal::kFloat),
   ob25(ANtpDefVal::kFloat),
   onbranch(ANtpDefVal::kInt),

   eeprob(ANtpDefVal::kFloat),
   oeprob(ANtpDefVal::kFloat),
   ealpha(ANtpDefVal::kFloat),
   oalpha(ANtpDefVal::kFloat),

   ebeta(ANtpDefVal::kFloat),
   obeta(ANtpDefVal::kFloat)
{
   eallw1=new float[2880];
   oallw1=new float[2880];
   eallm1=new float[2880];
   oallm1=new float[2880];

   for(int i=0;i<2880;i++){
      eallw1[i]=ANtpDefVal::kFloat;
      oallw1[i]=ANtpDefVal::kFloat;
      eallm1[i]=ANtpDefVal::kFloat;
      oallm1[i]=ANtpDefVal::kFloat;
   }
}

MSTCalc::MSTCalc(const MSTCalc &mst):
  TObject(),
  enmsg(mst.enmsg),
  ew1(mst.ew1),
  enn1(mst.enn1),
  esm1(mst.esm1),
  ewtot(mst.ewtot),
  enntot(mst.enntot),
  esmtot(mst.esmtot),
  e4w(mst.e4w),
  e4sm(mst.e4sm),
  e4nn(mst.e4nn),
  eb1(mst.eb1),
  eb25(mst.eb25),
  enbranch(mst.enbranch),
  
  onmsg(mst.onmsg),
  ow1(mst.ow1),
  onn1(mst.onn1),
  osm1(mst.osm1),
  owtot(mst.owtot),
  onntot(mst.onntot),
  osmtot(mst.osmtot),
  o4w(mst.o4w),
  o4sm(mst.o4sm),
  o4nn(mst.o4nn),
  ob1(mst.ob1),
  ob25(mst.ob25),
  onbranch(mst.onbranch),
  
  eeprob(mst.eeprob),
  oeprob(mst.oeprob),
  ealpha(mst.ealpha),
  oalpha(mst.oalpha),
  
  ebeta(mst.ebeta),
  obeta(mst.obeta)
{
   eallw1=new float[2880];
   oallw1=new float[2880];
   eallm1=new float[2880];
   oallm1=new float[2880];

   for(int i=0;i<2880;i++){
      if(i<mst.enn1){
	 eallw1[i]=mst.eallw1[i];
	 eallm1[i]=mst.eallm1[i];
      }
      else{
	 eallw1[i]=ANtpDefVal::kFloat;
	 eallm1[i]=ANtpDefVal::kFloat;
      }
      if(i<mst.onn1){
	 oallw1[i]=mst.oallw1[i];
	 oallm1[i]=mst.oallm1[i];
      }
      else{
	 oallw1[i]=ANtpDefVal::kFloat;
	 oallm1[i]=ANtpDefVal::kFloat;
      }
   }
}



MSTCalc::~MSTCalc()
{
   delete[] eallw1;
   delete[] oallw1;
   delete[] eallm1;
   delete[] oallm1;
}

void MSTCalc::Clear(Option_t* /* option */)
{
  Reset();
  delete[] eallw1;
  delete[] oallw1;
  delete[] eallm1;
  delete[] oallm1;
}

void MSTCalc::Reset()
{
   enmsg=ANtpDefVal::kInt;
   ew1=ANtpDefVal::kFloat;
   enn1=ANtpDefVal::kInt;
   esm1=ANtpDefVal::kFloat;
   ewtot=ANtpDefVal::kFloat;
   enntot=ANtpDefVal::kInt;
   esmtot=ANtpDefVal::kFloat;
   e4w=ANtpDefVal::kFloat;
   e4sm=ANtpDefVal::kFloat;
   e4nn=ANtpDefVal::kInt;
   eb1=ANtpDefVal::kFloat;
   eb25=ANtpDefVal::kFloat;
   enbranch=ANtpDefVal::kInt;
   
   onmsg=ANtpDefVal::kInt;
   ow1=ANtpDefVal::kFloat;
   onn1=ANtpDefVal::kInt;
   osm1=ANtpDefVal::kFloat;
   owtot=ANtpDefVal::kFloat;
   onntot=ANtpDefVal::kInt;
   osmtot=ANtpDefVal::kFloat;
   o4w=ANtpDefVal::kFloat;
   o4sm=ANtpDefVal::kFloat;
   o4nn=ANtpDefVal::kInt;
   ob1=ANtpDefVal::kFloat;
   ob25=ANtpDefVal::kFloat;
   onbranch=ANtpDefVal::kInt;

   eeprob=ANtpDefVal::kFloat;
   oeprob=ANtpDefVal::kFloat;
   ealpha=ANtpDefVal::kFloat;
   oalpha=ANtpDefVal::kFloat;

   ebeta=ANtpDefVal::kFloat;
   obeta=ANtpDefVal::kFloat;

   delete[] eallw1;
   delete[] oallw1;
   delete[] eallm1;
   delete[] oallm1;

   eallw1=new float[2880];
   oallw1=new float[2880];
   eallm1=new float[2880];
   oallm1=new float[2880];  
 

   for(int i=0;i<2880;i++){
     eallw1[i]=ANtpDefVal::kFloat;
     oallw1[i]=ANtpDefVal::kFloat;
     eallm1[i]=ANtpDefVal::kFloat;
     oallm1[i]=ANtpDefVal::kFloat;
   }
}

void MSTCalc::Print(Option_t* /*option*/) const
{
  std::cout<<"MSTCalcAna::Print******"<<std::endl
	   <<"NMSG even: "<<enmsg<<" odd: "<<onmsg<<std::endl
	   <<"W1  even: "<<ew1<<" odd: "<<ow1<<std::endl
	   <<"NN1 even: "<<enn1<<" odd: "<<onn1<<std::endl
	   <<"SM1 even: "<<esm1<<" odd: "<<osm1<<std::endl
	   <<"WTOT even: "<<ewtot<<" odd: "<<owtot<<std::endl
	   <<"NNTOT even: "<<enntot<<" odd: "<<onntot<<std::endl
	   <<"SMTOT even: "<<esmtot<<" odd: "<<osmtot<<std::endl
	   <<"4W even: "<<e4w<<" odd: "<<o4w<<std::endl
	   <<"4SM even: "<<e4sm<<" odd: "<<o4sm<<std::endl
	   <<"Ele Prob even: "<<eeprob<<" odd: "<<oeprob<<std::endl
	   <<"Alpha even: "<<ealpha<<" odd: "<<oalpha<<std::endl
	   <<"Ratio B1 even "<<eb1<<" odd: "<<ob1<<std::endl
	   <<"Ratio B25 even "<<eb25<<" odd: "<<ob25<<std::endl
	   <<"NBranches even "<<enbranch<<" odd: "<<onbranch<<std::endl;
}
