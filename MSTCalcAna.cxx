#include "TH1F.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TClonesArray.h"
#include "StandardNtuple/NtpStRecord.h"
#include "Conventions/PlaneView.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "MessageService/MsgService.h"
#include "CalDetDST/UberRecord.h"
#include "CalDetDST/UberHit.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/NueAnaTools/DCGraph.h"
#include "NueAna/NueAnaTools/DCVertex.h"
#include "NueAna/NueAnaTools/DCHit.h"
#include "NueAna/MSTCalcAna.h"
#include "NueAna/MSTCalc.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

TH1* eth=0;
TH1* bth=0;

CVSID("$Id: MSTCalcAna.cxx,v 1.13 2009/06/30 20:52:49 vahle Exp $");

Double_t AFit(Double_t *x, Double_t *par){ 
   if(eth==0){
      return 0.;
   }
   if(bth==0){
      return 0.;
   }
   double xo = (*x);
   int sbin = eth->FindBin(xo,0.,0.);
   int bbin = bth->FindBin(xo,0.,0.);
   double sig = 1.*eth->GetBinContent(sbin);
   double back = 1.*bth->GetBinContent(bbin);
   return par[0]*sig+(1-par[0])*back;
//   return par[0]*sig+par[1]*back;
}

MSTCalcAna::MSTCalcAna(MSTCalc &mst):
  fMSTCalc(mst),
  minsig(0.),
  maxmetric(0.),
  minfarsigcor(0.),
  maxmetriclowz(0.),
  even(),
  odd(),
  etemplate(0),
  btemplate(0),
  emtemplate(0),
  bmtemplate(0)
{}

MSTCalcAna::~MSTCalcAna()
{
   for(unsigned int i=0;i<even.size();i++){
      if(even[i]!=0){
	 delete even[i];
	 even[i]=0;
      }
   }
   for(unsigned int i=0;i<odd.size();i++){
      if(odd[i]!=0){
	 delete odd[i];
	 odd[i]=0;
      }
   }
}

//void MSTCalcAna::Analyze(int evtn, NtpSRRecord *srobj, 
//			 NtpMCRecord */*mc*/, NtpTHRecord */*th*/)
void MSTCalcAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  set<DCHit> evenhitset;
  set<DCHit> oddhitset;
  Reset();
  FillHitSets(evtn,srobj,evenhitset,oddhitset);
  FindAllGraphs(evenhitset,oddhitset);
  ComputeGraphQuantities();
  //  cout<<"done with ana"<<endl;
}

void MSTCalcAna::Analyze(UberRecord *ur)
{
  set<DCHit> evenhitset;
  set<DCHit> oddhitset;
  Reset();
  FillHitSets(ur,evenhitset,oddhitset);
  FindAllGraphs(evenhitset,oddhitset);
  ComputeGraphQuantities();
}


void MSTCalcAna::ComputeGraphQuantities()
{
  //  cout<<"In computegraph"<<endl;
  //compute quantities
  fMSTCalc.enmsg=even.size();
  fMSTCalc.onmsg=odd.size();

  fMSTCalc.esmtot = 0.0;
  fMSTCalc.ewtot  = 0.0;
  fMSTCalc.enntot = 0;
  fMSTCalc.osmtot = 0.0;
  fMSTCalc.owtot  = 0.0;
  fMSTCalc.onntot = 0;
    
  const int NYBINS = etemplate->GetNbinsY()+1;

  TH1F ew("ew","ew",NYBINS-1,etemplate->GetYaxis()->GetXbins()->GetArray());
  TH1F ow("ow","ow",NYBINS-1,etemplate->GetYaxis()->GetXbins()->GetArray());
  
  const int NYMBINS = emtemplate->GetNbinsY()+1;
  TH1F em("em","em",NYMBINS-1,emtemplate->GetYaxis()->GetXbins()->GetArray());
  TH1F om("om","om",NYMBINS-1,emtemplate->GetYaxis()->GetXbins()->GetArray());
  
  for(unsigned int i=0;i<even.size();i++){
    if(even[i]==0){
      continue;
    }
    fMSTCalc.esmtot+=even[i]->GetSummedZ();
    fMSTCalc.ewtot+=even[i]->GetSummedMetric();
    fMSTCalc.enntot+=even[i]->GetNVerts();
  }
  for(unsigned int i=0;i<odd.size();i++){
    if(odd[i]==0){
      continue;
    }
    fMSTCalc.osmtot+=odd[i]->GetSummedZ();
    fMSTCalc.owtot+=odd[i]->GetSummedMetric();
    fMSTCalc.onntot+=odd[i]->GetNVerts();
  }
  
  DCGraph<DCHit> *me = GetMaxGraph(0);
  if(me!=0){
    fMSTCalc.esm1=me->GetSummedZ();
    fMSTCalc.ew1=me->GetSummedMetric();
    fMSTCalc.enn1=me->GetNVerts();	    	    
      vector<float> eweights = me->GetAllWeights();
      int NW=eweights.size();
      for(int i=0;i<NW;i++){
	 if(NW<2880){
	   fMSTCalc.eallw1[i] = eweights[i];
	 }
	 if(eweights[i]!=0){
	    ew.Fill(eweights[i]);
	 }
      }

      int nabove=0;
      float aveph=fMSTCalc.esm1/fMSTCalc.enn1;
      vector<float> emips = me->GetAllZ();
      int NM=emips.size();
      for(int i=0;i<NM;i++){
	 if(NM<2880){
	   fMSTCalc.eallm1[i] = emips[i];
	 }
	 em.Fill(emips[i]);
	 if(emips[i]>aveph){
	    nabove++;
	 }
      }
      
      if(nabove>0){
	 DCGraph<DCHit> *emax4 = me->GraphMaxHits(nabove);
	 fMSTCalc.e4w = emax4->GetSummedMetric();
	 fMSTCalc.e4sm = emax4->GetSummedZ();
	 fMSTCalc.e4nn = emax4->GetNVerts();
	 delete emax4;
      }
      else{
	fMSTCalc.e4w = ANtpDefVal::kFloat;
	fMSTCalc.e4sm = ANtpDefVal::kFloat;
	fMSTCalc.e4nn = 0;
      }	 
  }
  GraphLoop(me,0);
  

  DCGraph<DCHit> *mo = GetMaxGraph(1);
  if(mo!=0){
    fMSTCalc.osm1=mo->GetSummedZ();
    fMSTCalc.ow1=mo->GetSummedMetric();
    fMSTCalc.onn1=mo->GetNVerts();	    	    
    vector<float> oweights = mo->GetAllWeights();
     int NW=oweights.size();
     for(int i=0;i<NW;i++){
       if(NW<2880){
	 fMSTCalc.oallw1[i] = oweights[i];
       }
       if(oweights[i]!=0){
	 ow.Fill(oweights[i]);
       }
     }

     int nabove=0;
     float aveph=fMSTCalc.osm1/fMSTCalc.onn1;
     vector<float> omips = mo->GetAllZ();
     int NM=omips.size();
     for(int i=0;i<NM;i++){
       if(NM<2880){
	 fMSTCalc.oallm1[i] = omips[i];
       }
       om.Fill(omips[i]);
       if(omips[i]>aveph){
	 nabove++;
       }
     }
     
     if(nabove>0){
       DCGraph<DCHit> *omax4 = mo->GraphMaxHits(nabove);
       fMSTCalc.o4w = omax4->GetSummedMetric();
       fMSTCalc.o4sm = omax4->GetSummedZ();
       fMSTCalc.o4nn = omax4->GetNVerts();
       delete omax4;
     }
     else{
       fMSTCalc.o4w=ANtpDefVal::kFloat;
       fMSTCalc.o4sm=ANtpDefVal::kFloat;
       fMSTCalc.o4nn=ANtpDefVal::kInt;
     }
  }
  GraphLoop(mo,1);
  
  if(etemplate==0||btemplate==0){
    MSG("MSTCalcAna",Msg::kError)<<"Template Histograms are not available"
				 <<" can not calculate likelihoods"
				 <<" or alphas"<<endl;
    fMSTCalc.eeprob=ANtpDefVal::kFloat;
     fMSTCalc.oeprob=ANtpDefVal::kFloat;
     fMSTCalc.ealpha=ANtpDefVal::kFloat;
     fMSTCalc.oalpha=ANtpDefVal::kFloat;
     fMSTCalc.ebeta=ANtpDefVal::kFloat;
     fMSTCalc.obeta=ANtpDefVal::kFloat;
     ew.Reset();
     ow.Reset();
     em.Reset();
     om.Reset();

     return;
   }

   TH1F *esigt = (TH1F *)FindProjection(etemplate,
					fMSTCalc.enn1,"esigt");
   TH1F *osigt = (TH1F *)FindProjection(etemplate,
					fMSTCalc.onn1,"osigt");
   TH1F *ebgt = (TH1F *)FindProjection(btemplate,
				       fMSTCalc.enn1,"ebgt");
   TH1F *obgt = (TH1F *)FindProjection(btemplate,
				       fMSTCalc.onn1,"obgt");

   int endof, ondof;
   fMSTCalc.eeprob = ComputeLLike(&ew,esigt,endof);
   fMSTCalc.oeprob = ComputeLLike(&ow,osigt,ondof);


   FindLikeAlpha(&ew,esigt,ebgt,fMSTCalc.ealpha,fMSTCalc.ebeta);
   FindLikeAlpha(&ow,osigt,obgt,fMSTCalc.oalpha,fMSTCalc.obeta);

   esigt->Reset();
   osigt->Reset();
   ebgt->Reset();
   obgt->Reset();
   
   ew.Reset();
   ow.Reset();
   em.Reset();
   om.Reset();
   //   cout<<"done with compute graph"<<endl;
}

void MSTCalcAna::FindAllGraphs(set<DCHit> &evenhitset, 
			       set<DCHit> &oddhitset)
{
  //  cout<<"In findalagraphs"<<endl;
  DCGraph<DCHit> geven;
  DCGraph<DCHit> godd;
   set<DCHit>::reverse_iterator eit=evenhitset.rbegin();
   while(eit!=evenhitset.rend()){
     DCHit th=*eit;
     geven.AddVertex(&th);
     eit++;
   }
   
   set<DCHit>::reverse_iterator oit=oddhitset.rbegin();
   while(oit!=oddhitset.rend()){
      DCHit th=*oit;
      godd.AddVertex(&th);
      oit++;
   }

  if(geven.GetNVerts()!=0){
    even=geven.FindAllMinSpan(maxmetric,minfarsigcor,maxmetriclowz);
  }
  if(godd.GetNVerts()!=0){
    odd=godd.FindAllMinSpan(maxmetric,minfarsigcor,maxmetriclowz);
  }
  //  cout<<"done with findallgraphs"<<endl;
}


//void MSTCalcAna::FillHitSets(int evtn, NtpSRRecord *srobj, 
//     set<DCHit> &evenhitset, set<DCHit> &oddhitset)
void MSTCalcAna::FillHitSets(int evtn, RecRecordImp<RecCandHeader> *srobj, 
			     set<DCHit> &evenhitset, set<DCHit> &oddhitset)
{

  //  cout<<"In fill hitsets"<<endl;
  if(srobj==0){
    return;
  }
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }


  NtpSREvent *event = SntpHelpers::GetEvent(evtn,srobj);
  if(!event){
    MSG("MSTCalcAna",Msg::kError)<<"Couldn't get event "<<evtn
				 <<" from Snarl "
				 <<srobj->GetHeader().GetSnarl()
				 <<endl;
    return;
  }
  
  //  cout<<"gotevent"<<endl;
  
  int nstrips=event->nstrip;
  //  cout<<"nstrips "<<nstrips<<endl;
  for(int i=0;i<nstrips;i++){
    Int_t index = SntpHelpers::GetStripIndex(i,event);
    NtpSRStrip *s = SntpHelpers::GetStrip(index,srobj);
    if(!s){
      continue;
    }
    if(!evtstp0mip){
       MSG("MSTCalcAna",Msg::kError)<<"No mip strip information"<<endl;
       continue;
    }


    float energy = evtstp0mip[index] + evtstp1mip[index];
    
    if(energy<minsig){
      continue;
    }
    
    double time = (evtstp0mip[index]*s->time0+
		   evtstp1mip[index]*s->time1);
    if(energy!=0){
      time/=(energy);
    }
    else{
      time=0.;
    }
    DCHit h(s->z*100.,s->tpos*100., energy,time);
    h.SetData("zpos","tpos","pulseheight");
    
    if(s->planeview==PlaneView::kU){
      evenhitset.insert(h);
    }
    else if(s->planeview==PlaneView::kV){
      oddhitset.insert(h);
    }
    else{
    }
  }
  //  cout<<"done fil hitsets"<<endl;
}


void MSTCalcAna::FillHitSets(UberRecord *ur,
			     set<DCHit> &evenhitset, 
			     set<DCHit> &oddhitset)
{

  const TClonesArray *hits = ur->GetHitList();
  for(int i=0;i<ur->nhitstrips;i++){
    UberHit *uh=static_cast<UberHit *>((*hits)[i]);
    if(uh->GetPosMIP()+uh->GetNegMIP()<minsig){
      continue;
    }
    DCHit h(uh);
    h.SetData("zpos","tpos","pulseheight");
    
    if(uh->GetPlane()%2==0){
      evenhitset.insert(h);
    }
    else{
      oddhitset.insert(h);
    }
  }
}

void MSTCalcAna::Reset()
{
   for(unsigned int i=0;i<even.size();i++){
      if(even[i]!=0){
	 delete even[i];
	 even[i]=0;
      }
   }
   for(unsigned int i=0;i<odd.size();i++){
      if(odd[i]!=0){
	 delete odd[i];
	 odd[i]=0;
      }
   }
   even.clear();
   odd.clear();

   fMSTCalc.Reset();
}


void MSTCalcAna::DrawMaxGraph(int view)
{
   float maxmip=0.;
   int index=0;
   vector<DCGraph<DCHit> *> *g;
   if(view==0){
      g=&even;
   }
   else{
      g=&odd;
   }
   for(unsigned int i=0;i<g->size();i++){
      if((*g)[i]->GetSummedZ()>maxmip){
	 maxmip=(*g)[i]->GetSummedZ();
	 index=i;
      }
   }
   if(g->size()>0){
      (*g)[index]->Draw();
   }
}

DCGraph<DCHit> *MSTCalcAna::GetMaxGraph(int view)
{
  //  cout<<"in getmaxgraph "<<view<<endl;
   float maxmip=0.;
   int index=0;
   vector<DCGraph<DCHit> *> *g;
   if(view==0){
      g=&even;
   }
   else{
      g=&odd;
   }
   if(g->size()==0){
      return 0;
   }
   for(unsigned int i=0;i<g->size();i++){
      if((*g)[i]->GetSummedZ()>maxmip){
	 maxmip=(*g)[i]->GetSummedZ();
	 index=i;
      }
   }
   return (*g)[index];
}

void MSTCalcAna::Print() const
{
  fMSTCalc.Print();
}

TH1D *MSTCalcAna::FindProjection(TProfile2D *tmplt, int n, const char* name)
{
//   string n2d=(string)(name)+"_2d";
//   TH2D *proj2d = tmplt->ProjectionXY(n2d.c_str(),"E");
   int bin = tmplt->GetXaxis()->FindFixBin(n);
   TH1D *proj = tmplt->ProjectionY(name,bin,bin,"E,goff");
//   delete proj2d;
//   proj->Sumw2();
   return proj;
}

double MSTCalcAna::ComputeChi2(TH1 *h1, TH1 *h2, 
			      double &chi2, int &ndof)
{
   chi2=0.;
   double NENTH1=h1->Integral();
   double NENTH2=h2->Integral();

   ndof = 0;
   if(NENTH2==0||NENTH1==0){
      chi2=-1.;
      ndof=0;
      return -1.;
   }
   
//   double weight1 = sqrt(NENTH2/NENTH1);
//   double weight2 = sqrt(NENTH1/NENTH2);


   double weight1 = 1.;
   double weight2 = 1.;
   for(int i=1;i<=h1->GetNbinsX();i++){
      double h1bc=h1->GetBinContent(i);
      double h2bc=h2->GetBinContent(i);
      double eh1bc=h1->GetBinError(i);
      double eh2bc=h2->GetBinError(i);
      
      double num = pow((weight1*h1bc-weight2*h2bc),2);
      double denom = eh1bc*eh1bc+eh2bc*eh2bc;

      if(denom!=0){
	 chi2+=num/denom;
	 ndof++;
      }
   }
   ndof--; //one constraint, number of entries must be equal to n
   Double_t prob = TMath::Prob(chi2/2.,Int_t(ndof/2.));

   return prob;
}

double MSTCalcAna::ComputeLLike(TH1 *data, TH1 *tmpl, int &ndof)
{

   if(data==0){
      return ANtpDefVal::kFloat;
   }
   if(tmpl==0){
      return ANtpDefVal::kFloat;
   }
   if(data->Integral()<2){
      return ANtpDefVal::kFloat;
   }
   if(data->GetNbinsX()!=tmpl->GetNbinsX()){
     MSG("MSTCalcAna",Msg::kError)<<"Data bins "
				  <<data->GetNbinsX()<<endl
				  <<"Template bins "
				  <<tmpl->GetNbinsX()<<endl
				  <<"Can not do likelihood"<<endl;
     return ANtpDefVal::kFloat;
   }
   
   ndof=0;
   double llike=0.;

   for(int i=1;i<=tmpl->GetNbinsX();i++){
      float x = data->GetBinContent(i);
      float mu = tmpl->GetBinContent(i);
      double l=0.;
      if(x>0&&mu>0){
	 l=x*(TMath::Log(mu)-TMath::Log(x)+1)-mu;
      }
      else if(x==0&&mu!=0){
	 l=-mu;
      }
      else if(x!=0&&mu==0){
	 l=0.;
      }
      else{
	 l=0.;
      }
	 
      if(x!=0){
	 ndof+=(int)(x);
      }
      llike+=l;	   

   }

   double chi2 = -2*llike;
   return chi2;
}

float MSTCalcAna::FindAlpha(TH1 *data, TH1 *ehist,TH1 *bhist)
{

   if(data==0){      
      MSG("MSTCalcAna",Msg::kError)<<"data==0"<<endl;
      return ANtpDefVal::kFloat;
   }
   if(ehist==0){
      MSG("MSTCalcAna",Msg::kError)<<"ehist==0"<<endl;
      return ANtpDefVal::kFloat;
   }
   if(bhist==0){
      MSG("MSTCalcAna",Msg::kError)<<"bhist==0"<<endl;
      return ANtpDefVal::kFloat;
   }

   double num=0.;
   double denom=0.;
   double error=0.;
   float ea=0.;
   if(data->GetNbinsX()==ehist->GetNbinsX()&&
      data->GetNbinsX()==bhist->GetNbinsX()){
      for(int i=1;i<=data->GetNbinsX();i++){
	 double n=1.*data->GetBinContent(i);
	 double c=1.*ehist->GetBinContent(i);
	 double d=1.*bhist->GetBinContent(i);
	 double en=1.*data->GetBinError(i);
	 double ec=1.*ehist->GetBinError(i);
	 double ed=1.*bhist->GetBinError(i);

	 error = en*en+ec*ec+ed*ed;
	 if(error!=0){
	    num+=(d-n)*(d-c)/error;
	    denom+=(d-c)*(d-c)/error;
	 }
      }
      if(denom>1.e-10){
	 ea=num/denom;
      }
      else{
	 ea=ANtpDefVal::kFloat;
      }
   }
   return ea;

}


void MSTCalcAna::FindLikeAlpha(TH1 *data, TH1 *ehist,TH1 *bhist, 
			      double &alpha, double &beta)
{

   if(data==0){      
      MSG("MSTCalcAna",Msg::kError)<<"data==0"<<endl;
      alpha=ANtpDefVal::kFloat;
      beta=ANtpDefVal::kFloat;
      return;
   }
   if(ehist==0){
      MSG("MSTCalcAna",Msg::kError)<<"ehist==0"<<endl;
      alpha=ANtpDefVal::kFloat;
      beta=ANtpDefVal::kFloat;
      return;
   }
   if(bhist==0){
      MSG("MSTCalcAna",Msg::kError)<<"bhist==0"<<endl;
      alpha=ANtpDefVal::kFloat;
      beta=ANtpDefVal::kFloat;
      return;
   }

   if(data->GetNbinsX()!=ehist->GetNbinsX()){
      MSG("MSTCalcAna",Msg::kError)<<"Data bins: "
				   <<data->GetNbinsX()<<endl
				   <<"etemplate bins: "
				   <<ehist->GetNbinsX()<<endl
				   <<"Cant do like alpha"<<endl;
      alpha=ANtpDefVal::kFloat;
      beta=ANtpDefVal::kFloat;
      return;
   }
   if(data->GetNbinsX()!=bhist->GetNbinsX()){
      MSG("MSTCalcAna",Msg::kError)<<"Data bins: "
				   <<data->GetNbinsX()<<endl
				   <<"btemplate bins: "
				   <<bhist->GetNbinsX()<<endl
				   <<"Cant do like alpha"<<endl;
      alpha=ANtpDefVal::kFloat;
      beta=ANtpDefVal::kFloat;
      return;
   }
   eth=ehist;
   bth=bhist;

   TF1 *histfit = new TF1("histfit",AFit,0,100,1);
   histfit->SetParName(0,"alpha");
   histfit->SetParameter(0,0.9);
   histfit->SetParLimits(0,0.,1.);

   data->Fit(histfit,"RLQNO");
//   data->Fit(histfit,"RIL");
   
   alpha = histfit->GetParameter(0);
   beta = histfit->GetParError(0);
      
   if(histfit!=0){
     delete histfit;
     histfit=0;
   }
   return;
}

void MSTCalcAna::GraphLoop(DCGraph<DCHit> *g, int view)
{
  //  cout<<"in graph loop "<<view<<endl;
  if(g==0){
    return;
  }

  float mipsum=0;
  float b1mipsum=0;
  float b25mipsum=0;
  int nbranch=0;
  DCVertex<DCHit> *v = g->GetFirst();
  while(v!=0){
    DCEdge<DCHit> *e = v->GetFirstEdge();
    float vx1,vy1,vz1;
    v->GetData()->GetData(vx1,vy1,vz1);
    int nedge=0;
    while(e!=0){
      float vx2=0,vy2=0,vz2=0;
      e->GetEnd()->GetData()->GetData(vx2,vy2,vz2);
      if(e->GetMetric()<5){
	b1mipsum+=vz1+vz2;
      }
      if(e->GetMetric()>25){
	b25mipsum+=vz1+vz2;
      }
      mipsum+=vz1+vz2;
      nedge++;
      e=e->GetNext();
    }
    if(nedge>1){ nbranch++; }
    v=v->GetNextVertex();
  }


  if(view==0){
    if(mipsum>0){
      fMSTCalc.eb1=b1mipsum/mipsum;
      fMSTCalc.eb25=b25mipsum/mipsum;
    }
    fMSTCalc.enbranch=nbranch;
  }
  else if(view==1){
    if(mipsum>0){
      fMSTCalc.ob1=b1mipsum/mipsum;
      fMSTCalc.ob25=b25mipsum/mipsum;
    }
    fMSTCalc.onbranch=nbranch;
  }
  //  cout<<"done with graph loop"<<endl;
}
