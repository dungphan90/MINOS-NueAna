#include "NueAna/ParticlePID/ParticleAna/ParticleCheck.h"
#include "TH3F.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"

#include "TFile.h"
#include "TMath.h"

#include "MCNtuple/NtpMCStdHep.h"

#include "NueAna/ParticlePID/ParticleAna/TruthCompareAna.h"
#include "MessageService/MsgService.h"

ParticleCheck::ParticleCheck(const char*fname)
{
        input=new TChain("PO");
        inputPA=new TChain("PA");
        input->Add(fname);
        inputPA->Add(fname);
        poh=new ParticleObjectHolder();
        prec=new PRecord();
        input->SetBranchAddress("ParticleObjectHolder",&poh);
        inputPA->SetBranchAddress("PRecord",&prec);
}

ParticleCheck::~ParticleCheck()
{}


void ParticleCheck::Process(int cut, int type, double weight)
{

        //do we have a single electron match?
        
        for(int i=0;i<(int)prec->truthcompare.reco_to_true_match.size();i++)
        {
                std::vector<int> truematch=prec->truthcompare.reco_to_true_match[i];
                Particle3D * p = &poh->particles3d[prec->truthcompare.reco_matched[i]];

                if(truematch.size()==1)
                {
                        //we have a single match...
                        //are both electons?
                        if(p->particletype!=Particle3D::electron)continue;
                        NtpMCStdHep *mc = &(prec->mctrue.stdhep[truematch[0]]); 
                        double trueE=sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass);       
                
                        if(TMath::Abs(mc->IdHEP)==11)
                        {
                                elec_single_match_eres[cut][type]->Fill((p->sum_e/25.-trueE)/trueE,trueE,weight);
                                elec_single_match_angle[cut][type]->Fill(prec->truthcompare.matchangle,trueE,weight);
                        }

                        if(TMath::Abs(mc->IdHEP)==111)
                        {
                                elec_single_pi0_eres[cut][type]->Fill((p->sum_e/25.-trueE)/trueE,trueE,weight);
                                elec_single_pi0_angle[cut][type]->Fill(prec->truthcompare.matchangle,trueE,weight);
                        }
                }
                        
                //here we can have multiple true matched particles
                if(p->particletype==Particle3D::electron)
                {
                        int le=-1;
                        double leE=0;
                        
                        int lp0=-1;
                        double lp0E=0;
                        
                        double trueEtotal=0;
                        for(int j=0;j<(int)truematch.size();j++)
                        {
                                NtpMCStdHep *mc = &(prec->mctrue.stdhep[truematch[j]]);
                                double trueE=sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass);
                                trueEtotal+=trueE;
                        }
                        
                        for(int j=0;j<(int)truematch.size();j++)
                        {
                                NtpMCStdHep *mc = &(prec->mctrue.stdhep[truematch[j]]);
                                double trueE=sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass);

                                elec_matched_type_vs_recoE[cut][type]->Fill(TMath::Abs(mc->IdHEP),p->sum_e/25.,weight);
                                elec_matched_typefrac_vs_recoE[cut][type]->Fill(TMath::Abs(mc->IdHEP),trueE/trueEtotal,weight);
                                
                                if(TMath::Abs(mc->IdHEP)==11)
                                {
                                        if(trueE>leE)
                                        {
                                                leE=trueE;
                                                le=j;
                                        }
                                }

                                if(TMath::Abs(mc->IdHEP)==111)
                                {
                                        if(trueE>lp0E)
                                        {
                                                lp0E=trueE;
                                                lp0=j;
                                        }
                                }

                        }
                        if(le>-1)
                        {
                                NtpMCStdHep *mc = &(prec->mctrue.stdhep[truematch[le]]);
                                double trueE=sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass);
                                elec_matched_elec_eres[cut][type]->Fill((p->sum_e/25.-trueE)/trueE,trueE,weight);
                        }else{
                                elec_matched_elec_eres[cut][type]->Fill(0.,0.,weight);
                        }
                        
                        if(lp0>-1)
                        {
                                NtpMCStdHep *mc = &(prec->mctrue.stdhep[truematch[lp0]]);
                                double trueE=sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass);
                                elec_matched_pi0_eres[cut][type]->Fill((p->sum_e/25.-trueE)/trueE,trueE,weight);
                        }else{
                                elec_matched_pi0_eres[cut][type]->Fill(0.,0.,weight);
                        }
                
                
                
                        
                }
                
        }



}

void ParticleCheck::MakeHistos()
{
        for(int i=0;i<nCut;i++)
        for(int j=0;j<nType;j++)
        {
                char tmp[200];
                
                sprintf(tmp,"elec_single_match_eres_%d_%d",i,j);
                elec_single_match_eres[i][j]=new TH2F(tmp,tmp,100,-2,2,100,0,10);
                
                sprintf(tmp,"elec_single_match_angle_%d_%d",i,j);
                elec_single_match_angle[i][j]=new TH2F(tmp,tmp,100,-2,2,100,0,10);

                sprintf(tmp,"elec_single_pi0_eres_%d_%d",i,j);
                elec_single_pi0_eres[i][j]=new TH2F(tmp,tmp,100,-2,2,100,0,10);
                
                sprintf(tmp,"elec_single_pi0_angle_%d_%d",i,j);
                elec_single_pi0_angle[i][j]=new TH2F(tmp,tmp,100,-2,2,100,0,10);
                        
                
                sprintf(tmp,"elec_matched_elec_eres_%d_%d",i,j);
                elec_matched_elec_eres[i][j]=new TH2F(tmp,tmp,100,-2,2,100,0,10);       
                
                sprintf(tmp,"elec_matched_pi0_eres_%d_%d",i,j);
                elec_matched_pi0_eres[i][j]=new TH2F(tmp,tmp,100,-2,2,100,0,10);        

                sprintf(tmp,"elec_matched_type_vs_recoE_%d_%d",i,j);
                elec_matched_type_vs_recoE[i][j]=new TH2F(tmp,tmp,2300,0,2300,100,0,10);                                                

                sprintf(tmp,"elec_matched_typefrac_vs_recoE_%d_%d",i,j);
                elec_matched_typefrac_vs_recoE[i][j]=new TH2F(tmp,tmp,2300,0,2300,100,0,10);            
                        
        }
}



void ParticleCheck::SaveHistos()
{
        TFile *f = TFile::Open("pcheck.root","RECREATE");
        f->cd();
        for(int i=0;i<nCut;i++)
        for(int j=0;j<nType;j++)
        {
                elec_single_match_eres[i][j]->Write();
                elec_single_match_angle[i][j]->Write();
                elec_single_pi0_eres[i][j]->Write();
                elec_single_pi0_angle[i][j]->Write();
                elec_matched_elec_eres[i][j]->Write();
                elec_matched_pi0_eres[i][j]->Write();
                elec_matched_type_vs_recoE[i][j]->Write();
                elec_matched_typefrac_vs_recoE[i][j]->Write();
        }
        f->Close();
}


int ParticleCheck::PassesPreselec()
{
        int pass=1;
        
        pass = pass && prec->event.inFiducial==1;
        
        pass = pass && prec->event.nstrips>10 && prec->event.nstrips<50;
        
        pass = pass && prec->particles.ntot>0;



        pass = pass &&  prec->event.vtx_z-prec->event.min_z<0.35;
        pass = pass &&  prec->event.nstrips>15 && prec->event.nstrips<45;
        pass = pass &&  prec->particles.rms_r<1.1;
        pass = pass &&  prec->particles.prim_par_b>0.1&&prec->particles.prim_par_b<1.2;
        pass = pass &&  prec->particles.ntot/prec->event.visenergy<0.14;
        pass = pass &&  prec->particles.mol_rad_r<0.0425;
        pass = pass &&  prec->particles.largest_particle_avg_rms<0.017;
        pass = pass &&  prec->particles.largest_particle_avg_rms>0.005;
        pass = pass &&  prec->particles.ntot<8 ;
        pass = pass &&  prec->particles.longest_z<0.8;
        pass = pass &&  prec->particles.primary_long_e<60;
        
        

        
/*      
        pass = pass && prec->particles.rough_primary_theta_z<0.7;
        pass = pass && prec->particles.primary_long_e<15;
        
        pass = pass && sqrt((prec->event.max_u-prec->event.min_u)*(prec->event.max_u-prec->event.min_u)
                        +(prec->event.max_v-prec->event.min_v)*(prec->event.max_v-prec->event.min_v))>0.3;
                        
        pass = pass && prec->particles.largest_particle_avg_rms>0.0045;
        pass = pass && prec->particles.largest_particle_e/prec->event.visenergy>0.85;

        pass = pass && prec->particles.total_long_e<25;
        pass = pass && prec->particles.mol_rad_r>0.001;
        pass = pass && prec->particles.mol_rad_r<0.035;

        pass = pass && prec->particles.nelec>0;
        pass = pass && prec->particles.ntot<4;
        pass = pass && prec->particles.prim_par_b>0.1;
        pass = pass && prec->particles.prim_par_b<1;

        pass = pass && prec->particles.longest_z<1;
        pass = pass && prec->particles.longest_z>0.1;
        pass = pass && prec->particles.longest_particle_type<14;
*/
        return pass;
}

void ParticleCheck::Run()
{

        MakeHistos();



        int ent =input->GetEntries();
        int entpa =inputPA->GetEntries();
        printf("%d %d total entries\n",ent,entpa);
        
        
        double typecnt[5];
        for(int i=0;i<5;i++)typecnt[i]=0;
        
        
        int saved=0;
        int ji=0;
        input->GetEntry(0);
        for(int i=0;i<entpa;i++)
        {
                inputPA->GetEntry(i);
                
                while(poh->GetHeader().GetRun() != prec->GetHeader().GetRun() || poh->GetHeader().GetSnarl()!=prec->GetHeader().GetSnarl()  || poh->GetHeader().GetEvent()!=prec->GetHeader().GetEvent())
                        input->GetEntry(ji++);
                        
        //      MsgService * ms = MsgService::Instance();
        //      ms->GetStream("TruthCompareAna")->SetLogLevel(Msg::kDebug);
                
                if(prec->truthcompare.stage<2)
                {
                        TruthCompareAna a;
                        a.ana(poh,&prec->truthcompare);
                }
                
                
        
        
        
                double weight = poh->mctrue.totbeamweight*poh->mctrue.oscprob*3.25/6.5/5;
                double paweight = prec->mctrue.totbeamweight*prec->mctrue.oscprob*3.25/6.5/5;
                
                
                
                if(TMath::Abs(weight-paweight)>0.001)printf("weight mismatch!\n");
                
                int type=poh->mctrue.type;
        
                if(type<0 || type>5)continue;
        
                if(type==0)weight/=3.;
                if(type==0)paweight/=3.;
        

                for(int nn=0;nn<nCut;nn++)
                {
                
                        if(nn==1 && !PassesPreselec())continue;
                
                        Process(nn,type,weight);        
                }
                
                typecnt[type]+=weight;



                if(i%10000==0)printf("Entry %d\n",i);
        }

        printf("used %d entries\n",saved);

        printf("values\n");
        for(int i=0;i<5;i++)
                printf("%d   %2.2f\n",i,typecnt[i]);
        
        printf("FOM %2.2f\n",typecnt[1]/sqrt(typecnt[0]+typecnt[2]+typecnt[3]+typecnt[4]));
        
        SaveHistos();

}
