#include "NueAna/ParticlePID/ParticleAna/TruthCompareAna.h"
#include "TMath.h"

#include "MCNtuple/NtpMCStdHep.h"
#include <map>
#include "MessageService/MsgService.h"

CVSID("$Id: TruthCompareAna.cxx,v 1.3 2009/06/24 23:06:46 gmieg Exp $");

TruthCompareAna::TruthCompareAna()
{
        Reset();
        
        meupergev = 25.;//should use something more precise here!
}

TruthCompareAna::~TruthCompareAna()
{
        Reset();
}

void TruthCompareAna::Reset()
{
        
        
        true_idx.clear();
        matched.clear();
        reco_matched_visible_e=0;
        reco_visible_e=0;

}

void TruthCompareAna::ana(ParticleObjectHolder * poh, TruthCompare * t)
{
        
        t->Reset();
        
        MSG("TruthCompareAna",Msg::kDebug)<<"----------------------------\n";
        MSG("TruthCompareAna",Msg::kDebug)<<"Run: "<<poh->GetHeader().GetRun()<<" Snarl: "<<poh->GetHeader().GetSnarl()<<" Event: "<<poh->GetHeader().GetEvent()<<"\n"; 
        ProcessTrueParticles(poh);
         MSG("TruthCompareAna",Msg::kDebug)<<matched.size()<<" true particles to be considered\n";


        //find lep energy from CC
        t->truelepE=0;
        t->truelepType=0;
        for(unsigned int i=0;i<poh->mctrue.stdhep.size();i++)
        {
                NtpMCStdHep *mc = &(poh->mctrue.stdhep[i]);
                
                
                if(mc->parent[0]!=0)continue;
                                
                if(TMath::Abs(mc->IdHEP)==11 || TMath::Abs(mc->IdHEP)==13 || TMath::Abs(mc->IdHEP)==15)
                {
                        
                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);

                
                        t->truelepE=trueE;
                        t->truelepType=mc->IdHEP;
                        break;
                }
        }
                
                

        //find em energy
        t->emenergy=0;
        for(unsigned int i=0;i<poh->mctrue.stdhep.size();i++)
        {
                NtpMCStdHep *mc = &(poh->mctrue.stdhep[i]);
                if(mc->IstHEP!=1)continue;
                
                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);

                

                if(TMath::Abs(mc->IdHEP)==11 || TMath::Abs(mc->IdHEP)==111 || TMath::Abs(mc->IdHEP)==22)
                {
                        t->emenergy+=trueE;
                }
        }
                
                



        
        std::map<double,int> particle_e_index;
        for(unsigned int i=0;i<poh->particles3d.size();i++)
        {
                Particle3D * p = &poh->particles3d[i];

                particle_e_index.insert(std::pair<double,int>(p->sum_e,i));
        
        }
        
        std::map<double,int>::reverse_iterator it;
        for(it=particle_e_index.rbegin();it!=particle_e_index.rend();it++)
        {
                Particle3D * p = &poh->particles3d[it->second];
                MSG("TruthCompareAna",Msg::kDebug)<<"looking for match for recoed particle type "<<p->particletype<<" with e "<<p->sum_e/meupergev<<" ("<<p->start_u<<", "<<p->start_v<<", "<<p->start_z<<")\n";
                std::vector<int> truematch = MatchRecoParticle(p,poh);
                
                reco_visible_e+=p->sum_e/meupergev;
                
                if(truematch.size()>0)
                {
                        t->reco_to_true_match.push_back(truematch);
                        t->reco_matched.push_back(it->second);
                        t->matchangle=matchangle;
                
                        MSG("TruthCompareAna",Msg::kDebug)<<"reco particle "<<p->particletype<<" with e "<<p->sum_e/meupergev<<" matched to:\n";
                        for(unsigned int i=0;i<truematch.size();i++)
                        {
                                NtpMCStdHep *mc = &(poh->mctrue.stdhep[truematch[i]]);

                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);



                                MSG("TruthCompareAna",Msg::kDebug)<<"\t true particle "<<mc->IdHEP<<" with e "<<trueE<<"\n";
                        }
                
                
                        //mark the true particles as matched...
                        for(unsigned int k=0;k<truematch.size();k++)
                                for(unsigned int p=0;p<matched.size();p++)
                                        if(true_idx[p]==truematch[k])matched[p]=truematch[k];
                        reco_matched_visible_e=p->sum_e/meupergev;
                        
                        int left=0;for(unsigned int i=0;i<matched.size();i++)if(!matched[i])left++;
                        MSG("TruthCompareAna",Msg::kDebug)<<"true particles left "<<left<<"\n";
                }else{
                        t->reco_not_matched.push_back(it->second);
                }
                
                
        }       
        

                for(unsigned int i=0;i<true_idx.size();i++)
                {
                        if(matched[i])continue;//already used!
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[true_idx[i]]);


                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);


                        MSG("TruthCompareAna",Msg::kDebug)<<"Did NOT match true type "<<mc->IdHEP<< " with e "<<trueE<<"\n";
                }                       

                
        
        
        //we've done the matching...
        //compute some statistical values about the matching now...
        ComputeStatistics(poh,t);
        
        
        t->stage=1;
        
        MSG("TruthCompareAna",Msg::kDebug)<<"\n";
}



void TruthCompareAna::ComputeStatistics(ParticleObjectHolder *poh, TruthCompare *t)
{

         t->frac_particles_found=0;

        for(unsigned int i=0;i<matched.size();i++)
                if(!matched[i])t->true_not_recoed.push_back(true_idx[i]);
                else t->frac_particles_found++;

        t->particles_matched_to_true=t->reco_to_true_match.size();
        t->possible_true_particles=matched.size();
        
        if(t->frac_particles_found)t->frac_particles_found  /= t->possible_true_particles;
        
        t->reco_matched_visible_e=reco_matched_visible_e;
        t->reco_visible_e=reco_visible_e;
        
        t->true_visible_e=poh->mctrue.visible_energy;

        if(t->true_visible_e)t->frac_e_found=t->reco_matched_visible_e/t->true_visible_e;


        for(unsigned int i=0;i<t->reco_matched.size();i++)
        {
                Particle3D * p = &poh->particles3d[t->reco_matched[i]];
                MSG("TruthCompareAna",Msg::kDebug)<<"Matched reco particle type"<<p->particletype<<" with e "<<p->sum_e/meupergev<<"\n";
                std::vector<int> trues=t->reco_to_true_match[i];
                for(unsigned int j=0;j<trues.size();j++)
                {
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[trues[j]]);

                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);


                        MSG("TruthCompareAna",Msg::kDebug)<<"\tto type: "<<mc->IdHEP<<" e "<<trueE<<"\n";
                }
                
        }


}


//prepare the lists of the true particles that we should attempt to match
void TruthCompareAna::ProcessTrueParticles(ParticleObjectHolder *poh)
{
                ////calculate expected particles.......
//              printf("\nstdhep....run %d snarl %d event %d\n",poh->GetHeader().GetRun() , poh->GetHeader().GetSnarl() , poh->GetHeader().GetEvent());


                for(unsigned int i=0;i<poh->mctrue.stdhep.size();i++)
                {
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[i]);
                        
                        if(mc->IstHEP!=1)continue;


                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);

                        
                        //remove particles too small to consider!
                        if(TMath::Abs(mc->IdHEP)==13 || TMath::Abs(mc->IdHEP)==211)
                        {
                                if(trueE<0.5)continue;
                        }else{
                                if(trueE<0.2)continue;
                        }
                        
                        if(TMath::Abs(mc->IdHEP)==12 || TMath::Abs(mc->IdHEP)==14 || TMath::Abs(mc->IdHEP)==16)continue;//can't match a neutrino!
                        
                        //we'll consider it!
                        true_idx.push_back(i);
                        matched.push_back(0);
                        
        //              MSG("TruthCompareAna",Msg::kDebug)<<"true type "<<mc->IdHEP<<" e "<<trueE<<" from "<<mc->dethit[0].vtx[0]<<" "<<mc->dethit[0].vtx[1]<<" "<<mc->dethit[0].vtx[2]<<" to "<<mc->dethit[1].vtx[0]<<" "<<mc->dethit[1].vtx[1]<<" "<<mc->dethit[1].vtx[2]<<"\n";
                }
//              printf("\n");
                
                
                ///////////////////////
                
        

}

std::vector<int> TruthCompareAna::MatchRecoElectron(Particle3D *p,ParticleObjectHolder *poh)
{
        std::vector<int>truematch;

                        //calculate some numbers about this particle....
                MSG("TruthCompareAna",Msg::kDebug)<<"attempting to match elec\n";       

                        double sum_z2=0;
                        double sum_z=0;
                        double sum_u=0;
                        double sum_zu=0;                        
                        double sum_v=0;
                        double sum_zv=0;
                        int n=0;        
                                                
                        for(int i=0;i<p->entries;i++)
                        {
                                sum_z2+=p->z[i]*p->z[i];
                                sum_z+=p->z[i];
                                sum_zu+=p->z[i]*p->u[i];
                                sum_u+=p->u[i];
                                sum_zv+=p->z[i]*p->v[i];
                                sum_v+=p->v[i]; 
                                n++;                    
                        }
                
                        if(n<2)return truematch;
        
                        double off_u=0, slope_u=0, off_v=0, slope_v=0;
                        if(n*sum_z2-sum_z*sum_z>1e-10)
                        {
                                off_u=  (sum_u*sum_z2 - sum_z*sum_zu)/(n*sum_z2-sum_z*sum_z);
                                slope_u=  (n*sum_zu - sum_z*sum_u)/(n*sum_z2-sum_z*sum_z);
                                                
                                off_v=  (sum_v*sum_z2 - sum_z*sum_zv)/(n*sum_z2-sum_z*sum_z);
                                slope_v=  (n*sum_zv - sum_z*sum_v)/(n*sum_z2-sum_z*sum_z);
                        }

                        ////////////////////
                        double dz=p->end_z-p->start_z;
                        double du=dz*slope_u;
                        double dv=dz*slope_v;                   
                        
                        double r2 = du*du+dv*dv+dz*dz;
                        r2=sqrt(r2);
                        double pu= r2 ? p->sum_e/meupergev * du/r2 : 0;
                        double pv= r2 ? p->sum_e/meupergev * dv/r2 : 0;
                        double pz= r2 ? p->sum_e/meupergev * dz/r2 : 0;
                        
        
        MSG("TruthCompareAna",Msg::kDebug)<<"recoed elec puvz "<<pu<<" "<<pv<<" "<<pz<<"\n";
                        
                        //look for large true electrons
                        std::map<double,int>true_elec;
                        std::map<double,int>true_pi0;
                        std::map<double,int>true_gamma;
                        std::map<double,int>true_proton;
                                                                        
                for(unsigned int i=0;i<true_idx.size();i++)
                {
                        if(matched[i])continue;//already used!
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[true_idx[i]]);

                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);


                        if(TMath::Abs(mc->IdHEP)==11)true_elec.insert(std::pair<double,int>(trueE,true_idx[i]));
                        if(TMath::Abs(mc->IdHEP)==111)true_pi0.insert(std::pair<double,int>(trueE,true_idx[i]));
                        if(TMath::Abs(mc->IdHEP)==22)true_gamma.insert(std::pair<double,int>(trueE,true_idx[i]));
                        if(TMath::Abs(mc->IdHEP)==2212)true_proton.insert(std::pair<double,int>(trueE,true_idx[i]));
                }                       

                MSG("TruthCompareAna",Msg::kDebug)<<"true elec: "<<true_elec.size()<<" pi0: "<<true_pi0.size()<<" gamma: "<<true_gamma.size()<<" prot "<<true_proton.size()<<"\n";
                
                std::map<double,int>::reverse_iterator it;
                
                double trueE_matched=0;
                
                //need to add up true vectors to get potential angle!
                double cur_pu=0;
                double cur_pv=0;
                double cur_pz=0;
                
                for(it=true_elec.rbegin();it!=true_elec.rend();it++)
                {               
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[it->second]);

                        double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0)+cur_pu;
                        double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0)+cur_pv;
                        double true_pz=mc->p4[2]+cur_pz;

        MSG("TruthCompareAna",Msg::kDebug)<<"current true puvz "<<true_pu<<" "<<true_pv<<" "<<true_pz<<"\n";
                        
                        double reco_r = sqrt(pu*pu+pv*pv+pz*pz);        
                        double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
                                        
                        double angle = reco_r ? acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) ) : 0;
                        
                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);


                        
                        MSG("TruthCompareAna",Msg::kDebug)<<"comparing with angle "<<angle <<" trueE "<<trueE<<" recoE "<<p->sum_e/meupergev<<"\n";

                        //require somewhat the same direction
                        if(angle<0.5)
                        {
                                //matched
                                truematch.push_back(it->second);
                                trueE_matched+=trueE;
                                MSG("TruthCompareAna",Msg::kDebug)<<" matched!\n";
                                cur_pu=true_pu;
                                cur_pv=true_pv;
                                cur_pz=true_pz;
                        }else{
                                MSG("TruthCompareAna",Msg::kDebug)<<"\n";
                        }
                        
                        if( trueE_matched> 0.9 *p->sum_e/meupergev)
                        {
                                matchangle=angle;
                                return truematch;
                        }
                }
                
                
                for(it=true_pi0.rbegin();it!=true_pi0.rend();it++)
                {               
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[it->second]);

                        double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0)+cur_pu;
                        double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0)+cur_pv;
                        double true_pz=mc->p4[2]+cur_pz;

        MSG("TruthCompareAna",Msg::kDebug)<<"current true puvz "<<true_pu<<" "<<true_pv<<" "<<true_pz<<"\n";
                        
                        double reco_r = sqrt(pu*pu+pv*pv+pz*pz);        
                        double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
                                        
                        double angle = reco_r ? acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) ) : 0;
                        
                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }

                        double trueE=sqrt(trueE2);

                        
                        MSG("TruthCompareAna",Msg::kDebug)<<"comparing with angle "<<angle <<" trueE "<<trueE<<" recoE "<<p->sum_e/meupergev<<"\n";

                        //require somewhat the same direction
                        if(angle<0.5)
                        {
                                //matched
                                truematch.push_back(it->second);
                                trueE_matched+=trueE;
                                MSG("TruthCompareAna",Msg::kDebug)<<" matched!\n";
                                cur_pu=true_pu;
                                cur_pv=true_pv;
                                cur_pz=true_pz;
                        }
                        
                        if( trueE_matched> 0.9 *p->sum_e/meupergev)
                        {
                                matchangle=angle;
                                return truematch;
                        }
                }
                                
                
                
                for(it=true_gamma.rbegin();it!=true_gamma.rend();it++)
                {               
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[it->second]);

                        double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0)+cur_pu;
                        double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0)+cur_pv;
                        double true_pz=mc->p4[2]+cur_pz;

        MSG("TruthCompareAna",Msg::kDebug)<<"current true puvz "<<true_pu<<" "<<true_pv<<" "<<true_pz<<"\n";
                        
                        double reco_r = sqrt(pu*pu+pv*pv+pz*pz);        
                        double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
                                        
                        double angle = reco_r ? acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) ) : 0;
                        

                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);
                
                        
                        MSG("TruthCompareAna",Msg::kDebug)<<"comparing with angle "<<angle <<" trueE "<<trueE<<" recoE "<<p->sum_e/meupergev<<"\n";

                        //require somewhat the same direction
                        if(angle<0.5)
                        {
                                //matched
                                truematch.push_back(it->second);
                                trueE_matched+=trueE;
                                MSG("TruthCompareAna",Msg::kDebug)<<" matched!\n";
                                cur_pu=true_pu;
                                cur_pv=true_pv;
                                cur_pz=true_pz;
                        }
                        
                        if( trueE_matched> 0.9 *p->sum_e/meupergev){
                                matchangle=angle;
                                return truematch;
                        }
                }
                                                

                
                for(it=true_proton.rbegin();it!=true_proton.rend();it++)
                {               
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[it->second]);

                        double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0)+cur_pu;
                        double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0)+cur_pv;
                        double true_pz=mc->p4[2]+cur_pz;

        MSG("TruthCompareAna",Msg::kDebug)<<"current true puvz "<<true_pu<<" "<<true_pv<<" "<<true_pz<<"\n";
                        
                        double reco_r = sqrt(pu*pu+pv*pv+pz*pz);        
                        double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
                                        
                        double angle = reco_r ? acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) ) : 0;
                        
                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);

                        
                        MSG("TruthCompareAna",Msg::kDebug)<<"comparing with angle "<<angle <<" trueE "<<trueE<<" recoE "<<p->sum_e/meupergev<<"\n";

                        //require somewhat the same direction
                        if(angle<0.5)
                        {
                                //matched
                                truematch.push_back(it->second);
                                trueE_matched+=trueE;
                                MSG("TruthCompareAna",Msg::kDebug)<<" matched!\n";
                                cur_pu=true_pu;
                                cur_pv=true_pv;
                                cur_pz=true_pz;
                        }
                        
                        if( trueE_matched> 0.9 *p->sum_e/meupergev){
                                matchangle=angle;
                                return truematch;
                        }
                }
                                        
                        
        return truematch;               
}



std::vector<int> TruthCompareAna::MatchRecoMuon(Particle3D *p,ParticleObjectHolder *poh)
{
        std::vector<int>truematch;

                        //calculate some numbers about this particle....
                MSG("TruthCompareAna",Msg::kDebug)<<"attempting to match muon\n";       

                        double sum_z2=0;
                        double sum_z=0;
                        double sum_u=0;
                        double sum_zu=0;                        
                        double sum_v=0;
                        double sum_zv=0;
                        int n=0;        
                                                
                        int ent=p->entries<10?p->entries:10;                    
                        for(int i=0;i<ent;i++)
                        {
                                sum_z2+=p->z[i]*p->z[i];
                                sum_z+=p->z[i];
                                sum_zu+=p->z[i]*p->u[i];
                                sum_u+=p->u[i];
                                sum_zv+=p->z[i]*p->v[i];
                                sum_v+=p->v[i]; 
                                n++;                    
                        }
                
                        if(n<2)return truematch;

        
                        double denom = (n*sum_z2-sum_z*sum_z);

                        double off_u=0, slope_u=0, off_v=0, slope_v =0;

                        if(denom)
                        {
                                off_u=  (sum_u*sum_z2 - sum_z*sum_zu)/(n*sum_z2-sum_z*sum_z);
                                slope_u=  (n*sum_zu - sum_z*sum_u)/(n*sum_z2-sum_z*sum_z);

                                off_v=  (sum_v*sum_z2 - sum_z*sum_zv)/(n*sum_z2-sum_z*sum_z);
                                slope_v=  (n*sum_zv - sum_z*sum_v)/(n*sum_z2-sum_z*sum_z);
                        }                       


                        ////////////////////
                        double dz=p->end_z-p->start_z;
                        double du=dz*slope_u;
                        double dv=dz*slope_v;                   
                        
                        double r2 = du*du+dv*dv+dz*dz;
                        r2=sqrt(r2);
                        double pu= r2 ? p->sum_e/meupergev * du/r2 : 0;
                        double pv= r2 ? p->sum_e/meupergev * dv/r2 : 0;
                        double pz= r2 ? p->sum_e/meupergev * dz/r2 : 0;
                        
                        
                        //look for large true electrons
                        std::map<double,int>true_muon;
                        std::map<double,int>true_charged_pion;
                                                                        
                for(unsigned int i=0;i<true_idx.size();i++)
                {
                        if(matched[i])continue;//already used!
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[true_idx[i]]);
                        

                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);


                        if(TMath::Abs(mc->IdHEP)==13)true_muon.insert(std::pair<double,int>(trueE,true_idx[i]));
                        if(TMath::Abs(mc->IdHEP)==211 || TMath::Abs(mc->IdHEP)==-211) 
                                true_charged_pion.insert(std::pair<double,int>(trueE,true_idx[i]));
                }


                std::map<double,int>::reverse_iterator it;
                
                double trueE_matched=0;
                
                for(it=true_muon.rbegin();it!=true_muon.rend();it++)
                {               
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[it->second]);

                        double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0);
                        double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0);
                        double true_pz=mc->p4[2];
                        
                        double reco_r = sqrt(pu*pu+pv*pv+pz*pz);        
                        double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
                                        
                        double angle = reco_r ? acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) ) : 0;
                        
                
                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);


                        
                        MSG("TruthCompareAna",Msg::kDebug)<<"comparing with angle "<<angle <<" trueE "<<trueE<<" recoE "<<p->sum_e/meupergev<<"\n";
                        MSG("TruthCompareAna",Msg::kDebug)<<"reco p "<<pu<<" "<<pv<<" "<<pz<<" true p "<<true_pu<<" "<<true_pv<<" "<<true_pz<<"\n";
                        
                        //require somewhat the same direction
                        if(angle<0.5)
                        {
                                //matched
                                truematch.push_back(it->second);
                                trueE_matched+=trueE;
                                MSG("TruthCompareAna",Msg::kDebug)<<" matched!\n";
                        }
                        
                        if( trueE_matched> 0.9 *p->sum_e/meupergev){
                                matchangle=angle;
                                return truematch;
                        }
                }
                
                for(it=true_charged_pion.rbegin();it!=true_charged_pion.rend();it++)
                {               
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[it->second]);

                        double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0);
                        double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0);
                        double true_pz=mc->p4[2];
                        
                        double reco_r = sqrt(pu*pu+pv*pv+pz*pz);        
                        double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
                                        
                        double angle = reco_r ? acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) ) : 0;
                        
                        
                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);
        
                
                        MSG("TruthCompareAna",Msg::kDebug)<<"comparing with angle "<<angle <<" trueE "<<trueE<<" recoE "<<p->sum_e/meupergev<<"\n";
                        MSG("TruthCompareAna",Msg::kDebug)<<"reco p "<<pu<<" "<<pv<<" "<<pz<<" true p "<<true_pu<<" "<<true_pv<<" "<<true_pz<<"\n";
                        
                        
                        //require somewhat the same direction
                        if(angle<0.5)
                        {
                                //matched
                                truematch.push_back(it->second);
                                trueE_matched+=trueE;
                                MSG("TruthCompareAna",Msg::kDebug)<<" matched!\n";
                        }
                        
                        if( trueE_matched> 0.9 *p->sum_e/meupergev){
                                matchangle=angle;
                                return truematch;
                        }
                }                               
                        
        return truematch;               
}



std::vector<int> TruthCompareAna::MatchRecoProton(Particle3D *p,ParticleObjectHolder *poh)
{
        std::vector<int>truematch;

                        //calculate some numbers about this particle....
                MSG("TruthCompareAna",Msg::kDebug)<<"attempting to match muon\n";       

                        double sum_z2=0;
                        double sum_z=0;
                        double sum_u=0;
                        double sum_zu=0;                        
                        double sum_v=0;
                        double sum_zv=0;
                        int n=0;        
                                                
                        int ent=p->entries<10?p->entries:10;                    
                        
                        double off_u=  0;
                        double slope_u=  0;
                        
                        
                        double off_v=  0;
                        double slope_v=0;
                        
                        if(ent>1)
                        {
                                for(int i=0;i<ent;i++)
                                {
                                        sum_z2+=p->z[i]*p->z[i];
                                        sum_z+=p->z[i];
                                        sum_zu+=p->z[i]*p->u[i];
                                        sum_u+=p->u[i];
                                        sum_zv+=p->z[i]*p->v[i];
                                        sum_v+=p->v[i]; 
                                        n++;                    
                                }
                        
                                if(n*sum_z2-sum_z*sum_z)
                                {
                                        off_u=  (sum_u*sum_z2 - sum_z*sum_zu)/(n*sum_z2-sum_z*sum_z);
                                        slope_u=  (n*sum_zu - sum_z*sum_u)/(n*sum_z2-sum_z*sum_z);
        
                                        off_v=  (sum_v*sum_z2 - sum_z*sum_zv)/(n*sum_z2-sum_z*sum_z);
                                        slope_v=  (n*sum_zv - sum_z*sum_v)/(n*sum_z2-sum_z*sum_z);
                                }
                        }else if(TMath::Abs(p->z[0]-poh->event.vtx_z)>1e-5){
                                slope_u=(p->u[0]-poh->event.vtx_u)/(p->z[0]-poh->event.vtx_z);
                                slope_v=(p->v[0]-poh->event.vtx_v)/(p->z[0]-poh->event.vtx_z);
                                
                                off_u=poh->event.vtx_u-slope_u*poh->event.vtx_z;
                                off_v=poh->event.vtx_v-slope_v*poh->event.vtx_z;
                                
                        }
                        
                        

                        ////////////////////
                        double dz=p->end_z-p->start_z;
                        double du=dz*slope_u;
                        double dv=dz*slope_v;                   
                        
                        double r2 = du*du+dv*dv+dz*dz;
                        r2=sqrt(r2);
                        double pu= r2? p->sum_e/meupergev * du/r2:0;
                        double pv= r2? p->sum_e/meupergev * dv/r2:0;
                        double pz= r2? p->sum_e/meupergev * dz/r2:0;
                        
                        
                        //look for large true electrons
                        std::map<double,int>true_proton;
                        std::map<double,int>true_neutron;
                                                                        
                for(unsigned int i=0;i<true_idx.size();i++)
                {
                        if(matched[i])continue;//already used!
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[true_idx[i]]);
                        
                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0   = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);

                        if(TMath::Abs(mc->IdHEP)==2212)true_proton.insert(std::pair<double,int>(trueE,true_idx[i]));
                        if(TMath::Abs(mc->IdHEP)==2112 ) true_neutron.insert(std::pair<double,int>(trueE,true_idx[i]));
                }


                std::map<double,int>::reverse_iterator it;
                
                double trueE_matched=0;
                
                for(it=true_proton.rbegin();it!=true_proton.rend();it++)
                {               
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[it->second]);

                        double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0);
                        double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0);
                        double true_pz=mc->p4[2];
                        
                        double reco_r = sqrt(pu*pu+pv*pv+pz*pz);        
                        double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
                                        
                        double angle = reco_r ? acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) ) : 0;
                

                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);
        

                        
                        MSG("TruthCompareAna",Msg::kDebug)<<"comparing with angle "<<angle <<" trueE "<<trueE<<" recoE "<<p->sum_e/meupergev<<"\n";
                        MSG("TruthCompareAna",Msg::kDebug)<<"reco p "<<pu<<" "<<pv<<" "<<pz<<" true p "<<true_pu<<" "<<true_pv<<" "<<true_pz<<"\n";
                        
                        //require somewhat the same direction
                        if(angle<0.5)
                        {
                                //matched
                                truematch.push_back(it->second);
                                trueE_matched+=trueE;
                                MSG("TruthCompareAna",Msg::kDebug)<<" matched!\n";
                        }
                        
                        if( trueE_matched> 0.9 *p->sum_e/meupergev){
                                matchangle=angle;
                                return truematch;
                        }
                }
                
                for(it=true_neutron.rbegin();it!=true_neutron.rend();it++)
                {               
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[it->second]);

                        double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0);
                        double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0);
                        double true_pz=mc->p4[2];
                        
                        double reco_r = sqrt(pu*pu+pv*pv+pz*pz);        
                        double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
                                        
                        double angle = reco_r ? acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) ) : 0;

                        double trueE2 = mc->p4[3]*mc->p4[3]-mc->mass*mc->mass;
                        if(trueE2<0)
                        {
                                trueE2=-trueE2;
                                MSG("TruthCompareAna",Msg::kError)<<"true E2 < 0 = "<<trueE2<<"\n";
                        }
                        double trueE=sqrt(trueE2);

                        

                        
                        MSG("TruthCompareAna",Msg::kDebug)<<"comparing with angle "<<angle <<" trueE "<<trueE<<" recoE "<<p->sum_e/meupergev<<"\n";
                        MSG("TruthCompareAna",Msg::kDebug)<<"reco p "<<pu<<" "<<pv<<" "<<pz<<" true p "<<true_pu<<" "<<true_pv<<" "<<true_pz<<"\n";
                        
                        
                        //require somewhat the same direction
                        if(angle<0.5)
                        {
                                //matched
                                truematch.push_back(it->second);
                                trueE_matched+=trueE;
                                MSG("TruthCompareAna",Msg::kDebug)<<" matched!\n";
                        }
                        
                        if( trueE_matched> 0.9 *p->sum_e/meupergev){
                                matchangle=angle;
                                return truematch;
                        }
                }                               
                        
        return truematch;               
}





std::vector<int> TruthCompareAna::MatchRecoParticle(Particle3D *p, ParticleObjectHolder *poh)
{




                //now iterate over true particles for consideration........
                //this may be done differently depending on the type of recoed particle!!!
                std::vector<int>truematch;

                        if(p->particletype==Particle3D::electron ||p->particletype==Particle3D::other)truematch= MatchRecoElectron(p,poh);

                        if(p->particletype==Particle3D::muon)truematch= MatchRecoMuon(p,poh);
                
                        if(p->particletype==Particle3D::proton || p->particletype==Particle3D::neutron) truematch= MatchRecoProton(p,poh);
                
                return truematch;
                
                
/*              
                for(unsigned int i=0;i<true_idx.size();i++)
                {
                        if(matched[i])continue;//already used!
                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[true_idx[i]]);
                        
                
                }
*/

/*
                        if(p->particletype==Particle3D::electron)
                        {
                                printf("largest elec %d with e %f  (%f %f %f) p3 (%f %f %f)\n",p->particletype,p->sum_e/meupergev,p->start_u,p->start_v,p->start_z,pu,pv,pz);
                                if(true_elec_idx>-1)
                                {
                                        NtpMCStdHep *mc = &(poh->mctrue.stdhep[true_elec_idx]);
                        
                                        printf("compare to true elec %d with e %f at %d  (%f %f %f)  p3 (%f %f %f)\n",mc->IdHEP,sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass),mc->IstHEP,(mc->vtx[0]+mc->vtx[1])/sqrt(2.0),(-mc->vtx[0]+mc->vtx[1])/sqrt(2.0),mc->vtx[2],(mc->p4[0]+mc->p4[1])/sqrt(2.0),(-mc->p4[0]+mc->p4[1])/sqrt(2.0),mc->p4[2]);             
                                        //calculate the angle!!!!
                                        
                                        //angle between p vectors....
                                        
                                        double true_pu=(mc->p4[0]+mc->p4[1])/sqrt(2.0);
                                        double true_pv=(-mc->p4[0]+mc->p4[1])/sqrt(2.0);
                                        double true_pz=mc->p4[2];
                                        
                                        double reco_r = sqrt(pu*pu+pv*pv+pz*pz);        
                                        double true_r = sqrt(true_pu*true_pu+true_pv*true_pv+true_pz*true_pz);
                                        
                                        double angle = acos( (pu*true_pu+pv*true_pv+pz*true_pz) / (reco_r*true_r) );
                                        printf("angle between: %f\n",angle);
                                        
                                //      hist_prim_angle_vs_prim_eres[type]->Fill(angle,(sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass)-p->sum_e/meupergev)/(sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass)),weight);          
                                //      hist_prim_angle_vs_prim_true[type]->Fill(angle,sqrt(mc->p4[3]*mc->p4[3]-mc->mass*mc->mass));    
                                        
                                        
                                }
                        }
*/                      
        
}
