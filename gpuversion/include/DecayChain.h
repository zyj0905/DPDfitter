//DecayChain class
#ifndef DECAYCHAIN_H 
#define DECAYCHAIN_H

#include "Particle.h"

class DecayChain{
    public:
        particle intermediate;
        int idx_isobar;
        int N1_LS,N2_LS;
        int L1[10]; int L2[10];
        double S1[10]; double S2[10];
        double rho1_coff[10]; double phi1_coff[10];
        double rho2_coff[10]; double phi2_coff[10];
        
        DecayChain(){
            N1_LS = 0; N2_LS = 0; idx_isobar = 0;
            for(int i=0;i<10;i++){
                L1[i] = 0; S1[i] = 0;
                L2[i] = 0; S2[i] = 0;
                rho1_coff[i] = 0; phi1_coff[i] = 0;
                rho2_coff[i] = 0; phi2_coff[i] = 0;
            }
        }

        DecayChain(particle m_intermediate, int m_idx_isobar){
            intermediate = m_intermediate;
            idx_isobar = m_idx_isobar;
            N1_LS = 0; N2_LS = 0;
            for(int i=0;i<10;i++){
                L1[i] = 0; S1[i] = 0;
                L2[i] = 0; S2[i] = 0;
                rho1_coff[i] = 0; phi1_coff[i] = 0;
                rho2_coff[i] = 0; phi2_coff[i] = 0;
            }
        }

        DecayChain& operator=(const DecayChain& other){
            if(this == &other) return *this;
            
            intermediate = other.intermediate;
            idx_isobar = other.idx_isobar;
            N1_LS = other.N1_LS;
            N2_LS = other.N2_LS;
            
            for(int i=0;i<10;i++){
                L1[i] = other.L1[i]; S1[i] = other.S1[i];
                L2[i] = other.L2[i]; S2[i] = other.S2[i];
                rho1_coff[i] = other.rho1_coff[i]; phi1_coff[i] = other.phi1_coff[i];
                rho2_coff[i] = other.rho2_coff[i]; phi2_coff[i] = other.phi2_coff[i];
            }
            return *this;
        }

        void Add_LS1(int l1, double s1){L1[N1_LS] = l1;S1[N1_LS] = s1;N1_LS++;}
        void Add_LS2(int l2, double s2){L2[N2_LS] = l2;S2[N2_LS] = s2;N2_LS++;}

        void Cal_LS_auto(particle particle_list[4]){
            if(intermediate.p_parity==0){cout<<"Intermediate has not been set!"<<endl;return;}
            std::cout<<intermediate.name<<":"<<endl;
            //The decay 0-> isobar + intermediate
            int m_N_LS1 = 0; double m_L1[10]; double m_S1[10];
            particle_list[0].Get_decay_LS(intermediate,particle_list[idx_isobar],m_L1,m_S1,m_N_LS1);
            //The decay intermedia -> A + B
            int m_N_LS2 = 0; double m_L2[10]; double m_S2[10];
            if(idx_isobar==1) intermediate.Get_decay_LS(particle_list[2],particle_list[3],m_L2,m_S2,m_N_LS2);
            if(idx_isobar==2) intermediate.Get_decay_LS(particle_list[1],particle_list[3],m_L2,m_S2,m_N_LS2);
            if(idx_isobar==3) intermediate.Get_decay_LS(particle_list[1],particle_list[2],m_L2,m_S2,m_N_LS2);
            //Add all and set coff to DeviceComplex(1.0,0.0)
            for(int idx_LS1=0;idx_LS1<m_N_LS1;idx_LS1++){
                Add_LS1(m_L1[idx_LS1],m_S1[idx_LS1]);rho1_coff[idx_LS1]=1;phi1_coff[idx_LS1]=0;
                printf("L1:%.1f, S1:%.1f\n",m_L1[idx_LS1],m_S1[idx_LS1]);
            }
            for(int idx_LS2=0;idx_LS2<m_N_LS2;idx_LS2++){
                Add_LS2(m_L2[idx_LS2],m_S2[idx_LS2]);rho2_coff[idx_LS2]=1;phi2_coff[idx_LS2]=0;
                printf("L2:%.1f, S2:%.1f\n",m_L2[idx_LS2],m_S2[idx_LS2]);
            }
        }

        __device__ DeviceComplex Get_LScoff(int type_ls, int idx_ls){
            double rho(0.),phi(0.);
            if(type_ls==1){rho = rho1_coff[idx_ls]; phi = phi1_coff[idx_ls];}
            if(type_ls==2){rho = rho2_coff[idx_ls]; phi = phi2_coff[idx_ls];}
            return DeviceComplex(rho*cos(phi),rho*sin(phi));
        }
};
#endif //DECAYCHAIN_H