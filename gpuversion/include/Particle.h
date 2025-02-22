//particle class
#ifndef PARTICLE_H 
#define PARTICLE_H

#include "Dynamic.h"

class particle{
    public:
        string name;
        double spin;
        int p_parity;
        int c_parity;
        Dynamic dynamic;

        particle(){
            name = "void";
            spin = 0;
            p_parity = 0;
            c_parity = 0;
        }

        particle(string m_name, double m_spin, double m_p_parity, double m_c_parity){
            name = m_name;
            spin = m_spin;
            p_parity = m_p_parity;
            c_parity = m_c_parity;
        }

        void Set_Dynamic(Dynamic m_dynamic){dynamic = m_dynamic;}

        void Update_Dynamic_par(int idx, double par_new){
            if(idx>=dynamic.num_par){
                cout<<"Update_Dynamic_par: index out of range"<<endl;
                return;
            }
            dynamic.pars[idx] = par_new;
        }

        int Get_Num_Dynamic_par(){return dynamic.num_par;}

        particle &operator=(const particle& other){
            if(this == &other)
                return *this;

            spin = other.spin;
            name = other.name;
            p_parity = other.p_parity;
            c_parity = other.c_parity;
            dynamic = other.dynamic;
            return *this;
        }

        void Get_decay_LS(particle dau_A, particle dau_B, double m_L[10], double m_S[10], int& nLS){
            
            nLS = 0;

            if(dau_A.c_parity * dau_B.c_parity != c_parity && dau_A.c_parity * dau_B.c_parity != 0) return;

            for(double S=fabs(dau_A.spin-dau_B.spin);S<=(dau_A.spin+dau_B.spin);S++){
                for(double L=fabs(spin-S);L<=(spin+S);L++){
                    //P-parity
                    if(dau_A.p_parity * dau_B.p_parity == p_parity){
                        if(int(L)%2!=0) continue;
                    }
                    else{
                        if(int(L)%2==0) continue;
                    }
                    if(spin>(L+S) || spin < fabs(L-S)) continue;

                    m_L[nLS] = L;
                    m_S[nLS] = S;
                    nLS++;
                    
                }
            }
        }
};
#endif //PARTICLE_H
