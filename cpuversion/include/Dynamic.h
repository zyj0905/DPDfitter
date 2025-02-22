//Dynamic class
#ifndef DYNAMIC_H 
#define DYNAMIC_H

#include "DeviceComplex.h"
#include "Event.h"
#include <iostream>
using namespace std;

class Dynamic{
    public:
        int type;
        int num_par;
        double pars[100];

        Dynamic(){type = 1; num_par = 0;}

        Dynamic(int m_type, int m_num_par, double* m_pars){
            type = m_type; num_par = m_num_par;
            for(int i=0;i<num_par;i++){pars[i] = m_pars[i];}
        }

        Dynamic &operator=(const Dynamic& other){
            if(this == &other) return *this;
            type = other.type;
            num_par = other.num_par;
            for(int i=0;i<num_par;i++){pars[i] = other.pars[i];}
            return *this;
        }

        DeviceComplex eval(double sqrt_s, Event* evt, int idx_isobar){
            if(type==0){return fun0(sqrt_s);} //constant width BW
            if(type==1){return fun1(sqrt_s);} //constant
            if(type==2){return fun2(sqrt_s,evt,idx_isobar);} //BW with running width
            if(type==3){return fun3(sqrt_s,evt,idx_isobar);} // Flatte-like formula S-wave
            printf("NO available type in eval!\n");
            return 0.0;
        }

        DeviceComplex fun0(double sqrt_s){
            double mass = pars[0];
            double width = pars[1];
            return 1.0/DeviceComplex(mass*mass-sqrt_s*sqrt_s,-mass*width);
        }

        DeviceComplex fun1(double sqrt_s){return 1.0;}


        double break_mom_zero(double M, double m1, double m2){
            double lambda = pow(M,4) + pow(m1,4) + pow(m2,4) - 2*pow(M,2)*pow(m1,2) - 2*pow(M,2)*pow(m2,2) - 2*pow(m2,2)*pow(m1,2);
            if(lambda<0){return 0.0;}
            return sqrt(fabs(lambda))/2.0/M;
        }

        double break_mom_abs(double M, double m1, double m2){
            double lambda = pow(M,4) + pow(m1,4) + pow(m2,4) - 2*pow(M,2)*pow(m1,2) - 2*pow(M,2)*pow(m2,2) - 2*pow(m2,2)*pow(m1,2);
            return sqrt(fabs(lambda))/2.0/M;
        }

        DeviceComplex break_mom_complex(double M, double m1, double m2){
            double lambda = pow(M,4) + pow(m1,4) + pow(m2,4) - 2*pow(M,2)*pow(m1,2) - 2*pow(M,2)*pow(m2,2) - 2*pow(m2,2)*pow(m1,2);
            DeviceComplex mysqrt(0,0);
            if(lambda<0){mysqrt = DeviceComplex(0,sqrt(fabs(lambda)));}
            else{mysqrt = DeviceComplex(sqrt(fabs(lambda)),0);}
            return mysqrt/2.0/M;
        }

        double Blatt_Weisskopf_factor(double Q, double Q0, int L){
            if(L==0){return 1.0;}
            if(L==1){return sqrt(2.0/(Q*Q+Q0*Q0));}
            if(L==2){return sqrt(13.0/(Q*Q*Q*Q+3.0*Q0*Q0*Q*Q+9.0*Q0*Q0*Q0*Q0));}
            if(L==3){return sqrt(277.0/(pow(Q,6)+6.0*pow(Q,4)*pow(Q0,2)+45.0*pow(Q,2)*pow(Q0,4)+225.0*pow(Q0,6)));}
            if(L==4){return sqrt(12746.0/(pow(Q,8)+10.0*pow(Q,6)*pow(Q0,2)+135.0*pow(Q,4)*pow(Q0,4)+1575.0*pow(Q,2)*pow(Q0,6)+11025.0*pow(Q0,8)));}
            if(L>4){printf("NO valid Blatt Weisskopf factor!\n");}
            return 0.0;
        }

        DeviceComplex fun2(double sqrt_s, Event* evt, int idx_isobar){
            double mass = pars[0];
            double width = pars[1];
            int L = pars[2];
            double Q0 = pars[3];
            double m1(0.),m2(0.);
            //if(idx_isobar==1){m1 = sqrt(evt->_mass2_dau2);m2 = sqrt(evt->_mass2_dau3);}
            //if(idx_isobar==2){m1 = sqrt(evt->_mass2_dau1);m2 = sqrt(evt->_mass2_dau3);}
            //if(idx_isobar==3){m1 = sqrt(evt->_mass2_dau1);m2 = sqrt(evt->_mass2_dau2);}
            m1 = pars[4];
            m2 = pars[5];
            
            double mom = break_mom_abs(sqrt_s,m1,m2);
            double mom0 = break_mom_abs(mass,m1,m2);
            double Gamma = width*pow(mom/mom0,2*L+1)*(mass/sqrt_s)*pow(Blatt_Weisskopf_factor(mom,Q0,L),2);
            return 1.0/DeviceComplex(mass*mass-sqrt_s*sqrt_s,-mass*Gamma);
        }

        //https://docbes3.ihep.ac.cn/DocDB/0002/000216/030/zc_memo.pdf
        //Eq.~11
        DeviceComplex fun3(double sqrt_s, Event* evt, int idx_isobar){
            double mass = pars[0];
            int Nchannel = (num_par-1)/3;
            DeviceComplex gterm(0,0);
            for(int i=0;i<Nchannel;i++){
                double m1,m2,g;
                g = pars[i*3+1];
                m1 = pars[i*3+2]; m2 = pars[i*3+3];
                DeviceComplex break_mom = break_mom_complex(sqrt_s,m1,m2);
                DeviceComplex phase_space = 2.0*break_mom/sqrt_s;
                gterm = gterm + g*phase_space;
            }
            DeviceComplex result = sqrt_s*sqrt_s - mass*mass + DeviceComplex(0,1)*gterm;
            return 1./result;
        }
};

#endif //DYNAMIC_H