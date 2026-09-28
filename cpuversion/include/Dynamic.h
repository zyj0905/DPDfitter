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
            if(type==1){return fun1(sqrt_s);} //constant no resonant
            if(type==2){return fun2(sqrt_s);} //BW with running width
            if(type==3){return fun3(sqrt_s);} //Flatte-like formula S-wave
            if(type==4){return fun4(sqrt_s);} //Zc lineshape with D*D (in S and D wave) contribution + the others
            if(type==5){return fun5(sqrt_s);} //Dpi S-wave based on the CubicSpline from PhysRevD.94.072001
            if(type==6){return fun6(sqrt_s, evt, idx_isobar);}//Including the triangle-plot into D1D amplitude, according to 2201.08253v2, Eq. 26
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

        DeviceComplex fun2(double sqrt_s){
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
        DeviceComplex fun3(double sqrt_s){
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

        //Zc lineshape with D*D (in S and D wave) contribution + the others
        DeviceComplex fun4(double sqrt_s){
            double mass = pars[0];
            DeviceComplex g01 = DeviceComplex(pars[1]*cos(pars[2]),pars[1]*sin(pars[2]));
            DeviceComplex g21 = DeviceComplex(pars[3]*cos(pars[4]),pars[3]*sin(pars[4]));
            const double Q0 = 0.197321/0.59;
            double mom = break_mom_abs(sqrt_s,1.870,2.010);
            double prefix = (1./8/3.14159)*(mom/pow(sqrt_s,2))*(1/3.);
            //LS coupling
            DeviceComplex h01 = g01*pow(mom,0)*Blatt_Weisskopf_factor(mom,Q0,0);
            DeviceComplex h21 = g21*pow(mom,2)*Blatt_Weisskopf_factor(mom,Q0,2);
            //helicity coupling
            DeviceComplex H_m1_0 = sqrt(1/3.)*h01 + sqrt(2/3.)*h21;
            DeviceComplex H_0_0 = sqrt(1/3.)*h01 - sqrt(2/5.)*h21;
            DeviceComplex H_p1_0 = H_m1_0;//P-parity conversation
            //Decay Xsec
            double Gamma_DstarD = prefix*(H_m1_0.rho2()+H_0_0.rho2()+H_p1_0.rho2());
            double Gamma_others = pars[5];
            DeviceComplex result = sqrt_s*sqrt_s - mass*mass + DeviceComplex(0,1)*sqrt_s*(Gamma_DstarD+Gamma_others);

            return 1./result;
        }

        //Dpi S-wave based on the CubicSpline from PhysRevD.94.072001
        DeviceComplex fun5(double sqrt_s){
            double mass_Dpi[200] = {2.0100, 2.0131, 2.0163, 2.0194, 2.0225, 2.0257, 2.0288, 2.0319, 2.0351, 2.0382, 2.0413, 2.0445, 2.0476, 2.0507, 2.0539, 2.0570, 2.0601, 2.0633, 2.0664, 2.0695, 2.0727, 2.0758, 2.0789, 2.0821, 2.0852, 2.0883, 2.0915, 2.0946, 2.0977, 2.1009, 2.1040, 2.1071, 2.1103, 2.1134, 2.1165, 2.1197, 2.1228, 2.1259, 2.1291, 2.1322, 2.1353, 2.1385, 2.1416, 2.1447, 2.1479, 2.1510, 2.1541, 2.1573, 2.1604, 2.1635, 2.1667, 2.1698, 2.1729, 2.1761, 2.1792, 2.1823, 2.1855, 2.1886, 2.1917, 2.1949, 2.1980, 2.2011, 2.2043, 2.2074, 2.2105, 2.2137, 2.2168, 2.2199, 2.2231, 2.2262, 2.2293, 2.2325, 2.2356, 2.2387, 2.2419, 2.2450, 2.2481, 2.2513, 2.2544, 2.2575, 2.2607, 2.2638, 2.2669, 2.2701, 2.2732, 2.2763, 2.2794, 2.2826, 2.2857, 2.2888, 2.2920, 2.2951, 2.2982, 2.3014, 2.3045, 2.3076, 2.3108, 2.3139, 2.3170, 2.3202, 2.3233, 2.3264, 2.3296, 2.3327, 2.3358, 2.3390, 2.3421, 2.3452, 2.3484, 2.3515, 2.3546, 2.3578, 2.3609, 2.3640, 2.3672, 2.3703, 2.3734, 2.3766, 2.3797, 2.3828, 2.3860, 2.3891, 2.3922, 2.3954, 2.3985, 2.4016, 2.4048, 2.4079, 2.4110, 2.4142, 2.4173, 2.4204, 2.4236, 2.4267, 2.4298, 2.4330, 2.4361, 2.4392, 2.4424, 2.4455, 2.4486, 2.4518, 2.4549, 2.4580, 2.4612, 2.4643, 2.4674, 2.4706, 2.4737, 2.4768, 2.4800, 2.4831, 2.4862, 2.4894, 2.4925, 2.4956, 2.4988, 2.5019, 2.5050, 2.5082, 2.5113, 2.5144, 2.5176, 2.5207, 2.5238, 2.5270, 2.5301, 2.5332, 2.5364, 2.5395, 2.5426, 2.5458, 2.5489, 2.5520, 2.5552, 2.5583, 2.5614, 2.5646, 2.5677, 2.5708, 2.5740, 2.5771, 2.5802, 2.5834, 2.5865, 2.5896, 2.5928, 2.5959, 2.5990, 2.6022, 2.6053, 2.6084, 2.6116, 2.6147, 2.6178, 2.6210, 2.6241, 2.6272, 2.6304, 2.6335};
            double Re_Amp[200] = {-0.1100, -0.1143, -0.1179, -0.1208, -0.1230, -0.1246, -0.1254, -0.1256, -0.1252, -0.1241, -0.1225, -0.1202, -0.1174, -0.1140, -0.1101, -0.1057, -0.1007, -0.0953, -0.0893, -0.0830, -0.0761, -0.0689, -0.0612, -0.0531, -0.0447, -0.0358, -0.0267, -0.0172, -0.0073, 0.0028, 0.0132, 0.0239, 0.0349, 0.0461, 0.0575, 0.0692, 0.0810, 0.0930, 0.1052, 0.1175, 0.1300, 0.1425, 0.1552, 0.1680, 0.1808, 0.1937, 0.2066, 0.2195, 0.2325, 0.2454, 0.2583, 0.2712, 0.2840, 0.2967, 0.3094, 0.3219, 0.3343, 0.3466, 0.3588, 0.3707, 0.3825, 0.3941, 0.4055, 0.4167, 0.4276, 0.4383, 0.4488, 0.4591, 0.4690, 0.4788, 0.4883, 0.4975, 0.5064, 0.5151, 0.5235, 0.5316, 0.5394, 0.5469, 0.5540, 0.5609, 0.5675, 0.5737, 0.5796, 0.5851, 0.5904, 0.5952, 0.5997, 0.6039, 0.6076, 0.6110, 0.6141, 0.6167, 0.6189, 0.6208, 0.6222, 0.6232, 0.6239, 0.6242, 0.6240, 0.6236, 0.6227, 0.6215, 0.6199, 0.6179, 0.6157, 0.6130, 0.6100, 0.6067, 0.6031, 0.5991, 0.5948, 0.5902, 0.5853, 0.5801, 0.5746, 0.5688, 0.5627, 0.5563, 0.5496, 0.5427, 0.5355, 0.5280, 0.5202, 0.5123, 0.5040, 0.4955, 0.4868, 0.4779, 0.4688, 0.4595, 0.4501, 0.4406, 0.4310, 0.4213, 0.4115, 0.4018, 0.3920, 0.3823, 0.3726, 0.3630, 0.3534, 0.3440, 0.3347, 0.3256, 0.3166, 0.3079, 0.2994, 0.2911, 0.2831, 0.2754, 0.2681, 0.2610, 0.2544, 0.2481, 0.2422, 0.2368, 0.2318, 0.2273, 0.2233, 0.2197, 0.2166, 0.2139, 0.2115, 0.2095, 0.2079, 0.2065, 0.2055, 0.2047, 0.2042, 0.2039, 0.2037, 0.2038, 0.2040, 0.2044, 0.2049, 0.2054, 0.2060, 0.2067, 0.2074, 0.2081, 0.2087, 0.2093, 0.2098, 0.2103, 0.2106, 0.2107, 0.2107, 0.2106, 0.2102, 0.2095, 0.2087, 0.2076, 0.2063, 0.2049, 0.2032, 0.2014, 0.1994, 0.1973, 0.1951, 0.1928};
            double Im_Amp[200] = {-0.0400, -0.0680, -0.0954, -0.1220, -0.1480, -0.1732, -0.1978, -0.2216, -0.2448, -0.2673, -0.2892, -0.3104, -0.3309, -0.3507, -0.3700, -0.3885, -0.4064, -0.4237, -0.4403, -0.4563, -0.4717, -0.4865, -0.5006, -0.5141, -0.5271, -0.5394, -0.5511, -0.5622, -0.5727, -0.5827, -0.5920, -0.6008, -0.6090, -0.6167, -0.6237, -0.6302, -0.6362, -0.6416, -0.6465, -0.6508, -0.6546, -0.6578, -0.6605, -0.6627, -0.6644, -0.6656, -0.6662, -0.6664, -0.6660, -0.6651, -0.6638, -0.6619, -0.6596, -0.6568, -0.6535, -0.6498, -0.6456, -0.6409, -0.6357, -0.6302, -0.6241, -0.6176, -0.6107, -0.6034, -0.5956, -0.5875, -0.5790, -0.5702, -0.5610, -0.5515, -0.5418, -0.5317, -0.5214, -0.5108, -0.5000, -0.4890, -0.4778, -0.4665, -0.4549, -0.4433, -0.4315, -0.4196, -0.4076, -0.3955, -0.3834, -0.3713, -0.3591, -0.3469, -0.3348, -0.3226, -0.3106, -0.2985, -0.2866, -0.2748, -0.2631, -0.2515, -0.2400, -0.2287, -0.2175, -0.2065, -0.1957, -0.1850, -0.1744, -0.1641, -0.1539, -0.1440, -0.1342, -0.1247, -0.1154, -0.1063, -0.0974, -0.0888, -0.0804, -0.0723, -0.0645, -0.0569, -0.0496, -0.0426, -0.0359, -0.0295, -0.0233, -0.0176, -0.0121, -0.0070, -0.0022, 0.0023, 0.0064, 0.0102, 0.0136, 0.0167, 0.0195, 0.0219, 0.0241, 0.0260, 0.0275, 0.0288, 0.0298, 0.0306, 0.0310, 0.0312, 0.0312, 0.0309, 0.0303, 0.0296, 0.0286, 0.0273, 0.0259, 0.0243, 0.0224, 0.0204, 0.0182, 0.0158, 0.0132, 0.0105, 0.0076, 0.0045, 0.0013, -0.0020, -0.0055, -0.0091, -0.0128, -0.0166, -0.0204, -0.0243, -0.0283, -0.0323, -0.0363, -0.0403, -0.0443, -0.0482, -0.0522, -0.0560, -0.0598, -0.0636, -0.0672, -0.0707, -0.0741, -0.0773, -0.0804, -0.0833, -0.0861, -0.0886, -0.0910, -0.0931, -0.0950, -0.0966, -0.0980, -0.0990, -0.0998, -0.1003, -0.1005, -0.1004, -0.1000, -0.0994, -0.0985, -0.0974, -0.0961, -0.0947, -0.0931, -0.0913};
            int idx_sqrt_s = int(( sqrt_s - 2.01 )/((5.14-2.01)/1000));
            if(idx_sqrt_s<0){idx_sqrt_s=0;}
            if(idx_sqrt_s>=200){idx_sqrt_s=199;}
            return DeviceComplex(Re_Amp[idx_sqrt_s],Im_Amp[idx_sqrt_s]);
        }

        //Including the triangle-plot into D1D amplitude, according to 2201.08253v2, Eq. 26
        DeviceComplex fun6(double sqrt_s, Event* evt, int idx_isobar){
            double mass = pars[0];
            double width = pars[1];
            DeviceComplex BW_D1 = 1.0/DeviceComplex(sqrt_s*sqrt_s-mass*mass,mass*width);
            //Some input paras
            double md1 = 2.420;
            double gd1 = 0.031;
            double mdstr = 2.010;
            double md = 1.869;
            double mpi = 0.139;
            double mjpsi = 3.097;
            double hc = 0.197327;
            double CZ = -0.177*pow(1/hc,2);
            double C12 = 0.005*pow(1/hc,2);
            double b = 0.0;
            double mu = 1.0;
            double pi = 3.1415926;
            double a2 = -3.0;

            //First cal_I_nom
            double s = evt->_mass2_mom0;
            double M = sqrt(s);
            double m23 = sqrt(evt->sigma3_func());//for D* D, hence the order of particles should be D* D pi!
            double eps = 1E-8;

            double m1 = md1;
            double m2 = md;
            double m3 = mdstr;
            double mk = mpi;

            double mu_12 = m1 * m2 / (m1 + m2);
            double mu_23 = m2 * m3 / (m2 + m3);

            double b12 = m1 + m2 - M;
            double c1 = 2 * mu_12 * b12;

            double ek = (s - m23*m23 + mk*mk)/(2*M);
            //double qk = sqrt( (s - (mk + m23)*(mk + m23)) * (s - (mk - m23)*(mk - m23)) )/(2*M);
            double qk = break_mom_abs(M,mk,m23);
            double c2 = 2*mu_23*(m2 + m3 + ek - M) + qk*qk * mu_23 / m3;
            double a = pow(mu_23*qk/m3,2);

            double fac = mu_12*m3/(2*pi*qk) * 1.0/(8*md1*md*mdstr);
            DeviceComplex arc1 = ( 0.5*(c2 - c1)/ (a*(c1 - DeviceComplex(0,1)*eps)).sqrt() ).atan();
            DeviceComplex arc2 = ( 0.5*(c2 - c1 - 2*a)/(a*(c2 - a - DeviceComplex(0,1)*eps)).sqrt() ).atan();
            DeviceComplex cal_I_nom = fac*(arc1 - arc2);

            //Then cal_G_nom
            double w = m23;
            double sprime = w*w;
            double delta = md*md- mdstr*mdstr;
            double qcm = break_mom_abs(w,md,mdstr);
            double prefix = (1.0/16.0/pi/pi);
            double term1 = a2+2*log(md/mu)+2*(mdstr*mdstr-md*md+sprime)/2.0/sprime*log(mdstr/md);
            DeviceComplex term2 = qcm/w*(
                (DeviceComplex(sprime-delta+2*qcm*w,0)).ln()
                +(DeviceComplex(sprime+delta+2*qcm*w,0)).ln()
                -(DeviceComplex(-sprime+delta+2*qcm*w,0)).ln()
                -(DeviceComplex(-sprime-delta+2*qcm*w,0)).ln()
            );
            DeviceComplex cal_G_nom = prefix*(term1+term2);

            //Then cal_G_jpsipi_nom
            double a2p = -2.77;
            delta = mjpsi*mjpsi - mpi*mpi;
            qcm = break_mom_abs(w,mjpsi,mpi);
            term1 = a2p+2*log(mjpsi/mu)+2*(mpi*mpi-mjpsi*mjpsi+sprime)/2.0/sprime*log(mpi/mjpsi);
            term2 = qcm/w*(
                (DeviceComplex(sprime-delta+2*qcm*w,0)).ln()
                +(DeviceComplex(sprime+delta+2*qcm*w,0)).ln()
                -(DeviceComplex(-sprime+delta+2*qcm*w,0)).ln()
                -(DeviceComplex(-sprime-delta+2*qcm*w,0)).ln()
            );
            DeviceComplex cal_G_jpsipi_nom = prefix*(term1+term2);

            //Then cal_T
            double V12 = C12*4*sqrt(mjpsi*mpi*md*mdstr);
            double V22_slaph = (CZ+b/2.0/(md+mdstr)*(s-pow(md+mdstr,2)));
            double V22 = V22_slaph * 4*sqrt(md*mdstr*md*mdstr);
            DeviceComplex g1 = cal_G_jpsipi_nom;
            DeviceComplex g2 = cal_G_nom;
            DeviceComplex den = 1 - g2*V22 - g1*g2*V12*V12;
            DeviceComplex num22 = V22 + g1*V12*V12;
            DeviceComplex cal_T = num22/den;

            DeviceComplex T22 = cal_T;
            DeviceComplex I = cal_I_nom;

            return BW_D1 + T22 * I;    
        }
};

#endif //DYNAMIC_H
