//The Amplitude of Dalitz Plot Decomposition (DPD)
#include "../include/DPD.h"
#include "../include/CGcoeff_HASH.h"

//The weigner small D function
//Cited from https://en.wikipedia.org/wiki/Wigner_D-matrix
//In the form of $d^{j}_{m^\prime,m}(\beta)$
__device__ double DPD::factorial(int x){

    if(x<0){
        //printf("X is %d in factorial!\n",x);
        return 0.0;
    }
    double result = 1;
    for(int i=1;i<=x;i++){
        result = result*i;
    }
    return result;
}

__device__ double DPD::Wigner_smallD_Jle4(int J, int _M, int _Mp, double beta){

    if( abs(_M)> J || abs(_Mp)> J ){return 0.0;}

    //symmetry
    int M = _M;
    int Mp = _Mp;
    int factor = 1;

    //d^j_{mp,m} = d^j_{-m,-mp}
    if((M+Mp)<0){
        int temp = M;
        M = -Mp;
        Mp = -temp;
    }
    //d^j_{mp,m} = pow(-1,m-mp)*d^j_{m,mp}
    if(M<Mp){
        factor = factor * powf(-1,(M-Mp)/2);
        int temp = M;
        M = Mp;
        Mp = temp;
    }

    
    if(J==0){return factor*1.0;}
    if(J==1){
        if(M==+1&&Mp==+1){return factor*cos(0.5*beta);}
        if(M==+1&&Mp==-1){return factor*-sin(0.5*beta);}
    }
    if(J==2){
        if(M==+2&&Mp==+2){return factor*(1+cos(beta))/2;}
        if(M==+2&&Mp==+0){return factor*(-sin(beta)/sqrt(2.));}
        if(M==+2&&Mp==-2){return factor*(1-cos(beta))/2;}

        if(M==+0&&Mp==+0){return factor*cos(beta);}
    }
    if(J==3){
        if(M==+3&&Mp==+3){return factor*power_double_int(cos(0.5*beta),3);}
        if(M==+3&&Mp==+1){return factor*-sqrt(3.)*sin(0.5*beta)*power_double_int(cos(0.5*beta),2);}
        if(M==+3&&Mp==-1){return factor*+sqrt(3.)*cos(0.5*beta)*power_double_int(sin(0.5*beta),2);}
        if(M==+3&&Mp==-3){return factor*-power_double_int(sin(0.5*beta),3);}

        if(M==+1&&Mp==+1){return factor*cos(0.5*beta)*(3*power_double_int(cos(0.5*beta),2)-2);}
        if(M==+1&&Mp==-1){return factor*sin(0.5*beta)*(3*power_double_int(sin(0.5*beta),2)-2);}
    }
    if(J==4){
        if(M==+4&&Mp==+4){return factor*power_double_int(1+cos(beta),2)/4;}
        if(M==+4&&Mp==+2){return factor*-sin(beta)*(1+cos(beta))/2;}
        if(M==+4&&Mp==+0){return factor*sqrt(3.0/2.0)/2*power_double_int(sin(beta),2);}
        if(M==+4&&Mp==-2){return factor*-sin(beta)*(1-cos(beta))/2;}
        if(M==+4&&Mp==-4){return factor*power_double_int(1-cos(beta),2)/4;}

        if(M==+2&&Mp==+2){return factor*(2*cos(beta)*cos(beta)+cos(beta)-1)/2;}
        if(M==+2&&Mp==+0){return factor*-sqrt(3.0/2.0)*sin(beta)*cos(beta);}
        if(M==+2&&Mp==-2){return factor*-(2*cos(beta)*cos(beta)-cos(beta)-1)/2;}

        if(M==+0&&Mp==+0){return factor*(3*cos(beta)*cos(beta)-1)/2;}
    }
    if(J==5){
        if(M==+5&&Mp==+5){return factor*power_double_int(cos(0.5*beta),5);}
        if(M==+5&&Mp==+3){return factor*-sqrt(5.)*sin(0.5*beta)*power_double_int(cos(0.5*beta),4);}
        if(M==+5&&Mp==+1){return factor*+sqrt(10.)*power_double_int(sin(0.5*beta),2)*power_double_int(cos(0.5*beta),3);}
        if(M==+5&&Mp==-1){return factor*-sqrt(10.)*power_double_int(sin(0.5*beta),3)*power_double_int(cos(0.5*beta),2);}
        if(M==+5&&Mp==-3){return factor*+sqrt(5.)*power_double_int(sin(0.5*beta),4)*cos(beta);}
        if(M==+5&&Mp==-5){return factor*-power_double_int(sin(0.5*beta),5);}

        if(M==+3&&Mp==+3){return factor*power_double_int(cos(0.5*beta),3)*(1-5*power_double_int(sin(0.5*beta),2));}
        if(M==+3&&Mp==+1){return factor*-sqrt(2.)*sin(0.5*beta)*power_double_int(cos(0.5*beta),2)*(2-5*power_double_int(sin(0.5*beta),2));}
        if(M==+3&&Mp==-1){return factor*-sqrt(2.)*cos(0.5*beta)*power_double_int(sin(0.5*beta),2)*(2-5*power_double_int(cos(0.5*beta),2));}
        if(M==+3&&Mp==-3){return factor*power_double_int(sin(0.5*beta),3)*(1-5*power_double_int(cos(0.5*beta),2));}

        if(M==+1&&Mp==+1){return factor*cos(0.5*beta)*(3-12*power_double_int(cos(0.5*beta),2)+10*power_double_int(cos(0.5*beta),4));}
        if(M==+1&&Mp==-1){return factor*-sin(0.5*beta)*(3-12*power_double_int(sin(0.5*beta),2)+10*power_double_int(sin(0.5*beta),4));}
    }
    if(J==6){
        if(M==+6&&Mp==+6){return factor*power_double_int(1+cos(beta),3)/8;}
        if(M==+6&&Mp==+4){return factor*-sqrt(6.)*power_double_int(1+cos(beta),2)*sin(beta)/8;}
        if(M==+6&&Mp==+2){return factor*+sqrt(15.)*power_double_int(sin(beta),2)*(1+cos(beta))/8;}
        if(M==+6&&Mp==+0){return factor*-sqrt(5.)*power_double_int(sin(beta),3)/4;}
        if(M==+6&&Mp==-2){return factor*+sqrt(15.)*power_double_int(sin(beta),2)*(1-cos(beta))/8;}
        if(M==+6&&Mp==-4){return factor*-sqrt(6.)*power_double_int(1-cos(beta),2)*sin(beta)/8;}
        if(M==+6&&Mp==-6){return factor*power_double_int(1-cos(beta),3)/8;}

        if(M==+4&&Mp==+4){return factor*-power_double_int(1+cos(beta),2)*(2-3*cos(beta))/4;}
        if(M==+4&&Mp==+2){return factor*+sqrt(10.)*sin(beta)*(1-2*cos(beta)-3*power_double_int(cos(beta),2))/8;}
        if(M==+4&&Mp==+0){return factor*sqrt(30.)*power_double_int(sin(beta),2)*cos(beta)/4;}
        if(M==+4&&Mp==-2){return factor*-sqrt(10.)*sin(beta)*(1+2*cos(beta)-3*power_double_int(cos(beta),2))/8;}
        if(M==+4&&Mp==-4){return factor*+power_double_int(1-cos(beta),2)*(2+3*cos(beta))/4;}

        if(M==+2&&Mp==+2){return factor*-(1+cos(beta))*(1+10*cos(beta)-15*power_double_int(cos(beta),2))/8;}
        if(M==+2&&Mp==+0){return factor*sqrt(3.)*sin(beta)*(1-5*power_double_int(cos(beta),2))/4;}
        if(M==+2&&Mp==-2){return factor*-(1-cos(beta))*(1-10*cos(beta)-15*power_double_int(cos(beta),2))/8;}

        if(M==+0&&Mp==+0){return factor*-cos(beta)*(3-5*power_double_int(cos(beta),2))/2;}
    }
    if(J==7){
        if(M==+7&&Mp==+7){return factor*power_double_int(cos(0.5*beta),7);}
        if(M==+7&&Mp==+5){return factor*-sqrt(7.)*power_double_int(cos(0.5*beta),6)*sin(0.5*beta);}
        if(M==+7&&Mp==+3){return factor*+sqrt(21.)*power_double_int(cos(0.5*beta),5)*power_double_int(sin(0.5*beta),2);}
        if(M==+7&&Mp==+1){return factor*-sqrt(35.)*power_double_int(cos(0.5*beta),4)*power_double_int(sin(0.5*beta),3);}
        if(M==+7&&Mp==-1){return factor*+sqrt(35.)*power_double_int(cos(0.5*beta),3)*power_double_int(sin(0.5*beta),4);}
        if(M==+7&&Mp==-3){return factor*-sqrt(21.)*power_double_int(cos(0.5*beta),2)*power_double_int(sin(0.5*beta),5);}
        if(M==+7&&Mp==-5){return factor*+sqrt(7.)*power_double_int(cos(0.5*beta),1)*power_double_int(sin(0.5*beta),6);}
        if(M==+7&&Mp==-7){return factor*-power_double_int(sin(0.5*beta),7);}

        if(M==+5&&Mp==+5){return factor*power_double_int(cos(0.5*beta),5)*(1-7*power_double_int(sin(0.5*beta),2));}
        if(M==+5&&Mp==+3){return factor*-sqrt(3.)*power_double_int(cos(0.5*beta),4)*sin(0.5*beta)*(2-7*power_double_int(sin(0.5*beta),2));}
        if(M==+5&&Mp==+1){return factor*+sqrt(5.)*power_double_int(cos(0.5*beta),3)*power_double_int(sin(0.5*beta),2)*(3-7*power_double_int(sin(0.5*beta),2));}
        if(M==+5&&Mp==-1){return factor*+sqrt(5.)*power_double_int(cos(0.5*beta),2)*power_double_int(sin(0.5*beta),3)*(3-7*power_double_int(cos(0.5*beta),2));}
        if(M==+5&&Mp==-3){return factor*-sqrt(3.)*power_double_int(cos(0.5*beta),1)*power_double_int(sin(0.5*beta),4)*(2-7*power_double_int(cos(0.5*beta),2));}
        if(M==+5&&Mp==-5){return factor*power_double_int(sin(0.5*beta),5)*(1-7*power_double_int(cos(0.5*beta),2));}

        if(M==+3&&Mp==+3){return factor*power_double_int(cos(0.5*beta),3)*(10-30*power_double_int(cos(0.5*beta),2)+21*power_double_int(cos(0.5*beta),4));}
        if(M==+3&&Mp==+1){return factor*-sqrt(15.)*power_double_int(cos(0.5*beta),2)*sin(0.5*beta)*(2-8*power_double_int(cos(0.5*beta),2)+7*power_double_int(cos(0.5*beta),4));}
        if(M==+3&&Mp==-1){return factor*+sqrt(15.)*power_double_int(sin(0.5*beta),2)*cos(0.5*beta)*(2-8*power_double_int(sin(0.5*beta),2)+7*power_double_int(sin(0.5*beta),4));}
        if(M==+3&&Mp==-3){return factor*-power_double_int(sin(0.5*beta),3)*(10-30*power_double_int(sin(0.5*beta),2)+21*power_double_int(sin(0.5*beta),4));}

        if(M==+1&&Mp==+1){return factor*-cos(0.5*beta)*(4-30*power_double_int(cos(0.5*beta),2)+60*power_double_int(cos(0.5*beta),4)-35*power_double_int(cos(0.5*beta),6));}
        if(M==+1&&Mp==-1){return factor*-sin(0.5*beta)*(4-30*power_double_int(sin(0.5*beta),2)+60*power_double_int(sin(0.5*beta),4)-35*power_double_int(sin(0.5*beta),6));}
    }

    if(J==8){
        if(M==+8&&Mp==+8){return factor*power_double_int(1+cos(beta),4)/16;}
        if(M==+8&&Mp==+6){return factor*-sqrt(2.)*sin(beta)*power_double_int(1+cos(beta),3)/8;}
        if(M==+8&&Mp==+4){return factor*+sqrt(7.)*power_double_int(sin(beta),2)*power_double_int(1+cos(beta),2)/8;}
        if(M==+8&&Mp==+2){return factor*-sqrt(14.)*power_double_int(sin(beta),3)*power_double_int(1+cos(beta),1)/8;}
        if(M==+8&&Mp==+0){return factor*sqrt(70.)*power_double_int(sin(beta),4)/16;}
        if(M==+8&&Mp==-2){return factor*-sqrt(14.)*power_double_int(sin(beta),3)*power_double_int(1-cos(beta),1)/8;}
        if(M==+8&&Mp==-4){return factor*+sqrt(7.)*power_double_int(sin(beta),2)*power_double_int(1-cos(beta),2)/8;}
        if(M==+8&&Mp==-6){return factor*-sqrt(2.)*sin(beta)*power_double_int(1-cos(beta),3)/8;}
        if(M==+8&&Mp==-8){return factor*power_double_int(1-cos(beta),4)/16;}

        if(M==+6&&Mp==+6){return factor*-power_double_int(1+cos(beta),3)*(3-4*cos(beta))/8;}
        if(M==+6&&Mp==+4){return factor*sqrt(14.)*sin(beta)*power_double_int(1+cos(beta),2)*(1-2*cos(beta))/8;}
        if(M==+6&&Mp==+2){return factor*-sqrt(7.)*power_double_int(sin(beta),2)*(1+cos(beta))*(1-4*cos(beta))/8;}
        if(M==+6&&Mp==+0){return factor*-sqrt(35.)*power_double_int(sin(beta),3)*cos(beta)/4;}
        if(M==+6&&Mp==-2){return factor*+sqrt(7.)*power_double_int(sin(beta),2)*(1-cos(beta))*(1+4*cos(beta))/8;}
        if(M==+6&&Mp==-4){return factor*-sqrt(14.)*sin(beta)*power_double_int(1-cos(beta),2)*(1+2*cos(beta))/8;}
        if(M==+6&&Mp==-6){return factor*+power_double_int(1-cos(beta),3)*(3+4*cos(beta))/8;}

        if(M==+4&&Mp==+4){return factor*power_double_int(1+cos(beta),2)*(1-7*cos(beta)+7*power_double_int(cos(beta),2))/4;}
        if(M==+4&&Mp==+2){return factor*sqrt(2.)*sin(beta)*(1+cos(beta))*(1+7*cos(beta)-14*power_double_int(cos(beta),2))/8;}
        if(M==+4&&Mp==+0){return factor*-sqrt(10.)*power_double_int(sin(beta),2)*(1-7*power_double_int(cos(beta),2))/8;}
        if(M==+4&&Mp==-2){return factor*sqrt(2.)*sin(beta)*(1-cos(beta))*(1-7*cos(beta)-14*power_double_int(cos(beta),2))/8;}
        if(M==+4&&Mp==-4){return factor*power_double_int(1-cos(beta),2)*(1+7*cos(beta)+7*power_double_int(cos(beta),2))/4;}

        if(M==+2&&Mp==+2){return factor*(1+cos(beta))*(3-6*cos(beta)-21*power_double_int(cos(beta),2)+28*power_double_int(cos(beta),3))/8;}
        if(M==+2&&Mp==+0){return factor*sqrt(5.)*sin(beta)*cos(beta)*(3-7*power_double_int(cos(beta),2))/4;}
        if(M==+2&&Mp==-2){return factor*-(1-cos(beta))*(3+6*cos(beta)-21*power_double_int(cos(beta),2)-28*power_double_int(cos(beta),3))/8;}

        if(M==+0&&Mp==+0){return factor*(3-30*power_double_int(cos(beta),2)+35*power_double_int(cos(beta),4))/8;}

    }

    return 0.0;
}

__device__ double DPD::Wigner_smallD(double j, double mp, double m, double beta){
    if(j<=4){return Wigner_smallD_Jle4(2*j,2*mp,2*m,beta);}

    if( fabs(m)>j || fabs(mp)>j ){return 0.0;}

    double prefix = sqrt(factorial(j+mp)*factorial(j-mp)*factorial(j+m)*factorial(j-m));
    //double smin = min(0,int(m-mp));
    double smin(0);
    if((m-mp)>=0){smin = (m-mp);}
    double smax(j-mp);
    if((j+m)<=(j-mp)){smax = (j+m);}

    double sum = 0.0;
    for(double s=smin;s<=smax;s++){
        double numerator = powf(-1,mp-m+s)*power_double_int(cos(0.5*beta),2*j+m-mp-2*s)*power_double_int(sin(0.5*beta),mp-m+2*s);
        double denominator = factorial(j+m-s)*factorial(s)*factorial(mp-m+s)*factorial(j-mp-s);
        sum = sum + numerator/denominator;
    }

    return prefix*sum;
}

__device__ DeviceComplex DPD::Wigner_bigD(double j, double mp, double m, double alpha, double beta, double gamma){

    DeviceComplex part1(cos(-mp*alpha),sin(-mp*alpha));
    double part2 = Wigner_smallD(j,mp,m,beta);
    DeviceComplex part3(cos(-m*gamma),sin(-m*gamma));

    return part1*part2*part3;
}

__device__ double DPD::CG_coeff_le224(double j1, double m1, double j2, double m2, double j3, double m3){
    
    int _j1(2*j1),_m1(2*m1),_j2(2*j2),_m2(2*m2),_j3(2*j3),_m3(2*m3);
    if((_m1+_m2)!=_m3) return 0;
    if(_j1+_j2<_j3) return 0;
    if(_j1+_j3<_j2) return 0;
    if(_j2+_j3<_j1) return 0;
    if((_j1+_j2+_j3)%2!=0) return 0;

    int myfactor = 1;
    if(_j1<_j2){
        myfactor = myfactor*powf(-1,(_j1+_j2-_j3)/2);
        int temp = _j1; _j1 =_j2; _j2 = temp;
        temp = _m1; _m1 = _m2; _m2 = temp;
    }

    if(_m1<0){
        myfactor = myfactor*powf(-1,(_j1+_j2-_j3)/2);
        _m1 = -_m1;
        _m2 = -_m2;
    }

    int idx = (10000*_m2 + 1000*_m1 + 100*_j1 + 10*_j2 + _j3)%733+733;
    double test = myfactor*mycg[idx];
    return myfactor*mycg[idx];
}

//https://en.wikipedia.org/wiki/Table_of_Clebsch%E2%80%93Gordan_coefficients#_j1_=_1,_j2_=_1
__device__ double DPD::CG_coeff(double j1, double m1, double j2, double m2, double j3, double m3){

    if(j1<=2 && j2<=2){return CG_coeff_le224(j1,m1,j2,m2,j3,m3);}

    double prefix1 = (2*j3 + 1)*factorial(j3+j1-j2)*factorial(j3-j1+j2)*factorial(j1+j2-j3);
    prefix1 = prefix1/factorial(j1+j2+j3+1);
    prefix1 = sqrt(prefix1);

    double prefix2 = factorial(j3+m3)*factorial(j3-m3)*factorial(j1-m1)*factorial(j1+m1)*factorial(j2-m2)*factorial(j2+m2);
    prefix2 = sqrt(prefix2);

    double result = 0;
    for(int k=0;k<=(j1+j2-j3);k++){
        double f1 = j1+j2-j3-k;
        double f2 = j1-m1-k;
        double f3 = j2+m2-k;
        double f4 = j3-j2+m1+k;
        double f5 = j3-j1-m2+k;
        double dom = factorial(k)*factorial(f1)*factorial(f2)*factorial(f3)*factorial(f4)*factorial(f5);
        if(dom==0) continue;
        result = result + powf(-1,k) / dom;
    }
    
    return result*prefix1*prefix2;
}

__device__ double DPD::Blatt_Weisskopf_factor(double Q, int L){
    //Here Q0 is a hadron "scale" parameter Q0 =0.197321/R GeV/c, where R is the radius of the centrifugal barrier in fm.
    if(L==0){return 1.0;}
    if(L==1){return sqrt(2.0/(Q*Q+Q0*Q0));}
    if(L==2){return sqrt(13.0/(Q*Q*Q*Q+3.0*Q0*Q0*Q*Q+9.0*Q0*Q0*Q0*Q0));}
    if(L==3){return sqrt(277.0/(power_double_int(Q,6)+6.0*power_double_int(Q,4)*power_double_int(Q0,2)+45.0*power_double_int(Q,2)*power_double_int(Q0,4)+225.0*power_double_int(Q0,6)));}
    if(L==4){return sqrt(12746.0/(power_double_int(Q,8)+10.0*power_double_int(Q,6)*power_double_int(Q0,2)+135.0*power_double_int(Q,4)*power_double_int(Q0,4)+1575.0*power_double_int(Q,2)*power_double_int(Q0,6)+11025.0*power_double_int(Q0,8)));}
    if(L>4){printf("NO valid Blatt Weisskopf factor!\n");}
    return 0.0;
}

//LS couping approach
__device__ double DPD::Helicity_LScoupling(double p, int L){
    return power_double_int(p,L)*Blatt_Weisskopf_factor(p,L);
}

//Helicity coupling approach
//Ref: (J,tau-lambda_k)->(s,tau)+(jk,lambda_k)
//I rewrite it as (J,lam_1-lam2)->(j1,lam_1)+(j2,lam_2)
__device__ DeviceComplex DPD::Helicity_HCcoupling(double J, double j1, double lam1, double j2, double lam2, double mom, DecayChain* mychain, int type_ls){
    
    int nLS(0);
    if(type_ls==1){ nLS = mychain->N1_LS; }
    if(type_ls==2){ nLS = mychain->N2_LS; }
    DeviceComplex sum(0.);
    for(int i=0;i<nLS;i++){
        double S(0),L(0);
        if(type_ls==1){S = mychain->S1[i]; L = mychain->L1[i];}
        if(type_ls==2){S = mychain->S2[i]; L = mychain->L2[i];} 
        double cgcoff1 = CG_coeff(j1,lam1,j2,-lam2,S,(lam1-lam2));
        double cgcoff2 = CG_coeff(L,0,S,(lam1-lam2),J,(lam1-lam2));
        sum = sum + mychain->Get_LScoff(type_ls,i) * Helicity_LScoupling(mom,L)*sqrt((2*L+1.0)/(2*J+1.0))*cgcoff1*cgcoff2;
        //printf("%d %d %d %f %f\n",type_ls,S,L,cgcoff1,cgcoff2);
    }
    
    return sum;

}

//Dalitz Plot function of 0->ijK
//Decay chain of 0->(ij)+k,(ij)->i+j
//The spin of 0 is 'J', spin of (ij) is 's'
//nu is the helicity of '0'
//lam[3] is the helicity of i,j,k
//spin[3] is the spin of i,j,k
__device__ DeviceComplex DPD::Dalita_plot_function(double J, double s, double nu, double lam[3], double spin[3], Event* evt, int k, DecayChain* mychain){
    double prefix = sqrt(2*J+1)*sqrt(2*s+1);

    //All angles HERE are (theta)
    double theta1_k = evt->alignment_angle_func(k);
    double theta2_ij = evt->scatter_angle_func(k);
    double mom1_k = evt->break_mom1_func(k);
    double mom2_k = evt->break_mom2_func(k);
    double xi1_k = evt->wrotation_angle_func(1,k);
    double xi2_k = evt->wrotation_angle_func(2,k);
    double xi3_k = evt->wrotation_angle_func(3,k);

    double spin_i(-1),spin_j(-1),spin_k(-1);
    if(k==1){spin_i = spin[1]; spin_j = spin[2]; spin_k = spin[0];}
    if(k==2){spin_i = spin[2]; spin_j = spin[0]; spin_k = spin[1];}
    if(k==3){spin_i = spin[0]; spin_j = spin[1]; spin_k = spin[2];}
    double lamp_1(-1),lamp_2(-1),lamp_3(-1);

    DeviceComplex sum = 0.0;
    //Loop all tau, tau is the helicity of (ij) with spin s
    for(double tau=-s;tau<=+s;tau++){
        //Loop all lamp[3]
        for(double lamp_i = -spin_i; lamp_i<=spin_i; lamp_i++){
            for(double lamp_j = -spin_j; lamp_j<=spin_j; lamp_j++){
                for(double lamp_k = -spin_k; lamp_k<=spin_k; lamp_k++){
                    
                    //First decay
                    double wigner_smallD1 = Wigner_smallD(J,nu,tau-lamp_k,theta1_k);
                    if(wigner_smallD1==0){continue;}
                    DeviceComplex Helicity_coup1 = pow(-1,spin_k-lamp_k)*Helicity_HCcoupling(J,s,tau,spin_k,lamp_k,mom1_k,mychain,1);
                    //The Breit-Wigner term is excluded
                    //Second decay
                    double wigner_smallD2 = Wigner_smallD(s,tau,lamp_i-lamp_j,theta2_ij);
                    if(wigner_smallD2==0){continue;}
                    DeviceComplex Helicity_coup2 = pow(-1,spin_j-lamp_j)*Helicity_HCcoupling(s,spin_i,lamp_i,spin_j,lamp_j,mom2_k,mychain,2);
                    //Wigner rotations
                    if(k==1){lamp_1 = lamp_k; lamp_2 = lamp_i; lamp_3 = lamp_j;}
                    if(k==2){lamp_1 = lamp_j; lamp_2 = lamp_k; lamp_3 = lamp_i;}
                    if(k==3){lamp_1 = lamp_i; lamp_2 = lamp_j; lamp_3 = lamp_k;}

                    double rotation0 = Wigner_smallD(spin[0],lamp_1,lam[0],xi1_k);
                    double rotation1 = Wigner_smallD(spin[1],lamp_2,lam[1],xi2_k);
                    double rotation2 = Wigner_smallD(spin[2],lamp_3,lam[2],xi3_k);

                    DeviceComplex rslt = wigner_smallD1*Helicity_coup1*wigner_smallD2*Helicity_coup2*rotation0*rotation1*rotation2;
                    sum = rslt + sum;
                }
            }
        }
    }

    return sum*prefix;
}
