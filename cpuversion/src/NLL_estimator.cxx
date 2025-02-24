//Load the BKG file
#include "../include/NLL_estimator.h"

void NLL_estimator::Load_file(int file_type, TString file_name, TString chain_name, TString p4_1_name, TString p4_2_name, TString p4_3_name){

    TChain *chain = new TChain(chain_name);
    chain->Add(file_name);

    int Ntot = chain->GetEntries();
    if(file_type==0){ evt_dt = Ntot; array_evt_dt = new Event[Ntot];}
    if(file_type==1){ evt_bg = Ntot; array_evt_bg = new Event[Ntot];}
    if(file_type==2){ evt_mc = Ntot; array_evt_mc = new Event[Ntot];}
    if(file_type==3){ evt_mcT = Ntot; array_evt_mcT = new Event[Ntot];}

    double p4_1[4],p4_2[4],p4_3[4];
    TLorentzVector T4_1,T4_2,T4_3;
    chain->SetBranchAddress(p4_1_name,&p4_1);
    chain->SetBranchAddress(p4_2_name,&p4_2);
    chain->SetBranchAddress(p4_3_name,&p4_3);

    double p4_sec1[4], p4_sec2[4];
    if(with_sec==true){
        chain->SetBranchAddress(p4_dau1_name_sec,&p4_sec1);
        chain->SetBranchAddress(p4_dau2_name_sec,&p4_sec2);
    }

    int counter = 0;

    for(int j=0;j<Ntot;j++ ){
        chain->GetEntry(j);

        T4_1.SetPxPyPzE(p4_1[0],p4_1[1],p4_1[2],p4_1[3]);
        T4_2.SetPxPyPzE(p4_2[0],p4_2[1],p4_2[2],p4_2[3]);
        T4_3.SetPxPyPzE(p4_3[0],p4_3[1],p4_3[2],p4_3[3]);

        double sigma1 = pow((T4_2+T4_3).M(),2.0);
        double sigma2 = pow((T4_1+T4_3).M(),2.0);

        double mass2_mom0 = pow((T4_1+T4_2+T4_3).M(),2.0);
        double mass2_dau1 = pow((T4_1).M(),2.0);
        double mass2_dau2 = pow((T4_2).M(),2.0);
        double mass2_dau3 = pow((T4_3).M(),2.0);

        TVector3 boostVector = -(T4_1+T4_2+T4_3).BoostVector();
        T4_1.Boost(boostVector);
        T4_2.Boost(boostVector);
        T4_3.Boost(boostVector);

        double alpha = (-T4_1).Phi();
        double beta = (-T4_1).Theta();
        TVector3 plane23 = T4_2.Vect().Cross(T4_3.Vect());
        //TVector3 Y_axis(0,1,0);
        //double gamma = Y_axis.Angle(plane23);
        TVector3 Z_axis(0,0,1);
        TVector3 plane1Z = T4_1.Vect().Cross(Z_axis);
        double gamma = (plane1Z).Angle(plane23);

        //Create Event
        Event evt;
        evt._alpha = alpha; evt._beta = beta; evt._gamma = gamma;
        evt._sigma1 = sigma1; evt._sigma2 = sigma2;
        evt._mass2_mom0 = mass2_mom0; evt._mass2_dau1 = mass2_dau1; evt._mass2_dau2 = mass2_dau2; evt._mass2_dau3 = mass2_dau3; 

        //secondary decay
        if(with_sec==true){
            double theta_sec(0),phi_sec(0);
            int idx_sec = amp_obj->idx_sec;

            TLorentzVector T4_123 = T4_1 + T4_2 + T4_3;
            TLorentzVector T4_sec_1, T4_sec_2;
            T4_sec_1.SetPxPyPzE(p4_sec1[0],p4_sec1[1],p4_sec1[2],p4_sec1[3]);
            T4_sec_2.SetPxPyPzE(p4_sec2[0],p4_sec2[1],p4_sec2[2],p4_sec2[3]);

            TVector3 decay_plane_sec = (T4_sec_2.Vect()).Cross(T4_sec_1.Vect());
            TVector3 decay_plane_mom = (T4_1.Vect()).Cross(T4_2.Vect());
            phi_sec = ((decay_plane_sec.Cross(decay_plane_mom)).Dot(T4_1.Vect())>0? 1.0 : -1.0) * decay_plane_sec.Angle(decay_plane_mom);

            TVector3 boostVector(0,0,0);
            if(idx_sec==1) boostVector = -T4_1.BoostVector();
            if(idx_sec==2) boostVector = -T4_2.BoostVector();
            if(idx_sec==3) boostVector = -T4_3.BoostVector();

            T4_123.Boost(boostVector);
            T4_sec_1.Boost(boostVector);
            theta_sec = (-T4_123.Vect()).Angle(T4_sec_1.Vect());
            evt._second_theta = theta_sec;
            evt._second_phi = phi_sec;
        }

        if(file_type==0){ array_evt_dt[counter] = evt;}
        if(file_type==1){ array_evt_bg[counter] = evt;}
        if(file_type==2){ array_evt_mc[counter] = evt;}
        if(file_type==3){ array_evt_mcT[counter] = evt;}

        counter++;
    }

    if(file_type==0){ std::cout<<"DATA sample loaded! Number:"<<counter<<endl;}
    if(file_type==1){ std::cout<<"BKG sample loaded! Number: "<<counter<<endl;}
    if(file_type==2){ std::cout<<"MC sample loaded! Number:  "<<counter<<endl;}
    if(file_type==3){ std::cout<<"MC Truth sample loaded! Number:  "<<counter<<endl;}

    delete chain;
}

LSCoeff_par* NLL_estimator::Search_LSCoeff_par(int type_ls, int m_idx_chain, int m_idx_LS, LSCoeff_par* m_LSCoeffpar_list){
    int n_Cpar(0);
    if(type_ls==1) n_Cpar = amp_obj->Get_total_LS1coeff_par();
    if(type_ls==2) n_Cpar = amp_obj->Get_total_LS2coeff_par();
    for(int idx_Cpar=0;idx_Cpar<n_Cpar;idx_Cpar++){
        int idx_chain = m_LSCoeffpar_list[idx_Cpar].idx_chain;
        int idx_LS = m_LSCoeffpar_list[idx_Cpar].idx_LS;
        if(m_idx_chain==idx_chain && m_idx_LS == idx_LS){return &m_LSCoeffpar_list[idx_Cpar];}
    }
    LSCoeff_par* myvoid(nullptr);
    return myvoid;
}

Reson_par* NLL_estimator::Search_Reson_par(int m_idx_chain, int m_idx_dynamic_par, Reson_par* m_Respar_list){
    int n_Rpar = amp_obj->Get_total_dynamic_par();
    for(int idx_Rpar=0;idx_Rpar<n_Rpar;idx_Rpar++){
        int idx_chain = m_Respar_list[idx_Rpar].idx_chain;
        int idx_dynamic_par = m_Respar_list[idx_Rpar].idx_dynamic_par;
        if(m_idx_chain==idx_chain && m_idx_dynamic_par == idx_dynamic_par){return &m_Respar_list[idx_Rpar];}
    }
    Reson_par* myvoid(nullptr);
    return myvoid;
}



Para* NLL_estimator::Search_idx_minuit(int m_idx_minuit, bool& isfind){
    for(int idx_par=0;idx_par<N_totpar;idx_par++){
        //Resonance parameter
        for(int idx_rpar=0;idx_rpar<N_Respar;idx_rpar++){
            if(Respar_list[idx_rpar].par.idx_minuit==m_idx_minuit){
                isfind = true;
                return &(Respar_list[idx_rpar].par);
            }
        }
        //Coupling coeff of LS1
        for(int idx_cpar=0;idx_cpar<N1_LSCoeff;idx_cpar++){
            if(LS1Coeffpar_list[idx_cpar].par_rho.idx_minuit==m_idx_minuit){
                isfind = true;
                return &(LS1Coeffpar_list[idx_cpar].par_rho);
            }
            if(LS1Coeffpar_list[idx_cpar].par_phi.idx_minuit==m_idx_minuit){
                isfind = true;
                return &(LS1Coeffpar_list[idx_cpar].par_phi);
            }
        }
        //Coupling coeff of LS2
        for(int idx_cpar=0;idx_cpar<N2_LSCoeff;idx_cpar++){
            if(LS2Coeffpar_list[idx_cpar].par_rho.idx_minuit==m_idx_minuit){
                isfind = true;
                return &(LS2Coeffpar_list[idx_cpar].par_rho);
            }
            if(LS2Coeffpar_list[idx_cpar].par_phi.idx_minuit==m_idx_minuit){
                isfind = true;
                return &(LS2Coeffpar_list[idx_cpar].par_phi);
            }
        }
    }
    Para* myvoid(nullptr);
    return myvoid;
}

Para* NLL_estimator::Search_par_name(string m_name, bool& isfind){
    for(int idx_par=0;idx_par<N_totpar;idx_par++){
        //Resonance parameter
        for(int idx_rpar=0;idx_rpar<N_Respar;idx_rpar++){
            if(Respar_list[idx_rpar].par.para_name==m_name){
                isfind = true;
                return &(Respar_list[idx_rpar].par);
            }
        }
        //Coupling coeff of LS1
        for(int idx_cpar=0;idx_cpar<N1_LSCoeff;idx_cpar++){
            if(LS1Coeffpar_list[idx_cpar].par_rho.para_name==m_name){
                isfind = true;
                return &(LS1Coeffpar_list[idx_cpar].par_rho);
            }
            if(LS1Coeffpar_list[idx_cpar].par_phi.para_name==m_name){
                isfind = true;
                return &(LS1Coeffpar_list[idx_cpar].par_phi);
            }
        }
        //Coupling coeff of LS2
        for(int idx_cpar=0;idx_cpar<N2_LSCoeff;idx_cpar++){
            if(LS2Coeffpar_list[idx_cpar].par_rho.para_name==m_name){
                isfind = true;
                return &(LS2Coeffpar_list[idx_cpar].par_rho);
            }
            if(LS2Coeffpar_list[idx_cpar].par_phi.para_name==m_name){
                isfind = true;
                return &(LS2Coeffpar_list[idx_cpar].par_phi);
            }
        }
    }
    Para* myvoid(nullptr);
    return myvoid;
}

void NLL_estimator::Update_Paras(){

    int n_Cpar1 = amp_obj->Get_total_LS1coeff_par();
    int n_Cpar2 = amp_obj->Get_total_LS2coeff_par();
    int n_Rpar = amp_obj->Get_total_dynamic_par();

    for(int idx_Cpar=0;idx_Cpar<n_Cpar1;idx_Cpar++){
        int idx_chain = LS1Coeffpar_list[idx_Cpar].idx_chain;
        int idx_LS = LS1Coeffpar_list[idx_Cpar].idx_LS;
        (amp_obj->array_chain[idx_chain]).rho1_coff[idx_LS] = LS1Coeffpar_list[idx_Cpar].par_rho.Get_Val();
        (amp_obj->array_chain[idx_chain]).phi1_coff[idx_LS] = LS1Coeffpar_list[idx_Cpar].par_phi.Get_Val();
    }

    for(int idx_Cpar=0;idx_Cpar<n_Cpar2;idx_Cpar++){
        int idx_chain = LS2Coeffpar_list[idx_Cpar].idx_chain;
        int idx_LS = LS2Coeffpar_list[idx_Cpar].idx_LS;
        (amp_obj->array_chain[idx_chain]).rho2_coff[idx_LS] = LS2Coeffpar_list[idx_Cpar].par_rho.Get_Val();
        (amp_obj->array_chain[idx_chain]).phi2_coff[idx_LS] = LS2Coeffpar_list[idx_Cpar].par_phi.Get_Val();
    }

    for(int idx_Rpar=0;idx_Rpar<n_Rpar;idx_Rpar++){
        int idx_chain = Respar_list[idx_Rpar].idx_chain;
        int idx_dynamic_par = Respar_list[idx_Rpar].idx_dynamic_par;
        (amp_obj->array_chain[idx_chain]).intermediate.dynamic.pars[idx_dynamic_par] = Respar_list[idx_Rpar].par.Get_Val();
    }
}

void NLL_estimator::Clear_LScoeff(){
    int n_Cpar1 = amp_obj->Get_total_LS1coeff_par();
    int n_Cpar2 = amp_obj->Get_total_LS2coeff_par();
    int n_Rpar = amp_obj->Get_total_dynamic_par();
    for(int i=0;i<n_Cpar1;i++){
        LS1Coeffpar_list[i].par_rho.Val = 0;
        LS1Coeffpar_list[i].par_phi.Val = 0;
        //Remove the fix
        LS1Coeffpar_list[i].par_rho.is_fixed_to = false;
        LS1Coeffpar_list[i].par_phi.is_fixed_to = false;
    }
    for(int i=0;i<n_Cpar2;i++){
        LS2Coeffpar_list[i].par_rho.Val = 0;
        LS2Coeffpar_list[i].par_phi.Val = 0;
        //Remove the fix
        LS2Coeffpar_list[i].par_rho.is_fixed_to = false;
        LS2Coeffpar_list[i].par_phi.is_fixed_to = false;
    }
    Update_Paras();
}

void NLL_estimator::CalPDF(Event* evt_arr, int num_evt, double* amp2){
    for(int idx_evt=0; idx_evt<num_evt;idx_evt++){
        amp2[idx_evt] = amp_obj->SumOverLam(&evt_arr[idx_evt]);
    }
}

void NLL_estimator::save_root(int file_type, TString file_name, int save_component){
    TFile* file = new TFile(file_name, "recreate");
    TTree* my_tree;
    if(file_type==2){my_tree = new TTree("MC","MC");}
    if(file_type==3){my_tree = new TTree("MCT","MCT");}

    double *PDF_MC_tot; double norm_factor(0.);

    if(file_type==2||file_type==3){
        if(file_type==2){PDF_MC_tot = new double[evt_mc]; CalPDF(array_evt_mc,evt_mc,PDF_MC_tot);}
        if(file_type==3){PDF_MC_tot = new double[evt_mcT]; CalPDF(array_evt_mcT,evt_mcT,PDF_MC_tot);}
        norm_factor = 1.0;
    }

    double ***PDF_MC_component;
    if((file_type==2||file_type==3) && save_component==1){
        int nchain = amp_obj->nchain;
        //create the PDF matrix
        PDF_MC_component = new double**[nchain];
        for(int idx_ch1=0;idx_ch1<nchain;idx_ch1++){
            PDF_MC_component[idx_ch1] = new double*[nchain];
            for(int idx_ch2=0;idx_ch2<nchain;idx_ch2++){
                if(file_type==2){PDF_MC_component[idx_ch1][idx_ch2] = new double[evt_mc];}
                if(file_type==3){PDF_MC_component[idx_ch1][idx_ch2] = new double[evt_mcT];}
            }
        }

        int n_Cpar1 = amp_obj->Get_total_LS1coeff_par();
        int n_Cpar2 = amp_obj->Get_total_LS2coeff_par();
        int n_Rpar = amp_obj->Get_total_dynamic_par();
        LSCoeff_par* LS1Coeffpar_list_backup = new LSCoeff_par[n_Cpar1];
        LSCoeff_par* LS2Coeffpar_list_backup = new LSCoeff_par[n_Cpar2];
        Reson_par* Respar_list_backup = new Reson_par[n_Rpar];
        for(int i=0;i<n_Cpar1;i++){LS1Coeffpar_list_backup[i] = LS1Coeffpar_list[i];}
        for(int i=0;i<n_Cpar2;i++){LS2Coeffpar_list_backup[i] = LS2Coeffpar_list[i];}
        for(int i=0;i<n_Rpar;i++){Respar_list_backup[i] = Respar_list[i];}

        for(int idx_ch1=0;idx_ch1<nchain;idx_ch1++){
            for(int idx_ch2=0;idx_ch2<nchain;idx_ch2++){
                Clear_LScoeff();

                //Decay chain 1
                int nLS1_ch1 = amp_obj->array_chain[idx_ch1].N1_LS;
                int nLS2_ch1 = amp_obj->array_chain[idx_ch1].N2_LS;
                for(int idx_ch1_ls=0;idx_ch1_ls<nLS1_ch1;idx_ch1_ls++){
                    Search_LSCoeff_par(1,idx_ch1,idx_ch1_ls,LS1Coeffpar_list)->par_rho.Val = Search_LSCoeff_par(1,idx_ch1,idx_ch1_ls,LS1Coeffpar_list_backup)->par_rho.Val;
                    Search_LSCoeff_par(1,idx_ch1,idx_ch1_ls,LS1Coeffpar_list)->par_phi.Val = Search_LSCoeff_par(1,idx_ch1,idx_ch1_ls,LS1Coeffpar_list_backup)->par_phi.Val;
                }
                for(int idx_ch1_ls=0;idx_ch1_ls<nLS2_ch1;idx_ch1_ls++){
                    Search_LSCoeff_par(2,idx_ch1,idx_ch1_ls,LS2Coeffpar_list)->par_rho.Val = Search_LSCoeff_par(2,idx_ch1,idx_ch1_ls,LS2Coeffpar_list_backup)->par_rho.Val;
                    Search_LSCoeff_par(2,idx_ch1,idx_ch1_ls,LS2Coeffpar_list)->par_phi.Val = Search_LSCoeff_par(2,idx_ch1,idx_ch1_ls,LS2Coeffpar_list_backup)->par_phi.Val;
                }

                //Decay chain 2
                int nLS1_ch2 = amp_obj->array_chain[idx_ch2].N1_LS;
                int nLS2_ch2 = amp_obj->array_chain[idx_ch2].N2_LS;
                for(int idx_ch2_ls=0;idx_ch2_ls<nLS1_ch2;idx_ch2_ls++){
                    Search_LSCoeff_par(1,idx_ch2,idx_ch2_ls,LS1Coeffpar_list)->par_rho.Val = Search_LSCoeff_par(1,idx_ch2,idx_ch2_ls,LS1Coeffpar_list_backup)->par_rho.Val;
                    Search_LSCoeff_par(1,idx_ch2,idx_ch2_ls,LS1Coeffpar_list)->par_phi.Val = Search_LSCoeff_par(1,idx_ch2,idx_ch2_ls,LS1Coeffpar_list_backup)->par_phi.Val;
                }
                for(int idx_ch2_ls=0;idx_ch2_ls<nLS2_ch2;idx_ch2_ls++){
                    Search_LSCoeff_par(2,idx_ch2,idx_ch2_ls,LS2Coeffpar_list)->par_rho.Val = Search_LSCoeff_par(2,idx_ch2,idx_ch2_ls,LS2Coeffpar_list_backup)->par_rho.Val;
                    Search_LSCoeff_par(2,idx_ch2,idx_ch2_ls,LS2Coeffpar_list)->par_phi.Val = Search_LSCoeff_par(2,idx_ch2,idx_ch2_ls,LS2Coeffpar_list_backup)->par_phi.Val;
                }
                
                Update_Paras();
                if(file_type==2){CalPDF(array_evt_mc,evt_mc,PDF_MC_component[idx_ch1][idx_ch2]);}
                if(file_type==3){CalPDF(array_evt_mcT,evt_mcT,PDF_MC_component[idx_ch1][idx_ch2]);}
           }
        }

        for(int i=0;i<n_Cpar1;i++){LS1Coeffpar_list[i] = LS1Coeffpar_list_backup[i];}
        for(int i=0;i<n_Cpar2;i++){LS2Coeffpar_list[i] = LS2Coeffpar_list_backup[i];}
        for(int i=0;i<n_Rpar;i++){Respar_list[i] = Respar_list_backup[i];}
        Update_Paras();
    }
    

    double m_alpha; double m_beta; double m_gamma;
    double m_sigma1; double m_sigma2; double m_sigma3;
    double m_scatter_angle1; double m_scatter_angle2; double m_scatter_angle3;
    double m_align_angle1; double m_align_angle2; double m_align_angle3;
    double m_second_theta; double m_second_phi;
    double m_weight_tot;
    const int max_nres = 20;
    double m_weight_component[max_nres][max_nres];

    my_tree->Branch("alpha",&m_alpha,"m_alpha/D");
    my_tree->Branch("beta",&m_beta,"m_beta/D");
    my_tree->Branch("gamma",&m_gamma,"m_gamma/D");
    my_tree->Branch("sigma1",&m_sigma1,"m_sigma1/D");
    my_tree->Branch("sigma2",&m_sigma2,"m_sigma2/D");
    my_tree->Branch("sigma3",&m_sigma3,"m_sigma3/D");
    my_tree->Branch("scatter_angle1",&m_scatter_angle1,"m_scatter_angle1/D");
    my_tree->Branch("scatter_angle2",&m_scatter_angle2,"m_scatter_angle2/D");
    my_tree->Branch("scatter_angle3",&m_scatter_angle3,"m_scatter_angle3/D");
    my_tree->Branch("align_angle1",&m_align_angle1,"m_align_angle1/D");
    my_tree->Branch("align_angle2",&m_align_angle2,"m_align_angle2/D");
    my_tree->Branch("align_angle3",&m_align_angle3,"m_align_angle3/D");
    my_tree->Branch("second_theta",&m_second_theta,"m_second_theta/D");
    my_tree->Branch("second_phi",&m_second_phi,"m_second_phi/D");
    if(file_type==2||file_type==3){
        my_tree->Branch("weight_tot",&m_weight_tot,"m_weight_tot/D");
        my_tree->Branch("weight_component", m_weight_component, Form("m_weight_component[%d][%d]/D",max_nres,max_nres));
    }

    int Ntot(0); Event* arr_evt;
    if(file_type==2){Ntot = evt_mc;arr_evt = array_evt_mc;}
    if(file_type==3){Ntot = evt_mcT;arr_evt = array_evt_mcT;}

    for (Int_t evt_loop=0; evt_loop<Ntot; evt_loop++)
    {
        m_alpha = arr_evt[evt_loop]._alpha;
        m_beta = arr_evt[evt_loop]._beta;
        m_gamma = arr_evt[evt_loop]._gamma;
        m_sigma1 = arr_evt[evt_loop]._sigma1;
        m_sigma2 = arr_evt[evt_loop]._sigma2;
        m_sigma3 = arr_evt[evt_loop].sigma3_func();
        m_scatter_angle1 = arr_evt[evt_loop].scatter_angle_func(1);
        m_scatter_angle2 = arr_evt[evt_loop].scatter_angle_func(2);
        m_scatter_angle3 = arr_evt[evt_loop].scatter_angle_func(3);
        m_align_angle1 = arr_evt[evt_loop].alignment_angle_func(1);
        m_align_angle2 = arr_evt[evt_loop].alignment_angle_func(2);
        m_align_angle3 = arr_evt[evt_loop].alignment_angle_func(3);
        m_second_phi = arr_evt[evt_loop]._second_phi;
        m_second_theta = arr_evt[evt_loop]._second_theta;
        if(file_type==2||file_type==3){
            m_weight_tot = PDF_MC_tot[evt_loop]*norm_factor;
            int nchain = amp_obj->nchain;
            if(save_component==1){
                for(int idx_ch1=0;idx_ch1<nchain;idx_ch1++){
                    for(int idx_ch2=0;idx_ch2<nchain;idx_ch2++){
                        m_weight_component[idx_ch1][idx_ch2] = PDF_MC_component[idx_ch1][idx_ch2][evt_loop]*norm_factor;
                    }
                }
            }
        }
        my_tree->Fill();
    }

    if((file_type==2||file_type==3) && save_component==1){
        delete PDF_MC_tot;
        int nchain = amp_obj->nchain;
        for(int idx_ch1=0;idx_ch1<nchain;idx_ch1++){
            for(int idx_ch2=0;idx_ch2<nchain;idx_ch2++){
                delete PDF_MC_component[idx_ch1][idx_ch2];
            }
            delete PDF_MC_component[idx_ch1];
        }
    }

    file->Write();
    file->Close();
}