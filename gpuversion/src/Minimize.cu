#include "../include/Minimize.h"
#include <Math/GSLRndmEngines.h>

void objective_function(Int_t &nrPar, Double_t *grad, Double_t &Result, Double_t *par, Int_t flag_type){
    Minimize* Minimize_obj = (Minimize*)gMinuit->GetObjectFit();

    int Npar_tot = Minimize_obj->npar_tot;
    double *arr_NLL = new double[Minimize_obj->num_NLL];
    
    Result = 0;

    for(int i=0;i<Npar_tot;i++){
        for(int idx_nll=0;idx_nll<Minimize_obj->num_NLL;idx_nll++){
            NLL_estimator* mynll = Minimize_obj->NLL_estimator_array[idx_nll];
            bool isfind = false;
            Para* mypar = mynll->Search_idx_minuit(i,isfind);
            if(isfind==true){mypar->Val = par[i];}
        }
    }

    for(int idx_nll=0;idx_nll<Minimize_obj->num_NLL;idx_nll++){
        Minimize_obj->NLL_estimator_array[idx_nll]->Update_Paras();
        arr_NLL[idx_nll] = Minimize_obj->NLL_estimator_array[idx_nll]->Cal_log_likelihood();
        Result = Result - arr_NLL[idx_nll];
    }


    cout<<Result<<" ";
    for(int idx_nll=0;idx_nll<Minimize_obj->num_NLL;idx_nll++){cout<<arr_NLL[idx_nll]<<" ";}
    cout<<endl;

    delete arr_NLL;

}

void Minimize::Initial_optimize(){

    Int_t ierflg = 0;
    double arglist[10];
    gMinuit = new TMinuit(npar_tot);
    gMinuit->SetObjectFit(this);
    gMinuit->SetFCN(objective_function);

    //3 maximum output, showing progress of minimizations.
    arglist[0] = 3.0;
    gMinuit->mnexcm("SET PRINT",arglist ,1,ierflg);

    //for chisquared fits UP=1, and for negative log likelihood, UP=0.5.
    arglist[0] = 0.5;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    //Set pars for gMinuit
    int par_count = 0;
    
    for(int idx_nll=0;idx_nll<num_NLL;idx_nll++){
        NLL_estimator* mynll = NLL_estimator_array[idx_nll];
        int n_par_tot = mynll->N_totpar;
        for(int idx_par=0;idx_par<n_par_tot;idx_par++){
            bool isfind = false;
            Para* mypar = mynll->Search_idx_minuit(par_count,isfind);
            if(isfind==true){
                TString rpar_name = mypar->para_name;
                double initial_val = mypar->Val;
                double step = mypar->Step;
                double lower = mypar->low_bound;
                double upper = mypar->up_bound;
                gMinuit->mnparm(par_count, rpar_name, initial_val, step, lower, upper, ierflg);
                par_count++;
            }
        }
    }
}

void Minimize::Do_SCAN(){
    //scanparlist=[parno , numpts]
    Int_t ierflg = 0;
    double arglist[10];
    arglist[0] = 0;//0 means all
    arglist[1] = 40;//the points
    gMinuit->mnexcm("SCAN", arglist,  2,   ierflg);
}

void Minimize::Do_MIGRAD(){
    //ARGLIS = [self.cal, self.tol]
    Int_t ierflg = 0;
    double arglist[10];
    arglist[0] = 300000;
    arglist[1] = 0.1;
    gMinuit->mnexcm("MIGRAD", arglist , 2 ,ierflg);
}

void Minimize::Do_HESSE(){
    //ARGLIS = [self.cal, self.tol]
    Int_t ierflg = 0;
    double arglist[10];
    arglist[0] = 300000;
    gMinuit->mnexcm("HESSE", arglist , 1 ,ierflg);
}


void Minimize::UpdateMypar(){
    
    par_fitted = new double[npar_tot];
    par_error = new double[npar_tot];
    for(int par_loop=0;par_loop<npar_tot;par_loop++){
        gMinuit->GetParameter(par_loop,par_fitted[par_loop],par_error[par_loop]);
    }

    for(int i=0;i<npar_tot;i++){
        //if( par_error[i] == 0 ) continue;
        for(int idx_nll=0;idx_nll<num_NLL;idx_nll++){
            NLL_estimator* mynll = NLL_estimator_array[idx_nll];
            bool isfind = false;
            Para* mypar = mynll->Search_idx_minuit(i,isfind);
            if(isfind==true){
                mypar->Val = par_fitted[i];
                mypar->Step = par_error[i];
                //Remove the fix
                if(mypar->is_fixed_to==true){
                    mypar->Val = mypar->Get_Val();
                }
            }
        }
    }

    for(int idx_nll=0;idx_nll<num_NLL;idx_nll++){NLL_estimator_array[idx_nll]->Update_Paras();}
}

Json::Value Minimize::WritePar(Para* mypar){
    Json::Value chain_json;
    string parname = mypar->para_name;
    if(mypar->is_fixed_to==false){
        double parval = mypar->Val; chain_json.append(parval);
        double parerr = mypar->Step; chain_json.append(parerr);
        double parlow = mypar->low_bound; chain_json.append(parlow);
        double parup = mypar->up_bound; chain_json.append(parup);
    }
    else{
        string fixname = mypar->fixed_to_this->para_name;  chain_json.append(fixname);
        double fixratio = mypar->fix_ratio;  chain_json.append(fixratio);
        double parval = mypar->Val;  chain_json.append(parval);
    }
    return Json::Value(chain_json);
}

void Minimize::SavePars(string filename, int my_idx_nll){
    NLL_estimator* mynll = NLL_estimator_array[my_idx_nll];
    Amplitude* myamp = mynll->amp_obj;
    int n_par_tot = mynll->N_totpar;
    Json::Value root;

    int num_chain = myamp->nchain;
    int respar_counter(0),LS1par_counter(0),LS2par_counter(0);
    for(int idx_chain=0;idx_chain<num_chain;idx_chain++){
        DecayChain* mychain = &(myamp->array_chain[idx_chain]);
        
        if(mychain->N1_LS==0 || mychain->N2_LS==0){continue;}
        int N_Respar = mychain->intermediate.Get_Num_Dynamic_par();
        int N1_LSpar = mychain->N1_LS;
        int N2_LSpar = mychain->N2_LS;

        for(int idx_respar=0;idx_respar<N_Respar;idx_respar++){
            Para* mypar = &(mynll->Respar_list[respar_counter].par);
            root[mypar->para_name] = WritePar(mypar);
            respar_counter++;
        }

        for(int idx_LSpar=0;idx_LSpar<N1_LSpar;idx_LSpar++){
            Para* mypar = &(mynll->LS1Coeffpar_list[LS1par_counter].par_rho);
            root[mypar->para_name] = WritePar(mypar);
            mypar = &(mynll->LS1Coeffpar_list[LS1par_counter].par_phi);
            root[mypar->para_name] = WritePar(mypar);
            LS1par_counter++;
        }

        for(int idx_LSpar=0;idx_LSpar<N2_LSpar;idx_LSpar++){
            Para* mypar = &(mynll->LS2Coeffpar_list[LS2par_counter].par_rho);
            root[mypar->para_name] = WritePar(mypar);
            mypar = &(mynll->LS2Coeffpar_list[LS2par_counter].par_phi);
            root[mypar->para_name] = WritePar(mypar);
            LS2par_counter++;
        }
    }

    Json::StreamWriterBuilder sw;
    ofstream os;
    os.open(filename);
    Json::StreamWriterBuilder builder;
    builder["commentStyle"] = "None";
    //builder["indentation"] = "";
    builder.settings_["precision"] = 6;	
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    writer->write(root, &os);
	os.close();
}


void Minimize::SaveCovMatrix(){

    Int_t npar_float = gMinuit->GetNumPars();
    //double** cova = new double*[npar_float];
    //for(int i=0;i<npar_float;i++){cova[i] = new double[npar_float];}
    double cova[npar_float][npar_float];
    cov_matrix = new double[npar_float*npar_float];
    gMinuit->mnemat(&cova[0][0],npar_float);

    //Save the Covariant Matrix
    ofstream outfile_Cova;
    outfile_Cova.open("./Cov_matrix.dat");
    for(int par_loop1=0;par_loop1<npar_float;par_loop1++){
        for(int par_loop2=0;par_loop2<npar_float;par_loop2++){
            outfile_Cova<<cova[par_loop1][par_loop2]<<" ";
            cov_matrix[par_loop1*npar_float+par_loop2] = cova[par_loop1][par_loop2];
        }
        outfile_Cova<<endl;
    }
    outfile_Cova.close();

}

void Minimize::GenerateAlterPars(int rndseed, double* genpars){

    Int_t npar_float = gMinuit->GetNumPars();
    double* par_float = new double[npar_float];
    int count = 0;
    for(int par_loop=0;par_loop<npar_tot;par_loop++){
        if(par_error[par_loop]!=0){par_float[count] = par_fitted[par_loop];count++;}
    }

    ROOT::Math::GSLRandomEngine rnd;
    rnd.Initialize();
    rnd.SetSeed(rndseed);
    rnd.GaussianND(npar_float, par_float, cov_matrix, genpars);

    delete par_float;
}

void Minimize::Cal_FitFraction(int my_idx_nll, double** Fit_fraction, double** Fit_fraction_std, TString save_file){
    NLL_estimator* mynll = NLL_estimator_array[my_idx_nll];
    mynll->GetFitFraction(Fit_fraction);
    
    int nchain = mynll->amp_obj->nchain;

    Int_t nstep = 100;
    Int_t npar_float = gMinuit->GetNumPars();
    double* par_float = new double[npar_float];

    double*** Fit_fraction_alter = new double**[nstep];
    for(int i=0;i<nstep;i++){
        Fit_fraction_alter[i] = new double*[nchain];
        for(int j=0;j<nchain;j++){
            Fit_fraction_alter[i][j] =  new double[nchain];
        } 
    }

    for(int idx_step=0;idx_step<nstep;idx_step++){

        GenerateAlterPars(idx_step+666,par_float);

        //Updated pars
        int count = 0;
        for(int i=0;i<npar_tot;i++){
            if(par_error[i]==0) continue;
            for(int idx_nll=0;idx_nll<num_NLL;idx_nll++){
                NLL_estimator* mynll = NLL_estimator_array[idx_nll];
                bool isfind = false;
                Para* mypar = mynll->Search_idx_minuit(i,isfind);
                if(isfind==true){
                    mypar->Val = par_float[count];
                    if(mypar->is_fixed_to==true){
                        mypar->Val = mypar->Get_Val();
                    }
                }
            }
            count++;
        }
        mynll->Update_Paras();

        mynll->GetFitFraction(Fit_fraction_alter[idx_step]);
    }

    //Save the obtained distribution of FFs
    TFile* file = new TFile(save_file, "recreate");
    TTree* my_tree = new TTree("FF_distribution","FF_distribution");
    const int max_nres = 20;
    double FF[max_nres][max_nres];
    my_tree->Branch("FF", FF, Form("FF[%d][%d]/D",max_nres,max_nres));
    for(int i=0;i<nstep;i++){
        for(int idx_ch1=0;idx_ch1<nchain;idx_ch1++){
            for(int idx_ch2=0;idx_ch2<nchain;idx_ch2++){
                FF[idx_ch1][idx_ch2] = Fit_fraction_alter[i][idx_ch1][idx_ch2];
            }
        }
        my_tree->Fill();
    }

    //Calculate the FFs' std
    for(int idx_ch1=0;idx_ch1<nchain;idx_ch1++){
        for(int idx_ch2=0;idx_ch2<nchain;idx_ch2++){

            //By calculating the standard deviation
            double average = 0.;double std = 0;
            for(int i=0;i<nstep;i++){average = average + Fit_fraction_alter[i][idx_ch1][idx_ch2]/nstep;}
            for(int i=0;i<nstep;i++){std = std + pow(Fit_fraction_alter[i][idx_ch1][idx_ch2]-average,2)/nstep;}
            std = sqrt(std);

            //By calculating the 68% area
            double std1 = 0.; int area = 0;
            while(area<(0.68*nstep) && (std1<2*std)){
                std1 = std1 + 0.02*std;
                area = 0;
                for(int i=0;i<nstep;i++){if(fabs(Fit_fraction_alter[i][idx_ch1][idx_ch2]-average)<std1){area++;}}
            }

            Fit_fraction_std[idx_ch1][idx_ch2] = std1;//use the later one
        }
    }

    //Save the obtained FFs and std of FFs
    TMatrixD FF_val(max_nres,max_nres);
    TMatrixD FF_err(max_nres,max_nres);
    for(int idx_ch1=0;idx_ch1<nchain;idx_ch1++){
        for(int idx_ch2=0;idx_ch2<nchain;idx_ch2++){
            FF_val(idx_ch1,idx_ch2) = Fit_fraction[idx_ch1][idx_ch2];
            FF_err(idx_ch1,idx_ch2) = Fit_fraction_std[idx_ch1][idx_ch2];
        }
    }
    FF_val.Write("FF_val");
    FF_err.Write("FF_err");

    file->Write();
    file->Close();

    for(int i=0;i<nchain;i++){
        for(int j=0;j<nchain;j++){
            delete Fit_fraction_alter[i][j];
        }
        delete Fit_fraction_alter[i];
    }
    delete Fit_fraction_alter;
    delete par_float;
    

}

void Minimize::Print_FitFraction(int my_idx_nll, TString save_file){
    std::cout<<"The Fit Fraction of the sample: "<<my_idx_nll<<endl;
    NLL_estimator* mynll = NLL_estimator_array[my_idx_nll];

    int nchain = mynll->amp_obj->nchain;
    double** Fit_fraction = new double*[nchain];
    for(int i=0;i<nchain;i++){Fit_fraction[i] = new double[nchain];}
    double** Fit_fraction_err = new double*[nchain];
    for(int i=0;i<nchain;i++){Fit_fraction_err[i] = new double[nchain];}

    Cal_FitFraction(my_idx_nll,Fit_fraction,Fit_fraction_err,save_file);

    for(int i=0;i<nchain;i++){
        for(int j=0;j<nchain;j++){
            if(i<j) continue;
            std::cout<<mynll->amp_obj->array_chain[i].intermediate.name<<" & "<<mynll->amp_obj->array_chain[j].intermediate.name<<" : ";
            std::cout<<std::fixed<<std::setprecision(4)<<Fit_fraction[i][j]<<" "<<Fit_fraction_err[i][j]<<endl;
        }
    }

    ofstream outfile;
    outfile.open(save_file+".dat");
    for(int i=0;i<nchain;i++){
        for(int j=0;j<nchain;j++){
            outfile<<std::fixed<<std::setprecision(4)<<Fit_fraction[i][j]<<" "<<Fit_fraction_err[i][j]<<" ";
        }
        outfile<<endl;
    }
    outfile.close();

    for(int i=0;i<nchain;i++){delete Fit_fraction[i];}
    delete Fit_fraction;
    for(int i=0;i<nchain;i++){delete Fit_fraction_err[i];}
    delete Fit_fraction_err;
}