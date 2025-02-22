#include "../include/Minimize.h"

void objective_function(Int_t &nrPar, Double_t *grad, Double_t &Result, Double_t *par, Int_t flag_type){
    Minimize* Minimize_obj = (Minimize*)gMinuit->GetObjectFit();

    int Npar_tot(0);
    double *arr_NLL = new double[Minimize_obj->num_NLL];
    for(int idx_nll=0;idx_nll<Minimize_obj->num_NLL;idx_nll++){Npar_tot = Npar_tot + Minimize_obj->NLL_estimator_array[idx_nll]->N_totpar;}

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

}

void Minimize::Initial_optimize(){
    //Count the par number
    int n_par_sum = 0;
    for(int idx_nll=0;idx_nll<num_NLL;idx_nll++){
        Amplitude* this_amp = NLL_estimator_array[idx_nll]->amp_obj;
        n_par_sum = n_par_sum + (this_amp->Get_total_dynamic_par()) + (this_amp->Get_total_LS1coeff_par()+this_amp->Get_total_LS2coeff_par()) * 2;
    }

    Int_t ierflg = 0;
    double arglist[10];
    gMinuit = new TMinuit(n_par_sum);
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

    int Ntot(0);
    for(int idx_nll=0;idx_nll<num_NLL;idx_nll++){
        NLL_estimator* mynll = NLL_estimator_array[idx_nll];
        int n_par_tot = mynll->N_totpar;
        Ntot = Ntot + n_par_tot;
    }
    
    double *par_fitted = new double[Ntot];
    double *par_error = new double[Ntot];
    for(int par_loop=0;par_loop<Ntot;par_loop++){
        gMinuit->GetParameter(par_loop,par_fitted[par_loop],par_error[par_loop]);
    }

    for(int i=0;i<Ntot;i++){
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

    delete par_fitted;
    delete par_error;
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
