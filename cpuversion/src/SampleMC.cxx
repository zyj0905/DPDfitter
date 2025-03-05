#include "../include/NLL_estimator.h"
#include "../include/json.h"

void Load_par(Para* mypar, Json::Value root_par, int Nsample, NLL_estimator** mynll){
    string name = mypar->para_name;
    //The first member is not string, not fixed to others
    if(root_par[name][0].isString()==false){
        double parval = root_par[name][0].asDouble(); mypar->Val = parval;
        double parerr = root_par[name][1].asDouble(); mypar->Step = parerr;
        double parlow = root_par[name][2].asDouble(); mypar->low_bound = parlow;
        double parup = root_par[name][3].asDouble(); mypar->up_bound = parup;
    }
    //The first member is string, fixed to others
    else{
        string fix_name = root_par[name][0].asString();
        double fix_ratio = root_par[name][1].asDouble();
        double parval = root_par[name][2].asDouble(); mypar->Val = parval;
    }
}

//arg[1]: 0 is not save component, 1 is save component projection
int main(int argc,char *argv[]){

    int Nsample = argc-2;
    const int max_num = 10;//at most 10 samples
    Amplitude* Amplitude_obj[max_num];
    particle* myres[max_num];
    Dynamic* mydynamic[max_num];
    DecayChain* mychain[max_num];
    LSCoeff_par* LS1par_list[max_num];
    LSCoeff_par* LS2par_list[max_num];
    Reson_par* Respar_list[max_num];
    NLL_estimator* mynll[max_num];
    int A_counter = 0;

    for(int idx_sample=0;idx_sample<Nsample;idx_sample++){
        std::string inputString(argv[idx_sample+2]);
        cout<<"NOW IS LOADING: "<<inputString<<endl;
        ifstream in(inputString, ios::binary);
        Json::Reader reader;
	    Json::Value root;
        reader.parse(in, root);
        //=================================
        //Read Spin-density matrix
        //=================================
        int sizeofsdm = root["SDM"].size();
        double* SDM = new double[sizeofsdm];
        for(int i=0;i<sizeofsdm;i++){SDM[i] = root["SDM"][i].asDouble();}
        //=================================
        //Read mom & daus information
        //=================================
        string mom_name = root["Mom"]["name"].asString();
        double mom_spin = root["Mom"]["spin"].asDouble();
        int mom_p_parity = root["Mom"]["p-parity"].asInt();
        particle Mom(mom_name,mom_spin,mom_p_parity,0);

        string dau1_name = root["Daus"]["Dau1"]["name"].asString();
        double dau1_spin = root["Daus"]["Dau1"]["spin"].asDouble();
        int dau1_p_parity = root["Daus"]["Dau1"]["p-parity"].asInt();
        particle dau1(dau1_name,dau1_spin,dau1_p_parity,0);

        string dau2_name = root["Daus"]["Dau2"]["name"].asString();
        double dau2_spin = root["Daus"]["Dau2"]["spin"].asDouble();
        int dau2_p_parity = root["Daus"]["Dau2"]["p-parity"].asInt();
        particle dau2(dau2_name,dau2_spin,dau2_p_parity,0);

        string dau3_name = root["Daus"]["Dau3"]["name"].asString();
        double dau3_spin = root["Daus"]["Dau3"]["spin"].asDouble();
        int dau3_p_parity = root["Daus"]["Dau3"]["p-parity"].asInt();
        particle dau3(dau3_name,dau3_spin,dau3_p_parity,0);

        particle particle_list[4] = {Mom,dau1,dau2,dau3};

        //=================================
        //Initial the Amplitude class
        //=================================
        Amplitude_obj[idx_sample] = new Amplitude(particle_list);
        Amplitude_obj[idx_sample]->SetSDM(SDM);

        bool is_set_sec = root["set_sec_decay"]["is_set_sec"].asBool();
        if(is_set_sec==true){
            int type_sec = root["set_sec_decay"]["type_sec"].asInt();
            int idx_sec = root["set_sec_decay"]["idx_sec"].asInt();
            Amplitude_obj[idx_sample]->Add_Secondary_Decay(type_sec,idx_sec);
        }

        int sizeofres = root["Res"].size();
        myres[idx_sample] = new particle[sizeofres];
        mydynamic[idx_sample] = new Dynamic[sizeofres];
        mychain[idx_sample] = new DecayChain[sizeofres];

        Json::Value::Members keysofRes = root["Res"].getMemberNames();

        for(int i=0;i<sizeofres;i++){
            string keyname = keysofRes[i];
            string name_res = root["Res"][keyname]["name"].asString();
            double spin = root["Res"][keyname]["spin"].asDouble();
            int p_parity = root["Res"][keyname]["p-parity"].asInt();
            myres[idx_sample][i] = particle(name_res,spin,p_parity,0);
            int sizeofrespar = root["Res"][keyname]["res_par"].size();
            double* res_par = new double[sizeofrespar];
            for(int j=0;j<sizeofrespar;j++){res_par[j] = root["Res"][keyname]["res_par"][j].asDouble();}
            int dynamic_type = root["Res"][keyname]["dynamic_type"].asInt();
            mydynamic[idx_sample][i] = Dynamic(dynamic_type,sizeofrespar,res_par);
            myres[idx_sample][i].Set_Dynamic(mydynamic[idx_sample][i]);
            int idx_isobar = root["Res"][keyname]["idx_isobar"].asInt();
            mychain[idx_sample][i] = DecayChain(myres[idx_sample][i],idx_isobar);
            delete res_par;
        }

        int num_chain = sizeofres;
        //=================================
        //Calculate potential LS
        //=================================
        int nCpar1 = 0;
        int nCpar2 = 0;
        int nRpar = 0;
        for(int idx_chain=0;idx_chain<num_chain;idx_chain++){
            mychain[idx_sample][idx_chain].Cal_LS_auto(particle_list);
            Amplitude_obj[idx_sample]->AddDecayChain(mychain[idx_sample][idx_chain]);
            nCpar1 = nCpar1 + mychain[idx_sample][idx_chain].N1_LS;
            nCpar2 = nCpar2 + mychain[idx_sample][idx_chain].N2_LS;
            nRpar = nRpar + mychain[idx_sample][idx_chain].intermediate.dynamic.num_par;
        }
        //=================================
        //Create corresponding par
        //=================================
        LS1par_list[idx_sample] = new LSCoeff_par[nCpar1];
        LS2par_list[idx_sample] = new LSCoeff_par[nCpar2];
        Respar_list[idx_sample] = new Reson_par[nRpar];
        int C1_counter = 0; int C2_counter = 0; int R_counter = 0;

        for(int idx_chain=0;idx_chain<num_chain;idx_chain++){
            if(mychain[idx_sample][idx_chain].N1_LS==0 || mychain[idx_sample][idx_chain].N2_LS==0){continue;}
            for(int idx_dynamic_par=0;idx_dynamic_par<mychain[idx_sample][idx_chain].intermediate.Get_Num_Dynamic_par();idx_dynamic_par++){
                char name_temp[100];
                sprintf(name_temp, "_Res_par%d", idx_dynamic_par);
                char name1_temp[100];
                sprintf(name1_temp, "s%d_", idx_sample);
                string name = string(name1_temp)+mychain[idx_sample][idx_chain].intermediate.name+string(name_temp);
                Para Respar(name,mychain[idx_sample][idx_chain].intermediate.dynamic.pars[idx_dynamic_par],0.000,0.0,10.0);
                Respar_list[idx_sample][R_counter].idx_chain = idx_chain;
                Respar_list[idx_sample][R_counter].idx_dynamic_par = idx_dynamic_par;
                Respar_list[idx_sample][R_counter].par = Respar;
                Respar_list[idx_sample][R_counter].par.idx_minuit = A_counter; A_counter++;
                R_counter++;
            }
            for(int idx_LS=0;idx_LS<mychain[idx_sample][idx_chain].N1_LS;idx_LS++){
                char name_temp[100];
                sprintf(name_temp, "_LS1coff_L%dS%.1f", mychain[idx_sample][idx_chain].L1[idx_LS], mychain[idx_sample][idx_chain].S1[idx_LS]);
                char name1_temp[100];
                sprintf(name1_temp, "s%d_", idx_sample);
                string name = string(name1_temp)+mychain[idx_sample][idx_chain].intermediate.name+string(name_temp);
                Para LSpar_rho(name+string("_rho"),1.0,0.001,0.0,10.0);
                Para LSpar_phi(name+string("_phi"),0.0,0.001,-4.0,+4.0);
                LS1par_list[idx_sample][C1_counter].idx_chain = idx_chain;
                LS1par_list[idx_sample][C1_counter].idx_LS = idx_LS;
                LS1par_list[idx_sample][C1_counter].par_rho = LSpar_rho;
                LS1par_list[idx_sample][C1_counter].par_rho.idx_minuit = A_counter; A_counter++;
                LS1par_list[idx_sample][C1_counter].par_phi = LSpar_phi;
                LS1par_list[idx_sample][C1_counter].par_phi.idx_minuit = A_counter; A_counter++;
                C1_counter++;
            }
            for(int idx_LS=0;idx_LS<mychain[idx_sample][idx_chain].N2_LS;idx_LS++){
                char name_temp[100];
                sprintf(name_temp, "_LS2coff_L%dS%.1f", mychain[idx_sample][idx_chain].L2[idx_LS], mychain[idx_sample][idx_chain].S2[idx_LS]);
                char name1_temp[100];
                sprintf(name1_temp, "s%d_", idx_sample);
                string name = string(name1_temp)+mychain[idx_sample][idx_chain].intermediate.name+string(name_temp);
                Para LSpar_rho(name+string("_rho"),1.0,0.001,0.0,10.0);
                Para LSpar_phi(name+string("_phi"),0.0,0.001,-4.0,+4.0);
                LS2par_list[idx_sample][C2_counter].idx_chain = idx_chain;
                LS2par_list[idx_sample][C2_counter].idx_LS = idx_LS;
                LS2par_list[idx_sample][C2_counter].par_rho = LSpar_rho;
                LS2par_list[idx_sample][C2_counter].par_rho.idx_minuit = A_counter; A_counter++;
                LS2par_list[idx_sample][C2_counter].par_phi = LSpar_phi;
                LS2par_list[idx_sample][C2_counter].par_phi.idx_minuit = A_counter; A_counter++;
                C2_counter++;
            }
        }

        mynll[idx_sample] = new NLL_estimator(Amplitude_obj[idx_sample],LS1par_list[idx_sample],LS2par_list[idx_sample],Respar_list[idx_sample]);
        //=================================
        //Set data load
        //=================================
        if(is_set_sec==true){
            string p4_sec1 = root["set_sec_decay"]["p4_sec"][0].asString();
            string p4_sec2 = root["set_sec_decay"]["p4_sec"][1].asString();
            mynll[idx_sample]->set_sec_decay(p4_sec1,p4_sec2);
        }

        string infile_name_mc = root["mc"]["filename"].asString();
        string chain_name_mc = root["mc"]["chainname"].asString();
        string p4_final1_mc = root["mc"]["p4_final"][0].asString();
        string p4_final2_mc = root["mc"]["p4_final"][1].asString();
        string p4_final3_mc = root["mc"]["p4_final"][2].asString();

        string infile_name_mcT = root["Truth"]["filename"].asString();
        string chain_name_mcT = root["Truth"]["chainname"].asString();
        string p4_final1_mcT = root["Truth"]["p4_final"][0].asString();
        string p4_final2_mcT = root["Truth"]["p4_final"][1].asString();
        string p4_final3_mcT = root["Truth"]["p4_final"][2].asString();

        mynll[idx_sample]->Load_file(2,infile_name_mc,chain_name_mc,p4_final1_mc,p4_final2_mc,p4_final3_mc);
        mynll[idx_sample]->Load_file(3,infile_name_mcT,chain_name_mcT,p4_final1_mcT,p4_final2_mcT,p4_final3_mcT);

        delete SDM;
    }

    //=================================
    //Load the par list from json file
    //=================================
    for(int idx_sample=0;idx_sample<Nsample;idx_sample++){

        cout<<"NOW IS LOADING Pars of Sample "<<idx_sample<<endl;
        std::string inputString(argv[idx_sample+2]);
        ifstream in(inputString, ios::binary);
        Json::Reader reader;
	    Json::Value root;
        reader.parse(in, root);

        string par_in = root["para_list"]["listin"].asString();
        ifstream in_par(par_in, ios::binary);
        if(in_par.is_open()){
            Json::Reader reader_par;
	        Json::Value root_par;
            reader_par.parse(in_par, root_par);

            //Resonance pars
            for(int idx_res=0;idx_res<mynll[idx_sample]->N_Respar;idx_res++){
                Load_par(&(Respar_list[idx_sample][idx_res].par),root_par,Nsample,mynll);
            }
            //LS1 coupling
            for(int idx_LS=0;idx_LS<mynll[idx_sample]->N1_LSCoeff;idx_LS++){
                Load_par(&(LS1par_list[idx_sample][idx_LS].par_rho),root_par,Nsample,mynll);
                Load_par(&(LS1par_list[idx_sample][idx_LS].par_phi),root_par,Nsample,mynll);
            }
            //LS2 coupling
            for(int idx_LS=0;idx_LS<mynll[idx_sample]->N2_LSCoeff;idx_LS++){
                Load_par(&(LS2par_list[idx_sample][idx_LS].par_rho),root_par,Nsample,mynll);
                Load_par(&(LS2par_list[idx_sample][idx_LS].par_phi),root_par,Nsample,mynll);
            }

        }
        mynll[idx_sample]->Update_Paras();
    }

    auto start = std::chrono::high_resolution_clock::now();

    for(int idx_sample=0;idx_sample<Nsample;idx_sample++){
        cout<<"NOW IS Saving Root of Sample "<<idx_sample<<endl;
        std::string inputString(argv[idx_sample+2]);
        ifstream in(inputString, ios::binary);
        Json::Reader reader;
	    Json::Value root;
        reader.parse(in, root);

        string outfile_mc = root["save_root"]["mc"].asString();
        string infile_name_mc = root["mc"]["filename"].asString();
        string chain_name_mc = root["mc"]["chainname"].asString();
        string p4_final1_mc = root["mc"]["p4_final"][0].asString();
        string p4_final2_mc = root["mc"]["p4_final"][1].asString();
        string p4_final3_mc = root["mc"]["p4_final"][2].asString();
        mynll[idx_sample]->save_root(2,outfile_mc,atoi(argv[1]),infile_name_mc,chain_name_mc,p4_final1_mc,p4_final2_mc,p4_final3_mc);


        string outfile_mcT = root["save_root"]["Truth"].asString();
        string infile_name_mcT = root["Truth"]["filename"].asString();
        string chain_name_mcT = root["Truth"]["chainname"].asString();
        string p4_final1_mcT = root["Truth"]["p4_final"][0].asString();
        string p4_final2_mcT = root["Truth"]["p4_final"][1].asString();
        string p4_final3_mcT = root["Truth"]["p4_final"][2].asString();
        mynll[idx_sample]->save_root(3,outfile_mcT,atoi(argv[1]),infile_name_mcT,chain_name_mcT,p4_final1_mcT,p4_final2_mcT,p4_final3_mcT);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout<<duration.count()<<endl;

    return 0;
}
