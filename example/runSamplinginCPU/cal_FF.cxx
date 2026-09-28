void cal_FF(){

    int Nres = 4;

    TChain MCT("MCT");
    MCT.Add("./out/output_truth_mc_*.root");
    double weight_component[20][20];
    double weight_tot;
    double weight_diag_sum = 0.;
    double FF[20][20] = {0.0};
    memset(FF, 0, sizeof(FF));
    MCT.SetBranchAddress("weight_component",weight_component);
    MCT.SetBranchAddress("weight_tot",&weight_tot);

    for(int i=0;i<MCT.GetEntries();i++){
        MCT.GetEntry(i);
        for(int a=0;a<Nres;a++){
            for(int b=a;b<Nres;b++){
                FF[a][b] = FF[a][b] + weight_component[a][b];
            }
        }
        weight_diag_sum = weight_diag_sum + weight_tot;
    }

    for(int a=0;a<Nres;a++){
        for(int b=0;b<Nres;b++){
            FF[a][b] = FF[a][b]/weight_diag_sum;
            if(FF[a][b]<1E-10){FF[a][b]=0;}
        }
    }

    ofstream outfile;
    outfile.open("Fulltable.dat");

    ofstream outfile1;
    outfile1.open("FF.dat");

    for(int a=0;a<Nres;a++){
        for(int b=0;b<Nres;b++){
            if(b<a){outfile<<" "<<Form("%.3f",0.00); continue;}      
            if(a==b){outfile<<" "<<Form("%.3f",FF[a][b]);outfile1<<Form("%.3f",FF[a][b])<<endl;}
            else{outfile<<" "<<Form("%.3f",FF[a][b]-FF[a][a]-FF[b][b]);}
        }
        outfile<<endl;
    }

}
