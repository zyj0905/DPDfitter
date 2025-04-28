#include "TGaxis.h"
#include "TStyle.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGraph.h"
#include "TROOT.h"
#include <iostream>

using namespace std;

const int nRes = 5;
TString ResName[nRes] = {
    "Zc_3900","D0_2300","D1_2420","D1_2430","D2_2460"
};
int color_array[nRes] = {46,12,kRed,27,42};

void Draw(TString chain_info[6], TString var_name, TString weight_name, TString weight_component_name,
    double low, double up, int Nbin, TString Xtitle, TString Ytitle, TString savename){
    
    double var;
    double weight; double weight_component[20][20];

    //Data
    TH1F *Hist_data = new TH1F("Hist_data","",Nbin,low,up);
    TChain* chain_dt = new TChain(chain_info[0]);
    chain_dt->Add(chain_info[1]);
    chain_dt->SetBranchAddress(var_name,&var);
    for(int i=0;i<chain_dt->GetEntries();i++){
        chain_dt->GetEntry(i);
        Hist_data->Fill(var);
    }

    //BKG
    TH1F *Hist_bg = new TH1F("Hist_bg","",Nbin,low,up);
    TChain* chain_bg = new TChain(chain_info[4]);
    chain_bg->Add(chain_info[5]);
    chain_bg->SetBranchAddress(var_name,&var);
    for(int i=0;i<chain_bg->GetEntries();i++){
        chain_bg->GetEntry(i);
        Hist_bg->Fill(var);
    }

    //MC
    TH1F* Hist_mc = new TH1F("Hist_mc","",Nbin,low,up);
    TH1F* Hist_mc_component[nRes];
    for(int i=0;i<nRes;i++){Hist_mc_component[i] = new TH1F(Form("Hist_mc_component%d",i),"",Nbin,low,up);}
    TChain* chain_mc = new TChain(chain_info[2]);
    chain_mc->Add(chain_info[3]);
    chain_mc->SetBranchAddress(var_name,&var);
    chain_mc->SetBranchAddress(weight_name,&weight);
    chain_mc->SetBranchAddress(weight_component_name,weight_component);

    double weight_sum = 0.0;
    for(int i=0;i<chain_mc->GetEntries();i++){chain_mc->GetEntry(i); weight_sum = weight_sum + weight;}
    double scale_factor = double(chain_dt->GetEntries() - chain_bg->GetEntries())/weight_sum;

    for(int i=0;i<chain_mc->GetEntries();i++){
        chain_mc->GetEntry(i);
        Hist_mc->Fill(var,weight*scale_factor);
        for(int j=0;j<nRes;j++){
            Hist_mc_component[j]->Fill(var,weight_component[j][j]*scale_factor);
        }
    }

    // draw Histograms
    TCanvas* c = new TCanvas("Mbc","Mbc",600,600);
    c->Update();

    TPad *c1 =new TPad("c1","",0.0,0.2,1.0,1);
    c1->Draw();
    c1->cd();
    c1->SetFillStyle(4000);
    c1->Range(0,0,1,1);
    c1->SetLeftMargin(0.10);
    c1->SetRightMargin(0.03);
    c1->SetTopMargin(0.10);
    c1->SetBottomMargin(0.10);
    c1->SetFrameFillColor(0);

    Hist_mc->Add(Hist_bg);
    Hist_data->Draw("ep");
    Hist_data->GetYaxis()->SetRangeUser(0.,2*Hist_data->GetMaximum());
    Hist_data->GetYaxis()->SetTitle(Ytitle);
    Hist_data->GetYaxis()->SetTitleSize(0.06);
    Hist_data->GetYaxis()->SetTitleOffset(0.8);
    Hist_data->GetYaxis()->CenterTitle(kTRUE);
    Hist_data->GetXaxis()->CenterTitle(kTRUE);
    Hist_data->GetYaxis()->SetNdivisions(5);
    Hist_data->GetYaxis()->SetMaxDigits(3);
    Hist_data->GetXaxis()->SetNdivisions(5);
    Hist_data->GetXaxis()->SetLabelSize(0.0);
    Hist_data->GetYaxis()->SetLabelSize(0.06);
    Hist_data->SetMarkerStyle(8);
    Hist_data->SetMarkerSize(0.5);
    Hist_data->SetLineWidth(2);
    Hist_mc->SetLineWidth(2);
    Hist_mc->SetLineColor(kBlue);

    Hist_bg->Draw("Histsame");
    Hist_bg->SetFillStyle(3004);
    Hist_bg->SetFillColor(kBlue);
    for(int j=0;j<nRes;j++){
        Hist_mc_component[j]->SetLineWidth(2);
        Hist_mc_component[j]->SetLineStyle(1);
        Hist_mc_component[j]->SetLineColor(color_array[j]);
        Hist_mc_component[j]->Draw("C HIST SAME");
    };

    Hist_mc->Draw("Histsame");
    Hist_data->Draw("epsame");

    double chisq = 0.;
    int N_bin = 0;

    double* pull_Y = new double[Nbin];
    double* pull_X = new double[Nbin];

    for(int i=0;i<Nbin;i++){
        if(Hist_data->GetBinContent(i+1)!=0&&Hist_mc->GetBinContent(i+1)!=0){
            chisq = chisq + pow(Hist_data->GetBinContent(i+1)-Hist_mc->GetBinContent(i+1),2.0)/(Hist_data->GetBinContent(i+1)+Hist_bg->GetBinContent(i+1)); 
            N_bin++;
            pull_Y[i] = (Hist_data->GetBinContent(i+1)-Hist_mc->GetBinContent(i+1))/sqrt(Hist_data->GetBinContent(i+1)+Hist_bg->GetBinContent(i+1));
            pull_X[i] = (up-low)/Nbin*(i+0.5) + low;
        }
    }

    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.045);
    lt.DrawLatex(0.25,0.83, Form("#chi^{2}/Nbin=%.2f/%d",chisq,N_bin));

    TLegend *leg = new TLegend(0.65,0.50,0.80,0.85);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->SetTextSize(0.045);
    leg->SetFillColor(0);
    leg->AddEntry(Hist_data,"Data","lep");
    leg->AddEntry(Hist_bg,"Background","f");
    leg->AddEntry(Hist_mc,"Total fit","l");
    for(int i=0;i<nRes;i++){
        leg->AddEntry(Hist_mc_component[i],ResName[i],"l");
    }
    leg->Draw();

    c->cd();

    TPad *c2 =new TPad("c2","",0.0,0.0,1.0,0.25);
    c2->Draw();
    c2->cd();
    c2->SetFillStyle(4000);
    c2->Range(0,0,1,1);
    c2->SetLeftMargin(0.10);
    c2->SetRightMargin(0.03);
    c2->SetTopMargin(0.10);
    c2->SetBottomMargin(0.50);
    c2->SetFrameFillColor(0);

    TGraph *gr_pull = new TGraph(Nbin,pull_X,pull_Y);
    gr_pull->Draw("AP");
    gr_pull->GetXaxis()->SetTitle(Xtitle);
    gr_pull->GetXaxis()->SetTitleSize(0.17);
    gr_pull->GetYaxis()->SetTitle("#chi");
    gr_pull->GetYaxis()->SetTitleSize(0.17);
    gr_pull->GetXaxis()->SetTitleOffset(1.3);
    gr_pull->GetYaxis()->SetTitleOffset(0.3);
    gr_pull->SetMarkerStyle(8);
    gr_pull->SetMarkerSize(0.5);
    gr_pull->GetXaxis()->SetLabelSize(0.16);
    gr_pull->GetYaxis()->SetLabelSize(0.16);
    gr_pull->GetYaxis()->CenterTitle(kTRUE);
    gr_pull->GetXaxis()->CenterTitle(kTRUE);
    gr_pull->GetYaxis()->SetNdivisions(5);
    gr_pull->GetXaxis()->SetNdivisions(5);
    gr_pull->GetXaxis()->SetLimits(low,up);
    gr_pull->GetYaxis()->SetRangeUser(-6,6);

    TLine *line = new TLine(low,0,up,0);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw("same");

    c->Print(savename);

    delete c1;
    delete c2;
    delete c;
    delete Hist_data;
    delete Hist_bg;
    for(int i=0;i<nRes;i++){delete Hist_mc_component[i];}
    delete Hist_mc;
    delete chain_dt; delete chain_bg; delete chain_mc;

}

void Draw_projection_pull(){

    TGaxis::SetMaxDigits(3);

    TStyle *bes3Style= new TStyle("bes3","bes3 style");

    Int_t icol=0;
    bes3Style->SetFrameBorderMode(icol);
    bes3Style->SetCanvasBorderMode(icol);
    bes3Style->SetPadBorderMode(icol);
    bes3Style->SetPadColor(icol);
    bes3Style->SetCanvasColor(icol);
    bes3Style->SetStatColor(icol);
    bes3Style->SetTitleFillColor(icol);
    bes3Style->SetPalette(1); // set a good color palette

    bes3Style->SetPaperSize(TStyle::kUSLetter);

    bes3Style->SetPadTopMargin(.12);
    bes3Style->SetPadLeftMargin(.15);
    bes3Style->SetPadRightMargin(.08);
    bes3Style->SetPadBottomMargin(.15);

    Int_t font=22;      //times new roman, reg
    Double_t tsize=0.04; //should be set between 0.03-0.05, is in units of "% of pad"

    bes3Style->SetTextFont(font);
    bes3Style->SetTextSize(tsize);
    bes3Style->SetLabelSize(0.04,"xyz");
    bes3Style->SetLabelOffset(0.01,"xyz");
    bes3Style->SetTitleFont(font,"xyz");
    bes3Style->SetLabelFont(font,"xyz");
    bes3Style->SetTitleSize(tsize,"xyz");
    bes3Style->SetTitleXOffset(1.0);
    bes3Style->SetTitleYOffset(1.2); //offset the title of y axis a bit
    bes3Style->SetTitleBorderSize(2.);

    bes3Style->SetMarkerStyle(0);
    bes3Style->SetMarkerSize(0.4);
    bes3Style->SetFrameBorderMode(0.);
    bes3Style->SetFrameLineWidth(1.0);
    bes3Style->SetLineWidth(2.0);
    bes3Style->SetHistLineWidth(1.);

    bes3Style->SetErrorX(0.001);
    bes3Style->SetOptTitle(0);     //no title box
    bes3Style->SetOptStat(0);    //no stat info
    bes3Style->SetOptDate(0);
    bes3Style->SetDateY(.98);
    bes3Style->SetStripDecimals(kFALSE);

    bes3Style->SetEndErrorSize(0.0); //make the end of error bar longer
    gROOT->SetStyle("bes3");
    gROOT->ForceStyle();
    gROOT->Reset();
    


    TString chain_info[6] = {"Data","./out/output_data_Dst*.root","MC","./out/output_mc_Dst*.root","BKG","./out/output_bg_Dst*.root"};
    
    Draw(chain_info,"cos1","weight_tot","weight_component",-1,1,50,"cos#theta_{1}","Events / (0.04)", "./plot/cos1.png");
    Draw(chain_info,"cos2","weight_tot","weight_component",-1,1,50,"cos#theta_{2}","Events / (0.04)", "./plot/cos2.png");
    Draw(chain_info,"cos3","weight_tot","weight_component",-1,1,50,"cos#theta_{3}","Events / (0.04)", "./plot/cos3.png");

    Draw(chain_info,"M23","weight_tot","weight_component",2.00,2.40,40,"M_{23}","Events / (10 (MeV/#font[12]{c}^{2}))", "./plot/M23.png");
    Draw(chain_info,"M13","weight_tot","weight_component",2.10,2.50,40,"M_{13}","Events / (10 (MeV/#font[12]{c}^{2}))", "./plot/M13.png");
    Draw(chain_info,"M12","weight_tot","weight_component",3.80,4.30,50,"M_{12}","Events / (10 (MeV/#font[12]{c}^{2}))", "./plot/M12.png");
    
    Draw(chain_info,"cos23","weight_tot","weight_component",-1,1,50,"cos#theta_{23}","Events / (0.04)", "./plot/cos23.png");
    Draw(chain_info,"cos13","weight_tot","weight_component",-1,1,50,"cos#theta_{13}","Events / (0.04)", "./plot/cos13.png");
    Draw(chain_info,"cos12","weight_tot","weight_component",-1,1,50,"cos#theta_{12}","Events / (0.04)", "./plot/cos12.png");

    Draw(chain_info,"phi23","weight_tot","weight_component",-M_PI,+M_PI,50,"#phi_{23}","Events / (0.12)", "./plot/phi23.png");
    Draw(chain_info,"phi13","weight_tot","weight_component",-M_PI,+M_PI,50,"#phi_{13}","Events / (0.12)", "./plot/phi13.png");
    Draw(chain_info,"phi12","weight_tot","weight_component",-M_PI,+M_PI,50,"#phi_{12}","Events / (0.12)", "./plot/phi12.png");

    Draw(chain_info,"second_cos_c1","weight_tot","weight_component",-1,1,50,"cos#theta^{c1}_{sec}","Events / (0.04)", "./plot/second_cos_c1.png");
    Draw(chain_info,"second_cos_c2","weight_tot","weight_component",-1,1,50,"cos#theta^{c2}_{sec}","Events / (0.04)", "./plot/second_cos_c2.png");
    Draw(chain_info,"second_cos_c3","weight_tot","weight_component",-1,1,50,"cos#theta^{c3}_{sec}","Events / (0.04)", "./plot/second_cos_c3.png");

    Draw(chain_info,"second_phi_c1","weight_tot","weight_component",-M_PI,+M_PI,50,"#phi^{c1}_{sec}","Events / (0.12)", "./plot/second_phi_c1.png");
    Draw(chain_info,"second_phi_c2","weight_tot","weight_component",-M_PI,+M_PI,50,"#phi^{c2}_{sec}","Events / (0.12)", "./plot/second_phi_c2.png");
    Draw(chain_info,"second_phi_c3","weight_tot","weight_component",-M_PI,+M_PI,50,"#phi^{c3}_{sec}","Events / (0.12)", "./plot/second_phi_c3.png");
    
}
