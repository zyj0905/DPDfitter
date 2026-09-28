# DPD fitter

DPD fitter is a c++ implementaiton of the three-body decay amplitudes following derivations in [Dalitz-Plot-Decomposition paper](https://inspirehep.net/literature/1758460).

# Checks on codes

We have made necessary checks on our codes, including numerical comparison with mature tool like TF-PWA for some specific decay chains, validation of behavior by chaning the particle orders, cross-checked by Chat-GPT 5.6 sol (Codes are developed humanly before AI being so powerful...Hopefully, everything should be correct.)

# How to use it

We have provided some examples (e^+e^-\to D*D\pi) about how to perform fit with your GPU-device, shown in ./example/runFitinGPU/; and how to sample MC events according to your fit results, shown in ./example/runSamplinginCPU/

# About input in json style

Input config files & parameter lists are loaded in json style.

For example in the example dir, config file "Dst0DmPip.json" has the following parts:

"SDM":[1,0,0,0,0,0,0,0,1]: spin density matrix, here it is gamma* produced in e+e- collision without longitudinal polarization

"Mom":{"name": "vpho","spin": 1,"p-parity": -1}: Mother particle informations

"Daus":{}: Daughter particle informations

"set_sec_decay":{"is_set_sec": true,"type_sec": 0,"idx_sec": 1,"p4_sec": ["p4_sub_D","p4_sub_Pi"]}: information about the secondary decay, in the example, we have D* further decaying into D pi: 
"type_sec" is the index used in ./gpuversion/src/Amplitude.cu, until now only V->PP secondary decay is considered, one needs to include their cases if not available; 
"idx_sec" is the index of Daughter particle who has the secondary decay, in our case is 1=D* (Index order: mom is 0, Daus is 1,2,3); 
"p4_sec" is the four momenta of D & pi from D* decaying, in a form of "SetPxPyPzE(p4_1[0],p4_1[1],p4_1[2],p4_1[3])" (see ./gpuversion/src/NLL_estimator.cu)

"set_fit_type":{"fit_type": "cFit","PDF_bg": "PDF_bg","bg_ratio": -1}: how the likelihood is constructed, "nFit" means the background contribution is directly subtracted from data, "cFit" means the background contribution is modelled by a PDF and use ln(|M|^2+PDF_bg) for fitting
If one use "cFit", the variable "PDF_bg", namely the branch in root file saving the PDF value of background for each data & MC events, should be given
If one use "cFit", the "bg_ratio" should be given. If "-1" is used, this ratio would be calculated based on the event numbers of input bkg & data root

"Res":{
        "A_Zc_3900":{
            "name": "A_Zc_3900",
            "res_par": [3.9015,0.075,0.135,3.096,2.0325,1.865,2.010],
            "spin": 1,
            "p-parity": 1,
            "dynamic_type": 3,
            "idx_isobar": 3
        },
        "...":{
        }
}
"res_par": is the resonance parameters, such as mass & width, or flatte-like parameters
"dynamic_type": is the type of lineshape, please refer ./gpuversion/include/Dynamic.h, one may add their lineshape by revising it
"idx_isobar": the bachelor particle, for example, gamma*->pi Zc, pi=3 is the bachelor particle

"data":{
        "filename": "./Dst0DmPip_data_withpdf.root",
        "chainname": "ana",
        "p4_final": ["p4_Dst0","p4_Dm","p4_Pip"]
}
"chainname": is the tree name
"p4_final": is the input four momentum, also in "SetPxPyPzE(p4_1[0],p4_1[1],p4_1[2],p4_1[3])", order should be consistent with "Daus"
The same for mc & bg

"Fit_stragety":{
        "is_scan": false,
        "is_migrad": true,
        "is_hesse": true
}
Choice for TMinuit

"para_list":{
        "random_number": -1,
        "listin": "./par_in/Dst0DmPip_par.json",
        "listout": "./par_out/Dst0DmPip_par.json"
},
input parameter list:
"random_number": -1 means not used random initial parameter, other means random seed
"list in": input par list, void "" case would generate one list.
"list out": fit result

"Cal_FitFraction": "./out/output_FFs_Dst0DmPip.root"
If used, the Fit fraction & interference term & corresponding Stat. Un. would be calculated, and saved in a root file as a matrix

"save_root":{
        "data": "./out/output_data_Dst0DmPip.root",
        "mc": "./out/output_mc_Dst0DmPip.root",
        "bg": "./out/output_bg_Dst0DmPip.root"
  }
The output root files, with the calculated kinematic variables & event weight for mc, used for plots

# Run jobs
./Fit.exe Dst0DmPip.json: fit according to input json
./Fit.exe Dst0DmPip.json DstmD0Pip.json: simultenous fit to two channels
The parameters are fixed by "s1_A_Zc_3900_LS1coff_L0S1.0_phi" : ["s0_A_Zc_3900_LS1coff_L0S1.0_phi",1,0]", which means ["Fixed to which par","ratio","fixed value"], for example, here s1_A_Zc_3900_LS1coff_L0S1.0_phi is enforced to be equal to 1*s0_A_Zc_3900_LS1coff_L0S1.0_phi, "0" does not work when input, and would be updated to fixed value, so that One can use it in MC sampling

./SampleMC.exe 1 Dst0DmPip.json: Sampling PHSP MC sample based on fit result. Please notice that the input root should be updated to your PHSP MC in reconstructed and truth-level










