import json
import numpy as np


with open('./par_out/Dst0DmPip_par.json', 'r', encoding='utf-8') as file:
    data = json.load(file)

#The overall factor

overall_factor_rho_LS1 = data[list(data.keys())[1]][0]
overall_factor_phi_LS1 = data[list(data.keys())[0]][0]

data[list(data.keys())[1]][1] = 0
data[list(data.keys())[0]][1] = 0

for key in data.keys():
    if("Res_par" in key):
        continue
    if("s" in str(data[key][0])):
        continue
    #rescale all chains by the overall factor in LS1
    if("_phi" in key):
        data[key][0] = data[key][0] - overall_factor_phi_LS1
    if("_rho" in key):
        data[key][0] = data[key][0]/overall_factor_rho_LS1
    #rescale the errors
    if(data[key][1]!=0):
        data[key][1] = 1.0


overall_factor_rho_LS2 = -8888
overall_factor_phi_LS2 = -8888

for key in data.keys():
    if("Res_par" in key):
        continue
    if("s" in str(data[key][0])):
        continue
    if("_LS2" in key and "rho" in key and overall_factor_rho_LS2==-8888):
        overall_factor_rho_LS2 = data[key][0]
    if("_LS2" in key and "phi" in key and overall_factor_phi_LS2==-8888):
        overall_factor_phi_LS2 = data[key][0]

for key in data.keys():
    if("Res_par" in key):
        continue
    if("s" in str(data[key][0])):
        continue
    if(("_LS2" in key)==False):
        continue
    #rescale all chains by the overall factor in LS2
    if("_phi" in key):
        data[key][0] = data[key][0] - overall_factor_phi_LS2
    if("_rho" in key):
        data[key][0] = data[key][0]/overall_factor_rho_LS2



list_keys = list(data.keys())
for key in data.keys():
    if("Res_par" in key):
        continue
    if("s" in str(data[key][0])):
        continue
    if("_phi" in key):
        index = list_keys.index(key)
        if("_rho" in list_keys[index + 1] and data[list_keys[index + 1]][0]==0):
            data[key][0] = 0.0

#The chain dependent factor

LS2_status_rho = False
LS2_status_phi = False
LS2_rho_factor = 1.0
LS2_phi_factor = 0.0
arr_rho = []
arr_phi = []

for key in data.keys():
    if("Res_par" in key):
        continue
    if("s" in str(data[key][0])):
        continue
    
    #rescale each chain in the LS2 level
    if("_LS2coff" in key):

        if("_rho" in key and data[key][0]==0 and data[key][1]==0):
            LS2_status_rho = False
            LS2_status_phi = False
            continue
        
        if("_phi" in key):
            index = list_keys.index(key)
            if("_rho" in list_keys[index + 1] and data[list_keys[index + 1]][0]==0):
                LS2_status_rho = False
                LS2_status_phi = False
                continue

        if("_phi" in key):
            if(LS2_status_phi==False):
                LS2_status_phi = True
                LS2_phi_factor = data[key][0]
                data[key][1] = 0.0
                arr_phi.append(LS2_phi_factor)
            data[key][0] = data[key][0] - LS2_phi_factor

        if("_rho" in key):
            if(LS2_status_rho==False):
                LS2_status_rho = True
                LS2_rho_factor = data[key][0]
                data[key][1] = 0.0
                arr_rho.append(LS2_rho_factor)
            data[key][0] = data[key][0] / LS2_rho_factor
    
    if("_LS1coff" in key):
        LS2_status_rho  = False
        LS2_status_phi = False
        LS2_rho_factor = 1.0
        LS2_phi_factor = 0.0

counter = 0
for key in data.keys():
    if("Res_par" in key):
        continue
    if("s" in str(data[key][0])):
        continue
    if(("_LS1" in key)==False):
        continue

    if("_rho" in key and data[key][0]==0 and data[key][1]==0):
        continue
        
    if("_phi" in key):
        index = list_keys.index(key)
        if("_rho" in list_keys[index + 1] and data[list_keys[index + 1]][0]==0):
            continue

    if("_phi" in key):
        data[key][0] = data[key][0] + arr_phi[counter]
    if("_rho" in key):
        data[key][0] = data[key][0] * arr_rho[counter]

    index = list_keys.index(key)
    if("_LS2" in list_keys[index + 1]):
        counter = counter + 1


    

#Control the output precision
for key in data.keys():
    if("Res_par" in key):
        continue
    if("s" in str(data[key][0])):
        continue
    if("_phi" in key and data[key][0]!=0):
        data[key][0] = np.mod(data[key][0],2*np.pi)
        if(data[key][0]>np.pi):
            data[key][0] = data[key][0] - 2*np.pi
        if(data[key][0]<-np.pi):
            data[key][0] = data[key][0] + 2*np.pi
    data[key][0] = round(data[key][0],4)
    data[key][1] = round(data[key][1],4)



json_string = json.dumps(data, ensure_ascii=False, indent=None, separators=(',', ':'))
formatted_json_string = json_string.replace('{', '{\n    ')
formatted_json_string = formatted_json_string.replace('],', '],\n    ')
formatted_json_string = formatted_json_string.replace('}', '\n}')

with open('output.json', 'w', encoding='utf-8') as json_file:
    json_file.write(formatted_json_string)
