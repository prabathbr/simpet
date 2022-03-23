#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from os.path import join, dirname, abspath, isdir, exists
import os, sys
import shutil
import datetime
import yaml
import simpet
from multiprocessing import Process
import numpy as np
import matplotlib.pyplot as plt 

def main():
    
    #INPUT ARGUMENTS
    # scanner_name: string "Discovery_ST", "Discovery_STE", "GE_Advance", "Siemens_mCT", "Vereos"
    # model_type: string "simple_pet", "cylindrical"
    # div: positive integer. Number of divisions (number of threads the user wants to execute)
    # simEnv: 0 ó 1. 0= Local computer. 1= Cesga
    # mode: string "SimSET" or "STIR" (at the moment, only simSET works)
    # log_file: string. Path to a log file.
    
    dir_path = dirname(abspath(__file__))
    
    scanner="Discovery_STE"
    model_type = "simple_pet"    
    div=8
    simEnv=0
    mode="SimSET"
    
    log_file=join(dir_path,"NEMA","NEMA_log_file_"+scanner+".txt")   
    
    ### Sensitivity test
    A_cal= 14.430 # MBq (activity for sensitivity test)
    t_cal = 60 # s (time of the sensitivity assay)
    real_sens_value = 8.8 # cps/kBq (value obtained for the Discovery ST)
    center_slice_Sens = 86
    ###
         
    ### Spatial resolution test
    dose_SR = 0.002 # mCi (activity used in the real test for the Discovery ST)
    length_SR = 1200 # s
    center_slices = [86, 130] # center and 1/4 FOV
    fwhm = {"0_0" : [4.6, 4.6, 4.3], "0_10":[5.7,4.7,6], "10_0":[5.7,4.7,6]} # published fwhm for Discovery_STE
    ### 
    
    ### Image Quality test
    dose_IQ = 1.25 #mCi
    length_IQ = 1800 #s ¿?
    center_slice_IQ = 69
    ###
    
    # simu_sens_value = sensitivity(dir_path, scanner, model_type, div, simEnv, mode, A_cal, t_cal, center_slice_Sens, log_file)
    # sens_factor = simu_sens_value/real_sens_value
    # print(simu_sens_value) # 57.420182893299916
    # print(sens_factor) # 6.525020783329535
    
    sens_factor = 6.525020783329535
    
    length_spat_res = np.round(length_SR/sens_factor,2)
    length_imag_qua = np.round(length_IQ/sens_factor,2)
    
    spat_res(dir_path, scanner, model_type, div, simEnv, mode, dose_SR, length_spat_res, center_slices, fwhm, log_file)  
    # image_quality(dir_path, scanner, model_type, div, simEnv, mode, dose_IQ, length_imag_qua, center_slice_IQ, log_file)
  
        
def spat_res(dir_path, scanner, model_type, divisions, simuEnvironment, mode, dose, length, center_slices, fwhm, log_file):
     
    # message="Starting Spatial Resolution measurements"
    # simpet.tools.log_message(log_file, message, 'info')
    
    # data_file_name="NEMA_spatRes"
    # data_path = join(dir_path,"Data",data_file_name)
    # if not exists(data_path):
    #     os.makedirs(data_path)            
      
    # message="Data folder created: "+data_path
    # simpet.tools.log_message(log_file, message, 'info')
    
    # params_file_path = join(dir_path,"NEMA","params.yml")
    # with open(params_file_path,'rb') as f:
    #     params_file = yaml.load(f.read(), Loader=yaml.FullLoader)
        
    # params_file[('simulation_environment')]=simuEnvironment
    # params_file[('model_type')]=model_type
    # params_file[('sim_type')]=mode
    # params_file[('divisions')]=divisions
    # params_file[('scanner')]=scanner
    
    # params_file[('recons_type')]="FBP2D"
    # params_file[('total_dose')]=dose #mCi
    # params_file[('simulation_time')]= float(length) #s
    # params_file[('patient_dirname')]=data_file_name
       
    # maps_path = join(dir_path,"NEMA","phantoms","spatialResolution")
    # for i in ["0_0", "10_0", "0_10"]:
    #     shutil.copy(join(maps_path,i+"_act.hdr"),data_path)
    #     shutil.copy(join(maps_path,i+"_act.img"),data_path)
    #     shutil.copy(join(maps_path,i+"_att.hdr"),data_path)
    #     shutil.copy(join(maps_path,i+"_att.img"),data_path)
        
    #     params_file[('act_map')]=i+"_act.hdr"
    #     params_file[('att_map')]=i+"_att.hdr"
    #     params_file[('output_dir')]="NEMA_SpatRes_"+i+"_C"
    #     params_file[('center_slice')]=center_slices[0] #center of the FOV (cm)
        
    #     new_params_file_path = join(data_path,"params_SpatRes_"+i+"_C.yml")
    #     with open(new_params_file_path,"w") as pf:
    #         yaml.dump(params_file,pf,sort_keys=False)          
    #     message="Params file created: "+ new_params_file_path
    #     simpet.tools.log_message(log_file, message, 'info')
        
    #     message="Starting simulation for Spatial Resolution"+i+"_center"
    #     simpet.tools.log_message(log_file, message, 'info')
        
    #     simu = simpet.SimPET(new_params_file_path)
    #     simu.run()
        
    #     params_file[('center_slice')]=center_slices[1] #1/4 FOV
        
    #     params_file[('act_map')]=i+"_act.hdr"
    #     params_file[('att_map')]=i+"_att.hdr"
    #     params_file[('output_dir')]="NEMA_SpatRes_"+i+"_OF"
        
    #     new_params_file_path = join(data_path,"params_SpatRes_"+i+"_OF.yml")
    #     with open(new_params_file_path,"w") as pf:
    #         yaml.dump(params_file,pf,sort_keys=False)  
            
    #     message="Params file created: "+ new_params_file_path
    #     simpet.tools.log_message(log_file, message, 'info')
        
    #     message="Starting simulation for Spatial Resolution: "+i+"_1/4 FOV"
    #     simpet.tools.log_message(log_file, message, 'info')
        
    #     simu = simpet.SimPET(new_params_file_path)
    #     simu.run()
        
    compute_spatRes(dir_path, scanner, mode, fwhm, log_file)

def compute_spatRes(dir_path, scanner, mode, fwhm, log_file):
    results_path = join(dir_path,"Results")
    image = join(results_path,"NEMA_SpatRes_0_0_C",mode+"_Sim_"+scanner,"FBP2D","rec_FBP2D.hv")
    half = np.zeros((6,3))
    tenth = np.zeros((6,3))
    
    with open(image) as f:
            lines = f.readlines()

    lines = [x.strip() for x in lines]    
    image_dims = {"size_x": int(lines[16].split()[4]) , "sf_x": float(lines[17].split()[5]) , "size_y": int(lines[19].split()[4]) , "sf_y": float(lines[20].split()[5]) , "size_z": int(lines[22].split()[4]) , "sf_z": float(lines[23].split()[5]) }
    
    cont=0
    for i in ["0_0", "0_10", "10_0"]:
        
        [p_X, p_Y, p_Z] = profiles(join(results_path,"NEMA_SpatRes_"+i+"_C",mode+"_Sim_"+scanner,"FBP2D","rec_FBP2D.hdr"),image_dims, fwhm.get(i), log_file)
    
        half[cont,0], tenth[cont,0]=compute_fw(p_Y,image_dims.get("sf_y"),log_file)
        half[cont,1], tenth[cont,1]=compute_fw(p_Z,image_dims.get("sf_z"),log_file)
        half[cont,2], tenth[cont,2]=compute_fw(p_X, image_dims.get("sf_x"),log_file)
        
        [p_X, p_Y, p_Z] = profiles(join(results_path,"NEMA_SpatRes_"+i+"_OF",mode+"_Sim_"+scanner,"FBP2D","rec_FBP2D.hdr"),image_dims, fwhm.get(i), log_file)
        half[cont+1,0], tenth[cont,0]=compute_fw(p_Y,image_dims.get("sf_y"),log_file)
        half[cont+1,1], tenth[cont,1]=compute_fw(p_Z,image_dims.get("sf_z"),log_file)
        half[cont+1,2], tenth[cont,2]=compute_fw(p_X, image_dims.get("sf_x"),log_file)
        
        cont=cont+2
    
    aux = (half[0,0]+half[0,2]+half[1,0]+half[1,2])/4
    res = {'1cm_H_trans': aux}
    aux = (half[0,1]+half[1,1])/2
    res.update({'1cm_H_axi': aux})
    aux = (tenth[0,0]+tenth[0,2]+tenth[1,0]+tenth[1,2])/4
    res.update({'1cm_T_trans' : aux})
    aux = (tenth[0,1]+tenth[1,1])/2
    res.update({'1cm_T_axi' : aux})
    aux = (half[2,0]+half[3,0]+half[4,2]+half[5,2])/4
    res.update({'10cm_H_rad' : aux})
    aux = (half[2,2]+half[3,2]+half[4,0]+half[5,0])/4
    res.update({'10cm_H_tang' : aux})
    aux = (half[2,1]+half[3,1]+half[4,1]+half[5,1])/4
    res.update({'10cm_H_axi' : aux})
    aux = (tenth[2,0]+tenth[3,0]+tenth[4,2]+tenth[5,2])/4
    res.update({'10cm_T_rad' : aux})
    aux = (tenth[2,2]+tenth[3,2]+tenth[4,0]+tenth[5,0])/4
    res.update({'10cm_T_tang' : aux})
    aux = (tenth[2,1]+tenth[3,1]+tenth[4,1]+tenth[5,1])/4
    res.update({'10cm_T_axi' : aux})
    
    message=("Spatial Resolution Results: \n 1 cm (Half) \n -Transverse: " + 
    str(res.get('1cm_H_trans')) +"\n -Axial: "+ str(res.get('1cm_H_axi'))+"\n"+
    " 1 cm (Tenth) \n -Transverse: "+str(res.get('1cm_T_trans'))+" \n -Axial: "+
    str(res.get('1cm_T_axi'))+"\n \n 10 cm (Half) \n -Radial: "+
    str(res.get('10cm_H_rad'))+"\n -Tangential: "+str(res.get('10cm_H_tang'))+
    "\n -Axial: "+str(res.get('10cm_H_axi')) +"\n 10 cm (Tenth) \n -Radial: "+
    str(res.get('10cm_T_rad'))+"\n -Tangential: "+str(res.get('10cm_T_tang'))+
    "\n -Axial: "+str(res.get('10cm_T_axi')))
    simpet.tools.log_message(log_file, message, 'info')
    print(message)
    
def compute_fw(y,sf,log_file):
    
    max_pos = y.argmax() #search for the position of the maximum in y
    x1 = [max_pos-1, max_pos, max_pos+1]
    
    y1 = y[0,x1]
    p = np.polyfit(x1,y1,2) #parabolic fit
    max_indx = -p[1]/(2*p[0]) # save the first coordinate of the maximum of the parable
    max_parab = np.polyval(p,max_indx)
    
    c_half = max_parab*0.5
    c_tenth = max_parab*0.1
    
    ind1_half, ind2_half = computeIndices(y,c_half,log_file)
    ind1_tenth, ind2_tenth = computeIndices(y,c_tenth,log_file)
    
    d_half = (ind2_half-ind1_half)*sf
    d_tenth = (ind2_tenth-ind1_tenth)*sf
    
    return d_half, d_tenth
    

def computeIndices(y,c,log_file):
    f = False
    i = 0
    while not f : 
        if y[0,i]>c:
            ind1_aux = i
            f=True
        i=i+1
    
    X_1 = ind1_aux-1
    X_2 = ind1_aux
    ind1 = (c-y[0,X_1-1])*(X_2-X_1)/(y[0,X_2-1]-y[0,X_1-1])+X_1
    
    f = False
    i = ind1_aux
    while not f : 
        if y[0,i]<c:
            ind2_aux = i
            f=True
        i=i+1
    
    X_1 = ind2_aux-1
    X_2 = ind2_aux
    ind2 = (c-y[0,X_1-1])*(X_2-X_1)/(y[0,X_2-1]-y[0,X_1-1])+X_1
    
    return ind1, ind2
    
def profiles(image_hdr, image_dims, fwhm, log_file):
    num_vox_x = int(2*np.round(fwhm[0]/image_dims.get("sf_x")))
    num_vox_y = int(2*np.round(fwhm[1]/image_dims.get("sf_y")))
    num_vox_z = int(2*np.round(fwhm[2]/image_dims.get("sf_z")))
    
    sum_Z = np.zeros((image_dims.get("size_x"),image_dims.get("size_y")))
    
    p_X = np.zeros((1,image_dims.get("size_x")))
    p_Y = np.zeros((1,image_dims.get("size_y")))
    p_Z = np.zeros((1,image_dims.get("size_z")))
    
    img, data = simpet.tools.nib_load(image_hdr)
   
    max_in_slices = 0.0
    max_slice = 0
    for z in range(0, image_dims.get("size_z")-1):        
        if max((data[:,:,z]).max(0)) > max_in_slices :
            max_in_slices = max((data[:,:,z]).max(0))
            max_slice = z
            
    
    max_num_vox_z = min(max_slice+num_vox_z, image_dims.get("size_z")-1)
    min_num_vox_z = max(max_slice-num_vox_z,0)
    
    data_aux = data[:,:,max_slice:max_num_vox_z]
    data_aux2 = data[:,:,min_num_vox_z:max_slice-1]
    print(data_aux2.shape)
    sum_Z=data_aux.sum(3)+data_aux2.sum(2) #sum along axis 2 (start in 0)
    #sum_Z = (data[:,:,max_slice:max_num_vox_z]).sum(2) + (data[:,:,min_num_vox_z:max_slice-1]).sum(2)
    print(sum_Z[65,65])
    ind_max_x = (sum_Z.max(0)).argmax()
    ind_max_y = (sum_Z.max(1)).argmax()

    max_num_vox_x = min(ind_max_x + num_vox_x, image_dims.get("size_x")-1)
    min_num_vox_x = max(ind_max_x - num_vox_x, 0) 
    max_num_vox_y = min(ind_max_y + num_vox_y, image_dims.get("size_y")-1)
    min_num_vox_y = max(ind_max_y - num_vox_y, 0)
    
    p_X = (sum_Z[min_num_vox_y:max_num_vox_y,:]).sum(0)
    p_X = p_X.reshape((1,image_dims.get("size_x")))
    p_Y = (sum_Z[:,min_num_vox_x:max_num_vox_x]).sum(1)
    p_Y = p_Y.reshape((1,image_dims.get("size_y")))
    
    ind_max_cols = ((data[:,:,max_slice]).max(0)).argmax()
    ind_max_rows = ((data[:,:,max_slice]).max(1)).argmax()
    min_col = max(ind_max_cols - num_vox_x, 0)
    max_col = min(ind_max_cols + num_vox_x, image_dims.get("size_x")-1)
    min_row = max(ind_max_rows - num_vox_y, 0)
    max_row = min(ind_max_rows + num_vox_y, image_dims.get("size_y")-1)
    
    for k1 in range(min_col,max_col+1,1):
        for k2 in range(min_row,max_row-1,1):
            p_Z = p_Z + (data[k2,k1,:]).reshape((1,image_dims.get("size_z")))
    
    return np.array(p_X), np.array(p_Y), np.array(p_Z)
    

def image_quality(dir_path, scanner, model_type, divisions, simuEnvironment, mode, dose, length, center_slice_IQ, log_file):
     
    message="Starting Image Quality measurements"
    simpet.tools.log_message(log_file, message, 'info')
    
    data_file_name="NEMA_imageQuality"
    data_path = join(dir_path,"Data",data_file_name)
    if not exists(data_path):
        os.makedirs(data_path)            
      
    message="Data folder created: "+data_path
    simpet.tools.log_message(log_file, message, 'info')
    
    params_file_path = join(dir_path,"NEMA","params.yml")
    with open(params_file_path,'rb') as f:
        params_file = yaml.load(f.read(), Loader=yaml.FullLoader)
        
    params_file[('simulation_environment')]=simuEnvironment
    params_file[('model_type')]=model_type
    params_file[('sim_type')]=mode
    params_file[('divisions')]=divisions
    params_file[('scanner')]=scanner
    
    params_file[('recons_type')]="OSEM3D"
    params_file[('total_dose')]= dose #mCi
    params_file[('simulation_time')]= float(length)  #s
    params_file[('patient_dirname')]=data_file_name
       
    maps_path = join(dir_path,"NEMA","phantoms","imageQuality")
    
    shutil.copy(join(maps_path,"IEC_NEMA_ACT.hdr"),data_path)
    shutil.copy(join(maps_path,"IEC_NEMA_ACT.img"),data_path)
    shutil.copy(join(maps_path,"IEC_NEMA_ATT.hdr"),data_path)
    shutil.copy(join(maps_path,"IEC_NEMA_ATT.img"),data_path)
        
    params_file[('act_map')]="IEC_NEMA_ACT.hdr"
    params_file[('att_map')]="IEC_NEMA_ATT.hdr"
    params_file[('output_dir')]="NEMA_imageQuality"
    params_file[('center_slice')]=center_slice_IQ
        
    new_params_file_path = join(data_path,"params_imageQuality.yml")
    with open(new_params_file_path,"w") as pf:
        yaml.dump(params_file,pf,sort_keys=False)          
    message="Params file created: "+ new_params_file_path
    simpet.tools.log_message(log_file, message, 'info')
        
    message="Starting simulation for Image Quality"
    simpet.tools.log_message(log_file, message, 'info')
        
    simu = simpet.SimPET(new_params_file_path)
    simu.run()
             
         
    
def sensitivity(dir_path, scanner, model_type, divisions, simuEnvironment, mode, A_cal, t_cal, center_slice_Sens, log_file):    
    message="Starting Sensitivity measurements"
    simpet.tools.log_message(log_file, message, 'info')
    
    data_file_name="NEMA_sensitivity"
    data_path = join(dir_path,"Data",data_file_name)
    if not exists(data_path):
        os.makedirs(data_path)            
      
    message="Data folder created: "+data_path
    simpet.tools.log_message(log_file, message, 'info')
    
    params_file_path = join(dir_path,"NEMA","params.yml")
    with open(params_file_path,'rb') as f:
        params_file = yaml.load(f.read(), Loader=yaml.FullLoader)
        
    params_file[('simulation_environment')]=simuEnvironment
    params_file[('model_type')]=model_type
    params_file[('sim_type')]=mode
    params_file[('divisions')]=divisions
    params_file[('scanner')]=scanner
    
    total_dose = np.round(A_cal*0.027027027,2)
    params_file[('total_dose')]=float(total_dose) #0.39 #mCi
    params_file[('simulation_time')]=t_cal #60  #s
    params_file[('center_slice')]=center_slice_Sens
    
    params_file[('do_reconstruction')]=0
    
    maps_path = join(dir_path,"NEMA","phantoms","sensitivity")
    shutil.copy(join(maps_path,"0_0_0_act.hdr"),data_path)
    shutil.copy(join(maps_path,"0_0_0_act.img"),data_path)
    
    params_file[('patient_dirname')]=data_file_name
    params_file[('act_map')]="0_0_0_act.hdr"
    
    for i in range(5):        
        shutil.copy(join(maps_path,"0_0_0_att_C"+str(i+1)+".hdr"),data_path)
        shutil.copy(join(maps_path,"0_0_0_att_C"+str(i+1)+".img"),data_path)        
        
        params_file[('att_map')]="0_0_0_att_C"+str(i+1)+".hdr"
        params_file[('output_dir')]="NEMA_Sensitivity_C"+str(i+1)
        
        new_params_file_path = join(data_path,"params_Sensitivity_C"+str(i+1)+".yml")
        with open(new_params_file_path,"w") as pf:
            yaml.dump(params_file,pf,sort_keys=False)  
            
        message="Params file created: "+ new_params_file_path
        simpet.tools.log_message(log_file, message, 'info')
        
        message="Starting simulation for Sensitivity: C"+str(i+1)
        simpet.tools.log_message(log_file, message, 'info')
        
        simu = simpet.SimPET(new_params_file_path)
        simu.run()        
    
        trues_file_hdr=join(dir_path,"Results",params_file.get('output_dir'),mode+"_Sim_"+scanner,"division_0","trues.hdr")
        
        if exists(trues_file_hdr):
            new_folder_SSRB= join(dir_path,"Results",params_file.get('output_dir'),mode+"_Sim_"+scanner,"SSRB")
            if not exists(new_folder_SSRB):
                os.makedirs(new_folder_SSRB)
            new_trues_file_path = join(new_folder_SSRB,"trues.hdr")
            simpet.tools.copy_analyze(trues_file_hdr,new_trues_file_path,new_folder_SSRB,log_file)
            
            from src.stir import stir_tools
            scannerParams = join(dir_path,"scanners",scanner+".yml")
            with open(scannerParams, 'rb') as f:
                scanner_file = yaml.load(f.read(), Loader=yaml.FullLoader)
            sinograms_stir = new_trues_file_path[0:-4]+"_stir.hs"
            stir_tools.create_stir_hs_from_detparams(scanner_file, sinograms_stir)                        
            simpet.tools.convert_simset_sino_to_stir(new_trues_file_path[0:-3]+"img")
            shutil.copy(new_trues_file_path[0:-4]+"_stir.img",sinograms_stir[0:-2]+"s")
            
            max_segment = scanner_file.get("max_segment")
            num_seg = max_segment*2 +1 
            num_views = 1
            do_norm = 0
            config = join(dir_path,"config.yml")
            with open(config,'rb') as f:
                config_file = yaml.load(f.read(), Loader=yaml.FullLoader)                
            sinos_stir_ssrb = stir_tools.SSRB(config_file, sinograms_stir, num_seg, num_views, do_norm, log_file)           
            
            
        else:
            message="Something wrong at the simulation of Sensitivity: C"+str(i+1)
            simpet.tools.log_message(log_file, message, 'error')
     
    sensitivity = compute_sens(dir_path, scannerParams, mode, os.path.basename(sinos_stir_ssrb), A_cal, t_cal, log_file)
    return sensitivity

def compute_sens(dir_path, scannerParams, mode, sinos_name, A_cal, t_cal, log_file):
    # t_acq=60 #seconds
    n_tubes=5
    X_j=np.array([1.25, 2.5, 3.75, 5, 6.25]) # mm
    # A_cal = 14.430 # MBq
    scanner_name=os.path.basename(scannerParams)
    scanner_name=scanner_name[0:-4] #remove extension ".yml"
    with open(scannerParams, 'rb') as f:
        scanner_file = yaml.load(f.read(), Loader=yaml.FullLoader) 
    results_path = join(dir_path,"Results") 
    views = scanner_file.get("num_aa_bins")
    bins = scanner_file.get("num_td_bins") 
    slices = scanner_file.get("max_segment")*2+1
    counts= np.zeros(slices)
    R = np.zeros((n_tubes,slices))
    for i in range(n_tubes):
        rel_path = join(results_path,"NEMA_Sensitivity_C"+str(i+1),mode+"_Sim_"+scanner_name,"SSRB")
        simpet.tools.create_analyze_from_imgdata(join(rel_path,sinos_name[0:-2]+"s"), join(rel_path,sinos_name[0:-2]+"hdr"), bins, slices, views,  1.0, 1.0, 1.0)
        img, data = simpet.tools.nib_load(join(rel_path,sinos_name[0:-2]+"hdr"))
        for j in range(slices):
            counts[j]=sum(sum(data[:,j,:]))
        R[i, :]=counts/t_cal
    
    R_j = sum(np.transpose(R)) 
    x = -2*X_j 

    coeff = np.polyfit(x,np.log(R_j),1)  
    R_corr0 = np.exp(coeff[1])
    S_Tot = R_corr0/A_cal #counts/sec/MBq    
    
    S_i = (R[0, :]/R_j[0])*S_Tot*0.001
    
    plt.figure(1)
    plt.plot(X_j,sum(np.transpose(R/A_cal)))
    plt.savefig(join(dir_path, "NEMA","fig1_sensitivity.png"))
    plt.figure(2)
    plt.plot(np.arange(slices),S_i)
    plt.savefig(join(dir_path, "NEMA","fig2_sensitivity.png"))
    
    return S_Tot*0.001 #counts/sec/kBq
    
if __name__== "__main__":
    main()
