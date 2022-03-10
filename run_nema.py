#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from os.path import join, dirname, abspath, isdir, exists
import os, sys
import shutil
import datetime
import yaml
import simpet
from multiprocessing import Process

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
    log_file=join(dir_path,"NEMA","NEMA_log_file.txt")            

    sensitivity(dir_path, scanner, model_type, div, simEnv, mode, log_file)
    #
    #spat_res(dir_path, scanner, model_type, div, simEnv, mode, log_file)  
    image_quality(dir_path, scanner, model_type, div, simEnv, mode, log_file)
    # return 
        
def spat_res(dir_path, scanner, model_type, divisions, simuEnvironment, mode, log_file):
     
    message="Starting Spatial Resolution measurements"
    simpet.tools.log_message(log_file, message, 'info')
    
    data_file_name="NEMA_spatRes"
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
    
    params_file[('recons_type')]="FBP2D"
    params_file[('total_dose')]=0.002 #mCi
    params_file[('simulation_time')]=255 #s
    params_file[('patient_dirname')]=data_file_name
       
    maps_path = join(dir_path,"NEMA","phantoms","spatialResolution")
    for i in ["0_0", "10_0", "0_10"]:
        shutil.copy(join(maps_path,i+"_act.hdr"),data_path)
        shutil.copy(join(maps_path,i+"_act.img"),data_path)
        shutil.copy(join(maps_path,i+"_att.hdr"),data_path)
        shutil.copy(join(maps_path,i+"_att.img"),data_path)
        
        params_file[('act_map')]=i+"_act.hdr"
        params_file[('att_map')]=i+"_att.hdr"
        params_file[('output_dir')]="NEMA_SpatRes_"+i+"_C"
        params_file[('center_slice')]=86 #center of the FOV (cm)
        
        new_params_file_path = join(data_path,"params_SpatRes_"+i+"_C.yml")
        with open(new_params_file_path,"w") as pf:
            yaml.dump(params_file,pf,sort_keys=False)          
        message="Params file created: "+ new_params_file_path
        simpet.tools.log_message(log_file, message, 'info')
        
        message="Starting simulation for Spatial Resolution"+i+"_center"
        simpet.tools.log_message(log_file, message, 'info')
        
        simu = simpet.SimPET(new_params_file_path)
        simu.run()
        
        params_file[('center_slice')]=130 #1/4 FOV
        
        params_file[('act_map')]=i+"_act.hdr"
        params_file[('att_map')]=i+"_att.hdr"
        params_file[('output_dir')]="NEMA_SpatRes_"+i+"_OF"
        
        new_params_file_path = join(data_path,"params_SpatRes_"+i+"_OF.yml")
        with open(new_params_file_path,"w") as pf:
            yaml.dump(params_file,pf,sort_keys=False)  
            
        message="Params file created: "+ new_params_file_path
        simpet.tools.log_message(log_file, message, 'info')
        
        message="Starting simulation for Spatial Resolution: "+i+"_1/4 FOV"
        simpet.tools.log_message(log_file, message, 'info')
        
        simu = simpet.SimPET(new_params_file_path)
        simu.run()
    
def image_quality(dir_path, scanner, model_type, divisions, simuEnvironment, mode, log_file):
     
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
    params_file[('total_dose')]= 1.25 #mCi
    params_file[('simulation_time')]=100  #s
    params_file[('patient_dirname')]=data_file_name
       
    maps_path = join(dir_path,"NEMA","phantoms","imageQuality")
    
    shutil.copy(join(maps_path,"IEC_NEMA_ACT.hdr"),data_path)
    shutil.copy(join(maps_path,"IEC_NEMA_ACT.img"),data_path)
    shutil.copy(join(maps_path,"IEC_NEMA_ATT.hdr"),data_path)
    shutil.copy(join(maps_path,"IEC_NEMA_ATT.img"),data_path)
        
    params_file[('act_map')]="IEC_NEMA_ACT.hdr"
    params_file[('att_map')]="IEC_NEMA_ATT.hdr"
    params_file[('output_dir')]="NEMA_imageQuality"
    params_file[('center_slice')]=69
        
    new_params_file_path = join(data_path,"params_imageQuality.yml")
    with open(new_params_file_path,"w") as pf:
        yaml.dump(params_file,pf,sort_keys=False)          
    message="Params file created: "+ new_params_file_path
    simpet.tools.log_message(log_file, message, 'info')
        
    message="Starting simulation for Image Quality"
    simpet.tools.log_message(log_file, message, 'info')
        
    simu = simpet.SimPET(new_params_file_path)
    simu.run()
             
         
    
def sensitivity(dir_path, scanner, model_type, divisions, simuEnvironment, mode, log_file):    
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
    
    params_file[('total_dose')]=0.39 #mCi
    params_file[('simulation_time')]=60  #s
    params_file[('center_slice')]=86
    
    params_file[('do_reconstruction')]=0
    
    maps_path = join(dir_path,"NEMA","phantoms","sensitivity")
    shutil.copy(join(maps_path,"0_0_0_act.hdr"),data_path)
    shutil.copy(join(maps_path,"0_0_0_act.img"),data_path)
    
    params_file[('patient_dirname')]=data_file_name
    params_file[('act_map')]="0_0_0_act.hdr"
    
    for i in range(6):        
        shutil.copy(join(maps_path,"0_0_0_att_C"+str(i)+".hdr"),data_path)
        shutil.copy(join(maps_path,"0_0_0_att_C"+str(i)+".img"),data_path)        
        
        params_file[('att_map')]="0_0_0_att_C"+str(i)+".hdr"
        params_file[('output_dir')]="NEMA_Sensitivity_C"+str(i)
        
        new_params_file_path = join(data_path,"params_Sensitivity_C"+str(i)+".yml")
        with open(new_params_file_path,"w") as pf:
            yaml.dump(params_file,pf,sort_keys=False)  
            
        message="Params file created: "+ new_params_file_path
        simpet.tools.log_message(log_file, message, 'info')
        
        message="Starting simulation for Sensitivity: C"+str(i)
        simpet.tools.log_message(log_file, message, 'info')
        
        simu = simpet.SimPET(new_params_file_path)
        simu.run()        
    
        trues_file_hdr=join(dir_path,"Results",params_file.get('output_dir'),"SimSET_Sim_"+scanner,"division_0","trues.hdr")
        
        if exists(trues_file_hdr):
            new_folder_SSRB= join(dir_path,"Results",params_file.get('output_dir'),"SimSET_Sim_"+scanner,"SSRB")
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
            message="Something wrong at the simulation of Sensitivity: C"+str(i)
            simpet.tools.log_message(log_file, message, 'error')
            
    sensitivity = compute_sens(dir_path, scannerParams, os.path.basename(sinos_stir_ssrb), log_file)
    return sensitivity

def compute_sens(dir_path, scannerParams, sinos_name, log_file):
    
    scanner_name=os.path.basename(scannerParams)
    with open(scannerParams, 'rb') as f:
        scanner_file = yaml.load(f.read(), Loader=yaml.FullLoader) 
    results_path = join(dir_path,"Results") 
    views = scanner_file.get("num_aa_bins")
    bins = scanner_file.get("num_td_bins") 
    slices = scanner_file.get("max_segment")*2+1
    
    for i in range(6):
        rel_path = join(results_path,"NEMA_Sensitivity_C"+str(i),"SimSET_Sim_"+scanner_name,"SSRB")
        simpet.tools.create_analyze_from_imgdata(join(rel_path,sinos_name[0:-2]+"s"), join(rel_path,sinos_name[0:-2]+"hdr"), bins, views, slices, 1.0, 1.0, 1.0)
        img, data = tools.nib_load(join(rel_path,sinos_name[0:-2]+"hdr")
        for j in range(slices):
            #
    
if __name__== "__main__":
    main()
