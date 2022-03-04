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

    spat_res(dir_path, scanner, model_type, div, simEnv, mode, log_file)      
    

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
    params_file[('total_dose')]=0.002
    params_file[('simulation_time')]=255 
       
    maps_path = join(dir_path,"NEMA","phantoms","spatialResolution")
    for i in ["0_0", "10_0", "0_10"]:
        shutil.copy(join(maps_path,i+"_act.hdr"),data_path)
        shutil.copy(join(maps_path,i+"_act.img"),data_path)
        shutil.copy(join(maps_path,i+"_att.hdr"),data_path)
        shutil.copy(join(maps_path,i+"_att.img"),data_path)
        
        params_file[('patient_dirname')]=data_file_name
        params_file[('act_map')]=i+"_act.hdr"
        params_file[('att_map')]=i+"_att.hdr"
        params_file[('output_dir')]="NEMA_SpatRes_"+i+"_C"
        params_file[('center_slice')]=86 #center of the FOV
        
        new_params_file_path = join(data_path,"params_SpatRes_"+i+"_C.yml")
        with open(new_params_file_path,"w") as pf:
            yaml.dump(params_file,pf,sort_keys=False)          
        message="Params file created: "+ new_params_file_path
        simpet.tools.log_message(log_file, message, 'info')
        
        message="Starting simulation for "+i+"_center"
        simpet.tools.log_message(log_file, message, 'info')
        
        simu = simpet.SimPET(new_params_file_path)
        simu.run()
        
        params_file[('center_slice')]=130 #1/4 FOV
        params_file[('patient_dirname')]=data_file_name
        params_file[('act_map')]=i+"_act.hdr"
        params_file[('att_map')]=i+"_att.hdr"
        params_file[('output_dir')]="NEMA_SpatRes_"+i+"_OF"
        
        new_params_file_path = join(data_path,"params_SpatRes_"+i+"_OF.yml")
        with open(new_params_file_path,"w") as pf:
            yaml.dump(params_file,pf,sort_keys=False)  
            
        message="Params file created: "+ new_params_file_path
        simpet.tools.log_message(log_file, message, 'info')
        
        message="Starting simulation for "+i+"_1/4 FOV"
        simpet.tools.log_message(log_file, message, 'info')
        
        simu = simpet.SimPET(new_params_file_path)
        simu.run()
    
            
         
    
    
    
    
    
    
    

# def sensitivity():
    
        
if __name__== "__main__":
    main()
