import os
import yaml
from typing import List,Tuple
from dataclasses import dataclass, field

from models.m3gnet import run_relax, predict_formation_energy,predict_bandgap
from models.screen import check_conductivity, check_stability
from models.dynamics import run_md , calculate_diffusion_coefficient , run_molecular_dynamics

from pymatgen.core import Structure


from data_prep.structure_builder import prepare_folders, substitute_materials
from utils.manage_files import  save_dataclass_list_to_json
from minio_lake.client import upload_files_to_remote
from minio_lake.client import download_folder_from_remote, upload_file

from flytekit import task, workflow


# Specify the path to the config.yaml file
config_path = os.path.join("..", "config", "config.yaml")

# Read the config.yaml file
with open(config_path, 'r') as file:
    config_data = yaml.safe_load(file)


@dataclass
class material_class:
    relaxed_cif_path: str
    unrelaxed_cif_path: str
    trjaectory_files: List[str] = None
    md_log_files: List[str] = None
    formation_energy: float = 0.0 
    bandgap: float = 0.0
    stable: bool = False
    insulator: bool = False
    diffusion_coefficient: float = 0.0
    fit: bool = False

@workflow(cache=False, container_image="docker.io/akshatvolta/m3gnet:new")
def workflow(cif_file : str,  config: dict ) -> material_class:


    # Perform geometry optimization (relaxation) on the unrelaxed structure
    local_cif, remote_cif = run_relax(cif_file, config['model']['m3gnet']['relaxer'], config['data'])

    
    # Perform formation energy prediction on the relaxed structure
    formation_energy = predict_formation_energy(local_cif, remote_cif, config['model']['m3gnet'])

    stability = check_stability(formation_energy, config['workflow']['screen_cutoffs']['formation_energy'])


    # Create instances of WorkflowResult data class
    screen_result = material_class(relaxed_cif_path=local_cif,
                                    unrelaxed_cif_path=cif_file,
                                    formation_energy=formation_energy,
                                    stable=stability,
                                    )
    
    #screening  - based on bandgap

    if screen_result.stable:
       
        # Perform bandgap prediction on the relaxed structure
        bandgap = predict_bandgap(local_cif, remote_cif, config['model']['m3gnet'])
        insulator = check_conductivity(bandgap, config['workflow']['screen_cutoffs']['bandgap'])

        # Create instances of WorkflowResult data class
        screen_result.bandgap = bandgap
        screen_result.insulator = insulator

    else:
        print("The material not stable.")
        return screen_result


    #screening  - based on diffusion
    
    if screen_result.insulator and screen_result.stable :
        print("The material fits the criteria.")
        screen_result.fit = True

        
        diffusion_coefficient_list , traj_filename_list , lof_file_list = run_molecular_dynamics(local_cif, remote_cif, config['data'] , config['workflow'] )
        upload_files_to_remote(traj_filename_list, config_data['data']['path']['md_traj_save_path'], config_data['data']['remote']['md_traj_save_path'] , config_data['data']['remote']['bucket_name'])
        upload_files_to_remote(lof_file_list, config_data['data']['path']['md_traj_save_path'], config_data['data']['remote']['md_traj_save_path'] , config_data['data']['remote']['bucket_name'])


        screen_result.diffusion_coefficient = diffusion_coefficient_list


        return screen_result 
    
        


    else:
        print("The material does not fit the criteria.")
        return screen_result

    


@workflow(cache=False, container_image="docker.io/akshatvolta/m3gnet:new")
def start():
    materials,substituted_materials = prepare_folders(config_data['data']['path'] , config_data['data']['substitution'] )

    raw_cif_paths , substituted_cif_paths = substitute_materials(materials, substituted_materials, config_data['data']['path'],  config_data['data']['substitution'])

    upload_files_to_remote(raw_cif_paths, config_data['data']['path']['raw_save_path'], config_data['data']['remote']['raw_save_path'] , config_data['data']['remote']['bucket_name'])
    upload_files_to_remote(substituted_cif_paths, config_data['data']['path']['processed_save_path'], config_data['data']['remote']['processed_save_path'] ,  config_data['data']['remote']['bucket_name'])
    
    download_folder_from_remote( config_data['data']['remote']['raw_save_path'] , config_data['data']['path']['raw_save_path'], config_data['data']['remote']['bucket_name'])
    download_folder_from_remote( config_data['data']['remote']['processed_save_path'] , config_data['data']['path']['processed_save_path'], config_data['data']['remote']['bucket_name'])
   


    results: List[material_class] = []
    for cif_file in substituted_cif_paths:

        workflow_result =  workflow(cif_file, config_data)
        results.append(workflow_result)

    save_dataclass_list_to_json(results, config_data['data']['path']['results'])



if __name__ == "__main__":
        
    materials,substituted_materials = prepare_folders(config_data['data']['path'] , config_data['data']['substitution'] )

    raw_cif_paths , substituted_cif_paths = substitute_materials(materials, substituted_materials, config_data['data']['path'],  config_data['data']['substitution'])

    upload_files_to_remote(raw_cif_paths, config_data['data']['path']['raw_save_path'], config_data['data']['remote']['raw_save_path'] , config_data['data']['remote']['bucket_name'])
    upload_files_to_remote(substituted_cif_paths, config_data['data']['path']['processed_save_path'], config_data['data']['remote']['processed_save_path'] ,  config_data['data']['remote']['bucket_name'])
    
    download_folder_from_remote( config_data['data']['remote']['raw_save_path'] , config_data['data']['path']['raw_save_path'], config_data['data']['remote']['bucket_name'])
    download_folder_from_remote( config_data['data']['remote']['processed_save_path'] , config_data['data']['path']['processed_save_path'], config_data['data']['remote']['bucket_name'])
   


    results: List[material_class] = []
    for cif_file in substituted_cif_paths:

        workflow_result =  workflow(cif_file, config_data)
        results.append(workflow_result)

    save_dataclass_list_to_json(results, config_data['data']['path']['results'])



 





