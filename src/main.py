import os
import yaml
from typing import List,Tuple
from dataclasses import dataclass, field



from data_prep.structure_builder import prepare_folders, substitute_materials
from utils.manage_files import  save_dataclass_list_to_json
from minio_lake.client import upload_files_to_remote
from minio_lake.client import download_folder_from_remote, upload_file , download_file

from flytekit import task, workflow, dynamic, map_task
from dataclasses_json import dataclass_json

from functools import partial


@dataclass_json
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


@task(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio") 
def pipeline(cif_file : str,  config: dict ) -> str:

    from models.m3gnet_ff import run_relax, predict_formation_energy,predict_bandgap
    from models.screen import check_conductivity, check_stability
    from models.dynamics import run_md , calculate_diffusion_coefficient , run_molecular_dynamics

    from pymatgen.core import Structure


    # Perform geometry optimization (relaxation) on the unrelaxed structure
    local_cif, remote_cif = run_relax(cif_file=cif_file, config_model=config['model']['m3gnet']['relaxer'], config_path=config['data'])
    
    # Perform formation energy prediction on the relaxed structure
    formation_energy = predict_formation_energy(local_cif=local_cif, remote_cif=remote_cif, config=config['model']['m3gnet'])

    print(formation_energy)
    print(formation_energy)
    print(formation_energy)

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
        bandgap = predict_bandgap(local_cif=local_cif, remote_cif=remote_cif, config=config['model']['m3gnet'])
        insulator = check_conductivity(bandgap, config['workflow']['screen_cutoffs']['bandgap'])

        # Create instances of WorkflowResult data class
        screen_result.bandgap = bandgap
        screen_result.insulator = insulator

    else:
        print("The material not stable.")


    #screening  - based on diffusion
    
    if screen_result.insulator and screen_result.stable :
        print("The material fits the criteria.")
        screen_result.fit = True

        
        diffusion_coefficient_list , traj_filename_list , lof_file_list = run_molecular_dynamics(local_cif=local_cif, remote_cif=remote_cif, config_path=config['data'] , config_workflow=config['workflow'] )
        upload_files_to_remote(traj_filename_list, config['data']['path']['md_traj_save_path'], config['data']['remote']['md_traj_save_path'] , config['data']['remote']['bucket_name'])
        upload_files_to_remote(lof_file_list, config['data']['path']['md_traj_save_path'], config['data']['remote']['md_traj_save_path'] , config['data']['remote']['bucket_name'])

        screen_result.diffusion_coefficient = diffusion_coefficient_list
        print(screen_result)
        return "hi"

    else:
        print("The material does not fit the criteria.")



@task(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio")  
def data_prep(minio_path:str , local_path: str) -> List[str]:

    download_file(minio_path,local_path)

    # Read the config.yaml file

    with open(local_path, 'r') as file:
        config_data = yaml.safe_load(file)

    print(config_data)

    materials,substituted_materials = prepare_folders(config_data['data'], config_data['data']['substitution'] )

    raw_cif_paths , substituted_cif_paths = substitute_materials(materials, substituted_materials, config_data['data']['path'],  config_data['data']['substitution'])

    upload_files_to_remote(raw_cif_paths, config_data['data']['path']['raw_save_path'], config_data['data']['remote']['raw_save_path'] , config_data['data']['remote']['bucket_name'])
    upload_files_to_remote(substituted_cif_paths, config_data['data']['path']['processed_save_path'], config_data['data']['remote']['processed_save_path'] ,  config_data['data']['remote']['bucket_name'])
    
    download_folder_from_remote( config_data['data']['remote']['raw_save_path'] , config_data['data']['path']['raw_save_path'], config_data['data']['remote']['bucket_name'])
    download_folder_from_remote( config_data['data']['remote']['processed_save_path'] , config_data['data']['path']['processed_save_path'], config_data['data']['remote']['bucket_name'])
   
    return substituted_cif_paths

@dynamic(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio") 
def parallel_workflow(substituted_cif_paths : List[str], config : dict)->str:
    print(substituted_cif_paths)
    result = list()
    for path in substituted_cif_paths:
        print(path)
        print(path)
        print(path)
        r = pipeline(cif_file=path, config=config)
        print(r) 
    return "hi"
    #return [pipeline(cif_file=path, config=config) for path in substituted_cif_paths]
    #return map_task(partial(pipeline, config=config) , cif_file=substituted_cif_paths)


@dynamic(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio") 
def start() -> None:

    minio_path = "data/config.yaml"
    local_path = "superionic_ai/src/data/config.yaml"

    download_file(object_name=minio_path,file_path=local_path)

    substituted_cif_paths = data_prep(minio_path=minio_path, local_path = local_path)


    # Read the config.yaml file
    with open(local_path, 'r') as file:
        config_data = yaml.safe_load(stream=file)

    parallel_workflow(substituted_cif_paths=substituted_cif_paths, config=config_data)
    
    #results: List[material_class] = []

    # for cif_file in substituted_cif_paths:

    #     workflow_result =  workflow(cif_file=cif_file, config=config_data)
    #     results.append(workflow_result)
    #print(results)
    #save_dataclass_list_to_json(results, config_data['data']['results'])


if __name__ == "__main__":

    start()
        
#     materials,substituted_materials = prepare_folders(config_data['data']['path'] , config_data['data']['substitution'] )

#     raw_cif_paths , substituted_cif_paths = substitute_materials(materials, substituted_materials, config_data['data']['path'],  config_data['data']['substitution'])

#     upload_files_to_remote(raw_cif_paths, config_data['data']['path']['raw_save_path'], config_data['data']['remote']['raw_save_path'] , config_data['data']['remote']['bucket_name'])
#     upload_files_to_remote(substituted_cif_paths, config_data['data']['path']['processed_save_path'], config_data['data']['remote']['processed_save_path'] ,  config_data['data']['remote']['bucket_name'])
    
#     download_folder_from_remote( config_data['data']['remote']['raw_save_path'] , config_data['data']['path']['raw_save_path'], config_data['data']['remote']['bucket_name'])
#     download_folder_from_remote( config_data['data']['remote']['processed_save_path'] , config_data['data']['path']['processed_save_path'], config_data['data']['remote']['bucket_name'])
   


#     results: List[material_class] = []
#     for cif_file in substituted_cif_paths:

#         workflow_result =  workflow(cif_file, config_data)
#         results.append(workflow_result)

#     save_dataclass_list_to_json(results, config_data['data']['results'])



 





