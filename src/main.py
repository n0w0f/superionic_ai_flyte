import os
import yaml
from typing import List,Tuple
from dataclasses import dataclass, field
from dataclasses import asdict



from data_prep.structure_builder import prepare_folders, substitute_materials
from utils.manage_files import  save_dataclass_list_to_json , catalyst_cif_paths
from minio_lake.client import upload_files_to_remote
from minio_lake.client import download_folder_from_remote, upload_file , download_file

from flytekit import task, workflow, dynamic, map_task, TaskMetadata, Resources
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
def pipeline_battery(cif_file : str,  minio_path : str, local_path :str ) -> int:


    from models.m3gnet_ff import run_relax, predict_formation_energy,predict_bandgap
    from models.screen import check_conductivity, check_stability
    from models.dynamics import run_md , calculate_diffusion_coefficient , run_molecular_dynamics

    from pymatgen.core import Structure


    download_file(minio_path,local_path)

    # Read the config.yaml file
    with open(local_path, 'r') as file:
        config = yaml.safe_load(file)


    # Perform geometry optimization (relaxation) on the unrelaxed structure
    local_cif, remote_cif = run_relax(cif_file=cif_file, config_model=config['model']['m3gnet']['relaxer'], config_path=config['data'])
    print(local_cif)



    # Perform formation energy prediction on the relaxed structure
    formation_energy = predict_formation_energy(local_cif=local_cif, remote_cif=remote_cif, config=config)
    print(formation_energy)


    stability = check_stability(formation_energy, config['workflow']['screen_cutoffs']['formation_energy'])
    print(stability)

    
    # Create instances of WorkflowResult data class
    screen_result = material_class(relaxed_cif_path=local_cif,
                                    unrelaxed_cif_path=cif_file,
                                    formation_energy=formation_energy,
                                    stable=stability,
                                    )

    #screening  - based on bandgap
    print(screen_result)


    if screen_result.stable:

        
        # Perform bandgap prediction on the relaxed structure
        bandgap = predict_bandgap(local_cif=local_cif, remote_cif=remote_cif, config=config)
        print(bandgap)


        insulator = check_conductivity(bandgap, config['workflow']['screen_cutoffs']['bandgap'])
        print(insulator)

        
        # Create instances of WorkflowResult data class
        screen_result.bandgap = bandgap
        screen_result.insulator = insulator

    else:
        print("The material not stable.")
        return 0


    #screening  - based on diffusion
    
    if screen_result.insulator and screen_result.stable :
        print("The material fits the criteria.")
        screen_result.fit = True

        
        diffusion_coefficient_list , traj_filename_list , lof_file_list = run_molecular_dynamics(local_cif=local_cif, remote_cif=remote_cif, config_path=config['data'] , config_workflow=config['workflow'] )
        upload_files_to_remote(traj_filename_list, config['data']['path']['md_traj_save_path'], config['data']['remote']['md_traj_save_path'] , config['data']['remote']['bucket_name'])
        upload_files_to_remote(lof_file_list, config['data']['path']['md_traj_save_path'], config['data']['remote']['md_traj_save_path'] , config['data']['remote']['bucket_name'])

        screen_result.diffusion_coefficient = diffusion_coefficient_list
        screen_result.md_log_files = lof_file_list
        screen_result.trjaectory_files = traj_filename_list
        print(screen_result)
        return 0

    else:
        print("The material does not fit the criteria.")
        return 0



@task(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio")  
def data_prep(minio_path:str , local_path: str) -> List[str]:

    download_file(minio_path,local_path)

    # Read the config.yaml file
    with open(local_path, 'r') as file:
        config_data = yaml.safe_load(file)

    download_file(minio_path,local_path)

    materials,substituted_materials = prepare_folders(config_data['data'], config_data['data']['substitution'] )

    raw_cif_paths , substituted_cif_paths = substitute_materials(materials, substituted_materials, config_data['data']['path'],  config_data['data']['substitution'])

    upload_files_to_remote(raw_cif_paths, config_data['data']['path']['raw_save_path'], config_data['data']['remote']['raw_save_path'] , config_data['data']['remote']['bucket_name'])
    upload_files_to_remote(substituted_cif_paths, config_data['data']['path']['processed_save_path'], config_data['data']['remote']['processed_save_path'] ,  config_data['data']['remote']['bucket_name'])
    
    download_folder_from_remote( config_data['data']['remote']['raw_save_path'] , config_data['data']['path']['raw_save_path'], config_data['data']['remote']['bucket_name'])
    download_folder_from_remote( config_data['data']['remote']['processed_save_path'] , config_data['data']['path']['processed_save_path'], config_data['data']['remote']['bucket_name'])
   
    
    return substituted_cif_paths


@dynamic(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio") 
def parallel_workflow_battery(substituted_cif_paths : List[str], minio_path : str, local_path: str)->int:
    
    print(substituted_cif_paths)
    #result = list()
    for path in substituted_cif_paths:
        print(path)

        r = pipeline_battery(cif_file=path, minio_path=minio_path, local_path = local_path)
        print(r) 
    return 0



@dynamic(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio") 
def start() -> None:

    minio_path = "data/config.yaml"
    local_path = "superionic_ai/src/data/config.yaml"

    download_file(object_name=minio_path,file_path=local_path)

    substituted_cif_paths = data_prep(minio_path=minio_path, local_path = local_path)

    parallel_workflow_battery(substituted_cif_paths=substituted_cif_paths, minio_path=minio_path, local_path = local_path)


'''
==========================================================
Catalyst work flow after here
==========================================================
'''


@task(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio")  
def data_prep_catalyst(minio_path:str , local_path: str) -> List[str]:
    import pickle
    from data_prep.catalyst_prep import create_data

    download_file(minio_path,local_path,bucket_name="catalyst")

    # Read the config.yaml file
    with open(local_path, 'r') as file:
        config_data = yaml.safe_load(file)

    print(config_data)

    
    with open('/home/nawaf/spectra/mp_data_dicts.pickle', 'rb') as file:
        mp_data_dicts = pickle.load(file)
    
    root_path = config_data['data']['path']['root_path']

    slab_cif_paths : List[str] = []
    adsorbate_cif_paths : List[str] = []

    for index, material in enumerate(mp_data_dicts):    

        bare_slab, ads_structure, molecule = prep_catalyst_workflow_structures(
        catalyst=material['structure'], molecule="CO2", miller_index=[1, 1, 1]
        )

        folder_name = material['material_id']
        slab_cif_path , adsorbate_cif_path = create_data( mpid=folder_name, catalyst=bare_slab, adsorbate=ads_structure , root_path=root_path )
        slab_cif_paths.append(slab_cif_path)
        adsorbate_cif_paths.append(adsorbate_cif_path)
    
    print(slab_cif_path)
    print(adsorbate_cif_path)
    upload_files_to_remote(slab_cif_paths, config_data['data']['path']['root_path'], config_data['data']['remote']['root_path'] , config_data['data']['remote']['bucket_name'])
    upload_files_to_remote(adsorbate_cif_paths, config_data['data']['path']['root_path'], config_data['data']['remote']['root_path'] ,  config_data['data']['remote']['bucket_name'])
    
    download_folder_from_remote( config_data['data']['remote']['root_path'] , config_data['data']['path']['root_path'], config_data['data']['remote']['bucket_name'])
    #download_folder_from_remote( config_data['data']['remote']['root_path'] , config_data['data']['path']['root_path'], config_data['data']['remote']['bucket_name'])

    return slab_cif_paths , adsorbate_cif_paths


@dataclass_json
@dataclass
class adsorption_class:
    slab_relaxed_cif_path: str
    adsorbate_relaxed_cif_path: str
    adsorption_energy: float = 0.0
 
import json

@task(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio", limits=Resources(mem="2Gi", cpu="2")) 
def pipeline_catalyst(slab_cif_file : str, adsorbate_cif_file : str, minio_path : str, local_path :str ) -> int:


    from models.m3gnet_ff import relax_catalyst
    from pymatgen.core import Structure


    download_file(minio_path,local_path)

    # Read the config.yaml file
    with open(local_path, 'r') as file:
        config = yaml.safe_load(file)


    # Perform geometry optimization (relaxation) on the unrelaxed structure
    slab_energy , slab_file_path_remote = relax_catalyst(cif_file=slab_cif_file, config_model=config['model']['m3gnet']['relaxer'], config_path=config['data'])
    print(slab_file_path_remote)


    adsorbate_energy , adsorbate_file_path_remote = relax_catalyst( cif_file=adsorbate_cif_file, config_model=config['model']['m3gnet']['relaxer'], config_path=config['data'])
    print(adsorbate_file_path_remote)

    result_json = slab_cif_file.replace("_slab.cif", "_result.json")

    screen_result = adsorption_class(slab_relaxed_cif_path=slab_file_path_remote,
                                adsorbate_relaxed_cif_path=adsorbate_file_path_remote,
                                adsorption_energy=float(adsorbate_energy - slab_energy),
                                )
    json_data = [asdict(screen_result)]

    with open(result_json, 'w') as json_file:
        json.dump(json_data, json_file, indent=4)

    
    file_path_remote = upload_file( result_json, config['data']['path']['root_path'], config['data']['remote']['root_path'] , config['data']['remote']['bucket_name'] )
    print(file_path_remote)

    return 0



@dynamic(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio") 
def parallel_workflow_catalyst(slab_cif_paths : List[str],  adsorbate_cif_paths : List[str],    minio_path : str, local_path: str)->int:
    
    print(slab_cif_paths)
    print(adsorbate_cif_paths)
    num_materials = 900

    # for slab, adsorbate in zip(slab_cif_paths[:num_materials], adsorbate_cif_paths[:num_materials]):
    
    #     r = pipeline_catalyst(slab_cif_file = slab, adsorbate_cif_file = adsorbate, minio_path=minio_path, local_path = local_path)
    #     print(r) 

    r = map_task( partial(pipeline_catalyst,  
                      minio_path=minio_path, 
                      local_path = local_path), 
                metadata=TaskMetadata(retries=1),
                concurrency=20, 
                min_success_ratio=0.75, ) (slab_cif_file = slab_cif_paths[:num_materials], adsorbate_cif_file = adsorbate_cif_paths[:num_materials] )

    return 0



@dynamic(cache=False, container_image="docker.io/aswanthkrshna/m3gnet:minio") 
def start_catalyst() -> None:

    minio_path = "config_catalyst.yaml"
    local_path = "catalyst/data/config_catalyst.yaml"

    download_file(object_name=minio_path,file_path=local_path)
    with open(local_path, 'r') as file:
        config_data = yaml.safe_load(file)

    print(config_data)

    download_folder_from_remote( config_data['data']['remote']['root_path'] , config_data['data']['path']['root_path'], config_data['data']['remote']['bucket_name'])

    slab_cif_paths , adsorbate_cif_paths = catalyst_cif_paths("catalyst/structures")

    parallel_workflow_catalyst(slab_cif_paths=slab_cif_paths, adsorbate_cif_paths=adsorbate_cif_paths, minio_path=minio_path, local_path = local_path)





from data_prep.catalyst_prep import prep_catalyst_workflow_structures

if __name__ == "__main__":

    minio_path = "config_catalyst.yaml"
    local_path = "catalyst/data/config_catalyst.yaml"

    #slab_cif_paths , adsorbate_cif_paths = data_prep_catalyst(minio_path=minio_path , local_path=local_path)

    slab_cif_paths , adsorbate_cif_paths = catalyst_cif_paths("catalyst/structures")

    print(len(slab_cif_paths) , len(adsorbate_cif_paths))

    parallel_workflow_catalyst(slab_cif_paths=slab_cif_paths, adsorbate_cif_paths=adsorbate_cif_paths, minio_path=minio_path, local_path = local_path)



        


 





