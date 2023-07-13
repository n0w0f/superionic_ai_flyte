import warnings

import os
from typing import Tuple
from m3gnet.models import M3GNet, Relaxer
from pymatgen.io.cif import CifParser
from pymatgen.core import Structure

from minio_lake.client import  upload_file, download_file , download_folder_from_remote
from flytekit import task

# for category in (UserWarning, DeprecationWarning):
#     warnings.filterwarnings("ignore", category=category, module="tensorflow")

# model = M3GNet.load()

from pymatgen.io.cif import CifWriter



def push_traj (traj : str , traj_save : str):
    import pickle
    from ase import Atoms
    from ase.io import Trajectory

    with open(traj, "rb") as f:
        traj = pickle.load(f)

    traj["atom_positions"]
    atoms_list = []
    for frame in range(len(traj['energy'])):
        atoms =  Atoms(symbols=traj["atomic_number"],
                    positions=traj["atom_positions"][frame],
                    cell=traj["cell"][frame],
        )
        atoms.info["energy"] = traj["energy"][frame]
        atoms.info["forces"] = traj["forces"][frame]
        atoms.info["stress"] = traj["stresses"][frame]
        atoms_list.append(atoms)

    trajectory = Trajectory(traj_save, "w")
    for atoms in atoms_list:
        trajectory.write(atoms)

    print("ASE trajectory saved successfully.")
    
    

#@task(cache=False, container_image = "docker.io/akshatvolta/m3gnet:new")
def run_relax(cif_file : str , config_model : dict , config_path : dict) -> Tuple[str, str] :

    updated_cif_filename = cif_file.replace("_updated.cif", "_relaxed.cif")
    relaxed_cif_filename = updated_cif_filename.replace(config_path['path']['processed_save_path'], config_path['path']['relaxed_save_path'])

    # Read the unrelaxed structure from the processed CIF
    unrelaxed_structure = Structure.from_file(cif_file)

    relaxer = Relaxer()  # This loads the default pre-trained model

    relax_results = relaxer.relax(
        unrelaxed_structure,    
        fmax = config_model['fmax'],
        steps = config_model['steps'],
        traj_file = config_model['traj_file'],
        interval = config_model['interval'],
        verbose = config_model['verbose'],
        )

    relaxed_structure = relax_results['final_structure']
    #Energy = float(relax_results['trajectory'].energies[-1]/len(structure))
    os.makedirs(os.path.dirname(relaxed_cif_filename), exist_ok=True)

    cif_writer = CifWriter(relaxed_structure)
    cif_writer.write_file(relaxed_cif_filename)

    print(relaxed_cif_filename)
    print(relaxed_cif_filename)
    print(relaxed_cif_filename)
    print(relaxed_cif_filename)


    file_path_remote = upload_file( relaxed_cif_filename, config_path['path']['relaxed_save_path'], config_path['remote']['relaxed_save_path'] , config_path['remote']['bucket_name'] )


    return relaxed_cif_filename, file_path_remote




#@task(cache=False, container_image = "docker.io/akshatvolta/m3gnet:new")
def predict_formation_energy(local_cif : str, remote_cif :str, config : dict)-> float:

    print("I am in formation energy")

    # download_file(remote_cif, local_cif )

    fe_checkpoints_local = "/superionic_ai/src/data/models/m3gnet_models/matbench_mp_e_form/0/m3gnet"
    folder_path = os.path.join("/superionic_ai/src/data/models/m3gnet_models/matbench_mp_e_form/0/m3gnet")
    os.makedirs(folder_path, exist_ok=True)

    fe_checkpoints_remote = "data/models/m3gnet_models/matbench_mp_e_form/0/m3gnet"
    download_folder_from_remote( fe_checkpoints_remote , fe_checkpoints_local, config['data']['remote']['bucket_name'])

    pymatgen_struct = Structure.from_file(local_cif)
    m3gnet_e_form = M3GNet.from_dir("/superionic_ai/src/data/models/m3gnet_models/matbench_mp_e_form/0/m3gnet")
    e_form_predict = m3gnet_e_form.predict_structure(pymatgen_struct)
    
    return e_form_predict.numpy().tolist()[0].pop()


#@task(cache=False, container_image = "docker.io/akshatvolta/m3gnet:new")
def predict_bandgap(local_cif : str, remote_cif :str, config : dict)-> float:

    # download_file(remote_cif, local_cif, )

    print("I am in bandgap prediction")

    pymatgen_struct = Structure.from_file(local_cif)

    bg_checkpoints_local = "/superionic_ai/src/data/models/m3gnet_models/matbench_mp_gap/0/m3gnet"
    folder_path = os.path.join("/superionic_ai/src/data/models/m3gnet_models/matbench_mp_gap/0/m3gnet")
    os.makedirs(folder_path, exist_ok=True)
    
    bg_checkpoints_remote =  "data/models/m3gnet_models/matbench_mp_gap/0/m3gnet"
    download_folder_from_remote( bg_checkpoints_remote , bg_checkpoints_local, config['data']['remote']['bucket_name'])
    

    m3gnet_bgap = M3GNet.from_dir("/superionic_ai/src/data/models/m3gnet_models/matbench_mp_gap/0/m3gnet")
    bgap_predict = m3gnet_bgap.predict_structure(pymatgen_struct)
    

    return bgap_predict.numpy().tolist()[0].pop()




def relax_catalyst(cif_file : str , config_model : dict , config_path : dict) -> Tuple[str, str] :

    relaxed_cif_filename = cif_file.replace(".cif", "_relaxed.cif")
    traj_filename = cif_file.replace(".cif", "_trajectory.traj")
    ase_traj_filename = cif_file.replace(".cif", "_ase_trajectory.traj")

    # Read the unrelaxed structure from the processed CIF
    unrelaxed_structure = Structure.from_file(cif_file)

    relaxer = Relaxer()  # This loads the default pre-trained model

    relax_results = relaxer.relax(
        unrelaxed_structure,    
        fmax = config_model['fmax'],
        steps = config_model['steps'],
        traj_file = traj_filename,
        interval = config_model['interval'],
        verbose = config_model['verbose'],
        )

    relaxed_structure = relax_results['final_structure']
    Energy = float(relax_results['trajectory'].energies[-1]/len(relaxed_structure))
    os.makedirs(os.path.dirname(relaxed_cif_filename), exist_ok=True)

    cif_writer = CifWriter(relaxed_structure)
    cif_writer.write_file(relaxed_cif_filename)

    print(relaxed_cif_filename)
    print(relaxed_cif_filename)
    print(relaxed_cif_filename)
    print(relaxed_cif_filename)

    push_traj(traj_filename,ase_traj_filename)


    file_path_remote = upload_file( relaxed_cif_filename, config_path['path']['root_path'], config_path['remote']['root_path'] , config_path['remote']['bucket_name'] )
    traj = upload_file( ase_traj_filename, config_path['path']['root_path'], config_path['remote']['root_path'] , config_path['remote']['bucket_name'] )
    print(traj)


    return Energy , file_path_remote
