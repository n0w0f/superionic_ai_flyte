import warnings

from typing import Tuple
from m3gnet.models import M3GNet, Relaxer
from pymatgen.io.cif import CifParser
from pymatgen.core import Structure

from minio_lake.client import  upload_file, download_file
from flytekit import task

# for category in (UserWarning, DeprecationWarning):
#     warnings.filterwarnings("ignore", category=category, module="tensorflow")

# model = M3GNet.load()

from pymatgen.io.cif import CifWriter

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

    cif_writer = CifWriter(relaxed_structure)
    cif_writer.write_file(relaxed_cif_filename)

    file_path_remote = upload_file( relaxed_cif_filename, config_path['path']['relaxed_save_path'], config_path['remote']['relaxed_save_path'] , config_path['remote']['bucket_name'] )


    return relaxed_cif_filename, file_path_remote




#@task(cache=False, container_image = "docker.io/akshatvolta/m3gnet:new")
def predict_formation_energy(local_cif : str, remote_cif :str, config : dict)-> float:

    download_file(remote_cif, local_cif )

    pymatgen_struct = Structure.from_file(local_cif)
    m3gnet_e_form = M3GNet.from_dir(config['checkpoints']['formation_energy_checkpoint'])
    e_form_predict = m3gnet_e_form.predict_structure(pymatgen_struct)
    
    return e_form_predict.numpy().tolist()[0].pop()


#@task(cache=False, container_image = "docker.io/akshatvolta/m3gnet:new")
def predict_bandgap(local_cif : str, remote_cif :str, config : dict)-> float:

    download_file(remote_cif, local_cif, )

    pymatgen_struct = Structure.from_file(local_cif)
    

    m3gnet_bgap = M3GNet.from_dir(config['checkpoints']['bandgap_checkpoint'])
    bgap_predict = m3gnet_bgap.predict_structure(pymatgen_struct)
    

    return bgap_predict.numpy().tolist()[0].pop()

