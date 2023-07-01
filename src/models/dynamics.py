from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core import Structure

from ase.io import Trajectory
from ase.data import atomic_numbers
from ase import units
from ase.md.analysis import DiffusionCoefficient


from m3gnet.models import MolecularDynamics
from minio_lake.client import download_file

from typing import List, Tuple

'''

from ase.constraints import FixAtoms
elements = [element.symbol for element in structure.composition.elements]
elements.remove('Na')
c1 = FixAtoms(indices=[atom.index for atom in sc_st if atom.symbol == 'P'])
c2 = FixAtoms(indices=[atom.index for atom in sc_st if atom.symbol == 'S'])

sc_st.set_constraint([c1, c2])

sc_st.set_pbc((True, True, True))

'''
from typing import List

def run_molecular_dynamics(local_cif : str , remote_cif : str , config_path : dict,  config_workflow : dict) -> Tuple[List[str] , List[str] , List[str]]:

    download_file(remote_cif, local_cif, )

    filename = local_cif.replace( "_relaxed.cif", "md.traj")   
    traj_filename = filename.replace(config_path['path']['relaxed_save_path'],config_path['path']['md_traj_save_path'])
    lof_file = traj_filename.replace( "md.traj", 'md.log')

    # traj_filename_remote = traj_filename.replace(config_path['path']['md_traj_save_path'],config_path['minio']['md_traj_save_path'])
    # lof_file_remote = traj_filename_remote.replace( "md.traj", 'md.log')


    traj_filename_list : List = []
    lof_file_list : List = []
    diffusion_coefficient_list : List = []


    for temperature in config_workflow['md_dynamics']['list_of_temp']:

        traj_ = traj_filename.replace('.traj', f'_T{temperature}.traj')
        lof_ = lof_file.replace('.log', f'_T{temperature}.log')


        config_workflow['md_dynamics']['temperature'] = temperature
        config_workflow['md_dynamics']['trajectory'] = traj_
        config_workflow['md_dynamics']['logfile'] = lof_

        run_md(config_workflow['md_dynamics'] , local_cif)
        traj_filename_list.append(traj_)
        lof_file_list.append(lof_)  
        diffusion_coefficient_list.append(calculate_diffusion_coefficient(traj_))

    return diffusion_coefficient_list , traj_filename_list , lof_file_list


def run_md(config: dict, path:str) -> None:

    _structure = Structure.from_file(path)
    structure = AseAtomsAdaptor.get_atoms(_structure)
    structure.repeat((2,2,2))
    structure.set_pbc((True, True, True))

    md = MolecularDynamics(
        atoms=structure,
        temperature=config['temperature'],
        ensemble=config['ensemble'],
        timestep=config['timestep'],
        trajectory=config['trajectory'],
        logfile=config['logfile'],
        loginterval=config['loginterval'],
        )
    md.run(steps=config['steps'])




def calculate_diffusion_coefficient(trajectory: str) -> List[float]:
    """
    Calculate the diffusion coefficient from a trajectory file using the pymatgen DiffusionAnalyzer.

    Args:
        trajectory (str): Path to the trajectory file.

    Returns:
        tuple: A tuple containing diffusion coefficients, printed data, and a plot of the diffusion analysis.
    """
    traj = Trajectory(trajectory)
    atoms = traj[0]

    Na_indices = [i for i, atom in enumerate(atoms) if atom.number == atomic_numbers['Na']]

    dc = DiffusionCoefficient( traj =  traj , timestep =  1 * units.fs * 100, atom_indices=Na_indices, molecule=False)
    diffusion_coefficients = dc.get_diffusion_coefficients()
    printed_data = dc.print_data()
    diffusion_plot = dc.plot()

    return diffusion_coefficients
    
    
    
    