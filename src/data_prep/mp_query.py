from constants import MP_API_KEY
from mp_api.client import MPRester
from pymatgen.core import structure

import ase
from pymatgen.core import structure, Molecule
from pymatgen.io.ase import AseAtomsAdaptor


def mp_query_id(composition: str = "ZnO", mp_api_key: str = MP_API_KEY) -> list:
    """Queries structures from Materials Project Database and returns a list of material project ids with different polymorphs
    Args:
        mp_api_key (str): materials project api key
        composition (str): composition of material to be queried

    Returns:
        list of mp-ids (list): Multiple materials project structure ids

    """
    with MPRester(mp_api_key) as mpr:
        mpr = MPRester(mp_api_key)
        mat_id = mpr.get_material_ids(composition)
    print(mat_id)
    return mat_id


def mp_query_structure(mp_id: int = 0, mp_api_key: str = MP_API_KEY) -> structure:
    """Queries mp structures from Materials Project Database and returns pymatgen structure
    Args:
        mp_api_key (str): materials project api key
        mp_id (str): materials project id of the structre
    Returns:
        structure  (structure):  materials project structure
    """
    with MPRester(mp_api_key) as mpr:
        structure = mpr.get_structure_by_material_id(mp_id)
    print(structure)
    return structure



#@task(cache=False, container_image = "docker.io/akshatvolta/m3gnet:new")
def mp_to_ase(structure: structure) -> ase.atoms:
    """convert ase atoms object to pymatgen structure object
    Args:
      structure (structure) : pymatgen structure object
    Returns:
      structure (ase.Atoms):  ASE atoms object
    """
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    return ase_atoms

#@task(cache=False, container_image = "docker.io/akshatvolta/m3gnet:new")
def ase_mp_mol(ase_mol):
    """convert ase molecule to pymatgen molecule"""
    pos = ase_mol.get_positions()
    symbols = ase_mol.get_chemical_symbols()
    return Molecule(symbols, pos)


