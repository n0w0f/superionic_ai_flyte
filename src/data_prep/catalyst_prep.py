from mp_api.client import MPRester
from pymatgen.core import structure


import ase

from pymatgen.analysis.adsorption import AdsorbateSiteFinder
from pymatgen.core.surface import generate_all_slabs
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import structure, Molecule
from ase.build import molecule

from dataclasses import dataclass
from typing import Tuple, List, NamedTuple

from flytekit import task, workflow


from .mp_query import mp_query_id , mp_query_structure , ase_mp_mol , mp_to_ase




def prep_catalyst_workflow_structures(
    mat_composition: str, molecule: str, miller_index: List[int]
) ->  Tuple:
    """Prepare structure, molecule for adsorption workflow
    Args:
        composition (str): composition of material to be queried
        molecule (Molecule) :  molecule that has to be adsorbed
        miller_index (tuple): orientaion of the plane where adsorption has to happen
    Returns:
        tuple (tuple) :  (bulkstructure,molecule,miller_index)

    """



    def ad_molecule(molecule_name="CO2") -> Molecule:
        """molecule structure in materials project format
        Args:
          molecule_name (str) : molecule's name
        Returns:
          mp_molecule (Molecule):  pymatgen molecule structure

        """
        ase_molecule = ase.build.molecule(molecule_name)
        mp_molecule = ase_mp_mol(ase_molecule)

        return mp_molecule

    def create_slab(structure: structure, miller_index: tuple = (1, 1, 1)) -> structure:
        """Create slab structures using pymatgen functions
        Args:
          structure (structure) : unit cell strucutre
          miller_index (tuple) : miller index, plane on which you want to adsorb
        Returns:
          bare_slab (structure): bare slab structure in pymatgen format
        """

        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
        slabs = generate_all_slabs(
            structure,
            max_index=1,
            min_slab_size=10.0,
            min_vacuum_size=10.0,
            center_slab=True,
        )
        bare_slab = [
            slab for slab in slabs if slab.miller_index == tuple(miller_index)
        ][0]

        return bare_slab

    def create_adstructure(structure, molecule, plane=(1, 1, 1)):
        """Create molecule adsorbed on slab structures using pymatgen functions
        Args:
          structure (structure) : slab strucutre
          molecule (Molecule) : MP molecule structure that has to be adsorbed
        Returns:
          ad_structs (structure):  pymatgen molecule + slab adsorbed structure
        """
        asf_slab = AdsorbateSiteFinder(structure)
        # asf = AdsorbateSiteFinder.from_bulk_and_miller(
        #     structure, miller_index=plane, undercoord_threshold=0.10
        # )
        ads_structs = asf_slab.generate_adsorption_structures(
            molecule, translate=True, repeat=[1, 1, 1]
        )

        return ads_structs[0]

    def prep_adsorption_structure(mat_composition: str, molecule: str) -> tuple:
        """Prepare structure, molecule for adsorption workflow
        Args:
            composition (str): composition of material to be queried
            molecule (Molecule) :  molecule that has to be adsorbed
        Returns:
            tuple (tuple) :  (bulkstructure,molecule)

        """

        ids = mp_query_id(mat_composition)
        # ids is a list of different polymorphs in mp
        bulkstructure = mp_query_structure(ids[0])
        mp_molecule = ad_molecule(molecule)

        return (bulkstructure, mp_molecule)

    bulkstructure, mp_mol = prep_adsorption_structure(mat_composition, molecule)
    bare_slab_mp = create_slab(bulkstructure, miller_index)
    ads_structure_mp = create_adstructure(bare_slab_mp, mp_mol)

    bare_slab = mp_to_ase(bare_slab_mp)
    ads_structure = mp_to_ase(ads_structure_mp)
    molecule = ase.build.molecule(molecule)

    return (bare_slab, ads_structure, molecule)

