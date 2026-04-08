from pymatgen.core import Lattice, Structure
from pymatgen.io.ase import AseAtomsAdaptor
import warnings
import matgl
from matgl.ext.ase import MolecularDynamics, PESCalculator, Relaxer

# To suppress warnings for clearer output
warnings.simplefilter("ignore")

pot = matgl.load_model("M3GNet-MP-2021.2.8-PES")
relaxer = Relaxer(potential=pot)
struct = Structure.from_file("mp-12403.cif")
struct.to_file("1.cif")
relax_results = relaxer.relax(struct, fmax=0.01)
# # extract results
final_structure = relax_results["final_structure"]
final_energy = relax_results["trajectory"].energies[-1]
# # print out the final relaxed structure and energy
#
print(final_structure)
print(f"The final energy is {float(final_energy):.3f} eV.")
