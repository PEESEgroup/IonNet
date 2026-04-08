from pymatgen.core import Structure, Composition
from math import gcd
import random
from fractions import Fraction
from functools import reduce

def read_structure(cif_path):
    structure = Structure.from_file(cif_path)
    print(f"Original formula: {structure.composition.reduced_formula}")
    return structure


def determine_minimal_supercell(target_formula):
    target_comp = Composition(target_formula)

    fractions = [Fraction(amount).limit_denominator() for amount in target_comp.values()]

    denominators = [frac.denominator for frac in fractions]

    min_supercell_factor = reduce(lambda x, y: x * y // gcd(x, y), denominators)
    print(min_supercell_factor)

    return min_supercell_factor

def scale_target_formula(target_comp, scale_factor):
    # expand the formula
    scaled_target_comp = {elem: target_comp[elem] * scale_factor for elem in target_comp}
    scaled_target_formula = Composition(scaled_target_comp)
    return scaled_target_formula

def generate_doped_structure(original_structure, target_comp):
    # calculate the minimum cell
    supercell_dimensions = determine_minimal_supercell(target_comp)
    expanded_structure = original_structure * [supercell_dimensions, 1, 1]


    # get the chemical formulas
    original_comp = scale_target_formula(original_structure.composition, supercell_dimensions)


    # build supercell
    scaled_target_comp = scale_target_formula(target_comp, supercell_dimensions)
    print(original_comp, scaled_target_comp)

    # replace
    new_structure = expanded_structure.copy()
    for elem in original_comp:
        if elem in scaled_target_comp:
            # check whether it has new element
            if scaled_target_comp[elem] < original_comp[elem]:
                print(scaled_target_comp[elem], original_comp[elem])
                # calculate the number of elements that need to be replaced
                replace_count = original_comp[elem] - scaled_target_comp[elem]
                replace_count = replace_count
                orig_indices = [i for i, site in enumerate(new_structure) if site.specie.symbol == str(elem)]

                # randomly select the atom
                if replace_count > len(orig_indices):
                    replace_count = len(orig_indices)

                selected_indices = random.sample(orig_indices, int(replace_count))

                for index in selected_indices:
                    # replace with the new elements
                    for new_elem in scaled_target_comp:
                        if new_elem != elem:
                            new_structure[index] = new_elem
                            break

    return new_structure

def main():
    cif_path = 'Li-structures\\mp-15999.cif'  
    original_structure = read_structure(cif_path)
    reduced_formula, factor = original_structure.composition.get_reduced_formula_and_factor()
    target_formula = "Li6AuAgS4" 
    target_comp = Composition(target_formula)
    target_comp = {elem: target_comp[elem] * factor for elem in target_comp}
    target_comp = Composition(target_comp)

    new_structure = generate_doped_structure(original_structure, target_comp)
    print(new_structure)
    filename = 'doped_structure.cif'
    new_structure.to(filename=filename)
    print(f"Saved as {filename}")

if __name__ == "__main__":
    main()

