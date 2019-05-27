#include <iostream>
#include "spear/Molecule.hpp"
#include "spear/FunctionalGroup.hpp"
#include "spear/Geometry.hpp"
#include "chemfiles.hpp"

int main(int argc, char **argv) {
    auto in_traj = chemfiles::Trajectory(argv[1]);
    auto ou_traj = chemfiles::Trajectory(argv[2], 'w');

    Spear::FunctionalGroup carboxylic_acid("C(=O)[OH1]");
    while (!in_traj.done()) {
        auto curr_frame = in_traj.read();
        auto my_name = curr_frame.get<chemfiles::Property::STRING>("CAS_NUMBER").value_or("NOCAS");
        curr_frame.set("name", my_name);

        auto curr_mol = Spear::Molecule(curr_frame);
        if (curr_mol.connected_components() != 1) {
            std::cout << "Skipping " << my_name << " as it has too many components." << std::endl;
            continue;
        }

        auto carboxylic_acids = find_functional_groups(curr_mol, carboxylic_acid);
        if (carboxylic_acids.size() == 0) {
            std::cout << "Skipping " << my_name << " as it is not a carboxylic acid." << std::endl;
            continue;
        }

        auto& carboxylic_atoms = *carboxylic_acids.begin();
        curr_frame[carboxylic_atoms[2]].set_charge(-1);

        ou_traj.write(curr_frame);
    }
}
