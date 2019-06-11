#include <iostream>

#include "spear/Molecule.hpp"
#include "spear/scoringfunctions/VinaScore.hpp"
#include "spear/atomtypes/VinaType.hpp"
#include "spear/Grid.hpp"
#include "spear/Geometry.hpp"
#include "chemfiles.hpp"

using sf_vector = std::vector<std::unique_ptr<Spear::ScoringFunction>>;

int main(int argc, char** argv) {
    auto prot = Spear::Molecule(chemfiles::Trajectory(argv[1]).read());
    auto grid = Spear::Grid(prot.positions());
    prot.add_atomtype<Spear::VinaType>();

    Spear::VinaScore scoring_func;

    std::cout << "name\tg1\tg2\trep\thydrogen\thydrophobic\tvina\n";

    auto ltraj = chemfiles::Trajectory(argv[2]);

    while (!ltraj.done()) {
        auto frame = ltraj.read();
        
        std::cout << frame.get<chemfiles::Property::STRING>("name").value_or("XXXX");
        std::cout << "\t";

        auto mol = Spear::Molecule(frame);
        mol.add_atomtype<Spear::VinaType>();

        auto thing = scoring_func.calculate_components(grid, prot, mol);
        std::cout << thing.g1 << "\t";
        std::cout << thing.g2 << "\t";
        std::cout << thing.rep << "\t";
        std::cout << thing.hydrogen  << "\t";
        std::cout << thing.hydrophobic << "\t";
        std::cout << scoring_func.score(grid, prot, mol) << "\n";
    }
}
