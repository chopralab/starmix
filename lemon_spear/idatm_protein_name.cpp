#include <iostream>
#include <memory>
#include <lemon/lemon.hpp>
#include <lemon/options.hpp>
#include <spear/Molecule.hpp>
#include "lemon/lemon.hpp"
#include "lemon/launch.hpp"
#include "lemon/geometry.hpp"
#include "spear/scoringfunctions/Bernard12.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "spear/Grid.hpp"

using Spear::Bernard12;
using Spear::IDATM;

int main(int argc, char** argv) {
    lemon::Options o;
    o.parse_command_line(argc, argv);

    auto worker = [](chemfiles::Frame entry,
                     const std::string& pdbid) {
        using idxset = std::unordered_set<size_t>;
        // Selection phase
        auto smallm = lemon::select::small_molecules(entry);

        // Pruning phase
        lemon::prune::identical_residues(entry, smallm);
        lemon::prune::cofactors(entry, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(entry, smallm, lemon::common_fatty_acids);

        if (smallm.empty()) {
            return std::string("");
        }

        Spear::Molecule mol(std::move(entry));
        Spear::IDATM idatm(mol, Spear::AtomType::GEOMETRY);
        auto alltypes = idatm.all_types();
        auto grid = Spear::Grid(mol.frame().positions());
        auto& topo = mol.frame().topology();

        // Output phase
        std::string result;
        for (auto smallm_id : smallm) {
            for (auto& smallm_atom : topo.residues()[smallm_id]) {
                auto neighbors = grid.neighbors(mol.frame().positions()[smallm_atom], 15.0);
                for (auto rec_atom : neighbors) {
                    auto rec_res = *topo.residue_for_atom(rec_atom);

                    const auto& comp_type = rec_res.get("composition_type")->as_string();
                    bool interaction_good = false;

                    if (comp_type == "PEPTIDE-LIKE") {
                        continue;
                    }

                    if (comp_type.find("PEPTIDE") != std::string::npos) {
                        interaction_good = true;
                    }

                    if (comp_type.find("DNA") != std::string::npos ||
                        comp_type.find("RNA") != std::string::npos) {
                        interaction_good = true;
                    }

                    if (lemon::common_cofactors.count(rec_res.name()) != 0) {
                        interaction_good = true;
                    }

                    if (rec_res.name() == "HOH") {
                        interaction_good = true;
                    }

                    if (rec_res.size() == 1 && topo[*rec_res.begin()].charge() == 2) {
                        interaction_good = true;
                    }

                    if (!interaction_good) {
                        continue;
                    }

                    auto dist = mol.frame().distance(smallm_atom, rec_atom);

                    if (dist > 15.0) {
                        continue;
                    }

                    result += pdbid + "\t" + rec_res.name() +
                              "_" + topo[rec_atom].name() + "\t" +
                              Spear::atomtype_name_for_id<Spear::IDATM>(alltypes[smallm_atom]) +
                              "\t" + std::to_string(dist) + "\n";
                }
            }
        }
        return result;
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
