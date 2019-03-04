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

typedef std::pair<std::string, size_t> DistanceBin;
typedef std::map<DistanceBin, size_t> DistanceCounts;

int main(int argc, char** argv) {
    lemon::Options o;
    o.parse_command_line(argc, argv);

    auto bin_size = 0.001;

    auto worker = [bin_size](chemfiles::Frame entry,
                     const std::string& pdbid) {
        DistanceCounts bins;

        // Selection phase
        auto smallm = lemon::select::small_molecules(entry);

        // Pruning phase
        lemon::prune::identical_residues(entry, smallm);
        lemon::prune::cofactors(entry, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(entry, smallm, lemon::common_fatty_acids);

        if (smallm.empty()) {
            return bins;
        }

        Spear::Molecule mol(std::move(entry));
        Spear::IDATM idatm(mol, Spear::AtomType::GEOMETRY);
        auto alltypes = idatm.all_types();
        auto grid = Spear::Grid(mol.frame().positions());
        auto& topo = mol.frame().topology();

        // Output phase
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

                    auto bin_name = rec_res.name() + "_" +
                              topo[rec_atom].name() + "\t" +
                              Spear::atomtype_name_for_id<Spear::IDATM>(alltypes[smallm_atom]);

                    auto dist_bin = static_cast<size_t>(std::floor(dist / bin_size));

                    DistanceBin sbin = {bin_name, dist_bin};
                    auto bin_iterator = bins.find(sbin);

                    if (bin_iterator == bins.end()) {
                        bins[sbin] = 1;
                        continue;
                    }

                    ++(bin_iterator->second);
                }
            }
        }
        return bins;
    };


    DistanceCounts total;
    auto collector = lemon::map_combine<DistanceCounts>(total);
    lemon::launch(o, worker, collector);

    for (const auto& i : total) {
        std::cout << i.first.first << "\t"
                  << static_cast<double>(i.first.second) * bin_size << "\t"
                  << i.second << "\n";
    }

    return 0;
}
