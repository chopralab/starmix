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

using Spear::IDATM;
using Spear::atomtype_name_for_id;
using Spear::van_der_waals;

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
                    auto dist = mol.frame().distance(smallm_atom, rec_atom);
                    auto vdw_sum = van_der_waals<IDATM>(alltypes[rec_atom]) +
                                   van_der_waals<IDATM>(alltypes[smallm_atom]);

                    if (dist > 15.0 || dist < vdw_sum) {
                        continue;
                    }

                    auto bin_name = atomtype_name_for_id<IDATM>(alltypes[rec_atom]) + "_" +
                                    atomtype_name_for_id<IDATM>(alltypes[smallm_atom]);

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
