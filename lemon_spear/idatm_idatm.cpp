#include <iostream>
#include "lemon/lemon.hpp"
#include "lemon/options.hpp"
#include "lemon/launch.hpp"
#include "spear/Molecule.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "spear/Grid.hpp"

using Spear::IDATM;
using Spear::atomtype_name_for_id;
using Spear::van_der_waals;

typedef std::pair<std::string, size_t> DistanceBin;
typedef std::map<DistanceBin, size_t> DistanceCounts;

int main(int argc, char** argv) {
    lemon::Options o;
    auto bin_size = 0.001;
    auto max_dist = 15.0;
    auto vdw_coef = 0.75;
    o.add_option("--bin_size,-b", bin_size,
                 "Bin size. Larger value is a coarser potential.");
    o.add_option("--max_dist,-r", max_dist,
                 "Maximum distance");
    o.add_option("--vdw_coef,-c", vdw_coef,
                 "Van der Waals scaling coeficient. Used to determine min distance");
    o.parse_command_line(argc, argv);

    auto worker = [bin_size,max_dist,vdw_coef](chemfiles::Frame entry,
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
        auto& alltypes = idatm.all_types();
        auto& positions = mol.frame().positions();
        auto& topo = mol.frame().topology();
        auto grid = Spear::Grid(positions);

        // Output phase
        for (auto smallm_id : smallm) {
            for (auto& smallm_atom : topo.residues()[smallm_id]) {
                auto neighbors = grid.neighbors(positions[smallm_atom], max_dist);
                for (auto rec_atom : neighbors) {
                    auto dist = mol.frame().distance(smallm_atom, rec_atom);
                    auto vdw_sum = van_der_waals<IDATM>(alltypes[rec_atom]) +
                                   van_der_waals<IDATM>(alltypes[smallm_atom]);

                    if (dist > max_dist || dist < vdw_sum * vdw_coef) {
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
