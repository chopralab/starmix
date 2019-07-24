// StarMix: A collection of programs which use the CANDIY suite
// Copyright (C) Purdue University -- BSD license

#include <iostream>
#include <memory>
#include "lemon/lemon.hpp"
#include "lemon/options.hpp"
#include "lemon/launch.hpp"
#include "spear/Molecule.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "spear/Grid.hpp"
#include "spear/Geometry.hpp"

using Spear::IDATM;
using Spear::atomtype_name_for_id;

typedef std::pair<std::string, size_t> DistanceBin;
typedef std::map<DistanceBin, size_t> DistanceCounts;

int main(int argc, char** argv) {
    lemon::Options o;
    auto bin_size = 0.001;
    auto max_dist = 15.0;
    o.add_option("--bin_size,-b", bin_size,
                 "Bin size. Larger value is a coarser potential.");
    o.add_option("--max_dist,-r", max_dist,
                 "Maximum distance");
    o.parse_command_line(argc, argv);

    auto worker = [bin_size,max_dist](chemfiles::Frame entry,
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

        Spear::Molecule mol(entry);
        Spear::IDATM idatm(mol, Spear::AtomType::GEOMETRY);
        auto& positions = mol.positions();
        auto& topo = mol.topology();
        auto grid = Spear::Grid(positions);

        // Output phase
        for (auto smallm_id : smallm) {
            for (auto& smallm_atom : topo.residues()[smallm_id]) {
                auto neighbors = grid.neighbors(positions[smallm_atom], max_dist);
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

                    // Metal ion
                    if (rec_res.size() == 1 &&
                        topo[*rec_res.begin()].charge() > 0) {
                        interaction_good = true;
                    }

                    if (!interaction_good) {
                        continue;
                    }

                    auto& smallm_atom_pos = mol[smallm_atom].position();
                    auto& rec_atom_pos = mol[rec_atom].position();
                    auto dist = Spear::distance(smallm_atom_pos, rec_atom_pos);

                    if (dist > max_dist) {
                        continue;
                    }

                    auto bin_name = rec_res.name() + "_" +
                              topo[rec_atom].name() + "\t" +
                              atomtype_name_for_id<IDATM>(idatm[smallm_atom]);

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
