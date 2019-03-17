#include <iostream>
#include <memory>
#include "lemon/lemon.hpp"
#include "lemon/options.hpp"
#include "lemon/launch.hpp"
#include "lemon/geometry.hpp"
#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/scoringfunctions/Bernard12.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "spear/Grid.hpp"

using Spear::Bernard12;
using Spear::IDATM;

int main(int argc, char** argv) {
    lemon::Options o;
    std::string distrib("data/csd_distributions.dat");
    o.add_option("--dist,-d", distrib, "Location of the distribution file.");
    o.parse_command_line(argc, argv);

    std::ifstream csd_disbrib(distrib);

    const Spear::AtomicDistributions atomic_distrib =
        Spear::read_atomic_distributions<IDATM>(csd_disbrib);

    const auto rmc = Bernard12::RADIAL | Bernard12::MEAN | Bernard12::COMPLETE;
    const auto rcc = Bernard12::RADIAL | Bernard12::CUMULATIVE | Bernard12::COMPLETE;
    const auto fmc = Bernard12::NORMALIZED_FREQUENCY | Bernard12::MEAN | Bernard12::COMPLETE;
    const auto fcc = Bernard12::NORMALIZED_FREQUENCY | Bernard12::CUMULATIVE | Bernard12::COMPLETE;

    std::vector<std::unique_ptr<Spear::ScoringFunction>> rmcs;
    std::vector<std::unique_ptr<Spear::ScoringFunction>> rccs;
    std::vector<std::unique_ptr<Spear::ScoringFunction>> fmcs;
    std::vector<std::unique_ptr<Spear::ScoringFunction>> fccs;

    for (auto r = 4.0; r <= 15.0; r += 1.0) {
        rmcs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(rmc),
                                                      r, atomic_distrib,
                                                      "IDATM_geometry"));

        rccs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(rcc),
                                                      r, atomic_distrib,
                                                      "IDATM_geometry"));

        fmcs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(fmc),
                                                      r, atomic_distrib,
                                                      "IDATM_geometry"));

        fccs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(fcc),
                                                      r, atomic_distrib,
                                                      "IDATM_geometry"));
    }

    auto worker = [&rmcs, &rccs, &fmcs, &fccs](
                    chemfiles::Frame entry,
                    const std::string& pdbid) {
        // Selection phase
        std::list<size_t> smallm;
        if (lemon::select::small_molecules(entry, smallm) == 0) {
            return std::string("");
        }

        // Pruning phase
        lemon::prune::identical_residues(entry, smallm);
        lemon::prune::cofactors(entry, smallm, lemon::common_cofactors);
        lemon::prune::cofactors(entry, smallm, lemon::common_fatty_acids);

        if (smallm.empty()) {
            return std::string("");
        }

        Spear::Molecule mol(std::move(entry));
        mol.add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);
        auto grid = Spear::Grid(mol.positions());

        // Output phase
        std::string result;
        for (auto smallm_id : smallm) {
            result += pdbid + '\t';
            result += mol.frame().topology().residues()[smallm_id].name() + '\t';
            for (auto& sf : rmcs) {
                result += std::to_string(sf->score(grid, mol, smallm_id)) + '\t';
            }
            for (auto& sf : rccs) {
                result += std::to_string(sf->score(grid, mol, smallm_id)) + '\t';
            }
            for (auto& sf : fmcs) {
                result += std::to_string(sf->score(grid, mol, smallm_id)) + '\t';
            }
            for (auto& sf : fccs) {
                result += std::to_string(sf->score(grid, mol, smallm_id)) + '\t';
            }
            result += '\n';
        }

        return result;
    };

    auto collector = lemon::print_combine(std::cout);
    return lemon::launch(o, worker, collector);
}
