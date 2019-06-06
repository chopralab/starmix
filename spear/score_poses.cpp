#include <iostream>
#include "spear/Molecule.hpp"
#include "spear/scoringfunctions/Bernard12.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "spear/Grid.hpp"
#include "spear/Geometry.hpp"
#include "chemfiles.hpp"

using Spear::Bernard12;
using Spear::IDATM;

using sf_vector = std::vector<std::unique_ptr<Spear::ScoringFunction>>;

int main(int argc, char** argv) {
    auto prot = Spear::Molecule(chemfiles::Trajectory(argv[1]).read());
    auto grid = Spear::Grid(prot.positions());
    auto idatm_name = prot.add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);
    auto types1 = prot.atomtype(idatm_name);

    auto lign = Spear::Molecule(chemfiles::Trajectory(argv[2]).read());
    auto types2 = lign.atomtype(lign.add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY));

    std::unordered_set<size_t> all_types;
    std::copy(types1->cbegin(), types1->cend(), std::inserter(all_types, all_types.begin()));
    all_types.erase(47);
    all_types.erase(48);

    std::ifstream csd_disbrib(argv[3]);

    const Spear::AtomicDistributions atomic_distrib =
        Spear::read_atomic_distributions<IDATM>(csd_disbrib);

    const auto rmr = Bernard12::RADIAL | Bernard12::MEAN | Bernard12::REDUCED;
    const auto rcr = Bernard12::RADIAL | Bernard12::CUMULATIVE | Bernard12::REDUCED;
    const auto fmr = Bernard12::NORMALIZED_FREQUENCY | Bernard12::MEAN | Bernard12::REDUCED;
    const auto fcr = Bernard12::NORMALIZED_FREQUENCY | Bernard12::CUMULATIVE | Bernard12::REDUCED;

    sf_vector rmrs;
    sf_vector rcrs;
    sf_vector fmrs;
    sf_vector fcrs;

    const auto rmc = Bernard12::RADIAL | Bernard12::MEAN | Bernard12::COMPLETE;
    const auto rcc = Bernard12::RADIAL | Bernard12::CUMULATIVE | Bernard12::COMPLETE;
    const auto fmc = Bernard12::NORMALIZED_FREQUENCY | Bernard12::MEAN | Bernard12::COMPLETE;
    const auto fcc = Bernard12::NORMALIZED_FREQUENCY | Bernard12::CUMULATIVE | Bernard12::COMPLETE;

    sf_vector rmcs;
    sf_vector rccs;
    sf_vector fmcs;
    sf_vector fccs;

    for (auto r = 4.0; r <= 15.0; r += 1.0) {

        rmrs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(rmr),
                                                      r, atomic_distrib,
                                                      idatm_name, all_types));

        rcrs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(rcr),
                                                      r, atomic_distrib,
                                                      idatm_name, all_types));

        fmrs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(fmr),
                                                      r, atomic_distrib,
                                                      idatm_name, all_types));

        fcrs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(fcr),
                                                      r, atomic_distrib,
                                                      idatm_name, all_types));


        rmcs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(rmc),
                                                      r, atomic_distrib,
                                                      idatm_name));

        rccs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(rcc),
                                                      r, atomic_distrib,
                                                      idatm_name));

        fmcs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(fmc),
                                                      r, atomic_distrib,
                                                      idatm_name));

        fccs.emplace_back(std::make_unique<Bernard12>(
                                                      static_cast<Bernard12::Options>(fcc),
                                                      r, atomic_distrib,
                                                      idatm_name));
    }

    std::cout << "name\t";

    auto print_names = [](std::string sf) {
        for (auto r = 4; r <= 15; r += 1) {
            std::cout << "\t" << sf << r;
        }
    };

    print_names("rmr");
    print_names("rcr");
    print_names("fmr");
    print_names("fcr");
    print_names("rmc");
    print_names("rcc");
    print_names("fmc");
    print_names("fcc");

    std::cout << "\tsize\n";

    auto ltraj = chemfiles::Trajectory(argv[2]);
    while (!ltraj.done()) {
        auto frame = ltraj.read();
        
        std::cout << frame.get<chemfiles::Property::STRING>("name").value_or("XXXX");

        auto mol = Spear::Molecule(frame);
        //mol.remove_hydrogens();
        mol.add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);

        for (auto& sf : rmrs) {
            std::cout << std::to_string(sf->score(grid, prot, mol)) << '\t';
        }
        for (auto& sf : rcrs) {
            std::cout << std::to_string(sf->score(grid, prot, mol)) << '\t';
        }
        for (auto& sf : fmrs) {
            std::cout << std::to_string(sf->score(grid, prot, mol)) << '\t';
        }
        for (auto& sf : fcrs) {
            std::cout << std::to_string(sf->score(grid, prot, mol)) << '\t';
        }
        for (auto& sf : rmcs) {
            std::cout << std::to_string(sf->score(grid, prot, mol)) << '\t';
        }
        for (auto& sf : rccs) {
            std::cout << std::to_string(sf->score(grid, prot, mol)) << '\t';
        }
        for (auto& sf : fmcs) {
            std::cout << std::to_string(sf->score(grid, prot, mol)) << '\t';
        }
        for (auto& sf : fccs) {
            std::cout << std::to_string(sf->score(grid, prot, mol)) << '\t';
        }
        std::cout << mol.size() << '\n';
    }
}
