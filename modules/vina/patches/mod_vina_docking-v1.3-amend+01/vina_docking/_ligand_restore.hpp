#pragma once

#include <vector>
#include <string>

#include "_hydrogen_optimize.hpp"

std::vector<std::string> __pybind_export__pdbqt_to_pdb(const std::string &original_pdb, const std::string &original_pdbqt, const std::vector<std::string> &docked_pdbqt);
std::vector<OpenBabel::OBMol*> pdbqt_to_pdb(OpenBabel::OBMol &original_pdb, OpenBabel::OBMol &original_pdbqt, std::vector<OpenBabel::OBMol*> &mols_to_patch);
