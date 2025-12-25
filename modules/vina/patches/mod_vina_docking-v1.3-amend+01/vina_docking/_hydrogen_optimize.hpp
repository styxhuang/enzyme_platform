#pragma once
#include <string>
#include <vector>

#include <openbabel/mol.h>

std::string __pybind_export__hydrogen_optimize(const std::string &pdb_content, int steps, double e_conv, const std::string & forcefield);
OpenBabel::OBMol*  hydrogen_optimize(OpenBabel::OBMol mol, int steps, double e_conv, const std::string & forcefield);
OpenBabel::OBMol*  hydrogen_optimize(OpenBabel::OBMol mol, int steps, double e_conv, const std::string & forcefield, const std::vector<unsigned int> & consider_atoms);
