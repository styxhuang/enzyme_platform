#include <cmath>
#include <iostream>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obconversion.h>
#include <openbabel/elements.h>
#include <openbabel/oberror.h>
#include <openbabel/forcefield.h>
#include <openbabel/obiter.h>

#include "_collision.hpp"

inline void
mask_error(void){
    OpenBabel::obErrorLog.SetOutputLevel(OpenBabel::obError);
}

double distance(double * r1, double * r2){
    double distance = 0;
    for (int i = 0; i < 3; i++){
        distance += (r1[i] - r2[i]) * (r1[i] - r2[i]);
    }
    return sqrt(distance);
}

std::vector<double> * get_coord_covalent_radii(const std::string & pdb_content){
    mask_error();

    auto mol = OpenBabel::OBMol();
    auto conv = OpenBabel::OBConversion();

    conv.SetInAndOutFormats("pdb", "pdb");
    conv.ReadString(&mol, pdb_content);

    std::size_t n_atoms = mol.NumAtoms();
    auto result = new std::vector<double>(n_atoms*4);
    std::size_t idx_atom = 0;

    for (OpenBabel::OBMolAtomIter atom(mol); atom; ++atom){
        (*result)[idx_atom * 4 + 0] = atom->GetX();
        (*result)[idx_atom * 4 + 1] = atom->GetY();
        (*result)[idx_atom * 4 + 2] = atom->GetZ();
        (*result)[idx_atom * 4 + 3] = OpenBabel::OBElements::GetCovalentRad(atom->GetAtomicNum());
        idx_atom++;
    }

    return result;
}

std::vector<double> get_atomic_gradient_norm(const std::string & pdb_content, const std::string & forcefield){
    mask_error();

    auto mol = OpenBabel::OBMol();
    auto conv = OpenBabel::OBConversion();

    conv.SetInAndOutFormats("pdb", "pdb");
    conv.ReadString(&mol, pdb_content);

    std::size_t n_atoms = mol.NumAtoms();
    std::vector<double> result(n_atoms, 1e8);

    auto ff = OpenBabel::OBForceField::FindForceField(forcefield);
    bool ff_ready = ff->Setup(mol);

    if (!ff_ready){
        std::cerr << "Fail to setup forcefield" << std::endl;
        return result;
    }

    auto energy = ff->Energy(true);
    
    for (auto atom_idx = 1; atom_idx <= n_atoms; atom_idx++){
        auto grad = ff->GetGradient(mol.GetAtom(atom_idx));
        auto grad_norm = grad.length();
        result[atom_idx - 1] = grad_norm;
    }

    return result;

}

std::vector<collisions> * geometry_aspect(std::vector<std::string *> &pdb_contents, COLLISION_SCOPE scope, float flexible = 1.05){
    std::vector<std::vector<double> *> coord_covalent_radii;

    for (auto content = pdb_contents.begin(); content != pdb_contents.end(); ++content){
        coord_covalent_radii.push_back(
            get_coord_covalent_radii(**content)
        );
    }

    auto result = new std::vector<collisions>;

    // means at least we want inter
    if (scope != COLLISION_SCOPE::intra){
        std::size_t idx_pdb = 0;

        for (auto coord_radii = coord_covalent_radii.begin(); coord_radii != coord_covalent_radii.end(); ++coord_radii){
            auto n_atoms = (*coord_radii)->size() / 4;

            // generate collision list
            for (unsigned long idx_atom_i = 0; idx_atom_i < n_atoms; idx_atom_i++){
                auto atom_i = (*coord_radii)->data() + (idx_atom_i * 4);
                for (unsigned long idx_atom_j = idx_atom_i + 1; idx_atom_j < n_atoms; idx_atom_j++){
                    auto atom_j = (*coord_radii)->data() + (idx_atom_j * 4);
                    double dist_i_j = distance(atom_i, atom_j);
                    double threshold = (atom_i[3] + atom_j[3]) * (2 - flexible);
                    if (dist_i_j <= threshold){
                        result->push_back(std::make_tuple(idx_pdb, idx_atom_i + 1, idx_pdb, idx_atom_j + 1));
                    }
                }
            }

            idx_pdb++;
        }
    }

    if (scope != COLLISION_SCOPE::inter){
        std::size_t n_pdbs = pdb_contents.size();

        // iterate pdb
        for (std::size_t idx_pdb_i = 0; idx_pdb_i < n_pdbs; idx_pdb_i++){
            for (std::size_t idx_pdb_j = idx_pdb_i + 1; idx_pdb_j < n_pdbs; idx_pdb_j++){

                auto n_atoms_i = coord_covalent_radii[idx_pdb_i]->size() / 4;
                auto n_atoms_j = coord_covalent_radii[idx_pdb_j]->size() / 4;

                // iterate atom
                for (unsigned long idx_atom_i = 0; idx_atom_i < n_atoms_i; idx_atom_i++){
                    for (unsigned long idx_atom_j = 0; idx_atom_j < n_atoms_j; idx_atom_j++){

                        auto atom_i = coord_covalent_radii[idx_pdb_i]->data() + (idx_atom_i * 4);
                        auto atom_j = coord_covalent_radii[idx_pdb_j]->data() + (idx_atom_j * 4);
                        double dist_i_j = distance(atom_i, atom_j);
                        double threshold = (atom_i[3] + atom_j[3]) * (2 - flexible);
                        if (dist_i_j <= threshold){
                            result->push_back(std::make_tuple(idx_pdb_i, idx_atom_i + 1, idx_pdb_j, idx_atom_j + 1));
                        }
                    }
                }
            }
        }
    }

    return result;
}


std::vector<collisions> * energy_aspect(std::vector<std::string *> &pdb_contents, float maximal_gradient_norm = 1e3, const std::string & forcefield = "UFF"){
    std::vector<std::vector<double>> atomic_grad_norm;

    for (auto content = pdb_contents.begin(); content != pdb_contents.end(); ++content){
        atomic_grad_norm.push_back(
            get_atomic_gradient_norm(**content, forcefield)
        );
    }

    auto result = new std::vector<collisions>;
    std::size_t idx_pdb = 0;

    for (auto grad = atomic_grad_norm.begin(); grad != atomic_grad_norm.end(); ++grad){
        std::size_t n_atoms = grad->size();
        auto grad_norm = grad->data();

        for (std::size_t atom_idx = 0; atom_idx < n_atoms; atom_idx++){
            if (grad_norm[atom_idx] >= maximal_gradient_norm){
                result->push_back(std::make_tuple(idx_pdb, atom_idx + 1, 0, 0));
            }
        }
    }

    return result;
}


std::vector<collisions> __pybind_export__detect(std::vector<std::string> &pdb_contents, COLLISION_SCOPE scope){
    std::vector<std::string *> pdb_contents_ptr;
    for (auto pdb_content=pdb_contents.begin(); pdb_content != pdb_contents.end(); ++pdb_content){
        pdb_contents_ptr.push_back(&(*pdb_content));
    }
    
    if (scope == COLLISION_SCOPE::inter){
        // first try energy aspect
        auto collisions = energy_aspect(pdb_contents_ptr);
        return *collisions;
    }

    std::cerr << "[_collision.so] Energy aspect diagnosis not implemented, fallback to geometry aspect." << std::endl;

    // fallback condition, try geometry aspect
    auto collisions = geometry_aspect(pdb_contents_ptr, scope);
    return *collisions;
}


#if INTERACTIVE_DEBUG
int main(int argc, char *argv[]) {

    std::ifstream pdb1;
    pdb1.open("/home/huzheyang/repos/cadd-modules/tests/assets/simple_protein/docked_complex000.pdb");

    std::stringstream strbuffer;
    strbuffer << pdb1.rdbuf();
    auto pdb1_content = strbuffer.str();

    std::vector<std::string *> buf;
    buf.push_back(&pdb1_content);
    std::vector<collisions> * result;

    std::cout << "Test on geometry aspect for simple docked structure ...";

    result = geometry_aspect(buf, COLLISION_SCOPE::inter);

    for (auto &r: *result){
        // it should not be here
        std::cout << "Wired" << std::endl;
        std::cout << "Atom " << std::get<1>(r) << " in " << std::get<0>(r) << std::endl;
        std::cout << "collision with" << std::endl;
        std::cout << "Atom " << std::get<3>(r) << " in " << std::get<2>(r) << std::endl;
    }

    std::cout << "ok" << std::endl;

    std::cout << "Test on energy aspect for simple docked structure ...";
    result = energy_aspect(buf);

    for (auto &r: *result){
        // it should not be here
        std::cout << "Wired" << std::endl;
        std::cout << "Atom " << std::get<1>(r) << " in " << std::get<0>(r) << std::endl;
        std::cout << "collision with" << std::endl;
        std::cout << "Atom " << std::get<3>(r) << " in " << std::get<2>(r) << std::endl;
    }

    std::cout << "ok" << std::endl;
};
#endif