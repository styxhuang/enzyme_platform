#include <Eigen/Core>
#include <LBFGS.h>
#include <vector>
#include <random>
#include <exception>
#include <cmath>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obconversion.h>
#include <openbabel/elements.h>
#include <openbabel/oberror.h>
#include <openbabel/forcefield.h>
#include <openbabel/obiter.h>

#include "_hydrogen_optimize.hpp"

#if DEBUG
#define PRINT_MORE true
#else
#define PRINT_MORE false
#endif

double rand01()
{
    return rand() / (RAND_MAX + 1.);
}

inline void queue_hydrogen_atoms(const OpenBabel::OBMol &mol, std::vector<unsigned int> &consider_atoms)
{
    auto n_atoms = mol.NumAtoms();
    for (unsigned int atom_idx = 0; atom_idx < n_atoms; atom_idx++)
    {
        switch (mol.GetAtom(atom_idx + 1)->GetAtomicNum())
        {
        case 1:
            consider_atoms.push_back(atom_idx);
            #ifdef DEBUG
            std::cout << "Considering atom at Array idx (" << atom_idx << ")" << std::endl;
            #endif
            break;
        
        default:
            break;
        }
    }
};

void radius_rand1(double &x, double &y, double &z)
{
    x = rand01() - .5;
    y = rand01() - .5;
    z = rand01() - .5;

    while (x == 0 && y == 0 && z == 0)
    {
        x = rand01() - .5;
        y = rand01() - .5;
        z = rand01() - .5;
    }

    auto rho = sqrt(x * x + y * y + z * z);
    x /= rho;
    y /= rho;
    z /= rho; 
}

void hydrogen_placement(const std::vector<unsigned int>& consider_atoms, OpenBabel::OBMol& mol)
{
    // set hydrogen near its anchor closely
    for (auto consider_atom: consider_atoms)
    {
        auto hydrogen_atom = mol.GetAtom(consider_atom + 1);
        OpenBabel::OBAtomAtomIter anchor_atom(hydrogen_atom);
        auto anchor_pos = (*anchor_atom).GetVector();
        #ifdef DEBUG
        std::cout << "placing atom at Array idx (" << consider_atom << ") from ";
        std::cout << (*hydrogen_atom).GetX() << ", " << (*hydrogen_atom).GetY() << ", "<< (*hydrogen_atom).GetZ();
        #endif
        auto ax = anchor_pos.GetX();
        auto ay = anchor_pos.GetY();
        auto az = anchor_pos.GetZ();
        radius_rand1(ax, ay, az);
        OpenBabel::vector3 hydrogen_pos = {ax, ay, az};
        hydrogen_atom->SetVector(hydrogen_pos);
        #ifdef DEBUG
        std::cout << " to ";
        std::cout << (*hydrogen_atom).GetX() << ", " << (*hydrogen_atom).GetY() << ", "<< (*hydrogen_atom).GetZ();
        std::cout << std::endl;
        #endif
    }
}

class SinglePoint{
    private:
        OpenBabel::OBForceField * ff;
        std::vector<unsigned int> consider_atoms;

        unsigned int n_atoms;
    public:
        OpenBabel::OBMol * mol;

        SinglePoint(OpenBabel::OBForceField * const ff, OpenBabel::OBMol * const mol, const std::vector<unsigned int> &consider_atoms);
        double operator()(const Eigen::VectorXd &coords, Eigen::VectorXd &grad);
};

SinglePoint::SinglePoint(OpenBabel::OBForceField * const ff, OpenBabel::OBMol * const mol, const std::vector<unsigned int> &consider_atoms)
{
    this->ff = ff;
    this->mol = new OpenBabel::OBMol(*mol);
    this->consider_atoms.assign(consider_atoms.begin(), consider_atoms.end());

    this->n_atoms = mol->NumAtoms();
};

double SinglePoint::operator()(const Eigen::VectorXd &coords, Eigen::VectorXd &grad){
    // evalute energy and update gradient


    if (this->n_atoms * 3 != coords.size()){
        throw "Cols of coordinations not match atom numbers";
    }

    if (this->n_atoms * 3 != grad.size()){
        throw "Cols of gradients not match atom numbers";
    }

    // reporting buffer
    // current XYZ (3) displacement XYZ (3) gradients XYZ (3)
#ifdef DEBUG
    auto buffer = new double[this->n_atoms * 9];
    memset(buffer, 0, sizeof(double) * this->n_atoms * 9);
#endif

    // update coordinate
    for (unsigned int atom_idx = 0; atom_idx < this->n_atoms; atom_idx++){
        auto atom = this->mol->GetAtom(atom_idx + 1);
        auto atom_vec_original = atom->GetVector();

        atom->SetVector(coords(atom_idx * 3 + 0), coords(atom_idx * 3 + 1),coords(atom_idx * 3 + 2));
#ifdef DEBUG
        buffer[atom_idx * 9 + 0] = coords(atom_idx * 3 + 0);
        buffer[atom_idx * 9 + 1] = coords(atom_idx * 3 + 1);
        buffer[atom_idx * 9 + 2] = coords(atom_idx * 3 + 2);

        buffer[atom_idx * 9 + 3] = coords(atom_idx * 3 + 0) - atom_vec_original.x();
        buffer[atom_idx * 9 + 4] = coords(atom_idx * 3 + 1) - atom_vec_original.y();
        buffer[atom_idx * 9 + 5] = coords(atom_idx * 3 + 2) - atom_vec_original.z();
#endif

    }
    auto update_complete = this->ff->SetConformers(*this->mol);
    if (!update_complete){
        throw "Update coordinate failed";
    }

    double energy = this->ff->Energy();

    // update gradients
    grad.setZero(grad.size());
    double pos_epslion = 1.0;
    for (auto &atom_idx : this->consider_atoms){
        auto atom = this->mol->GetAtom(atom_idx + 1);
        auto atom_grads = this->ff->GetGradient(atom);

        // gradient might be nan, look at the first element might be enough
        if (!std::isfinite(atom_grads.x())){
            // calculate 1st difference of bond energy
            #if DEBUG
            std::cerr << "Start to numerical 1st difference" << std::endl;
            #endif

            auto origin_pos = atom->GetVector();
            auto energy_bond = this->ff->E_Bond(false);
            
            // dx
            atom->SetVector(origin_pos.x() + pos_epslion, origin_pos.y(), origin_pos.z());
            this->ff->SetConformers(*this->mol);
            auto dEdx = (this->ff->E_Bond(false) - energy_bond) / pos_epslion;
            // dy
            atom->SetVector(origin_pos.x(), origin_pos.y() + pos_epslion, origin_pos.z());
            this->ff->SetConformers(*this->mol);
            auto dEdy = (this->ff->E_Bond(false) - energy_bond) / pos_epslion;
            // dz
            atom->SetVector(origin_pos.x(), origin_pos.y(), origin_pos.z() + pos_epslion);
            this->ff->SetConformers(*this->mol);
            auto dEdz = (this->ff->E_Bond(false) - energy_bond) / pos_epslion;

            auto gradx = -dEdx + rand01();
            auto grady = -dEdy + rand01();
            auto gradz = -dEdz + rand01();

            grad(atom_idx * 3 + 0) = gradx;
            grad(atom_idx * 3 + 1) = grady;
            grad(atom_idx * 3 + 2) = gradz;

#ifdef DEBUG
            buffer[atom_idx * 9 + 6] = gradx;
            buffer[atom_idx * 9 + 7] = grady;
            buffer[atom_idx * 9 + 8] = gradz;
#endif
            
        } else {
            // analytical gradient is ok
            auto gradx = atom_grads.x();
            auto grady = atom_grads.y();
            auto gradz = atom_grads.z();

            grad(atom_idx * 3 + 0) = gradx;
            grad(atom_idx * 3 + 1) = grady;
            grad(atom_idx * 3 + 2) = gradz;

#ifdef DEBUG
            buffer[atom_idx * 9 + 6] = gradx;
            buffer[atom_idx * 9 + 7] = grady;
            buffer[atom_idx * 9 + 8] = gradz;
#endif
        }
    }

    // print diagnosis information
#ifdef DEBUG
    std::cerr << "Energy: " << energy << std::endl;
    std::cerr << "Coordinate|Displacement|Gradients" << std::endl;
    for (unsigned int idx = 0; idx < this->n_atoms; idx++){
        std::cerr << "(" << buffer[idx * 9 + 0] << ", " << buffer[idx * 9 + 1] << ", " << buffer[idx * 9 + 2] << ") ";
        std::cerr << "(" << buffer[idx * 9 + 3] << ", " << buffer[idx * 9 + 4] << ", " << buffer[idx * 9 + 5] << ") ";
        std::cerr << "(" << buffer[idx * 9 + 6] << ", " << buffer[idx * 9 + 7] << ", " << buffer[idx * 9 + 8] << ") ";
        std::cerr << std::endl;
        }

    delete buffer;
#endif

    return energy;
};

std::string __pybind_export__hydrogen_optimize(const std::string &pdb_content, int steps = 5000, double e_conv = 1e-4, const std::string & forcefield = "UFF"){
    OpenBabel::OBMol mol;
    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("pdb", "pdb");
    conv.ReadString(&mol, pdb_content);

    std::string result = conv.WriteString(hydrogen_optimize(mol, steps, e_conv, forcefield));

    return result;
};

OpenBabel::OBMol* hydrogen_optimize(OpenBabel::OBMol mol, int steps, double e_conv, const std::string & forcefield = "UFF")
{
    std::vector<unsigned int> consider_atoms;
    queue_hydrogen_atoms(mol, consider_atoms);

    return hydrogen_optimize(mol, steps, e_conv, forcefield, consider_atoms);
};

OpenBabel::OBMol* hydrogen_optimize(OpenBabel::OBMol mol, int steps, double e_conv, const std::string & forcefield, const std::vector<unsigned int> & consider_atoms)
{
    unsigned int n_atoms = mol.NumAtoms();

    // Assign forcefield
    auto ff = OpenBabel::OBForceField::FindForceField(forcefield);
    bool ff_ready = ff->Setup(mol);

    if (!ff_ready){
        std::cerr << "Fail to setup forcefield" << std::endl;
        return nullptr;
    }

    #ifdef DEBUG
    // inspect Id of Atom vs. index from 0 to num_atom
    for (unsigned int atom_idx = 0; atom_idx<n_atoms; atom_idx++){
        auto atom = mol.GetAtom(atom_idx + 1);
        auto atomid = atom->GetId();
        auto atomidx = atom->GetIdx();
        std::cout << "Array idx: " << atom_idx << " -- Atom ID: " << atomid << " -- Atom Idx: " << atomidx << std::endl;
    }
    #endif

    // set hydrogen near its anchor closely
    for (auto consider_atom: consider_atoms){
        auto hydrogen_atom = mol.GetAtom(consider_atom + 1);
        OpenBabel::OBAtomAtomIter anchor_atom(hydrogen_atom);
        auto anchor_pos = (*anchor_atom).GetVector();
        #ifdef DEBUG
        std::cout << "placing atom at Array idx (" << consider_atom << ") from ";
        std::cout << (*hydrogen_atom).GetX() << ", " << (*hydrogen_atom).GetY() << ", "<< (*hydrogen_atom).GetZ();
        #endif
        anchor_pos.SetX(anchor_pos.GetX() + rand01());
        anchor_pos.SetY(anchor_pos.GetY() + rand01());
        anchor_pos.SetZ(anchor_pos.GetZ() + rand01());
        hydrogen_atom->SetVector(anchor_pos);
        #ifdef DEBUG
        std::cout << " to ";
        std::cout << (*hydrogen_atom).GetX() << ", " << (*hydrogen_atom).GetY() << ", "<< (*hydrogen_atom).GetZ();
        std::cout << std::endl;
        #endif
    }

    auto sp = SinglePoint(ff, &mol, consider_atoms);
    auto coords = Eigen::VectorXd(n_atoms * 3);
    double energy;

    // optimize configuration
    LBFGSpp::LBFGSParam<double> bfgs_param;
    bfgs_param.m = 20;
    bfgs_param.epsilon = e_conv;
    bfgs_param.max_iterations = steps;
    bfgs_param.max_step = 2.0;

    bool optimized = false;
    int ttl = 3;

    while ((!optimized) && (ttl > 0))
    {
        LBFGSpp::LBFGSSolver<double> bfgs_solver(bfgs_param);

        hydrogen_placement(consider_atoms, mol);

        for (auto atom_idx = 0; atom_idx < n_atoms; atom_idx++){
            auto atom = mol.GetAtom(atom_idx + 1);
            auto atom_coord = atom->GetVector();

            coords(atom_idx * 3 + 0) = atom_coord.x();
            coords(atom_idx * 3 + 1) = atom_coord.y();
            coords(atom_idx * 3 + 2) = atom_coord.z();
        }

        try
        {
            bfgs_solver.minimize(sp, coords, energy);
            optimized = true;
        }
        catch (std::runtime_error bfgs_err)
        {
            if (0 == strcmp(bfgs_err.what(), "the line search step became smaller than the minimum value allowed")) {
                // step is too small, be permissive for now
                // try it again
                std::cerr << "BFGS search step failure, restart calculation." << std::endl;
            }
            else
            {
                std::cerr << "BFGS failure: " << bfgs_err.what() << std::endl;
            }
        };

        ttl--;
    }

    // output final result
    return sp.mol;
}


#if INTERACTIVE_DEBUG
int main(){
    std::ifstream pdb_istream;
    pdb_istream.open("/home/huzheyang/repos/mod_vina_docking/vina_docking/test.pdb");

    std::stringstream strbuffer;
    strbuffer << pdb_istream.rdbuf();
    auto pdb_content = strbuffer.str();

    hydrogen_optimize(pdb_content, 1000, 1e-5, "UFF");

    return 0;
}
#endif