#include <exception>
#include <map>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obconversion.h>
#include <openbabel/elements.h>
#include <openbabel/oberror.h>
#include <openbabel/forcefield.h>
#include <openbabel/obiter.h>

#include "_ligand_restore.hpp"

#if defined(__GNUC__) || defined(__GNUG__)
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define unlikely(x) x
#endif


std::vector<std::string> __pybind_export__pdbqt_to_pdb(const std::string &original_pdb, const std::string &original_pdbqt, const std::vector<std::string> &docked_pdbqt){
    // read files
    OpenBabel::OBMol origin_full_mol, origin_reduced_mol;
    OpenBabel::OBConversion obconv;
    obconv.SetInAndOutFormats("pdb", "pdb");
    if (!obconv.ReadString(&origin_full_mol, original_pdb)){
        throw std::runtime_error("Cannot parse original PDB file");
    };
    obconv.SetInAndOutFormats("pdbqt", "pdb");
    if (!obconv.ReadString(&origin_reduced_mol, original_pdbqt)){
        throw std::runtime_error("Cannot parse original PDBQT file");
    };

    std::vector<OpenBabel::OBMol*> mols_to_patch;
    for (auto docked_content: docked_pdbqt){
        auto docked_mol = new OpenBabel::OBMol;
        if (!obconv.ReadString(docked_mol, docked_content)){
            throw std::runtime_error("Cannot parse docked PDBQT file");
        }
        mols_to_patch.push_back(docked_mol);
    }

    auto opt_reuslts = pdbqt_to_pdb(origin_full_mol, origin_reduced_mol, mols_to_patch);
    
    std::vector<std::string> result;
    for (auto opt_result: opt_reuslts){
        if (opt_result == nullptr){
            result.push_back("");
            continue;
        };
        auto final_opt_result = obconv.WriteString(opt_result);
        result.push_back(final_opt_result);
    }

    for (auto dirty_ptr: mols_to_patch){
        if (dirty_ptr != nullptr){
            delete dirty_ptr;
        }
    };
    for (auto dirty_ptr: opt_reuslts){
        if (dirty_ptr != nullptr){
            delete dirty_ptr;
        }
    };

    return result;
};

std::vector<OpenBabel::OBMol*> pdbqt_to_pdb(OpenBabel::OBMol &original_full, OpenBabel::OBMol &original_reduced, std::vector<OpenBabel::OBMol*> &mols_to_patch){
    // build up atom map from pdbqt to pdb
    std::map<unsigned long, unsigned long> atom_mapper;
    for (auto fatom = OpenBabel::OBMolAtomIter(original_full); fatom; fatom++){
        bool map_resolved = false;
        for (auto ratom = OpenBabel::OBMolAtomIter(original_reduced); ratom; ratom++){
            if ((*fatom).GetVector().IsApprox((*ratom).GetVector(), .01)){
                atom_mapper.emplace((*ratom).GetId(), (*fatom).GetId());
                map_resolved = true;
                break;
            };
        };
        if (!map_resolved){
            if ((*fatom).GetAtomicNum() != 1){
                throw std::runtime_error("There is at least one heavy atom missing in the PDBQT");
            };
        };
    };

    std::vector<OpenBabel::OBMol*> result;
    for (auto mol_to_patch: mols_to_patch){
        auto dummy_mol = OpenBabel::OBMol(original_full);
        bool ready_to_optimize = true;

        for (auto atom_to_patch = OpenBabel::OBMolAtomIter(mol_to_patch); atom_to_patch; atom_to_patch++)
        {
            auto id_atom_to_patch = atom_to_patch->GetId();
            auto vec_atom_to_patch = atom_to_patch->GetVector();

            // atom to patch might not exist, say the flexible docking
            // atoms in AA will included in docking result
            if (auto search = atom_mapper.find(id_atom_to_patch); search == atom_mapper.end())
            {
                // the atom is likely to be in the AA
                // ignore such atoms is fine so far
                ;
            }
            else
            {
                // the atom is found in dummy, resettle that atom
                // would be fine
                auto dummy_atom = dummy_mol.GetAtomById(atom_mapper.at(id_atom_to_patch));
                if (unlikely(dummy_atom == nullptr))
                {
                ready_to_optimize = false;
                break;
            }
            dummy_atom->SetVector((*atom_to_patch).GetVector());
        };
        };
        if (!ready_to_optimize){
            result.push_back(nullptr);
            continue;
        };

        auto opt_mol = hydrogen_optimize(dummy_mol, 1000, 1e-5, "UFF");
        result.push_back(opt_mol);
    };
    return result;
};


#if INTERACTIVE_DEBUG
int main(){
    std::ifstream pdb_istream;
    std::stringstream strbuffer;

    pdb_istream.open("/home/huzheyang/下载/1698650085798_1aadd/intermediate_files/partition_0000000/phes_1699328454_Untitled-3/00000.pdb");
    strbuffer.str("");
    strbuffer << pdb_istream.rdbuf();
    pdb_istream.close();
    auto pdb_content = strbuffer.str();

    pdb_istream.open("/home/huzheyang/下载/1698650085798_1aadd/intermediate_files/partition_0000000/phes_1699328454_Untitled-3/000000002.pdbqt");
    strbuffer.str("");
    strbuffer << pdb_istream.rdbuf();
    pdb_istream.close();
    auto pdbqt_content = strbuffer.str();

    pdb_istream.open("/home/huzheyang/下载/1698650085798_1aadd/intermediate_files/partition_0000000/phes_1699328454_Untitled-3/ligand_out_model6.pdbqt");
    strbuffer.str("");
    strbuffer << pdb_istream.rdbuf();
    pdb_istream.close();
    auto patching = strbuffer.str();


    __pybind_export__pdbqt_to_pdb(pdb_content, pdbqt_content, {patching});

    return 0;
}
#endif