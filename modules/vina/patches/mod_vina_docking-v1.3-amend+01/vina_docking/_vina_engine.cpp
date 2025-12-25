#include "_collision.hpp"
#include "_hydrogen_optimize.hpp"
#include "_ligand_restore.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_vina_engine, mod){
    mod.def("hydrogen_optimize", &__pybind_export__hydrogen_optimize, "pdb_content", "steps", "e_conv", "forcefield");
    mod.def("restore_ligand", &__pybind_export__pdbqt_to_pdb, "original_pdb", "original_pdbqt", "docked_pdbqt");

    py::enum_<COLLISION_SCOPE>(mod, "COLLISION_SCOPE")
        .value("inter", inter)
        .value("intra", intra)
        .value("all", all)
        .export_values();

    mod.def("detect", &__pybind_export__detect, "pdb_content", "scope");
}
