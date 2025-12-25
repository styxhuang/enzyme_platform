#pragma once
#include <vector>
#include <string>

enum COLLISION_SCOPE {inter, intra, all};

typedef std::tuple<std::size_t, unsigned long, std::size_t, unsigned long> collisions;

std::vector<collisions> __pybind_export__detect(std::vector<std::string> &pdb_contents, COLLISION_SCOPE scope);
