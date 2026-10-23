#include "SpeciesTable.hh"

#include <fstream>
#include <sstream>
#include <vector>

#include "G4Exception.hh"

namespace {

std::vector<std::vector<std::string>> ReadCsv(const std::string& path)
{
  std::ifstream ifs(path);
  if (!ifs) {
    G4Exception("SpeciesTable", "File001", FatalException,
                ("cannot open " + path).c_str());
  }
  std::vector<std::vector<std::string>> rows;
  std::string line;
  std::getline(ifs, line);  // header
  while (std::getline(ifs, line)) {
    if (line.empty()) continue;
    std::vector<std::string> cols;
    std::stringstream ss(line);
    std::string col;
    while (std::getline(ss, col, ',')) cols.push_back(col);
    rows.push_back(cols);
  }
  return rows;
}

}  // namespace

void SpeciesTable::LoadSpecies(const std::string& path)
{
  // name,A,Z,n_lambda,mass_mev,estimated
  for (const auto& c : ReadCsv(path)) {
    SpeciesInfo s{c[0], std::stoi(c[1]), std::stoi(c[2]),
                  std::stoi(c[3]), std::stod(c[4])};
    fSpecies[s.name] = s;
  }
}

void SpeciesTable::LoadTargets(const std::string& path)
{
  // name,A,Z,b_xi_mev,m_init_mev
  for (const auto& c : ReadCsv(path)) {
    SpeciesInfo s{c[0], std::stoi(c[1]), std::stoi(c[2]), 0,
                  std::stod(c[4])};
    fTargets[s.name] = s;
  }
}

const SpeciesInfo& SpeciesTable::Species(const std::string& name) const
{
  auto it = fSpecies.find(name);
  if (it == fSpecies.end()) {
    G4Exception("SpeciesTable", "Species001", FatalException,
                ("unknown species " + name).c_str());
  }
  return it->second;
}

double SpeciesTable::InitialMass(const std::string& target) const
{
  auto it = fTargets.find(target);
  if (it == fTargets.end()) {
    G4Exception("SpeciesTable", "Target001", FatalException,
                ("unknown target " + target).c_str());
  }
  return it->second.mass;
}

void SpeciesTable::InitialAZ(const std::string& target, int& a,
                             int& z) const
{
  InitialMass(target);  // existence check
  const auto& t = fTargets.at(target);
  a = t.a + 1;
  z = t.z - 1;
}
