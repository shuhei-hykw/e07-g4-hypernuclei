#include "SpeciesTable.hh"

#include <fstream>
#include <sstream>

#include "G4Exception.hh"
#include "G4IonTable.hh"
#include "G4Lambda.hh"
#include "G4Neutron.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"

namespace {

std::vector<std::string> Split(const std::string& s, char delim)
{
  std::vector<std::string> cols;
  std::stringstream ss(s);
  std::string col;
  while (std::getline(ss, col, delim)) cols.push_back(col);
  return cols;
}

// Minimal CSV reader handling double-quoted fields
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
    std::vector<std::string> cols(1);
    bool quoted = false;
    for (char c : line) {
      if (c == '"') {
        quoted = !quoted;
      } else if (c == ',' && !quoted) {
        cols.emplace_back();
      } else {
        cols.back() += c;
      }
    }
    rows.push_back(cols);
  }
  return rows;
}

}  // namespace

SpeciesTable& SpeciesTable::Instance()
{
  static SpeciesTable instance;
  return instance;
}

SpeciesTable::SpeciesTable()
{
  // Pions appear only as weak-decay products
  fSpecies["pi-"] = {"pi-", 0, -1, 0,
                     G4PionMinus::Definition()->GetPDGMass() / MeV};
  fSpecies["pi0"] = {"pi0", 0, 0, 0,
                     G4PionZero::Definition()->GetPDGMass() / MeV};
}

void SpeciesTable::LoadSpecies(const std::string& path)
{
  // name,A,Z,n_lambda,mass_mev,estimated
  for (const auto& c : ReadCsv(path)) {
    SpeciesInfo s{c[0], std::stoi(c[1]), std::stoi(c[2]),
                  std::stoi(c[3]), std::stod(c[4])};
    fSpecies[s.name] = s;
    fByZAL[{s.z, s.a, s.nLambda}] = s.name;
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

void SpeciesTable::LoadDecays(const std::string& path)
{
  // parent,mode,channel,q_mev,n_bodies,n_charged,uses_estimate
  for (const auto& c : ReadCsv(path)) {
    fDecays[c[0]].push_back({c[1], c[2], Parse(c[2]), std::stod(c[3])});
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

const SpeciesInfo* SpeciesTable::Find(int z, int a, int nLambda) const
{
  auto it = fByZAL.find({z, a, nLambda});
  return it == fByZAL.end() ? nullptr : &fSpecies.at(it->second);
}

std::vector<SpeciesInfo> SpeciesTable::Parse(const std::string& list) const
{
  std::vector<SpeciesInfo> products;
  for (const auto& name : Split(list, ',')) {
    products.push_back(Species(name));
  }
  return products;
}

const std::vector<DecayChannel>& SpeciesTable::Decays(
  const std::string& parent) const
{
  static const std::vector<DecayChannel> empty;
  auto it = fDecays.find(parent);
  return it == fDecays.end() ? empty : it->second;
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

G4ParticleDefinition* SpeciesTable::Definition(const SpeciesInfo& s)
{
  if (s.a == 0) {
    if (s.z == -1) return G4PionMinus::Definition();
    return G4PionZero::Definition();
  }
  if (s.a == 1) {
    if (s.nLambda == 1) return G4Lambda::Definition();
    if (s.z == 1) return G4Proton::Definition();
    return G4Neutron::Definition();
  }
  auto* table = G4IonTable::GetIonTable();
  if (s.nLambda == 0) return table->GetIon(s.z, s.a);
  return table->GetIon(s.z, s.a, s.nLambda, 0);
}
