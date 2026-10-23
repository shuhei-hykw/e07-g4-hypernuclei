#pragma once

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "G4ThreeVector.hh"

class G4ParticleDefinition;

// Masses of nuclei / hypernuclei, initial states and weak-decay channels,
// read from the CSV files produced by the scripts in module/.
struct SpeciesInfo {
  std::string name;
  int a = 0;
  int z = 0;
  int nLambda = 0;
  double mass = 0.;  // MeV
};

// Products with fixed momenta (MeV/c) reproducing an observed vertex
struct FixedVertex {
  std::vector<SpeciesInfo> products;
  std::vector<G4ThreeVector> momenta;
};

struct DecayChannel {
  std::string mode;  // "pi-", "pi0" or "nm"
  std::string label; // comma separated products
  std::vector<SpeciesInfo> products;
  double q = 0.;     // MeV
};

class SpeciesTable {
public:
  static SpeciesTable& Instance();

  void LoadSpecies(const std::string& path);
  void LoadTargets(const std::string& path);
  void LoadDecays(const std::string& path);
  // vertex,name,px_mev,py_mev,pz_mev (vertex: capture or parent name)
  void LoadReplay(const std::string& path);
  // nullptr unless a replay file defines the vertex
  const FixedVertex* Replay(const std::string& vertex) const;

  const SpeciesInfo& Species(const std::string& name) const;
  // nullptr if (z, a, nLambda) is not in the table
  const SpeciesInfo* Find(int z, int a, int nLambda) const;
  // Products of a comma separated list, e.g. "LL6He,4He,t"
  std::vector<SpeciesInfo> Parse(const std::string& list) const;
  const std::vector<DecayChannel>& Decays(const std::string& parent) const;

  // Initial mass M(target) + M(Xi-) - B_Xi in MeV
  double InitialMass(const std::string& target) const;
  // Baryon number and charge of the Xi- + target system
  void InitialAZ(const std::string& target, int& a, int& z) const;

  static G4ParticleDefinition* Definition(const SpeciesInfo& s);

private:
  SpeciesTable();

  std::map<std::string, SpeciesInfo> fSpecies;
  std::map<std::tuple<int, int, int>, std::string> fByZAL;
  std::map<std::string, SpeciesInfo> fTargets;
  std::map<std::string, std::vector<DecayChannel>> fDecays;
  std::map<std::string, FixedVertex> fReplay;
};
