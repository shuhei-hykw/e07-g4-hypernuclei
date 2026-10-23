#pragma once

#include <map>
#include <string>

// Masses of nuclei / hypernuclei and initial states, read from the CSV
// files produced by channels/enumerate_channels.py.
struct SpeciesInfo {
  std::string name;
  int a = 0;
  int z = 0;
  int nLambda = 0;
  double mass = 0.;  // MeV
};

class SpeciesTable {
public:
  void LoadSpecies(const std::string& path);
  void LoadTargets(const std::string& path);
  const SpeciesInfo& Species(const std::string& name) const;
  // Initial mass M(target) + M(Xi-) - B_Xi in MeV
  double InitialMass(const std::string& target) const;
  // Baryon number and charge of the Xi- + target system
  void InitialAZ(const std::string& target, int& a, int& z) const;

private:
  std::map<std::string, SpeciesInfo> fSpecies;
  std::map<std::string, SpeciesInfo> fTargets;
};
