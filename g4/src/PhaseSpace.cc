#include "PhaseSpace.hh"

#include <algorithm>
#include <map>

#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TRandom.h"

namespace {

// Upper bound of trials in accept-reject sampling of phase space
constexpr int kMaxTrials = 1000000;
// GetWtMax() can exceed the true maximum weight by orders of magnitude
// for many-body decays with very different masses; the maximum is then
// estimated from samples (with a safety factor) and cached per decay.
constexpr int kWeightSamples = 20000;
constexpr double kWeightSafety = 1.2;

double MaxWeight(TGenPhaseSpace& phaseSpace,
                 const std::vector<double>& key)
{
  static std::map<std::vector<double>, double> cache;
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;
  double wMax = 0.;
  for (int i = 0; i < kWeightSamples; ++i) {
    wMax = std::max(wMax, phaseSpace.Generate());
  }
  wMax = std::min(wMax * kWeightSafety, phaseSpace.GetWtMax());
  cache[key] = wMax;
  return wMax;
}

}  // namespace

std::vector<G4ThreeVector> GeneratePhaseSpace(
  double parentMass, const std::vector<SpeciesInfo>& products)
{
  // TGenPhaseSpace works in GeV
  TLorentzVector initial(0., 0., 0., parentMass * MeV / GeV);
  std::vector<double> masses;
  for (const auto& s : products) masses.push_back(s.mass * MeV / GeV);
  TGenPhaseSpace phaseSpace;
  if (!phaseSpace.SetDecay(initial, masses.size(), masses.data())) {
    G4Exception("GeneratePhaseSpace", "PS001", FatalException,
                "decay is kinematically forbidden");
  }

  // Accept-reject to obtain unweighted events
  std::vector<double> key(masses);
  key.push_back(initial.M());
  const double wtMax = MaxWeight(phaseSpace, key);
  int trial = 0;
  while (phaseSpace.Generate() < gRandom->Uniform(wtMax)) {
    if (++trial > kMaxTrials) {
      G4Exception("GeneratePhaseSpace", "PS002", FatalException,
                  "phase space sampling failed");
    }
  }

  std::vector<G4ThreeVector> momenta;
  for (std::size_t i = 0; i < products.size(); ++i) {
    const auto* p4 = phaseSpace.GetDecay(i);
    momenta.emplace_back(p4->Px() * GeV / MeV, p4->Py() * GeV / MeV,
                         p4->Pz() * GeV / MeV);
  }
  return momenta;
}
