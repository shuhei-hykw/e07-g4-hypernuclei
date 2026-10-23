#include "PhaseSpace.hh"

#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TRandom.h"

namespace {

// Upper bound of trials in accept-reject sampling of phase space
constexpr int kMaxTrials = 1000000;

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
  const double wtMax = phaseSpace.GetWtMax();
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
