#include "PrimaryGeneratorAction.hh"

#include <sstream>

#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4Lambda.hh"
#include "G4Neutron.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TRandom.h"

namespace {

constexpr int kDefaultSeed = 12345;
constexpr int kLambdaPerXi = 2;
// Upper bound of trials in accept-reject sampling of phase space
constexpr int kMaxTrials = 1000000;

}  // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : fSpeciesFile("../../output/species.csv"),
    fTargetsFile("../../output/targets.csv"),
    fTarget("12C"),
    fChannel("LL6He,4He,t"),
    fSeed(kDefaultSeed)
{
  fMessenger = std::make_unique<G4GenericMessenger>(
    this, "/xi/", "Xi- capture at rest");
  fMessenger->DeclareProperty("speciesFile", fSpeciesFile);
  fMessenger->DeclareProperty("targetsFile", fTargetsFile);
  fMessenger->DeclareProperty("target", fTarget, "12C, 14N or 16O");
  fMessenger->DeclareProperty("channel", fChannel,
                              "Products separated by commas");
  fMessenger->DeclareProperty("seed", fSeed, "Seed of TGenPhaseSpace");
}

void PrimaryGeneratorAction::Initialize()
{
  fTable.LoadSpecies(fSpeciesFile);
  fTable.LoadTargets(fTargetsFile);
  fInitialMass = fTable.InitialMass(fTarget);

  int aSum = 0, zSum = 0, lSum = 0;
  double mSum = 0.;
  std::stringstream ss(fChannel);
  std::string name;
  while (std::getline(ss, name, ',')) {
    const auto& s = fTable.Species(name);
    fProducts.push_back(s);
    aSum += s.a;
    zSum += s.z;
    lSum += s.nLambda;
    mSum += s.mass;
  }
  int aInit = 0, zInit = 0;
  fTable.InitialAZ(fTarget, aInit, zInit);
  if (aSum != aInit || zSum != zInit || lSum != kLambdaPerXi) {
    G4Exception("PrimaryGeneratorAction", "Gen001", FatalException,
                ("A, Z or S not conserved in " + fChannel).c_str());
  }
  fQValue = fInitialMass - mSum;
  if (fQValue <= 0.) {
    G4Exception("PrimaryGeneratorAction", "Gen002", FatalException,
                ("channel closed: " + fChannel).c_str());
  }
  gRandom->SetSeed(fSeed);
  G4cout << "Xi- + " << fTarget << " -> " << fChannel
         << "  Q = " << fQValue << " MeV" << G4endl;
  fInitialized = true;
}

G4ParticleDefinition* PrimaryGeneratorAction::Definition(
  const SpeciesInfo& s) const
{
  if (s.a == 1) {
    if (s.nLambda == 1) return G4Lambda::Definition();
    if (s.z == 1) return G4Proton::Definition();
    return G4Neutron::Definition();
  }
  auto* table = G4IonTable::GetIonTable();
  if (s.nLambda == 0) return table->GetIon(s.z, s.a);
  return table->GetIon(s.z, s.a, s.nLambda, 0);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if (!fInitialized) Initialize();

  // TGenPhaseSpace works in GeV
  TLorentzVector initial(0., 0., 0., fInitialMass / GeV);
  std::vector<double> masses;
  for (const auto& s : fProducts) masses.push_back(s.mass / GeV);
  TGenPhaseSpace phaseSpace;
  phaseSpace.SetDecay(initial, masses.size(), masses.data());

  // Accept-reject to obtain unweighted events
  const double wtMax = phaseSpace.GetWtMax();
  int trial = 0;
  while (phaseSpace.Generate() < gRandom->Uniform(wtMax)) {
    if (++trial > kMaxTrials) {
      G4Exception("PrimaryGeneratorAction", "Gen003", FatalException,
                  "phase space sampling failed");
    }
  }

  auto* vertex = new G4PrimaryVertex(G4ThreeVector(), 0.);
  for (std::size_t i = 0; i < fProducts.size(); ++i) {
    const auto* p4 = phaseSpace.GetDecay(i);
    auto* particle = new G4PrimaryParticle(Definition(fProducts[i]));
    // Mass of our table must be set before the momentum
    particle->SetMass(fProducts[i].mass * MeV);
    particle->SetMomentum(p4->Px() * GeV, p4->Py() * GeV, p4->Pz() * GeV);
    vertex->SetPrimary(particle);
  }
  event->AddPrimaryVertex(vertex);
}
