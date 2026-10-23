#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "PhaseSpace.hh"
#include "TRandom.h"

namespace {

constexpr int kDefaultSeed = 12345;
constexpr int kLambdaPerXi = 2;

}  // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : fSpeciesFile("../../output/species.csv"),
    fTargetsFile("../../output/targets.csv"),
    fDecaysFile("../../output/decays.csv"),
    fTarget("12C"),
    fChannel("LL6He,4He,t"),
    fSeed(kDefaultSeed)
{
  fMessenger = std::make_unique<G4GenericMessenger>(
    this, "/xi/", "Xi- capture at rest");
  fMessenger->DeclareProperty("speciesFile", fSpeciesFile);
  fMessenger->DeclareProperty("targetsFile", fTargetsFile);
  fMessenger->DeclareProperty("decaysFile", fDecaysFile);
  fMessenger->DeclareProperty("target", fTarget, "12C, 14N or 16O");
  fMessenger->DeclareProperty("channel", fChannel,
                              "Products separated by commas");
  fMessenger->DeclareProperty("seed", fSeed, "Seed of TGenPhaseSpace");
}

void PrimaryGeneratorAction::Initialize()
{
  auto& table = SpeciesTable::Instance();
  table.LoadSpecies(fSpeciesFile);
  table.LoadTargets(fTargetsFile);
  table.LoadDecays(fDecaysFile);
  fInitialMass = table.InitialMass(fTarget);
  fProducts = table.Parse(fChannel);

  int aSum = 0, zSum = 0, lSum = 0;
  double mSum = 0.;
  for (const auto& s : fProducts) {
    aSum += s.a;
    zSum += s.z;
    lSum += s.nLambda;
    mSum += s.mass;
  }
  int aInit = 0, zInit = 0;
  table.InitialAZ(fTarget, aInit, zInit);
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

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if (!fInitialized) Initialize();

  const auto momenta = GeneratePhaseSpace(fInitialMass, fProducts);
  auto* vertex = new G4PrimaryVertex(G4ThreeVector(), 0.);
  for (std::size_t i = 0; i < fProducts.size(); ++i) {
    auto* particle = new G4PrimaryParticle(
      SpeciesTable::Definition(fProducts[i]));
    // Mass of our table must be set before the momentum
    particle->SetMass(fProducts[i].mass * MeV);
    particle->SetMomentum(momenta[i].x() * MeV, momenta[i].y() * MeV,
                          momenta[i].z() * MeV);
    vertex->SetPrimary(particle);
  }
  event->AddPrimaryVertex(vertex);
}
