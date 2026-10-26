#include "PrimaryGeneratorAction.hh"

#include <fstream>
#include <sstream>

#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "PhaseSpace.hh"
#include "Randomize.hh"
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
    fSeed(kDefaultSeed),
    fVertexMode("random")
{
  fMessenger = std::make_unique<G4GenericMessenger>(
    this, "/xi/", "Xi- capture at rest");
  fMessenger->DeclareProperty("speciesFile", fSpeciesFile);
  fMessenger->DeclareProperty("targetsFile", fTargetsFile);
  fMessenger->DeclareProperty("decaysFile", fDecaysFile);
  fMessenger->DeclareProperty(
    "replayFile", fReplayFile,
    "Fixed momenta of the capture and decays (overrides /xi/channel)");
  fMessenger->DeclareProperty(
    "gunFile", fGunFile,
    "Single tracks per event (overrides the capture generation)");
  fMessenger->DeclareProperty("target", fTarget, "12C, 14N or 16O");
  fMessenger->DeclareProperty("channel", fChannel,
                              "Products separated by commas");
  fMessenger->DeclareProperty("seed", fSeed, "Seed of TGenPhaseSpace");
  auto& mode = fMessenger->DeclareProperty(
    "vertexMode", fVertexMode,
    "random: uniform in thick emulsion layers, fixed: /xi/vertex");
  mode.SetCandidates("random fixed");
  fMessenger->DeclarePropertyWithUnit("vertex", "um", fVertex,
                                      "Capture point for the fixed mode");
}

void PrimaryGeneratorAction::LoadGun()
{
  std::ifstream in(fGunFile);
  if (!in) {
    G4Exception("PrimaryGeneratorAction", "Gen005", FatalException,
                ("cannot open " + fGunFile).c_str());
  }
  std::string line;
  std::getline(in, line);  // header
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::vector<std::string> c;
    for (std::string f; std::getline(ss, f, ',');) c.push_back(f);
    constexpr std::size_t kColumns = 8;
    if (c.size() < kColumns) continue;
    const int ev = std::stoi(c[0]);
    const G4ThreeVector pos(std::stod(c[2]) * um, std::stod(c[3]) * um,
                            std::stod(c[4]) * um);
    if (c[1] == "center") {
      fGunCenter[ev] = pos;
      continue;
    }
    fGun[ev].push_back({c[1], pos,
                        {std::stod(c[5]), std::stod(c[6]),
                         std::stod(c[7])}});
  }
  G4cout << "Gun: " << fGun.size() << " events from " << fGunFile
         << G4endl;
}

void PrimaryGeneratorAction::GenerateGun(G4Event* event)
{
  const int ev = event->GetEventID();
  fVertex = fGunCenter.count(ev) ? fGunCenter[ev] : G4ThreeVector();
  auto& table = SpeciesTable::Instance();
  for (const auto& t : fGun[ev]) {
    const auto info = table.Parse(t.name).at(0);
    auto* vertex = new G4PrimaryVertex(t.position, 0.);
    auto* particle = new G4PrimaryParticle(SpeciesTable::Definition(info));
    particle->SetMass(info.mass * MeV);
    particle->SetMomentum(t.momentum.x() * MeV, t.momentum.y() * MeV,
                          t.momentum.z() * MeV);
    vertex->SetPrimary(particle);
    event->AddPrimaryVertex(vertex);
  }
}

void PrimaryGeneratorAction::Initialize()
{
  auto& table = SpeciesTable::Instance();
  table.LoadSpecies(fSpeciesFile);
  table.LoadTargets(fTargetsFile);
  table.LoadDecays(fDecaysFile);
  if (!fGunFile.empty()) {
    LoadGun();
    fChannel = "gun";
    G4Random::setTheSeed(fSeed);
    fInitialized = true;
    return;
  }
  fInitialMass = table.InitialMass(fTarget);
  if (!fReplayFile.empty()) {
    table.LoadReplay(fReplayFile);
    fReplay = table.Replay("capture");
    if (!fReplay) {
      G4Exception("PrimaryGeneratorAction", "Gen004", FatalException,
                  "replay file has no capture vertex");
    }
    fChannel.clear();
    for (const auto& s : fReplay->products) {
      fChannel += (fChannel.empty() ? "" : ",") + s.name;
    }
  }
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
  // Capture points are sampled with the Geant4 engine: seed it as well
  G4Random::setTheSeed(fSeed);
  G4cout << "Xi- + " << fTarget << " -> " << fChannel
         << "  Q = " << fQValue << " MeV" << G4endl;
  fInitialized = true;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if (!fInitialized) Initialize();
  if (!fGunFile.empty()) {
    GenerateGun(event);
    return;
  }

  if (fVertexMode == "random") {
    const auto* detector = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fVertex = detector->SampleVertex();
  }
  const auto momenta = fReplay
    ? fReplay->momenta : GeneratePhaseSpace(fInitialMass, fProducts);
  auto* vertex = new G4PrimaryVertex(fVertex, 0.);
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
