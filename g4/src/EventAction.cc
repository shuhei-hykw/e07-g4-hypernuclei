#include "EventAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorAction.hh"
#include "TTree.h"

namespace {

// Lighter particles (e+-, gamma, neutrinos) are not recorded
constexpr double kMinRecordedMass = 100. * MeV;

}  // namespace

void EventAction::Book(TTree* tree)
{
  fTree = tree;
  tree->Branch("event", &fEventId);
  tree->Branch("channel", &fChannel);
  tree->Branch("q", &fQValue);
  tree->Branch("id", &fTrackId);
  tree->Branch("parent", &fParentId);
  tree->Branch("process", &fProcess);
  tree->Branch("name", &fName);
  tree->Branch("z", &fZ);
  tree->Branch("a", &fA);
  tree->Branch("nl", &fNLambda);
  tree->Branch("ke", &fKE);
  tree->Branch("px", &fPx);
  tree->Branch("py", &fPy);
  tree->Branch("pz", &fPz);
  tree->Branch("range", &fRange);
  tree->Branch("length", &fLength);
  tree->Branch("x", &fX);
  tree->Branch("y", &fY);
  tree->Branch("zpos", &fZpos);
  tree->Branch("edep", &fEdep);
}

void EventAction::BeginOfEventAction(const G4Event* event)
{
  fEventId = event->GetEventID();
  fIndex.clear();
  for (auto* v : {&fName, &fProcess}) v->clear();
  for (auto* v : {&fTrackId, &fParentId, &fZ, &fA, &fNLambda}) v->clear();
  for (auto* v : {&fKE, &fPx, &fPy, &fPz, &fRange, &fLength}) v->clear();
  for (auto* v : {&fX, &fY, &fZpos, &fEdep}) v->clear();
}

void EventAction::RecordStep(const G4Step* step)
{
  const auto* track = step->GetTrack();
  const auto* def = track->GetDefinition();
  if (def->GetPDGMass() < kMinRecordedMass) return;  // e+-, gamma, nu

  const int id = track->GetTrackID();
  auto it = fIndex.find(id);
  if (it == fIndex.end()) {
    const auto* pre = step->GetPreStepPoint();
    const auto* creator = track->GetCreatorProcess();
    const auto mom = pre->GetMomentum();
    fIndex[id] = fName.size();
    fTrackId.push_back(id);
    fParentId.push_back(track->GetParentID());
    fProcess.push_back(creator ? creator->GetProcessName() : "primary");
    fName.push_back(def->GetParticleName());
    fZ.push_back(static_cast<int>(def->GetPDGCharge() / eplus));
    fA.push_back(def->GetBaryonNumber());
    fNLambda.push_back(def->GetNumberOfLambdasInHypernucleus());
    fKE.push_back(pre->GetKineticEnergy() / MeV);
    fPx.push_back(mom.x() / MeV);
    fPy.push_back(mom.y() / MeV);
    fPz.push_back(mom.z() / MeV);
    fRange.push_back(0.);
    fLength.push_back(0.);
    const auto pos = pre->GetPosition();
    fX.push_back({float(pos.x() / um)});
    fY.push_back({float(pos.y() / um)});
    fZpos.push_back({float(pos.z() / um)});
    fEdep.push_back({0.f});
    it = fIndex.find(id);
  }
  const std::size_t i = it->second;
  const auto pos = step->GetPostStepPoint()->GetPosition();
  const auto start = track->GetVertexPosition();
  fX[i].push_back(pos.x() / um);
  fY[i].push_back(pos.y() / um);
  fZpos[i].push_back(pos.z() / um);
  fEdep[i].push_back(step->GetTotalEnergyDeposit() / MeV);
  fLength[i] += step->GetStepLength() / um;
  fRange[i] = (pos - start).mag() / um;
}

void EventAction::EndOfEventAction(const G4Event*)
{
  const auto* gen = static_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  fChannel = gen->Target() + " -> " + gen->Channel();
  fQValue = gen->QValue();
  if (fTree) fTree->Fill();
}
