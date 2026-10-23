#include "Actions.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "TFile.h"
#include "TTree.h"

RunAction::RunAction(EventAction* eventAction) : fEventAction(eventAction)
{
  fMessenger = std::make_unique<G4GenericMessenger>(
    this, "/out/", "Output settings");
  fMessenger->DeclareProperty("file", fFileName, "Output ROOT file");
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  fFile = TFile::Open(fFileName.c_str(), "RECREATE");
  if (!fFile || fFile->IsZombie()) {
    G4Exception("RunAction", "Out001", FatalException,
                ("cannot open " + fFileName).c_str());
  }
  fTree = new TTree("tree", "Xi- capture at rest in emulsion");
  fEventAction->Book(fTree);
  WriteGeometry();
}

void RunAction::WriteGeometry() const
{
  const auto* detector = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  int sheet = 0;
  bool emulsion = false, thick = false;
  double zMin = 0., zMax = 0.;
  double halfX = detector->HalfX() / um, halfY = detector->HalfY() / um;
  TTree layers("layers", "Emulsion module layers (um)");
  layers.Branch("sheet", &sheet);
  layers.Branch("emulsion", &emulsion);
  layers.Branch("thick", &thick);
  layers.Branch("zmin", &zMin);
  layers.Branch("zmax", &zMax);
  layers.Branch("halfx", &halfX);
  layers.Branch("halfy", &halfY);
  for (const auto& l : detector->Layers()) {
    sheet = l.sheet;
    emulsion = l.emulsion;
    thick = l.thick;
    zMin = l.zMin / um;
    zMax = l.zMax / um;
    layers.Fill();
  }
  layers.Write();
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  fFile->cd();
  fTree->Write();
  fFile->Close();
  delete fFile;
  fFile = nullptr;
  G4cout << run->GetNumberOfEvent() << " events written to " << fFileName
         << G4endl;
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  fEventAction->RecordStep(step);
  // A track entering the world volume (or leaving it) has escaped from
  // the emulsion module.
  const auto* next = step->GetPostStepPoint()->GetPhysicalVolume();
  if (next == nullptr || next->GetMotherLogical() == nullptr) {
    auto* track = step->GetTrack();
    fEventAction->MarkEscaped(track->GetTrackID());
    track->SetTrackStatus(fStopAndKill);
  }
}

void ActionInitialization::Build() const
{
  auto* eventAction = new EventAction;
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(eventAction);
  SetUserAction(new RunAction(eventAction));
  SetUserAction(new SteppingAction(eventAction));
}
