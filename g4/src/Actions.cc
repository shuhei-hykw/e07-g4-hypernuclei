#include "Actions.hh"

#include "EventAction.hh"
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
}

void ActionInitialization::Build() const
{
  auto* eventAction = new EventAction;
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(eventAction);
  SetUserAction(new RunAction(eventAction));
  SetUserAction(new SteppingAction(eventAction));
}
