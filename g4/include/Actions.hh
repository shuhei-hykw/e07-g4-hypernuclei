#pragma once

#include <memory>
#include <string>

#include "G4GenericMessenger.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"

class EventAction;
class TFile;
class TTree;

class RunAction : public G4UserRunAction {
public:
  explicit RunAction(EventAction* eventAction);
  void BeginOfRunAction(const G4Run* run) override;
  void EndOfRunAction(const G4Run* run) override;

private:
  void WriteGeometry() const;

  EventAction* fEventAction;
  std::string fFileName = "emulsim.root";
  TFile* fFile = nullptr;
  TTree* fTree = nullptr;
  std::unique_ptr<G4GenericMessenger> fMessenger;
};

class SteppingAction : public G4UserSteppingAction {
public:
  explicit SteppingAction(EventAction* eventAction)
    : fEventAction(eventAction) {}
  void UserSteppingAction(const G4Step* step) override;

private:
  EventAction* fEventAction;
};

class ActionInitialization : public G4VUserActionInitialization {
public:
  void Build() const override;
};
