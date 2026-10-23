// Track simulation of Xi- capture products in nuclear emulsion.
//   ./emulsim run.mac

#include "Actions.hh"
#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4PrimaryTransformer.hh"
#include "G4RunManagerFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "HypernucleusDecay.hh"

namespace {

constexpr double kProductionCut = 1. * um;
// Our hypernuclear masses differ from those of G4HyperNucleiProperties
constexpr double kMassTolerance = 100. * MeV;

}  // namespace

int main(int argc, char** argv)
{
  G4UIExecutive* ui = (argc == 1) ? new G4UIExecutive(argc, argv) : nullptr;

  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetUserInitialization(new DetectorConstruction);

  auto* physics = new FTFP_BERT;
  physics->ReplacePhysics(new G4EmStandardPhysics_option4);
  physics->RegisterPhysics(new G4StepLimiterPhysics);
  physics->RegisterPhysics(new HypernucleusDecayPhysics);
  physics->SetDefaultCutValue(kProductionCut);
  runManager->SetUserInitialization(physics);
  runManager->SetUserInitialization(new ActionInitialization);

  G4PrimaryTransformer::SetKETolerance(kMassTolerance, IgnoreTheIssue);

  auto* uiManager = G4UImanager::GetUIpointer();
  if (ui) {
    ui->SessionStart();
    delete ui;
  } else {
    uiManager->ApplyCommand(G4String("/control/execute ") + argv[1]);
  }
  delete runManager;
  return 0;
}
