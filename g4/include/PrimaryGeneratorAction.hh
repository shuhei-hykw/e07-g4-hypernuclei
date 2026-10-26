#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "G4GenericMessenger.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "SpeciesTable.hh"

// Xi- capture at rest: the final state of the chosen channel is generated
// isotropically according to n-body phase space (TGenPhaseSpace).
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction();
  void GeneratePrimaries(G4Event* event) override;

  const std::string& Channel() const { return fChannel; }
  const std::string& Target() const { return fTarget; }
  double QValue() const { return fQValue; }
  const G4ThreeVector& Vertex() const { return fVertex; }

private:
  struct GunTrack {
    std::string name;
    G4ThreeVector position;  // um
    G4ThreeVector momentum;  // MeV/c
  };

  void Initialize();
  void LoadGun();
  void GenerateGun(G4Event* event);

  std::string fSpeciesFile;
  std::string fTargetsFile;
  std::string fDecaysFile;
  std::string fReplayFile;  // fixed kinematics of an observed event
  // Single-track gun: CSV event,name,x_um,y_um,z_um,px_mev,py_mev,pz_mev;
  // the row named "center" of an event only sets its reference point
  std::string fGunFile;
  std::map<int, std::vector<GunTrack>> fGun;
  std::map<int, G4ThreeVector> fGunCenter;
  std::string fTarget;
  std::string fChannel;  // comma separated, e.g. "LL6He,4He,t"
  int fSeed;
  std::string fVertexMode;  // "random" or "fixed"
  G4ThreeVector fVertex;

  bool fInitialized = false;
  std::vector<SpeciesInfo> fProducts;
  const FixedVertex* fReplay = nullptr;
  double fInitialMass = 0.;
  double fQValue = 0.;
  std::unique_ptr<G4GenericMessenger> fMessenger;
};
