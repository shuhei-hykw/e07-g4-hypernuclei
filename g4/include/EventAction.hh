#pragma once

#include <map>
#include <string>
#include <vector>

#include "G4UserEventAction.hh"

class G4Step;
class TTree;

// Collects the trajectories of all tracks except e+-, gamma and neutrinos
// (primaries, hypernuclear decay products, pi- capture products, ...)
// and fills the TTree.
class EventAction : public G4UserEventAction {
public:
  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;

  void Book(TTree* tree);
  void RecordStep(const G4Step* step);

private:
  TTree* fTree = nullptr;
  std::map<int, std::size_t> fIndex;  // track ID -> record index

  int fEventId = 0;
  std::string fChannel;
  double fQValue = 0.;
  std::vector<int> fTrackId, fParentId;
  std::vector<std::string> fName, fProcess;
  std::vector<int> fZ, fA, fNLambda;
  std::vector<double> fKE, fPx, fPy, fPz;     // MeV, MeV/c
  std::vector<double> fRange, fLength;         // um
  std::vector<std::vector<float>> fX, fY, fZpos, fEdep;  // um, MeV
};
