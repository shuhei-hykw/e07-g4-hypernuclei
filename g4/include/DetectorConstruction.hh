#pragma once

#include <memory>
#include <vector>

#include "G4GenericMessenger.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"

class G4Material;

// One layer of the emulsion stack along the beam axis (z)
struct Layer {
  int sheet;      // sheet index from upstream
  bool emulsion;  // emulsion or base film
  bool thick;     // belongs to a thick-type sheet
  double zMin;
  double zMax;
};

// E07 emulsion module: thin-type sheets at both ends of a stack of
// thick-type sheets, each sheet being emulsion | base | emulsion, stacked
// along z without gaps. Tracks leaving the module are flagged and killed.
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  G4VPhysicalVolume* Construct() override;

  const std::vector<Layer>& Layers() const { return fLayers; }
  double HalfX() const { return fSheetX / 2.; }
  double HalfY() const { return fSheetY / 2.; }
  // Uniform point inside the emulsion layers of the thick-type sheets
  G4ThreeVector SampleVertex() const;

private:
  G4Material* BuildEmulsion() const;
  void BuildLayers();

  double fDensity;       // g/cm3
  double fSheetX;
  double fSheetY;
  double fThickEmulsion; // per side
  double fThickBase;
  double fThinEmulsion;  // per side
  double fThinBase;
  int fNumThick;
  int fNumThinUp;
  int fNumThinDown;
  double fMinStep;       // lower bound of the step limit
  double fMaxStep;       // upper bound of the step limit
  double fRangeFraction; // step limit relative to the residual range

  std::vector<Layer> fLayers;
  std::unique_ptr<G4GenericMessenger> fMessenger;
};
