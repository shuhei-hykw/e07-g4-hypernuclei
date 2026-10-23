#pragma once

#include <memory>

#include "G4GenericMessenger.hh"
#include "G4VUserDetectorConstruction.hh"

class G4Material;

// A uniform block of nuclear emulsion with the capture point at its center.
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  G4VPhysicalVolume* Construct() override;

private:
  G4Material* BuildEmulsion() const;

  double fDensity;     // g/cm3
  double fHalfSize;    // half width of the emulsion block
  double fMinStep;       // lower bound of the step limit
  double fMaxStep;       // upper bound of the step limit
  double fRangeFraction; // step limit relative to the residual range
  std::unique_ptr<G4GenericMessenger> fMessenger;
};
