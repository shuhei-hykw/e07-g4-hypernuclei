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

  double fDensity;   // g/cm3
  double fMaxStep;   // mm, step limit for fine track sampling
  std::unique_ptr<G4GenericMessenger> fMessenger;
};
