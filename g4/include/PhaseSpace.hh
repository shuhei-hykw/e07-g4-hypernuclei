#pragma once

#include <vector>

#include "G4ThreeVector.hh"
#include "SpeciesTable.hh"

// Momenta (MeV/c) of an isotropic, unweighted n-body phase-space decay of
// a system of mass parentMass (MeV) at rest, sampled with TGenPhaseSpace.
std::vector<G4ThreeVector> GeneratePhaseSpace(
  double parentMass, const std::vector<SpeciesInfo>& products);
