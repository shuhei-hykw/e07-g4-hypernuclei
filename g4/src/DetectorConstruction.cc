#include "DetectorConstruction.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UserLimits.hh"

namespace {

// Nominal density of E07 emulsion (g/cm3); calibrate with alpha ranges
constexpr double kDefaultDensity = 3.53;
constexpr double kDefaultHalfSize = 2. * cm;
constexpr double kDefaultMinStep = 0.1 * um;
constexpr double kDefaultMaxStep = 10. * um;
constexpr double kDefaultRangeFraction = 0.01;
constexpr double kWorldMargin = 1.1;
// Bragg-Kleeman estimate of the proton range in emulsion,
// R = coef * E^exponent (R in um, E in MeV), used only for step limits
constexpr double kProtonRangeCoef = 9.8 * um;
constexpr double kRangeExponent = 1.75;

// Approximate mass fractions of AgBr + gelatin emulsion
struct Fraction {
  const char* symbol;
  double weight;
};
constexpr Fraction kComposition[] = {
  {"Ag", 0.4428}, {"Br", 0.3249}, {"I", 0.0119}, {"C", 0.1027},
  {"H", 0.0156}, {"N", 0.0356}, {"O", 0.0638}, {"S", 0.0027},
};

// Rough residual range scaled from a proton of the same velocity.
// G4LossTableManager::GetRange() is not used because calling it here
// alters the effective charge state of generic ions.
double EstimateRange(const G4Track& track)
{
  const double z = track.GetDynamicParticle()->GetCharge() / eplus;
  if (z == 0.) return DBL_MAX;
  const double massRatio = track.GetDynamicParticle()->GetMass()
                           / CLHEP::proton_mass_c2;
  const double eProton = track.GetKineticEnergy() / massRatio / MeV;
  return kProtonRangeCoef * std::pow(eProton, kRangeExponent) * massRatio
         / (z * z);
}

// Step limit proportional to the residual range of charged particles, so
// that short tracks and track ends are sampled finely
class RangeDependentLimits : public G4UserLimits {
public:
  RangeDependentLimits(double minStep, double maxStep, double fraction)
    : fMin(minStep), fMax(maxStep), fFraction(fraction) {}
  G4double GetMaxAllowedStep(const G4Track& track) override
  {
    return std::clamp(fFraction * EstimateRange(track), fMin, fMax);
  }

private:
  double fMin;
  double fMax;
  double fFraction;
};

}  // namespace

DetectorConstruction::DetectorConstruction()
  : fDensity(kDefaultDensity),
    fHalfSize(kDefaultHalfSize),
    fMinStep(kDefaultMinStep),
    fMaxStep(kDefaultMaxStep),
    fRangeFraction(kDefaultRangeFraction)
{
  fMessenger = std::make_unique<G4GenericMessenger>(
    this, "/emul/", "Emulsion settings");
  fMessenger->DeclareProperty("density", fDensity,
                              "Emulsion density in g/cm3");
  fMessenger->DeclarePropertyWithUnit("halfSize", "mm", fHalfSize,
                                      "Half width of the emulsion block");
  fMessenger->DeclarePropertyWithUnit("minStep", "um", fMinStep,
                                      "Lower bound of the step limit");
  fMessenger->DeclarePropertyWithUnit("maxStep", "um", fMaxStep,
                                      "Upper bound of the step limit");
  fMessenger->DeclareProperty("rangeFraction", fRangeFraction,
                              "Step limit / residual range");
}

G4Material* DetectorConstruction::BuildEmulsion() const
{
  auto* nist = G4NistManager::Instance();
  const int nElements = sizeof(kComposition) / sizeof(kComposition[0]);
  auto* mat = new G4Material("Emulsion", fDensity * g / cm3, nElements);
  for (const auto& f : kComposition) {
    mat->AddElement(nist->FindOrBuildElement(f.symbol), f.weight);
  }
  return mat;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  auto* vacuum = G4NistManager::Instance()->FindOrBuildMaterial(
    "G4_Galactic");
  const double worldHalf = fHalfSize * kWorldMargin;
  auto* worldLV = new G4LogicalVolume(
    new G4Box("World", worldHalf, worldHalf, worldHalf), vacuum, "World");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World",
                                    nullptr, false, 0);

  auto* emulLV = new G4LogicalVolume(
    new G4Box("Emulsion", fHalfSize, fHalfSize, fHalfSize),
    BuildEmulsion(), "Emulsion");
  emulLV->SetUserLimits(
    new RangeDependentLimits(fMinStep, fMaxStep, fRangeFraction));
  new G4PVPlacement(nullptr, {}, emulLV, "Emulsion", worldLV, false, 0);
  return worldPV;
}
