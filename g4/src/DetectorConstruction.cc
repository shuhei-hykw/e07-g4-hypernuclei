#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"

namespace {

// Nominal density of E07 emulsion (g/cm3); calibrate with alpha ranges
constexpr double kDefaultDensity = 3.53;
constexpr double kDefaultMaxStepUm = 0.1;
constexpr double kEmulsionHalfSize = 1. * mm;
constexpr double kWorldMargin = 1.1;

// Approximate mass fractions of AgBr + gelatin emulsion
struct Fraction {
  const char* symbol;
  double weight;
};
constexpr Fraction kComposition[] = {
  {"Ag", 0.4428}, {"Br", 0.3249}, {"I", 0.0119}, {"C", 0.1027},
  {"H", 0.0156}, {"N", 0.0356}, {"O", 0.0638}, {"S", 0.0027},
};

}  // namespace

DetectorConstruction::DetectorConstruction()
  : fDensity(kDefaultDensity), fMaxStep(kDefaultMaxStepUm * um)
{
  fMessenger = std::make_unique<G4GenericMessenger>(
    this, "/emul/", "Emulsion settings");
  fMessenger->DeclareProperty("density", fDensity,
                              "Emulsion density in g/cm3");
  fMessenger->DeclarePropertyWithUnit("maxStep", "um", fMaxStep,
                                      "Maximum step length");
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
  const double worldHalf = kEmulsionHalfSize * kWorldMargin;
  auto* worldLV = new G4LogicalVolume(
    new G4Box("World", worldHalf, worldHalf, worldHalf), vacuum, "World");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World",
                                    nullptr, false, 0);

  auto* emulLV = new G4LogicalVolume(
    new G4Box("Emulsion", kEmulsionHalfSize, kEmulsionHalfSize,
              kEmulsionHalfSize),
    BuildEmulsion(), "Emulsion");
  emulLV->SetUserLimits(new G4UserLimits(fMaxStep));
  new G4PVPlacement(nullptr, {}, emulLV, "Emulsion", worldLV, false, 0);
  return worldPV;
}
