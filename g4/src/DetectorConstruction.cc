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
#include "Randomize.hh"

namespace {

// Nominal density of E07 emulsion (g/cm3); calibrate with alpha ranges
constexpr double kDefaultDensity = 3.53;
// E07 module (nominal values before development)
constexpr double kDefaultSheetX = 345. * mm;
constexpr double kDefaultSheetY = 350. * mm;
constexpr double kDefaultThickEmulsion = 450. * um;
constexpr double kDefaultThickBase = 40. * um;
constexpr double kDefaultThinEmulsion = 100. * um;
constexpr double kDefaultThinBase = 180. * um;
constexpr int kDefaultNumThick = 11;
constexpr int kDefaultNumThin = 1;
constexpr double kDefaultMinStep = 0.1 * um;
constexpr double kDefaultMaxStep = 10. * um;
constexpr double kDefaultRangeFraction = 0.01;
constexpr double kWorldMargin = 1. * cm;
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
    fSheetX(kDefaultSheetX),
    fSheetY(kDefaultSheetY),
    fThickEmulsion(kDefaultThickEmulsion),
    fThickBase(kDefaultThickBase),
    fThinEmulsion(kDefaultThinEmulsion),
    fThinBase(kDefaultThinBase),
    fNumThick(kDefaultNumThick),
    fNumThinUp(kDefaultNumThin),
    fNumThinDown(kDefaultNumThin),
    fMinStep(kDefaultMinStep),
    fMaxStep(kDefaultMaxStep),
    fRangeFraction(kDefaultRangeFraction)
{
  fMessenger = std::make_unique<G4GenericMessenger>(
    this, "/emul/", "Emulsion module settings");
  auto& m = *fMessenger;
  m.DeclareProperty("density", fDensity, "Emulsion density in g/cm3");
  m.DeclarePropertyWithUnit("sheetX", "mm", fSheetX, "Sheet width (x)");
  m.DeclarePropertyWithUnit("sheetY", "mm", fSheetY, "Sheet height (y)");
  m.DeclarePropertyWithUnit("thickEmulsion", "um", fThickEmulsion,
                            "Emulsion layer of thick sheets (per side)");
  m.DeclarePropertyWithUnit("thickBase", "um", fThickBase,
                            "Base film of thick sheets");
  m.DeclarePropertyWithUnit("thinEmulsion", "um", fThinEmulsion,
                            "Emulsion layer of thin sheets (per side)");
  m.DeclarePropertyWithUnit("thinBase", "um", fThinBase,
                            "Base film of thin sheets");
  m.DeclareProperty("numThick", fNumThick, "Number of thick sheets");
  m.DeclareProperty("numThinUp", fNumThinUp, "Thin sheets upstream");
  m.DeclareProperty("numThinDown", fNumThinDown, "Thin sheets downstream");
  m.DeclarePropertyWithUnit("minStep", "um", fMinStep,
                            "Lower bound of the step limit");
  m.DeclarePropertyWithUnit("maxStep", "um", fMaxStep,
                            "Upper bound of the step limit");
  m.DeclareProperty("rangeFraction", fRangeFraction,
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

void DetectorConstruction::BuildLayers()
{
  fLayers.clear();
  std::vector<bool> sheets(fNumThinUp, false);
  sheets.insert(sheets.end(), fNumThick, true);
  sheets.insert(sheets.end(), fNumThinDown, false);

  double total = 0.;
  for (bool thick : sheets) {
    total += thick ? 2. * fThickEmulsion + fThickBase
                   : 2. * fThinEmulsion + fThinBase;
  }
  // The module is centered at z = 0
  double z = -total / 2.;
  for (std::size_t i = 0; i < sheets.size(); ++i) {
    const bool thick = sheets[i];
    const double emul = thick ? fThickEmulsion : fThinEmulsion;
    const double base = thick ? fThickBase : fThinBase;
    for (auto [isEmul, dz] : {std::pair{true, emul}, std::pair{false, base},
                              std::pair{true, emul}}) {
      fLayers.push_back({int(i), isEmul, thick, z, z + dz});
      z += dz;
    }
  }
}

G4ThreeVector DetectorConstruction::SampleVertex() const
{
  std::vector<const Layer*> targets;
  for (const auto& l : fLayers) {
    if (l.emulsion && l.thick) targets.push_back(&l);
  }
  // All thick emulsion layers have the same thickness
  const auto* l = targets[G4int(G4UniformRand() * targets.size())];
  return {(G4UniformRand() - 0.5) * fSheetX,
          (G4UniformRand() - 0.5) * fSheetY,
          l->zMin + G4UniformRand() * (l->zMax - l->zMin)};
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  BuildLayers();
  auto* nist = G4NistManager::Instance();
  auto* emulsion = BuildEmulsion();
  auto* base = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
  auto* air = nist->FindOrBuildMaterial("G4_AIR");

  const double halfZ = std::max(-fLayers.front().zMin, fLayers.back().zMax);
  auto* worldLV = new G4LogicalVolume(
    new G4Box("World", HalfX() + kWorldMargin, HalfY() + kWorldMargin,
              halfZ + kWorldMargin),
    air, "World");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World",
                                    nullptr, false, 0);

  // Mother volume of the sheets: a track leaving it has escaped. Using a
  // mother avoids false escapes at the shared faces of adjacent layers.
  auto* moduleLV = new G4LogicalVolume(
    new G4Box("Module", HalfX(), HalfY(), halfZ), air, "Module");
  new G4PVPlacement(nullptr, {}, moduleLV, "Module", worldLV, false, 0);

  auto* limits = new RangeDependentLimits(fMinStep, fMaxStep,
                                          fRangeFraction);
  for (std::size_t i = 0; i < fLayers.size(); ++i) {
    const auto& l = fLayers[i];
    const char* name = l.emulsion ? "Emulsion" : "Base";
    const double halfDz = (l.zMax - l.zMin) / 2.;
    auto* lv = new G4LogicalVolume(
      new G4Box(name, HalfX(), HalfY(), halfDz),
      l.emulsion ? emulsion : base, name);
    lv->SetUserLimits(limits);
    new G4PVPlacement(nullptr, {0., 0., l.zMin + halfDz}, lv, name,
                      moduleLV, false, int(i));
  }
  return worldPV;
}
