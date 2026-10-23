#include "HypernucleusDecay.hh"

#include <set>

#include "G4DynamicParticle.hh"
#include "G4GenericIon.hh"
#include "G4NuclearStopping.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"
#include "PhaseSpace.hh"
#include "Randomize.hh"
#include "SpeciesTable.hh"

namespace {

const std::string kUniform = "uniform";
const std::string kUniformPiMinus = "uniformPiMinus";
const std::string kPiMinus = "pi-";

}  // namespace

// ---------------------------------------------------------------------------
DecaySelector::DecaySelector()
{
  fMessenger = std::make_unique<G4GenericMessenger>(
    this, "/decay/", "Weak decay of hypernuclei");
  auto& policy = fMessenger->DeclareProperty(
    "policy", fPolicy, "simple, uniform or uniformPiMinus");
  policy.SetCandidates("simple uniform uniformPiMinus");
  fMessenger->DeclareMethod("channel", &DecaySelector::SetChannel,
                            "<parent> <comma separated products>");
}

void DecaySelector::SetChannel(const G4String& parent,
                               const G4String& products)
{
  fFixed[parent] = products;
}

const DecayChannel* DecaySelector::Select(const std::string& parent) const
{
  const auto& channels = SpeciesTable::Instance().Decays(parent);
  if (channels.empty()) return nullptr;

  auto fixed = fFixed.find(parent);
  if (fixed != fFixed.end()) {
    for (const auto& ch : channels) {
      if (ch.label == fixed->second) return &ch;
    }
    G4Exception("DecaySelector", "Decay001", FatalException,
                ("no channel " + parent + " -> " + fixed->second).c_str());
  }

  std::vector<const DecayChannel*> all, piMinus;
  for (const auto& ch : channels) {
    all.push_back(&ch);
    if (ch.mode == kPiMinus) piMinus.push_back(&ch);
  }
  if (fPolicy == kUniform) return all[G4int(G4UniformRand() * all.size())];

  // Restrict to pi- modes when the parent has any
  const auto& pool = piMinus.empty() ? all : piMinus;
  if (fPolicy == kUniformPiMinus) {
    return pool[G4int(G4UniformRand() * pool.size())];
  }

  // simple: fewest bodies, then highest Q
  const DecayChannel* best = nullptr;
  for (const auto* ch : pool) {
    if (!best || ch->products.size() < best->products.size()
        || (ch->products.size() == best->products.size()
            && ch->q > best->q)) {
      best = ch;
    }
  }
  return best;
}

// ---------------------------------------------------------------------------
HypernucleusDecay::HypernucleusDecay(const DecaySelector* selector)
  : G4VRestProcess("HypernucleusDecay", fUserDefined), fSelector(selector)
{}

G4bool HypernucleusDecay::IsApplicable(const G4ParticleDefinition& def)
{
  return &def == G4GenericIon::Definition()
         || def.GetNumberOfLambdasInHypernucleus() > 0;
}

G4double HypernucleusDecay::GetMeanLifeTime(const G4Track&,
                                            G4ForceCondition*)
{
  // Decay immediately once stopped; the time is irrelevant here
  return 0.;
}

G4VParticleChange* HypernucleusDecay::AtRestDoIt(const G4Track& track,
                                                 const G4Step&)
{
  aParticleChange.Initialize(track);
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  const auto* def = track.GetDefinition();
  const int nLambda = def->GetNumberOfLambdasInHypernucleus();
  if (nLambda <= 0) return &aParticleChange;

  const auto& table = SpeciesTable::Instance();
  const int z = static_cast<int>(def->GetPDGCharge() / eplus + 0.5);
  const auto* parent = table.Find(z, def->GetBaryonNumber(), nLambda);
  const auto* fixed = parent ? table.Replay(parent->name) : nullptr;
  const auto* channel = parent && !fixed
    ? fSelector->Select(parent->name) : nullptr;
  if (!fixed && !channel) {
    static std::set<std::string> warned;
    if (warned.insert(def->GetParticleName()).second) {
      G4Exception("HypernucleusDecay", "Decay002", JustWarning,
                  ("no decay channel for " + def->GetParticleName()
                   + "; killed").c_str());
    }
    return &aParticleChange;
  }

  // Fixed momenta of an observed decay, or phase space
  const auto& products = fixed ? fixed->products : channel->products;
  const auto momenta = fixed
    ? fixed->momenta : GeneratePhaseSpace(parent->mass, products);
  aParticleChange.SetNumberOfSecondaries(momenta.size());
  for (std::size_t i = 0; i < momenta.size(); ++i) {
    const auto& s = products[i];
    const G4ThreeVector p = momenta[i] * MeV;
    const double m = s.mass * MeV;
    const G4LorentzVector p4(p, std::sqrt(p.mag2() + m * m));
    aParticleChange.AddSecondary(
      new G4DynamicParticle(SpeciesTable::Definition(s), p4));
  }
  return &aParticleChange;
}

// ---------------------------------------------------------------------------
HypernucleusDecayPhysics::HypernucleusDecayPhysics()
  : G4VPhysicsConstructor("HypernucleusDecay")
{}

void HypernucleusDecayPhysics::ConstructProcess()
{
  auto* process = new HypernucleusDecay(&fSelector);
  auto* genericManager = G4GenericIon::Definition()->GetProcessManager();
  genericManager->AddRestProcess(process);

  // Light hypernuclei predefined by Geant4 (hypertriton, hyperHe5, ...)
  // have their own process managers: replace their decay and make sure
  // they lose energy like generic ions.
  auto* helper = G4PhysicsListHelper::GetPhysicsListHelper();
  auto* it = G4ParticleTable::GetParticleTable()->GetIterator();
  it->reset();
  while ((*it)()) {
    auto* particle = it->value();
    auto* manager = particle->GetProcessManager();
    if (particle->GetNumberOfLambdasInHypernucleus() <= 0 || !manager
        || manager == genericManager) {
      continue;
    }
    if (auto* decay = manager->GetProcess("Decay")) {
      manager->RemoveProcess(decay);
    }
    if (particle->GetPDGCharge() != 0. && !manager->GetProcess("ionIoni")
        && !manager->GetProcess("hIoni")) {
      helper->RegisterProcess(new G4hMultipleScattering("ionmsc"),
                              particle);
      helper->RegisterProcess(new G4ionIonisation, particle);
      helper->RegisterProcess(new G4NuclearStopping, particle);
    }
    manager->AddRestProcess(process);
  }
}
