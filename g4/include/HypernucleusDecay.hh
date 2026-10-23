#pragma once

#include <map>
#include <memory>
#include <string>

#include "G4GenericMessenger.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VRestProcess.hh"

struct DecayChannel;

// Chooses a weak-decay channel of a hypernucleus from decays.csv.
//   simple         : pi- mode with the fewest bodies, highest Q among them
//   uniform        : any listed channel with equal probability
//   uniformPiMinus : any pi- channel with equal probability
// /decay/channel <parent> <products> fixes the channel of one parent.
class DecaySelector {
public:
  DecaySelector();
  const DecayChannel* Select(const std::string& parent) const;

private:
  void SetChannel(const G4String& parent, const G4String& products);

  std::string fPolicy = "simple";
  std::map<std::string, std::string> fFixed;
  std::unique_ptr<G4GenericMessenger> fMessenger;
};

// Weak decay at rest of a stopped hypernucleus into the chosen channel,
// with the products sampled from n-body phase space. Stopped ordinary
// generic ions are simply killed, as without this process.
class HypernucleusDecay : public G4VRestProcess {
public:
  explicit HypernucleusDecay(const DecaySelector* selector);
  G4bool IsApplicable(const G4ParticleDefinition& def) override;
  G4VParticleChange* AtRestDoIt(const G4Track& track,
                                const G4Step& step) override;

protected:
  G4double GetMeanLifeTime(const G4Track&, G4ForceCondition*) override;

private:
  const DecaySelector* fSelector;
};

class HypernucleusDecayPhysics : public G4VPhysicsConstructor {
public:
  HypernucleusDecayPhysics();
  void ConstructParticle() override {}
  void ConstructProcess() override;

private:
  DecaySelector fSelector;
};
