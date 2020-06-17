// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

#include <TH1F.h>
#include <cmath>
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;

// Simple access to collision
struct ATask {
  OutputObj<TH1F> vertex{TH1F("vertex", "vertex", 100, -10, 10)};

  void process(aod::McCollision const& mcCollision)
  {
    LOGF(info, "MC. vtx-z = %f", mcCollision.posZ());
    vertex->Fill(mcCollision.posZ());
  }
};

// Grouping between MC particles and collisions
struct BTask {
  OutputObj<TH1F> phiH{TH1F("phi", "phi", 100, 0., 2. * M_PI)};
  OutputObj<TH1F> etaH{TH1F("eta", "eta", 102, -2.01, 2.01)};

  bool isStable(Int_t pdg) const
  {
    //
    // Decide whether particle (pdg) is stable
    //

    // All ions/nucleons are considered as stable
    // Nuclear code is 10LZZZAAAI
    if (pdg > 1000000000)
      return kTRUE;

    constexpr Int_t kNstable = 18;
    Int_t pdgStable[kNstable] = {
      kGamma,      // Photon
      kElectron,   // Electron
      kMuonPlus,   // Muon
      kPiPlus,     // Pion
      kKPlus,      // Kaon
      kK0Short,    // K0s
      kK0Long,     // K0l
      kProton,     // Proton
      kNeutron,    // Neutron
      kLambda0,    // Lambda_0
      kSigmaMinus, // Sigma Minus
      kSigmaPlus,  // Sigma Plus
      3312,        // Xsi Minus
      3322,        // Xsi
      3334,        // Omega
      kNuE,        // Electron Neutrino
      kNuMu,       // Muon Neutrino
      kNuTau       // Tau Neutrino
    };

    bool isStable = false;
    for (int i = 0; i < kNstable; i++) {
      if (pdg == TMath::Abs(pdgStable[i])) {
        isStable = true;
        break;
      }
    }

    return isStable;
  }

  template <typename T>
  bool isPhysicalPrimary(aod::McParticles& mcParticles, T const& particle)
  {
    //
    // Test if a particle is a physical primary according to the following definition:
    // Particles produced in the collision including products of strong and
    // electromagnetic decay and excluding feed-down from weak decays of strange
    // particles.
    //

    const int ist = particle.statusCode();
    const int pdg = TMath::Abs(particle.pdgCode());

    //
    // Initial state particle
    // Solution for K0L decayed by Pythia6
    // ->
    if ((ist > 1) && (pdg != 130) && particle.producedByGenerator())
      return false;
    if ((ist > 1) && !particle.producedByGenerator())
      return false;
    // <-

    if (!isStable(pdg))
      return false;
    if (particle.producedByGenerator()) {
      //
      // Particle produced by generator
      // Solution for K0L decayed by Pythia6
      // ->
      if (particle.mother()[0] != -1) {
        auto mother = mcParticles.iteratorAt(particle.mother()[0]);
        if (std::abs(mother.pdgCode()) == 130)
          return false;
      }
      // <-
      // check for direct photon in parton shower
      // ->
      if (pdg == 22 && particle.daughter()[0] != -1) {
        LOGF(info, "D %d", particle.daughter()[0]);
        //         auto daughter = particle;
        //         daughter[particle.daughter()[0]];
        auto daughter = mcParticles.iteratorAt(particle.daughter()[0]);
        if (daughter.pdgCode() == 22)
          return false;
      }
      // <-
      return true;
    } else {
      //
      // Particle produced during transport
      //

      auto mother = mcParticles.iteratorAt(particle.mother()[0]);
      int mpdg = std::abs(mother.pdgCode());

      // Check for Sigma0
      if ((mpdg == 3212) && mother.producedByGenerator())
        return true;
      //
      // Check if it comes from a pi0 decay
      //
      if ((mpdg == kPi0) && mother.producedByGenerator())
        return true;

      // Check if this is a heavy flavor decay product
      int mfl = int(mpdg / TMath::Power(10, int(TMath::Log10(mother.pdgCode()))));
      //
      // Light hadron
      if (mfl < 4)
        return false;

      // Heavy flavor hadron produced by generator
      if (mother.producedByGenerator())
        return true;

      // To be sure that heavy flavor has not been produced in a secondary interaction
      // Loop back to the generated mother
      LOGF(info, "M0 %d %d", mother.producedByGenerator(), mother.mother()[0]);
      // TODO should be while
      if (mother.mother()[0] != -1 && !mother.producedByGenerator()) {
        auto mid = mother.mother()[0];
        mother = mcParticles.iteratorAt(mid);
        LOGF(info, "M+ %d %d", mother.producedByGenerator(), mother.mother()[0]);
        mpdg = std::abs(mother.pdgCode());
        mfl = int(mpdg / TMath::Power(10, int(TMath::Log10(mpdg))));
      }

      if (mfl < 4)
        return false;
      else
        return true;
    }
  }

  void process(aod::McCollision const& mcCollision, aod::McParticles& mcParticles)
  {
    LOGF(info, "MC. vtx-z = %f", mcCollision.posZ());
    for (auto& mcParticle : mcParticles) {
      if (isPhysicalPrimary(mcParticles, mcParticle)) {
        phiH->Fill(mcParticle.phi());
        etaH->Fill(mcParticle.eta());
      }
    }
  }
};

// Access from tracks to MC particle
struct CTask {
  OutputObj<TH1F> etaDiff{TH1F("etaDiff", ";eta_{MC} - eta_{Rec}", 100, -2, 2)};
  OutputObj<TH1F> phiDiff{TH1F("phiDiff", ";phi_{MC} - phi_{Rec}", 100, -M_PI, M_PI)};

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {
    LOGF(info, "vtx-z (data) = %f | vtx-z (MC) = %f", collision.posZ(), collision.label().posZ());
    for (auto& track : tracks) {
      //if (track.trackType() != 0)
      //  continue;
      //if (track.labelMask() != 0)
      //  continue;
      etaDiff->Fill(track.label().eta() - track.eta());
      auto delta = track.label().phi() - track.phi();
      if (delta > M_PI)
        delta -= 2 * M_PI;
      if (delta < -M_PI)
        delta += 2 * M_PI;
      phiDiff->Fill(delta);
      //LOGF(info, "eta: %.2f %.2f \t phi: %.2f %.2f | %d", track.label().eta(), track.eta(), track.label().phi(), track.phi(), track.label().index());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<ATask>("vertex-histogram"),
    adaptAnalysisTask<BTask>("etaphi-histogram"),
    adaptAnalysisTask<CTask>("eta-resolution-histogram"),
  };
}
