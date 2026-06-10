//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B2/B2a/include/EventAction.hh
/// \brief Definition of the B2::EventAction class

#ifndef MyEventAction_h
#define MyEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4String.hh"
#include "analysis.hh"

#include <map>
#include <vector>

class G4Event;
class G4Track;


/// Event action class

class MyEventAction : public G4UserEventAction
{
  public:
    MyEventAction() : fRunAction(nullptr) { 
      ResetEventData(); 
      for (int i = 0; i < 5; ++i) fHitCollID[i] = -1;
    }
    MyEventAction(MyRunAction* runaction ) : fRunAction(runaction) { 
      ResetEventData(); 
      for (int i = 0; i < 5; ++i) fHitCollID[i] = -1;
    }
    ~MyEventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    int get_counter_Cerenkov();
    void increment_counter_Cerenkov();

    int get_counter_Scintillation();
    void increment_counter_Scintillation();

    G4bool IsEMComponentTrack(const G4Track* track);
    void RegisterSecondaryTrack(const G4Track* parent, const G4Track* secondary);
    void AddCrystalTruthEdep(const G4Track* track, G4double edep);
    void AddLeakageEnergy(G4double energy);

    void ResetEventData(); 

    private :
    G4bool IsEMParticle(const G4String& particleName) const;
    G4bool IsNeutralMeson(const G4String& particleName) const;

    MyRunAction* fRunAction; 
    G4int fHitCollID[5]; 
    G4int eventID;
    G4String particle;
    G4double MCtruth_energy;
    G4double MCtruth_dir_x;
    G4double MCtruth_dir_y;
    G4double MCtruth_dir_z;
    G4double MCtruth_pos_x;
    G4double MCtruth_pos_y;
    G4double MCtruth_pos_z;
    G4double EdepCrystal;
    G4double EdepFiberCore;
    G4double EdepFiberClad;
    G4double EdepCarbonFrame;
    G4int Nph_Cherenkov;
    G4int Nph_Scint;

    G4double truthEdepCrystalTotal;
    G4double truthEdepCrystalEM;
    G4double leakageEnergy;
    G4int leakageNParticles;
    std::map<G4int, G4bool> trackIsEM;
    std::map<G4int, G4String> trackParticleName;

    int counter_Cerenkov;
    int counter_Scintillation ;
};


#endif
