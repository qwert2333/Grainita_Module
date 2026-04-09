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

#include <vector>

class G4Event;


/// Event action class

class MyEventAction : public G4UserEventAction
{
  public:
    MyEventAction() : fRunAction(nullptr) { 
      ResetEventData(); 
      for (int i = 0; i < 4; ++i) fHitCollID[i] = -1;
    }
    MyEventAction(MyRunAction* runaction ) : fRunAction(runaction) { 
      ResetEventData(); 
      for (int i = 0; i < 4; ++i) fHitCollID[i] = -1;
    }
    ~MyEventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    int get_counter_Cerenkov();
    void increment_counter_Cerenkov();

    int get_counter_Scintillation();
    void increment_counter_Scintillation();

    void ResetEventData(); 

    private :
    MyRunAction* fRunAction; 
    G4int fHitCollID[4]; 
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

    int counter_Cerenkov;
    int counter_Scintillation ;
};


#endif
