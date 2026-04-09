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

#include "MySteppingAction.hh"

//#include "construction.hh"
#include "MyEventAction.hh"

#include "G4Step.hh"
#include "globals.hh"

#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"

#include <vector>
#include <string>



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


MySteppingAction::MySteppingAction(MyEventAction* eventAction) : fEventAction(eventAction) {}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MySteppingAction::UserSteppingAction(const G4Step* step)
{
        // Collect energy and track length step by step

        //std::cout << " ---- ENTER in MySteppingAction::UserSteppingAction " << std::endl;


        G4Track *aTrack = step->GetTrack(); 
        G4int    trackID          = aTrack->GetTrackID();
        G4int parentID  = aTrack->GetParentID();


        const G4ParticleDefinition* particleDefinition = aTrack->GetParticleDefinition();
        G4String                    particleName       = particleDefinition->GetParticleName();

        const G4double              kineticEnergy      = aTrack->GetKineticEnergy() /MeV;

        const G4ThreeVector trackMomentumDir   = aTrack->GetMomentumDirection() ;
        double ax = trackMomentumDir.x();
        double ay = trackMomentumDir.y();
        double az = trackMomentumDir.z();

        G4ThreeVector pos           = step->GetPreStepPoint()->GetPosition() /mm;
        G4double preStepPositionX   = pos.x();
        G4double preStepPositionY   = pos.y();
        G4double preStepPositionZ   = pos.z();
        const G4ThreeVector preMomentumDir = step->GetPreStepPoint()->GetMomentumDirection();

        G4int eventID   = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
        G4double stepEnergyDeposit  = step->GetTotalEnergyDeposit() /MeV;
        G4double edep = stepEnergyDeposit;
        G4double nonIonizingEnergyDeposit  = step->GetNonIonizingEnergyDeposit() /MeV ;

        G4double stepLength         = step->GetStepLength()         /mm;
        G4double dEdx               = stepEnergyDeposit / stepLength;

        G4double trackLength = aTrack -> GetTrackLength() / mm;

        const G4ThreeVector trackVertex = aTrack->GetVertexPosition();
        double xvtx = trackVertex.x() ;
        double yvtx = trackVertex.y() ;
        double zvtx = trackVertex.z() ;

        const G4StepPoint* pPost = step->GetPostStepPoint();
        const G4VProcess* proc = pPost->GetProcessDefinedStep();
        G4String procName = (proc) ? proc->GetProcessName() : "Unknown";

        const    G4VProcess*  theprocess = aTrack->GetCreatorProcess();
        G4String CreatorprocessName = "None" ;
        if (theprocess != 0) {
                CreatorprocessName          = theprocess->GetProcessName(); // Get the name of the processG4VProcess creatorProcess = aTrack->GetCreatorProcess();
        }


        G4ThreeVector posPost           = step->GetPostStepPoint()->GetPosition() /mm;
        G4double postStepPositionX   = posPost.x();
        G4double postStepPositionY   = posPost.y();
        G4double postStepPositionZ   = posPost.z();


        G4double deltaEnergy = step ->GetDeltaEnergy() / MeV ;
        G4double time       = step->GetPreStepPoint()->GetGlobalTime() /ns;

        G4String  preVol = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
        G4String postVol = "None";
        G4VPhysicalVolume* volume =  step->GetPostStepPoint()->GetPhysicalVolume();
        if ( volume != 0 ) postVol = volume ->GetName();


        if ( particleName == "opticalphoton") {


                if ( CreatorprocessName == "Cerenkov") {
                        fEventAction -> increment_counter_Cerenkov() ;
                }
                else {
                        fEventAction -> increment_counter_Scintillation() ;
                }


                // For scintillation photons only:

                if ( CreatorprocessName != "Cerenkov") { 

                        // Photons in the first 4 X0. Set X0 = 1.7 cm
                        //if ( zvtx < -200 + 4 * 17 ) {
                        //        fEventAction -> add_nPhotons_4X0( 1. ) ;
                        //}

                }


        }

        int counter_Cerenkov = fEventAction ->get_counter_Cerenkov() ;
        int counter_Scintillation = fEventAction ->get_counter_Scintillation() ;



        // --------------------------------------------------------------------------
        //
        // --- Secondaries in this step
        //

        const std::vector< const G4Track * > *   secondaries = step -> GetSecondaryInCurrentStep () ;

        int nSecondaries = secondaries->size();
        // total energy carried by these secondaries :
        double EnergySecondaries = 0;

        double EleadingSec = 0;
        double EkinleadingSec = 0;
        G4String leadingSecName="None";
        for (int isec=0; isec < nSecondaries; isec++) {
                const G4Track* atrack = secondaries -> at(isec);
                G4double esec = atrack -> GetTotalEnergy() / MeV;
                EnergySecondaries += esec ;
                if (esec > EleadingSec) {
                        EleadingSec = esec ;
                        EkinleadingSec = atrack -> GetKineticEnergy () / MeV;
                        leadingSecName = atrack -> GetParticleDefinition()->GetParticleName();
                }
        }


        // --------------------------------------------------------------------------
        //
        // --- Do not track the optical photons (save time when the geometry os not trivial)

        //if (aTrack->GetDefinition() == G4OpticalPhoton::Definition())
        if (particleName == "opticalphoton")
        {
                aTrack->SetTrackStatus(fStopAndKill);
        }


        // --------------------------------------------------------------------------




        // --------------------------------------------------------------------------
        //
        // ----   Fill step-wise ntuple
        //

        bool doFillStepwise = false;
        //bool doFillStepwise =  true;       // set it to true for tests only, over a small number of events...


        if ( doFillStepwise) {

                G4AnalysisManager *man5 = G4AnalysisManager::Instance();
                int idx = 0;

                G4int ntupleID = 1;

                man5->FillNtupleIColumn( ntupleID, idx, eventID);  idx ++;
                man5->FillNtupleIColumn( ntupleID, idx, trackID);  idx ++ ;
                man5->FillNtupleSColumn( ntupleID, idx, particleName);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, kineticEnergy);     idx ++;
                man5->FillNtupleSColumn( ntupleID, idx, procName);     idx ++;
                man5->FillNtupleSColumn( ntupleID, idx, CreatorprocessName);   idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, stepEnergyDeposit);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, preStepPositionX);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, preStepPositionY);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, preStepPositionZ);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, postStepPositionX);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, postStepPositionY);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, postStepPositionZ);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, deltaEnergy);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, time);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, stepLength );     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, trackLength );     idx ++;

                man5->FillNtupleDColumn( ntupleID, idx, dEdx );     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, nonIonizingEnergyDeposit);     idx ++;
                man5->FillNtupleIColumn( ntupleID, idx, parentID);  idx ++;
                man5->FillNtupleIColumn( ntupleID, idx, nSecondaries);     idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, EnergySecondaries);     idx ++;
                man5->FillNtupleIColumn( ntupleID, idx, counter_Cerenkov); idx ++;
                man5->FillNtupleIColumn( ntupleID, idx, counter_Scintillation);  idx ++;
                man5->FillNtupleSColumn( ntupleID, idx, preVol);     idx ++;
                man5->FillNtupleSColumn( ntupleID, idx, postVol);     idx ++;

                man5->FillNtupleDColumn( ntupleID, idx, xvtx ); idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, yvtx ); idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, zvtx ); idx ++;

                // track momentum :
                man5->FillNtupleDColumn( ntupleID, idx, ax ); idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, ay ); idx ++;
                man5->FillNtupleDColumn( ntupleID, idx, az ); idx ++;

                man5->AddNtupleRow( ntupleID);

        }




        //std::cout << " ----  EXIT MySteppingAction::UserSteppingAction " << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

