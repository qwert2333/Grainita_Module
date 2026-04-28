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

#include "MyEventAction.hh"
#include "generator.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4ios.hh"


#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"


// -- for MSC debug
#include "G4WentzelOKandVIxSection.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4AntiProton.hh"
#include "construction.hh"
#include "G4EmParameters.hh"
#include "G4SDManager.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyEventAction::BeginOfEventAction(const G4Event* event) {

  ResetEventData();
  const MyDetectorConstruction *detector = static_cast<const MyDetectorConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  if ( fHitCollID[0]<0 )
  {
      auto sdManager = G4SDManager::GetSDMpointer();
      fHitCollID[0] = sdManager->GetCollectionID("CrystalModule_hits");
      fHitCollID[1] = sdManager->GetCollectionID("FiberCore_hits");
      fHitCollID[2] = sdManager->GetCollectionID("FiberCladding_hits");
      fHitCollID[3] = sdManager->GetCollectionID("CarbonFrame_hits");
  }

  if(!fRunAction)
    fRunAction = static_cast<MyRunAction *>(const_cast<G4UserRunAction *>(G4RunManager::GetRunManager()->GetUserRunAction()));
  fRunAction->ResetEventData();

}

int MyEventAction::get_counter_Cerenkov() {
    return counter_Cerenkov;
}

void MyEventAction::increment_counter_Cerenkov() {
    counter_Cerenkov ++;
}

int MyEventAction::get_counter_Scintillation() {
    return counter_Scintillation;
}

void MyEventAction::increment_counter_Scintillation() {
    counter_Scintillation ++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyEventAction::EndOfEventAction(const G4Event* event)
{


  // get number of stored trajectories
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  std::size_t n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // periodic printing

  G4int eventID = event->GetEventID();
  if (eventID < 100 || eventID % 100 == 0) {
    G4cout << ">>> Event: " << eventID << G4endl;
  }

  // G4cout << fHitCollID[0] << '\t' << fHitCollID[1] << '\t' << fHitCollID[2] << '\t' << fHitCollID[3]<<G4endl;
  // Extract MC info
  MyPrimaryGenerator* fPrimaryGen = fRunAction->GetGenerator();
  if(fPrimaryGen){
    particle = fPrimaryGen->GetParticleGun()->GetParticleDefinition()->GetParticleName();
    MCtruth_energy = fPrimaryGen->GetParticleGun()->GetParticleEnergy();
    MCtruth_dir_x = fPrimaryGen->GetParticleGun()->GetParticleMomentumDirection().x();
    MCtruth_dir_y = fPrimaryGen->GetParticleGun()->GetParticleMomentumDirection().y();
    MCtruth_dir_z = fPrimaryGen->GetParticleGun()->GetParticleMomentumDirection().z();
    MCtruth_pos_x = fPrimaryGen->GetParticleGun()->GetParticlePosition().x();
    MCtruth_pos_y = fPrimaryGen->GetParticleGun()->GetParticlePosition().y();
    MCtruth_pos_z = fPrimaryGen->GetParticleGun()->GetParticlePosition().z();
  }

  //Get hit info
  G4HCofThisEvent *hitsCE = event->GetHCofThisEvent();
  EcalHitsCollection *CrystalHitCollection = nullptr;
  EcalHitsCollection *FiberCoreHitCollection = nullptr;
  EcalHitsCollection *FiberCladHitCollection = nullptr;
  EcalHitsCollection *CarbonFrameHitCollection = nullptr;
  if (hitsCE){
    CrystalHitCollection = (EcalHitsCollection *)(hitsCE->GetHC(fHitCollID[0]));
    FiberCoreHitCollection = (EcalHitsCollection *)(hitsCE->GetHC(fHitCollID[1]));
    FiberCladHitCollection = (EcalHitsCollection *)(hitsCE->GetHC(fHitCollID[2]));
    CarbonFrameHitCollection = (EcalHitsCollection *)(hitsCE->GetHC(fHitCollID[3]));
  }

  if(CrystalHitCollection){
    size_t n_hit = CrystalHitCollection->entries();
    //G4cout << " Read Crystal hits: " << n_hit << G4endl;
    for(size_t i=0; i<n_hit; i++){
      G4int cellID = (*CrystalHitCollection)[i]->GetCellID();
      EdepCrystal += (*CrystalHitCollection)[i]->GetEdep();
      Nph_Cherenkov += (*CrystalHitCollection)[i]->GetNphChren();
      Nph_Scint += (*CrystalHitCollection)[i]->GetNphScint();

      fRunAction->Fill_vecCellID( (*CrystalHitCollection)[i]->GetCellID() );
      fRunAction->Fill_vecEdep(   (*CrystalHitCollection)[i]->GetEdep() );
      fRunAction->Fill_vecNChren( (*CrystalHitCollection)[i]->GetNphChren() );
      fRunAction->Fill_vecNScint( (*CrystalHitCollection)[i]->GetNphScint() );
    }
    //for(size_t i=0; i<n_hit; i++){
    //  for(int istep = 0; istep<(*CrystalHitCollection)[i]->GetStepSize(); istep++){
    //    G4ThreeVector pos = (*CrystalHitCollection)[i]->GetStepPos(istep);
    //    G4double Edep = (*CrystalHitCollection)[i]->GetStepE(istep);
    //    fRunAction->Fill_stepPos( pos.x(), pos.y(), pos.z() );
    //    fRunAction->Fill_stepEn( Edep );
    //  }
    //}

  }

  if(FiberCoreHitCollection){
    size_t n_hit = FiberCoreHitCollection->entries();
    //G4cout << " Read Fiber Core Hits: " << n_hit << G4endl;
    for(size_t i=0; i<n_hit; i++){
      EdepFiberCore += (*FiberCoreHitCollection)[i]->GetEdep();
    }
  }

  if(FiberCladHitCollection){
    size_t n_hit = FiberCladHitCollection->entries();
    //G4cout << " Read Fiber Cladding Hits: " << n_hit << G4endl;
    for(size_t i=0; i<n_hit; i++){
      EdepFiberClad += (*FiberCladHitCollection)[i]->GetEdep();
    }
  }

  if(CarbonFrameHitCollection){
    size_t n_hit = CarbonFrameHitCollection->entries();
    //G4cout << " Read Carbon frame Hits: " << n_hit << G4endl;
    for(size_t i=0; i<n_hit; i++){
      EdepCarbonFrame += (*CarbonFrameHitCollection)[i]->GetEdep();
    }
  }
  //G4cout<<" Deposite energy: "<<EdepCrystal<<'\t'<<EdepFiberCore<<'\t'<<EdepFiberClad<<'\t'<<EdepCarbonFrame<<G4endl;

   // fill the ntuple 
   G4AnalysisManager *man5 = G4AnalysisManager::Instance();
   G4int ntupleID = 0;
   int idx = 0;

   //std::cout << " ==== FILL ntuple  event =  " << eventID << " " << captured << " " << backscattered << " " << traversed << std::endl;
   //G4cout << " Event #" << eventID << ": Scintillating photon " << Nph_Scint << ", Cherenkov photon " << Nph_Cherenkov << G4endl;
   man5->FillNtupleIColumn( ntupleID, idx, eventID);  idx ++;
   man5->FillNtupleSColumn( ntupleID, idx, particle);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, MCtruth_energy);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, MCtruth_dir_x);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, MCtruth_dir_y);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, MCtruth_dir_z);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, MCtruth_pos_x);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, MCtruth_pos_y);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, MCtruth_pos_z);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, EdepCrystal);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, EdepFiberCore);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, EdepFiberClad);  idx ++;
   man5->FillNtupleDColumn( ntupleID, idx, EdepCarbonFrame);  idx ++;
   man5->FillNtupleIColumn( ntupleID, idx, Nph_Cherenkov); idx ++;
   man5->FillNtupleIColumn( ntupleID, idx, Nph_Scint);  idx ++;
   man5->FillNtupleIColumn( ntupleID, idx, counter_Cerenkov); idx ++;
   man5->FillNtupleIColumn( ntupleID, idx, counter_Scintillation);  idx ++;

   man5->AddNtupleRow( ntupleID);

   //std::cout << " ------ EXIT MyEventAction::EndOfEventAction ---- " << std::endl;
}

void MyEventAction::ResetEventData()
{
  eventID = 0;
  particle = "";
  MCtruth_energy = 0;
  MCtruth_dir_x = 0;
  MCtruth_dir_y = 0;
  MCtruth_dir_z = 0;
  MCtruth_pos_x = 0;
  MCtruth_pos_y = 0;
  MCtruth_pos_z = 0;
  EdepCrystal = 0.;
  EdepFiberCore = 0.;
  EdepFiberClad = 0.;
  EdepCarbonFrame = 0.;
  Nph_Cherenkov = 0;
  Nph_Scint = 0;

  counter_Scintillation = 0;
  counter_Cerenkov = 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

