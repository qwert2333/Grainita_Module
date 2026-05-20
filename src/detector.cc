#include "detector.hh"

// Local parameters.
// Same treatment for crystals, fiber cores and fiber claddings. 
bool stopAndKillTracks = false;



MySensitiveDetector::MySensitiveDetector(G4String name) : G4VSensitiveDetector(name), fHitCollection(nullptr), fHitIndex(0), fHitCID(-1), fAttLength(1e10), fpitch(7.) {
  collectionName.insert(name + "_hits");
  fAttLength = 3.4; // in mm.
}

MySensitiveDetector::~MySensitiveDetector(){ }

void MySensitiveDetector::Initialize(G4HCofThisEvent *hitsCE){
  fHitCollection = new EcalHitsCollection();
  if (fHitCID < 0)
  {
      fHitCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  hitsCE->AddHitsCollection(fHitCID, fHitCollection);
  fHitIndex = 0;
  cellIDCol.clear();

  //std::cout<<"In MySD: SD name " << this->GetName()<<", fHitCID "<<fHitCID<<", collection name "<<collectionName[0]<<std::endl;
}


G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist)
{

    G4int eventID   = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    G4double edep = aStep->GetTotalEnergyDeposit(); 
    G4ThreeVector pos = (aStep->GetPreStepPoint()->GetPosition()/mm + aStep->GetPostStepPoint()->GetPosition()/mm ) * 0.5;
    G4Track *aTrack = aStep->GetTrack(); 
    G4String particleName = aTrack->GetParticleDefinition()->GetParticleName();

    const    G4VProcess*  theprocess = aTrack->GetCreatorProcess();
    G4String CreatorprocessName = "None" ;
    if (theprocess != 0) 
      CreatorprocessName = theprocess->GetProcessName();


    if (stopAndKillTracks == true){
        aTrack->SetTrackStatus(fStopAndKill);
        aTrack->SetTrackStatus(fKillTrackAndSecondaries);
    }

    G4int trackID = aTrack->GetTrackID(); 
    G4int parentID = aTrack->GetParentID();    

    G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    G4ThreeVector volPos = touchable->GetHistory()->GetTopTransform().Inverse().NetTranslation();

    //std::cout << "---- enter in MySensitiveDetector::ProcessHits  event = " << eventID << " trackID = " << trackID << std::endl;
    //std::cout << "  Step particle name: " << particleName << ", process name " << CreatorprocessName << std::endl;
    //std::cout << "  Step position: ("<<pos.x()<<", "<<pos.y()<<", "<<pos.z()<<"), stepE "<<edep<<std::endl;
    //std::cout << "  Cell volume position: ("<<volPos.x()<<", "<<volPos.y()<<", "<<volPos.z()<<") "<<std::endl;

    G4int cellID = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
    //G4int boxID = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
    //G4int globalID = 1e8*boxID + cellID;
    G4int globalID = cellID;
    // std::cout << "  Step E: "<<edep<<", cell ID: "<<cellID << ", boxID "<<boxID<<", global ID "<<globalID<<std::endl;

    // Mannually define the hit and cellID from global position 
    // WARNING: hard-coded to 1 mm x-y segmentation here!
    //G4int idx = (int)(pos.x()+500.);
    //G4int idy = (int)(pos.y()+500.);
    //G4int globalID = 1000*idy + idx; 


    // For attenuation & transverse light cross talk effect: each step create 5*5 hits. 
    // WARNING: hard-coded cellid coding. 
    G4int id_max = 24;

    G4int idx = cellID%100;
    G4int idy = (cellID / 100) % 100;
    G4int idz = cellID / 10000;
    //std::cout<<"  Decoded cellID: "<<idx<<", "<<idy<<", "<<idz<<std::endl;
    //std::cout<<"  Make lookup table for neighbor cells "<<std::endl;

    std::vector<G4int> cellIDvec;
    std::vector<G4ThreeVector> volPosvec; 
    cellIDvec.push_back(cellID);
    volPosvec.push_back(volPos);
    for(int i=-2; i<=2; i++ ){
      if(idx+i<1 || idx+i>id_max) continue;
      for(int j=-2; j<=2; j++ ){
        if(idy+j<1 || idy+j>id_max) continue;
        if(i==0 && j==0) continue;

        G4int neighborID = idz*10000 + (idy+j)*100 + idx+i;
        cellIDvec.push_back(neighborID);
        G4ThreeVector neighborVolPos = volPos + G4ThreeVector(i*fpitch, j*fpitch, 0 );
        volPosvec.push_back(neighborVolPos);

        //std::cout<<"    Neighbor ("<<i<<", "<<j<<"): cellID "<<neighborID<<", position ("
        //<<neighborVolPos.x()<<", "<<neighborVolPos.y()<<", "<<neighborVolPos.z()<<") "<<std::endl;
      }
    }

    //std::cout<<"  Create hits "<<std::endl;
    for(int ihit=0; ihit<cellIDvec.size(); ihit++){
      G4int cellID_local = cellIDvec[ihit];
      G4ThreeVector relPos = pos - volPosvec[ihit];
      G4double distance = relPos.perp();
      G4double edep_att = edep*exp(-distance/fAttLength);
      //std::cout<<"  Cell ID "<<cellID_local<<": distance to step "<<distance<<", Effective En "<<edep_att<<std::endl;

      EcalHit* hit = nullptr;
      if(cellIDCol.find(cellID_local) != cellIDCol.end()){
        hit = (*fHitCollection)[cellIDCol[cellID_local]];
        hit->addEdep(edep_att);
        if(particleName=="opticalphoton"){
          if(CreatorprocessName=="Cerenkov") hit->addNphChren(1); 
          else{ hit->addNphScint(1);}
        }
        hit->addStep(pos.x(), pos.y(), pos.z(), edep_att);
      }
      else{
        hit = new EcalHit(cellID_local);
        hit->setEdep(edep_att);
        if(particleName=="opticalphoton"){
          if(CreatorprocessName=="Cerenkov") hit->addNphChren(1);
          else{ hit->addNphScint(1);}
        }      
        hit->addStep(pos.x(), pos.y(), pos.z(), edep_att);

        fHitCollection->insert(hit);
        cellIDCol.insert(std::make_pair(cellID_local, fHitIndex));
        fHitIndex++;
      }

    }
    //std::cout <<"  Hit created. Current hit size: " << fHitCollection->entries() << std::endl;
    //std::cout << " exit from MySensitiveDetector::ProcessHits " << std::endl;
    return true;
}

