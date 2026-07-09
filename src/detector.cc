#include "detector.hh"
#include "construction.hh"

#include <cmath>
#include <numeric>

// Local parameters.
// Same treatment for crystals, fiber cores and fiber claddings. 
bool stopAndKillTracks = false;



MySensitiveDetector::MySensitiveDetector(G4String name)
    : G4VSensitiveDetector(name),
      fHitCollection(nullptr),
      fRawHitCollection(nullptr),
      fHitCID(-1),
      fRawHitCID(-1),
      fHitIndex(0),
      fRawHitIndex(0),
      fStoreRawHits(name == "CrystalModule"),
      fAttLength(3.4),
      fpitch(7.),
      fIdMax(24),
      fCellIDBase(100),
      fResponseX0(0.856),
      fResponseSlope(0.93),
      fResponseIntercept(0.206),
      fResponseNorm(1. / std::exp(-fResponseX0 / fAttLength)),
      fReflectCoeff(0.9) {
  collectionName.insert(name + "_hits");
  if (fStoreRawHits) {
    collectionName.insert(name + "_raw_hits");
  }
}

MySensitiveDetector::~MySensitiveDetector(){ }

G4double MySensitiveDetector::ResponseFunction(G4double distance) const
{
    if (distance < fResponseX0) {
        return fResponseSlope * distance + fResponseIntercept;
    }
    return fResponseNorm * std::exp(-distance / fAttLength);
}

G4double MySensitiveDetector::AttenuatedResponse(G4double stepX, G4double stepY, G4double hitX, G4double hitY) const
{
    // const G4double xMin = -0.5 * fIdMax * fpitch;
    // const G4double xMax =  0.5 * fIdMax * fpitch;
    // const G4double yMin = -0.5 * fIdMax * fpitch;
    // const G4double yMax =  0.5 * fIdMax * fpitch;
    // const G4double doubleReflectCoeff = fReflectCoeff * fReflectCoeff;

    // const G4double imageX[9] = {
    //     stepX,
    //     2. * xMin - stepX,
    //     2. * xMax - stepX,
    //     stepX,
    //     stepX,
    //     2. * xMin - stepX,
    //     2. * xMin - stepX,
    //     2. * xMax - stepX,
    //     2. * xMax - stepX
    // };
    // const G4double imageY[9] = {
    //     stepY,
    //     stepY,
    //     stepY,
    //     2. * yMin - stepY,
    //     2. * yMax - stepY,
    //     2. * yMin - stepY,
    //     2. * yMax - stepY,
    //     2. * yMin - stepY,
    //     2. * yMax - stepY
    // };
    // const G4double imageWeight[9] = {
    //     1.0,
    //     fReflectCoeff,
    //     fReflectCoeff,
    //     fReflectCoeff,
    //     fReflectCoeff,
    //     doubleReflectCoeff,
    //     doubleReflectCoeff,
    //     doubleReflectCoeff,
    //     doubleReflectCoeff
    // };

    // G4double response = 0.;
    // for (G4int i = 0; i < 9; ++i) {
    //     const G4double dx = hitX - imageX[i];
    //     const G4double dy = hitY - imageY[i];
    //     response += imageWeight[i] * ResponseFunction(std::sqrt(dx * dx + dy * dy));
    // }
    // return response;
    G4double distance = std::sqrt((hitX - stepX) * (hitX - stepX) + (hitY - stepY) * (hitY - stepY));
    G4double response = ResponseFunction(distance);
    // std::cout<<"  Hit position: ("<<hitX<<", "<<hitY<<"), Step position: ("<<stepX<<", "<<stepY<<"), distance: "<<distance<<", response: "<<response<<std::endl;
    return response;
}

void MySensitiveDetector::Initialize(G4HCofThisEvent *hitsCE){
  fHitCollection = new EcalHitsCollection();
  if (fHitCID < 0)
  {
      fHitCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  hitsCE->AddHitsCollection(fHitCID, fHitCollection);
  fHitIndex = 0;
  cellIDCol.clear();

  if (fStoreRawHits) {
      fRawHitCollection = new EcalHitsCollection();
      if (fRawHitCID < 0)
      {
          fRawHitCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[1]);
      }
      hitsCE->AddHitsCollection(fRawHitCID, fRawHitCollection);
      fRawHitIndex = 0;
      rawCellIDCol.clear();
  }

  const auto *detector = static_cast<const MyDetectorConstruction *>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  if (detector) {
      fIdMax = detector->GetFiberNum();
      fCellIDBase = detector->GetCellIDBase();
      fpitch = detector->GetPitchSize() / mm;
      fResponseX0 = detector->GetResponseX0() / mm;
      fAttLength = detector->GetAttLength() / mm;
      fResponseSlope = detector->GetResponseSlope();
      fResponseIntercept = detector->GetResponseIntercept();
      fReflectCoeff = detector->GetReflectCoeff();
      fResponseNorm = 1. / std::exp(-fResponseX0 / fAttLength);
  }

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
    if(edep <= 0.) return true;

    // std::cout << "---- enter in MySensitiveDetector::ProcessHits  event = " << eventID << " trackID = " << trackID << std::endl;
    // std::cout << "  Step particle name: " << particleName << ", process name " << CreatorprocessName << std::endl;
    // std::cout << "  Step position: ("<<pos.x()<<", "<<pos.y()<<", "<<pos.z()<<"), stepE "<<edep<<std::endl;
    // std::cout << "  Cell volume position: ("<<volPos.x()<<", "<<volPos.y()<<", "<<volPos.z()<<") "<<std::endl;

    G4int cellID = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);
    //G4int boxID = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1);
    //G4int globalID = 1e8*boxID + cellID;
    G4int globalID = cellID;

    if (fStoreRawHits) {
      EcalHit* rawHit = nullptr;
      if(rawCellIDCol.find(cellID) != rawCellIDCol.end()){
        rawHit = (*fRawHitCollection)[rawCellIDCol[cellID]];
        rawHit->addEdep(edep);
        if(particleName=="opticalphoton"){
          if(CreatorprocessName=="Cerenkov") rawHit->addNphChren(1);
          else{ rawHit->addNphScint(1);}
        }
      }
      else{
        rawHit = new EcalHit(cellID);
        rawHit->setEdep(edep);
        rawHit->setPosition(volPos);
        if(particleName=="opticalphoton"){
          if(CreatorprocessName=="Cerenkov") rawHit->addNphChren(1);
          else{ rawHit->addNphScint(1);}
        }
        fRawHitCollection->insert(rawHit);
        rawCellIDCol.insert(std::make_pair(cellID, fRawHitIndex));
        fRawHitIndex++;
      }
    }

    // Mannually define the hit and cellID from global position 
    // WARNING: hard-coded to 1 mm x-y segmentation here!
    //G4int idx = (int)(pos.x()+500.);
    //G4int idy = (int)(pos.y()+500.);
    //G4int globalID = 1000*idy + idx; 


    // For attenuation & transverse light cross talk effect: each step create 5*5 hits.
    G4int idx = cellID % fCellIDBase;
    G4int idy = (cellID / fCellIDBase) % fCellIDBase;
    G4int idz = cellID / (fCellIDBase * fCellIDBase);
    //std::cout<<"  Raw cellID: "<<cellID<<", decoded cellID: "<<idx<<", "<<idy<<", "<<idz<<std::endl;
    //std::cout<<"  Make lookup table for neighbor cells "<<std::endl;

    // std::cout<<" Light response function parameters: "<<std::endl;
    // std::cout<<"  ResponseX0:  "<<fResponseX0<<std::endl;
    // std::cout<<"  Attenuation length: "<<fAttLength<<std::endl;
    // std::cout<<"  Response slope: "<<fResponseSlope<<std::endl;
    // std::cout<<"  Response intercept: "<<fResponseIntercept<<std::endl;
    // std::cout<<"  Reflect coeff: "<<fReflectCoeff<<std::endl;


    std::vector<G4int> cellIDvec;
    std::vector<G4double> responseVec; 
    std::vector<G4ThreeVector> neighborPosVec;
    cellIDvec.push_back(cellID);
    responseVec.push_back(AttenuatedResponse(pos.x(), pos.y(), volPos.x(), volPos.y()));
    for(int i=-2; i<=2; i++ ){
      if(idx+i<1 || idx+i>fIdMax) continue;
      for(int j=-2; j<=2; j++ ){
        if(idy+j<1 || idy+j>fIdMax) continue;
        if(i==0 && j==0) continue;

        G4int neighborID = idz * fCellIDBase * fCellIDBase + (idy+j) * fCellIDBase + idx+i;
        cellIDvec.push_back(neighborID);
        G4double neighborX = volPos.x() + i * fpitch;
        G4double neighborY = volPos.y() + j * fpitch;
        G4double neighborZ = volPos.z();
        neighborPosVec.push_back(G4ThreeVector(neighborX, neighborY, neighborZ));
        G4double response = AttenuatedResponse(pos.x(), pos.y(), neighborX, neighborY); 
        responseVec.push_back(response);

        // std::cout<<"    Neighbor ("<<i<<", "<<j<<"): cellID "<<neighborID<<", position ("
        // <<neighborX<<", "<<neighborY<<", "<<"), response (non-uniformed): "<<response<<std::endl;
      }
    }

    //Normalize the response vector to the central cell response.
    // G4double responseSum = std::accumulate(responseVec.begin(), responseVec.end(), 0.0);
    // if (responseSum > 0.) {
    //   for (auto &r : responseVec) {
    //     r /= responseSum;
    //   }
    // }



    //std::cout<<"  Create hits "<<std::endl;
    for(size_t ihit = 0; ihit < cellIDvec.size(); ihit++){
      G4int cellID_local = cellIDvec[ihit];
      G4double edep_att = edep * responseVec[ihit];
      // std::cout<<"  Cell ID "<<cellID_local<<", Effective En "<<edep_att<<std::endl;

      //G4int cellID_local = cellID;
      //G4double edep_att = edep;
      EcalHit* hit = nullptr;
      if(cellIDCol.find(cellID_local) != cellIDCol.end()){
        hit = (*fHitCollection)[cellIDCol[cellID_local]];
        hit->addEdep(edep_att);
        if(ihit==0 && particleName=="opticalphoton"){
          if(CreatorprocessName=="Cerenkov") hit->addNphChren(1); 
          else{ hit->addNphScint(1);}
        }
        // hit->addStep(pos.x(), pos.y(), pos.z(), edep_att);

      }
      else{
        hit = new EcalHit(cellID_local);
        hit->setEdep(edep_att);
        hit->setPosition(neighborPosVec[ihit]);

        if(ihit==0 && particleName=="opticalphoton"){
          if(CreatorprocessName=="Cerenkov") hit->addNphChren(1);
          else{ hit->addNphScint(1);}
        }      
        // hit->addStep(pos.x(), pos.y(), pos.z(), edep_att);

        fHitCollection->insert(hit);
        cellIDCol.insert(std::make_pair(cellID_local, fHitIndex));
        fHitIndex++;
      }

    }
    // std::cout <<"  Hit created. Current hit size: " << fHitCollection->entries() << std::endl;
    // std::cout << " exit from MySensitiveDetector::ProcessHits " << std::endl;
    return true;
}
