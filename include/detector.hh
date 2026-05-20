#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
//#include "g4root.hh"
#include "G4AnalysisManager.hh" 
#include "G4SystemOfUnits.hh"
#include "EcalHit.hh"
#include "G4SDManager.hh"

class MySensitiveDetector : public G4VSensitiveDetector
{
public:
    MySensitiveDetector(G4String);
    ~MySensitiveDetector();
    void Initialize(G4HCofThisEvent *hitsCE) override;
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *) override;
    
private:
    EcalHitsCollection *fHitCollection; 
    G4int fHitCID; 
    G4int fHitIndex;
    G4double fAttLength; 
    G4double fpitch; 
    std::map<G4int, G4int> cellIDCol; // <cellID, hit index>
    
};

#endif

