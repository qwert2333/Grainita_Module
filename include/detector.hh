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
    G4double ResponseFunction(G4double distance) const;
    G4double AttenuatedResponse(G4double stepX, G4double stepY, G4double hitX, G4double hitY) const;

    EcalHitsCollection *fHitCollection; 
    G4int fHitCID; 
    G4int fHitIndex;
    G4double fAttLength; 
    G4double fpitch; 
    G4int fIdMax;
    G4double fResponseX0;
    G4double fResponseSlope;
    G4double fResponseIntercept;
    G4double fResponseNorm;
    G4double fReflectCoeff;
    std::map<G4int, G4int> cellIDCol; // <cellID, hit index>
    
};

#endif
