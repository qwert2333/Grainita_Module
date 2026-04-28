#ifndef RUN_HH
#define RUN_HH

#include "G4UserRunAction.hh"
#include "G4Run.hh"
//#include "g4root.hh"
#include "G4AnalysisManager.hh" 
#include "generator.hh"
#include "TFile.h"
#include "TTree.h"

class MyRunAction : public G4UserRunAction
{
public:
    MyRunAction(MyPrimaryGenerator *PG);
    MyRunAction();
    ~MyRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    void ResetEventData();
    MyPrimaryGenerator* GetGenerator() const {return fPrimaryGen;}

    void Fill_Edep_Layer(G4double data) { Edep_Layer.push_back(data); }
    void Fill_Nph_Cherenkov_Layer(G4int data) { Nph_Cherenkov_Layer.push_back(data); }
    void Fill_Nph_Scint_Layer(G4int data) { Nph_Scint_Layer.push_back(data); }
    void Fill_vecCellID(G4int data) { vecCellID.push_back(data); }
    void Fill_vecEdep(G4double data) { vecEdep.push_back(data); }
    void Fill_vecNChren(G4int data) { vecNChren.push_back(data); }
    void Fill_vecNScint(G4int data) { vecNScint.push_back(data); }
    void Fill_stepPos(G4double _x, G4double _y, G4double _z) { 
      stepPosx.push_back(_x);
      stepPosy.push_back(_y);
      stepPosz.push_back(_z); 
    }
    void Fill_stepEn(G4double data) { stepEn.push_back(data); }

private: 
    MyPrimaryGenerator *fPrimaryGen;
    //TFile *fileRun;
    //TTree *treeEvt;

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
    
    std::vector<G4double> Edep_Layer;
    std::vector<G4int> Nph_Cherenkov_Layer;
    std::vector<G4int> Nph_Scint_Layer;
    std::vector<G4int> vecCellID;
    std::vector<G4double> vecEdep;
    std::vector<G4int> vecNChren;
    std::vector<G4int> vecNScint;
    std::vector<G4double> stepPosx;
    std::vector<G4double> stepPosy;
    std::vector<G4double> stepPosz;
    std::vector<G4double> stepEn;

};

#endif
