#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GenericMessenger.hh"
#include "detector.hh"
#include "G4UserLimits.hh"
#include <string>  // ✅ Ensure std::string is included
#include "G4GenericMessenger.hh"  // ✅ Include Geant4 Messenger
#include "G4ThreeVector.hh"




class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();
    
    virtual G4VPhysicalVolume *Construct();

    //void SetTargetMaterial(G4String materialName); // ✅ Allows dynamic material selection

    G4Material* MakeZnWO4() ;
    G4Material* MakeBGO() ;
    G4Material* MakeFastFloat();
    G4Material* MakeMixture();

    void DefineMaterials();
    void DefineDim(G4GenericMessenger* fMessenger);
    void SetCMOSPosition(G4ThreeVector position);

    void SetRingPosition( G4double position );
    void SetPitchSize( G4double val ) { fiber_pitch = val; }
    void SetFiberNum( G4int num ) { fiber_nx = num; fiber_ny = num; }
    void SetZSeg( G4int num ) { nseg_z = num; }

    G4double GetPitchSize() const { return fiber_pitch; }
    G4int GetFiberNum() const {return fiber_nx; }

private:
    //std::string targetMaterial;  // ✅ Holds the selected material
    //G4String targetMaterialName;  // ✅ Stores the selected material
    G4GenericMessenger *fMessenger;  // ✅ Messenger for UI commands
    //G4GenericMessenger *fMessengerCMOS;
    //G4ThreeVector fCMOSPosition;
    //G4double gpositionRing ;
    //void DefineCommands();

    //Dimensions
    G4double worldSize;
    G4double xtal_x, xtal_y, xtal_half_length;
    G4double fiber_Radius, fiber_lendth, fiber_pitch;
    G4int fiber_nx, fiber_ny, nseg_z;
    G4double cladding_thick;
    bool fiber_fullfill = false;

    G4double pmtRadius;
    G4double pmtThickness;

    G4Box       *solidWorld;
    //G4Tubs    *solidSwYMCPsurface;
    //G4Box	*solidCMOSDetector;  

    G4LogicalVolume    *logicWorld;
    G4LogicalVolume   *logicCrystal;
    G4LogicalVolume   *logicCore;
    G4LogicalVolume   *logicClad;
    G4LogicalVolume   *logicPMT;


    G4VPhysicalVolume  *physWorld;
    G4VPhysicalVolume  **physCrystal;
    G4VPhysicalVolume  **physFiberClad;
    G4VPhysicalVolume  **physFiberCore;
    G4VPhysicalVolume  *physPMT;

    //G4VPhysicalVolume    *physSensor1 ;
    //G4VPhysicalVolume    *physSensor2 ;
    //G4VPhysicalVolume    *physSensor3 ;
    //G4VPhysicalVolume    *physSensor4 ;
    //G4VPhysicalVolume    *physSensor5 ;
    //G4VPhysicalVolume    *physSensor6 ;
    //G4VPhysicalVolume    *physSensor7 ;
    //G4VPhysicalVolume    *physSensor8 ;
    //G4VPhysicalVolume    *physSensor9 ;
    //G4VPhysicalVolume    *physSensor10 ;

    

    G4Material *matWorld ;   //  , *matMCPsurface, *matMCP;
    G4Material *matVacuum;  // ✅ Add this to class member variables
    //G4Material *matCMOSDetector;  // ✅ Add this to class member variables
    //G4Material *matSUS;  // ✅ Add this to class member variables

    G4Material *ZnWO4 ;
    G4Material *matCrystal ;
    
        // Vacuum tube/flange solids

//    G4Sphere *solidScoringTracker;
//    G4LogicalVolume *logicScoringTracker;
//    G4VPhysicalVolume *physScoringTracker;


    virtual void ConstructSDandField();
//    void CreateCommands();  // ✅ Function to create UI commands



};

#endif

