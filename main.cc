#include <iostream>
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "construction.hh"
#include "FTFP_BERT.hh"
#include "G4VModularPhysicsList.hh"
#include "G4StepLimiterPhysics.hh"
#include "physics.hh"
#include "action.hh"
#include "G4SystemOfUnits.hh"
#include <ctime>
#include <string>

int main(int argc, char** argv) {
    G4long seed = time(NULL);
    G4Random::setTheSeed(seed);

#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    // Setup initialization
    runManager->SetUserInitialization(new MyDetectorConstruction());
    G4VModularPhysicsList* physicsList = new FTFP_BERT();
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);
    runManager->SetUserInitialization(new MyActionInitialization());

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    G4String macroFile = "run.mac";
    if (argc >= 2) macroFile = argv[1];

    // Execute run command
    UImanager->ApplyCommand("/control/execute " + macroFile);

    // Visualization option
    bool visualization = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-vis") visualization = true;
    }
    if (visualization) {
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        G4VisManager* visManager = new G4VisExecutive();
        visManager->Initialize();
        ui->SessionStart();
        delete visManager;
        delete ui;
    }

    delete runManager;
    return 0;
}


