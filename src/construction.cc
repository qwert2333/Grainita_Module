#include "construction.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "CADMesh.hh"
#include "G4TessellatedSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UIcommand.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "detector.hh"
#include "G4UserLimits.hh"

#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4MultiUnion.hh"

#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

// Local parameters

MyDetectorConstruction::MyDetectorConstruction() : G4VUserDetectorConstruction()
{
    //Define a UI Messenger to allow runtime material selection
    fMessenger = new G4GenericMessenger(this, "/detector/", "Detector Control");

    DefineDim(fMessenger);
    DefineMaterials();

}



MyDetectorConstruction::~MyDetectorConstruction(){}



// =================================================================================
//
// ======                   Define dimensions
//
// =================================================================================


void MyDetectorConstruction::DefineDim(G4GenericMessenger* fMessenger)
{

    worldSize = 200 * cm;
    xtal_x = 20*cm;
    xtal_y = 20*cm;
    xtal_half_length = 21.04*cm;
    nbox_x = 4; 
    nbox_y = 4;

    nfiber_x = 14; 
    nfiber_y = 14;
    nseg_z = 25;
    fiber_Radius = 0.5*mm;
    fiber_lendth = 2.*xtal_half_length + 1.*cm; //1 cm interface
    cladding_thick = 0.3*mm;
    carbonframe_thick = 0.5*mm;

    pmtRadius = 80.*mm;
    pmtThickness = 1.*mm;

    fMessenger->DeclareMethodWithUnit("ModuleSize", "mm", &MyDetectorConstruction::SetModuleSize)
        .SetGuidance("Set module transverse size")
        .SetParameterName("moduleSize", true)
        .SetStates(G4State_PreInit, G4State_Idle);

    fMessenger->DeclareMethodWithUnit("ModuleDepth", "mm", &MyDetectorConstruction::SetModuleDepth)
        .SetGuidance("Set module depth")
        .SetParameterName("moduleDepth", true)
        .SetStates(G4State_PreInit, G4State_Idle);

    fMessenger->DeclareMethod("BoxNum", &MyDetectorConstruction::SetBoxNum)
        .SetGuidance("Set number of box in the module")
        .SetParameterName("boxNum", true)
        .SetStates(G4State_PreInit, G4State_Idle);

    fMessenger->DeclareMethod("fiberNum", &MyDetectorConstruction::SetFiberNum)
        .SetGuidance("Set number of fibers in one box")
        .SetParameterName("fiberNum", true)
        .SetStates(G4State_PreInit, G4State_Idle);

    fMessenger->DeclareMethod("ZSegNum", &MyDetectorConstruction::SetZSeg)
        .SetGuidance("Set Z direction segmentation")
        .SetParameterName("ZSegNum", true)
        .SetStates(G4State_PreInit, G4State_Idle);


}



// =================================================================================
//
// ======                   Define materials
//
// =================================================================================
 
G4Material* MyDetectorConstruction::MakeBGO() {
  // --- elements ---
  G4Element* elBi = new G4Element("Bismuth", "Bi", 83., 208.9804*g/mole);
  G4Element* elGe = new G4Element("Germanium", "Ge", 32., 72.63*g/mole);
  G4Element* elO  = new G4Element("Oxygen", "O", 8., 15.999*g/mole);

  // --- material: Bi4Ge3O12 ---
  G4double density = 7.13*g/cm3;   // typical
  G4Material* BGO = new G4Material("BGO", density, 3);
  BGO->AddElement(elBi, 4);
  BGO->AddElement(elGe, 3);
  BGO->AddElement(elO, 12);

  // --- optical / scintillation properties ---
  const G4int nEntries = 5;
  // Photon energy grid (eV) covering ~350 nm (3.543 eV) to ~650 nm (1.91 eV)
  G4double photonEnergy[nEntries] = {
      1.91*eV,   // ~650 nm
      2.33*eV,   // ~532 nm
      2.58*eV,   // ~480 nm (emission peak)
      3.10*eV,   // ~400 nm
      3.54*eV    // ~350 nm
  };

  // Refractive index: typical value ~2.15 at emission max; slight dispersion
  G4double rindex[nEntries] = {2.12, 2.14, 2.15, 2.16, 2.18};

  // Absorption length: BGO is relatively transparent near emission.
  // Use large values (tens to hundreds of cm). Replace with measured ABSLENGTH(λ) if available.
  G4double absLength[nEntries] = {
      200.*cm,  // 650 nm
      150.*cm,  // 532 nm
      100.*cm,  // 480 nm (conservative)
      80.*cm,   // 400 nm
      50.*cm    // 350 nm (UV absorption edge begins)
  };

  // Scintillation emission spectrum (relative) - simple approx peaked at 480 nm
  G4double scintSpectrum[nEntries] = {
      0.2,  // 650 nm tail
      0.6,  // 532 nm
      1.0,  // 480 nm peak
      0.4,  // 400 nm
      0.05  // 350 nm
  };

  auto mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", photonEnergy, rindex, nEntries);
  mpt->AddProperty("ABSLENGTH", photonEnergy, absLength, nEntries);
  mpt->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintSpectrum, nEntries);

  // Absolute light yield (typical reported values ~8000–10000 ph/MeV)
  //mpt->AddConstProperty("SCINTILLATIONYIELD", 8500.0/MeV); // choose 8500 ph/MeV as baseline
  mpt->AddConstProperty("SCINTILLATIONYIELD", 8.5/MeV); // choose 8500 ph/MeV as baseline
  mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);

  // Decay constants: main slow ~300 ns; a faster component (~60 ns) is often reported.
  mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 300.*ns); // main
  mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 60.*ns);  // fast (if present)
  // Relative contributions (normalize as you prefer; here ~80% slow, 20% fast)
  mpt->AddConstProperty("SCINTILLATIONYIELD1", 0.80);
  mpt->AddConstProperty("SCINTILLATIONYIELD2", 0.20);

  // Optionally add WLS (if you have a wavelength-shifter) or set Birks constant
  // BGO->GetIonisation()->SetBirksConstant(0.01*mm/MeV);

  BGO->SetMaterialPropertiesTable(mpt);

  return BGO;
}



G4Material* MyDetectorConstruction::MakeZnWO4() {

// --- define ZnWO4 crystal

  // Elements
  G4Element* O  = new G4Element("Oxygen", "O", 8., 15.999*g/mole);
  G4Element* Zn = new G4Element("Zinc", "Zn", 30., 65.38*g/mole);
  G4Element* W  = new G4Element("Tungsten", "W", 74., 183.84*g/mole);

  // Material
  ZnWO4 = new G4Material("ZnWO4", 7.87*g/cm3, 3);
  ZnWO4->AddElement(Zn, 1);
  ZnWO4->AddElement(W, 1);
  ZnWO4->AddElement(O, 4);

  // Optional optical properties
  G4MaterialPropertiesTable* mptZnWO4 = new G4MaterialPropertiesTable();

  // --- optical / scintillation properties ---
  const G4int nEntries = 3;
  G4double photonEnergy[nEntries] = {2.0*eV, 2.5*eV, 3.5*eV}; // 620 nm, 496 nm, 354 nm

  // refractive index (flat-ish approx in visible; use exact table if available)
   G4double rindex[nEntries] = {2.20, 2.20, 2.21};  // ~2.1-2.3 measured in literature


  // hack to count the Cerenkov for the refractive index of the heavy liquid
  //G4double rindex[nEntries] = {1.85, 1.85, 1.851};

  // emission spectrum (relative, normalized shape) -- peak near ~2.5 eV (~496 nm)
  G4double scintSpectrum[nEntries] = {0.2, 1.0, 0.05}; // approximate relative shape

  auto *mpt = new G4MaterialPropertiesTable();
  mpt->AddProperty("RINDEX", photonEnergy, rindex, nEntries);
  mpt->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintSpectrum, nEntries);
  // If you want multiple spectral components, add SCINTILLATIONCOMPONENT2 etc.

  // -----------------------------------------------------------------------------------------------

  // Light yield of ZnWO4 :

  // The true yield of ZnWO4 :
  // mpt->AddConstProperty("SCINTILLATIONYIELD", 9300.0/MeV); // Danevich et al. ~9300 ph/MeV

  // To switch off scintillation (for Cerenkov studies) : 
  // mpt->AddConstProperty("SCINTILLATIONYIELD", 0.);     // to switch off scintillation and have only cerenkov photons

  // trick: pass the actual yield divided by 1000. This way, the number of produced Scintillation
  // photons actually corresponds to the number of detected photo-electrons ( 10^4 pe per GeV, hence 0.1% of the
  // number of produced scintillation photons). Moreover, a smaller yield makes the simulation faster.

   mpt->AddConstProperty("SCINTILLATIONYIELD", 9.3/MeV);  

   mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);


  // -----------------------------------------------------------------------------------------------
  // --
  // -- The value that should be passed below, for  the absorption length in ZnWO4,  is not completely clear.
  // -- Anyway, as long as we do not really use the optical simulation, we do not care.

  //G4double absLenArray[nEntries] = { 0.2*cm, 0.2*cm, 0.2*cm } ;  // absorption
  G4double absLenArray[nEntries] = { 200*cm, 200*cm, 200*cm } ;  // ~ no absorption
  mpt->AddProperty("ABSLENGTH", photonEnergy, absLenArray, nEntries);



  // -----------------------------------------------------------------------------------------------
  //
  // ---------------   Pulse shapes 

  // -- Careful, what comes below has not been checked. In particular, I'm not sure if the "amplitudes"
  // -- that are passed should be the coefficient in front of the exponential, or the fraction
  // -- of the light that is emitted with thi time constant.
  // -- That said, as long as we do not really use the optical simulation, we do not care.

  // Multi-exponential decay constants measured by Shibata et al. (room T, example):
  //  tau1 = 24.6 μs, tau2 = 4.46 μs, tau3 = 0.34 μs (for gamma excitation)
  // We'll add three time constants and approximate relative yields (normalized).
  mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 24.6*us);
  mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 4.46*us);
  mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT3", 0.34*us);

  // Relative intensities: use the relative coefficients measured (Shibata table) and normalize them:
  // (these are approximate normalized fractions — replace with your crystal's measured fractions)
  G4double amp1 = 0.347; // long component fraction (normalized)
  G4double amp2 = 0.151; // mid
  G4double amp3 = 0.502; // fast
  mpt->AddConstProperty("SCINTILLATIONYIELD1", amp1);
  mpt->AddConstProperty("SCINTILLATIONYIELD2", amp2);
  mpt->AddConstProperty("SCINTILLATIONYIELD3", amp3);

  // -----------------------------------------------------------------------------------------------


  // tell Geant4 that scintillation is isotropic and produces photons per energy deposit
  ZnWO4->SetMaterialPropertiesTable(mpt);

  // set Birks constant to model quenching - same as Mike in his talk from end of January...
  ZnWO4->GetIonisation()->SetBirksConstant(0.01181*mm/MeV); // tune as needed


 return ZnWO4;

}


// --- description of the "FastFloat" heavy liquid

G4Material* MyDetectorConstruction::MakeFastFloat() {

    // --------------------------------------------------------------------
    // Elements
    // --------------------------------------------------------------------
    G4Element* H  = new G4Element("Hydrogen", "H",  1.,   1.008*g/mole);
    G4Element* O  = new G4Element("Oxygen",   "O",  8.,  16.00*g/mole);
    G4Element* Na = new G4Element("Sodium",   "Na", 11., 22.99*g/mole);
    G4Element* W  = new G4Element("Tungsten", "W",  74., 183.84*g/mole);

    // --------------------------------------------------------------------
    // Water (H2O)
    // --------------------------------------------------------------------
    G4Material* Water = new G4Material(
      "Water", 1.0*g/cm3, 2, kStateLiquid, 293.15*kelvin
    );
    Water->AddElement(H, 2);
    Water->AddElement(O, 1);

    // --------------------------------------------------------------------
    // Liquid mixture
    // --------------------------------------------------------------------
    G4double density = 2.8*g/cm3;

    G4Material* LiquidMixture = new G4Material(
      "W_Na_O_Water_Liquid",
      density,
      4,                // number of components
      kStateLiquid,
      293.15*kelvin
    );

    // Weight fractions
    LiquidMixture->AddMaterial(Water, 0.25); // 25 wt% water
    LiquidMixture->AddElement(W,      0.56); // 56 wt% tungsten
    LiquidMixture->AddElement(Na,     0.02); // 2 wt% sodium
    LiquidMixture->AddElement(O,      0.16); // 16 wt% additional oxygen


    const G4int NUM = 4;

    G4double photonEnergy[NUM] = {
      2.0*eV,  // ~620 nm
      2.5*eV,
      3.0*eV,
      3.5*eV   // ~355 nm
    };

    G4double refractiveIndex[NUM] = {
      1.451,
      1.512,
      1.548,
      1.55
    };

    G4double absorption[NUM] = {
      40.*cm,
      40.*cm,
      40.*cm,
      40.*cm,
    };

    G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

    MPT->AddProperty("RINDEX", photonEnergy, refractiveIndex, NUM);
    MPT->AddProperty("ABSLENGTH", photonEnergy, absorption, NUM);

    LiquidMixture ->SetMaterialPropertiesTable(MPT);

    return LiquidMixture;

}



// ----  for a mixture of ZnWO4 and FastFloat liquid :

G4Material* MyDetectorConstruction::MakeMixture()
{

    G4Material* grains = MakeZnWO4() ;
    G4Material* liquid = MakeFastFloat() ;

    G4Material* Mixture = new G4Material("Mixture", 5.72*g/cm3, 2);
    // the fractions corresponding to the prototype :
    Mixture -> AddMaterial( grains, 0.7923 );
    Mixture -> AddMaterial( liquid, 0.2077 );

    const G4int NUM = 4;

    G4double photonEnergy[NUM] = {
      2.0*eV,  // ~620 nm
      2.5*eV,
      3.0*eV,
      3.5*eV   // ~355 nm
    };

    G4double refractiveIndex[NUM] = {
      1.451,
      1.512,
      1.548,
      1.55
    };

    G4double absorption[NUM] = {    // dummy values here
      40.*cm,
      40.*cm,
      40.*cm,
      40.*cm,
    };

    G4double scintSpectrum[NUM] = {0.2, 1.0, 0.05, 0.0}; // approximate relative shape

    G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

    MPT->AddProperty("RINDEX", photonEnergy, refractiveIndex, NUM);
    MPT->AddProperty("ABSLENGTH", photonEnergy, absorption, NUM);

    //Scintillation properties: copy the main component from ZnWO4. 
    MPT->AddProperty("SCINTILLATIONCOMPONENT1", photonEnergy, scintSpectrum, NUM);
    MPT->AddConstProperty("SCINTILLATIONYIELD", 10./MeV); // Effective light yield from cosmic ray test
    MPT->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 24.6*us);
    MPT->AddConstProperty("SCINTILLATIONYIELD1", 1);
    MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);

    Mixture ->SetMaterialPropertiesTable(MPT);

    return Mixture;

}


void MyDetectorConstruction::DefineMaterials()
{
    G4NistManager *nist = G4NistManager::Instance();

    auto air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

    // Define vacuum for the world volume
    matVacuum = nist->FindOrBuildMaterial("G4_Galactic");
    
    //matWorld = matVacuum;
    matWorld = air ;


// --- define the material of the block :

    //matCrystal = MakeZnWO4() ;     // pure ZnWO4 crystal
    //matCrystal = MakeFastFloat() ;     // pure liquid
    matCrystal = MakeMixture() ;      // mixture ZNWO4 + heavy liquid

// for info : print the X0 of the mixed material :
    G4double X0 =  matCrystal ->GetRadlen();
    std::cout << " ----  X0 of material  = " <<  X0 / cm << " cm" << std::endl;


}



// Construct detector geometry
G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  static G4RotationMatrix* rotX = [](){
    G4RotationMatrix* m = new G4RotationMatrix();
    m->rotateX(90.*deg);
    return m;}();

    physFiberClad = new G4VPhysicalVolume*[nfiber_x*nbox_x*nfiber_y*nbox_y];
    physFiberCore = new G4VPhysicalVolume*[nfiber_x*nbox_x*nfiber_y*nbox_y];
    //physCrystal = new G4VPhysicalVolume*[nseg_z];
    for (G4int i = 0; i < nfiber_x*nbox_x*nfiber_y*nbox_y; i++) {
        physFiberClad[i] = nullptr;
        physFiberCore[i] = nullptr;
    }
    //for(G4int i=0; i<nseg_z; i++){
    //  physCrystal[i] = nullptr;
    //}

    //DefineMaterials();
    // Clean up previous geometries
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    G4NistManager* nist = G4NistManager::Instance();
    // WLS fiber core: polystyrene-like
    G4Material* polystyrene = new G4Material("Polystyrene", 1.05*g/cm3, 2);
    G4Element* elC = nist->FindOrBuildElement("C");
    G4Element* elH = nist->FindOrBuildElement("H");
    polystyrene->AddElement(elC, 8);
    polystyrene->AddElement(elH, 8);

    // Cladding: PMMA
    G4Material* PMMA = nist->FindOrBuildMaterial("G4_PLEXIGLASS"); // polymethylmethacrylate

    //Carbon fiber frame
    G4Material* carbonfiber = new G4Material("CarbonFiber", 1.466*g/cm3, 3);
    G4Element* elO = nist->FindOrBuildElement("O");    
    carbonfiber->AddElement(elC, 0.7941);
    carbonfiber->AddElement(elH, 0.0584);
    carbonfiber->AddElement(elO, 0.1475);

    // Optical properties for fiber core (WLS) - absorption in blue, emission in green
    // Copied here, but not needed for now. 
    /*
    const G4int nw = 6;
    G4double peW[nw] = {1.91*eV, 2.15*eV, 2.33*eV, 2.58*eV, 2.80*eV, 3.54*eV};
    // core refractive index ~1.59 (polystyrene), cladding ~1.49 -> total internal reflection
    G4double rIndexCore[nw] = {1.59,1.59,1.59,1.59,1.59,1.59};
    // ABSLENGTH of WLS dye: short in blue (absorb), long in emission region
    //G4double wlsAbs[nw] = {10.*cm, 5.*cm, 2.*cm, 50.*cm, 200.*cm, 300.*cm}; // example
    G4double wlsAbs[nw] = {300.*cm, 300.*cm, 300.*cm, 300.*cm, 300.*cm, 300.*cm}; // example NO ABSORPTION !!
    // WLS emission spectrum (normalized relative)
    G4double wlsEmit[nw] = {0.05, 0.2, 0.7, 1.0, 0.5, 0.05};
    auto mptCore = new G4MaterialPropertiesTable();
    mptCore->AddProperty("RINDEX", peW, rIndexCore, nw);
    mptCore->AddProperty("ABSLENGTH", peW, wlsAbs, nw);
    // WLS properties:
    mptCore->AddProperty("WLSABSLENGTH", peW, wlsAbs, nw);
    mptCore->AddProperty("WLSCOMPONENT", peW, wlsEmit, nw); // emission spectrum
    mptCore->AddConstProperty("WLSTIMECONSTANT", 8.*ns);   // re-emission time
    mptCore->AddConstProperty("WLSMEANNUMBERPHOTONS", 1.0); // typical: ~1 photon per absorbed
    polystyrene->SetMaterialPropertiesTable(mptCore);

    // Cladding refractive index slightly lower
    const G4int nc = nw;
    G4double rIndexClad[nc] = {1.49,1.49,1.49,1.49,1.49,1.49};
    auto mptClad = new G4MaterialPropertiesTable();
    mptClad->AddProperty("RINDEX", peW, rIndexClad, nc);
    PMMA->SetMaterialPropertiesTable(mptClad);
    */


    // Define World
    solidWorld = new G4Box("solidWorld", worldSize/2, worldSize/2, worldSize/2);
    logicWorld = new G4LogicalVolume(solidWorld, matWorld, "logicWorld");
    physWorld  = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);


    std::cout<<"--- Create the module: "<<xtal_x*2<<" * "<<xtal_y*2 << " * " << xtal_half_length*2 <<std::endl;
    // ------ Crystal box ------
    //auto solidCrystal = new G4Box("Crystal", xtal_x, xtal_y, xtal_half_length);
    G4double box_halfx = xtal_x/nbox_x; 
    G4double box_halfy = xtal_y/nbox_y; 
    G4double box_halfz = xtal_half_length; ; 
    G4Box* solidBoxOuter = new G4Box("OuterBox", box_halfx, box_halfy, box_halfz);
    G4Box* solidBoxInner = new G4Box("InnerBox", box_halfx-carbonframe_thick, box_halfy-carbonframe_thick, box_halfz);
    G4SubtractionSolid* solidBoxFrame = new G4SubtractionSolid("CarbonFrame", solidBoxOuter, solidBoxInner);
    logicBox = new G4LogicalVolume(solidBoxOuter, matWorld, "logicBox");
    logicFrame = new G4LogicalVolume(solidBoxFrame, carbonfiber, "logicFrame");
    logicFrame->SetVisAttributes(new G4VisAttributes(G4Colour(0.4, 0.4, 0.4)));
    physFrame = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicFrame, "physFrame", logicBox, false, 0);
    std::cout<<"---   Box size: "<<box_halfx*2<<" * "<<box_halfy*2<<" * "<<box_halfz*2<<std::endl;
    std::cout<<"---   Box number: "<<nbox_x<<" * "<<nbox_y<<std::endl;


    G4double cell_halfx = (box_halfx-carbonframe_thick)/nfiber_x;
    G4double cell_halfy = (box_halfy-carbonframe_thick)/nfiber_y;;
    G4double cell_halfL = xtal_half_length / nseg_z; 
    std::vector<G4ThreeVector> holePositions;

    G4Box* solidCrysCell = new G4Box("cell", cell_halfx, cell_halfy, cell_halfL);
    G4Tubs* solidHole = new G4Tubs("Hole", 0., fiber_Radius+cladding_thick, cell_halfL, 0., 360.*deg);
    // Subtract the hole
    G4SubtractionSolid* solidWithHole = new G4SubtractionSolid("cellWithHole", solidCrysCell, solidHole, 0, G4ThreeVector(0,0,0)); 

    // Create the logic volume cell
    logicCrystal = new G4LogicalVolume(solidWithHole, matCrystal, "cellWithHole");
    logicCrystal->SetVisAttributes(new G4VisAttributes(G4Colour(0, 0.75, 1.0)));


    // Loop place the cubes
    G4double fiber_pitch_x = 2*(box_halfx - carbonframe_thick)/nfiber_x;
    G4double fiber_pitch_y = 2*(box_halfy - carbonframe_thick)/nfiber_y;
    G4double x0 = - ( (nfiber_x-1)*fiber_pitch_x ) / 2.;
    G4double y0 = - ( (nfiber_y-1)*fiber_pitch_y ) / 2.;
    G4double z0 = - ( (nseg_z-1)*cell_halfL*2 ) / 2.;
    std::cout<<"---   In each box: fiber (cell) number "<<nfiber_x<<" * "<<nfiber_y<<std::endl;

    auto cellMax = std::to_string(std::max(std::max(nfiber_x, nfiber_y), nseg_z));
    G4int cellMaxCount = pow(10, cellMax.length());
    G4int copyNo = 0; 

    for(int iz=0; iz<nseg_z; iz++){
      G4double zPos = z0 + iz * cell_halfL*2;
      for(int iy=0; iy<nfiber_y; iy++){
        G4double yPos = y0 + iy * cell_halfy*2;
        copyNo = (iz + 1) * cellMaxCount * cellMaxCount + (iy + 1) * cellMaxCount;
        for(int ix=0; ix<nfiber_x; ix++){
          copyNo++;
          G4double xPos = x0 + ix * cell_halfx*2;
          G4ThreeVector position(xPos, yPos, zPos);
//std::cout<<"      cell position: "<<xPos<<", "<<yPos<<", "<<zPos<<", copyNo "<<copyNo<<std::endl;
          new G4PVPlacement(nullptr, position, logicCrystal,
                            "CrystalVoxel", logicBox, false, copyNo);
          if(iz==0){
            G4ThreeVector position_cent(xPos, yPos, 0);
            holePositions.push_back(position_cent);
          }
        }
      }
    }

    // ------ Implement the fibers in holes ------
    auto solidCore = new G4Tubs("Core", 0., fiber_Radius, fiber_lendth/2., 0.*deg, 360.*deg);
    auto solidClad = new G4Tubs("Clad", fiber_Radius, fiber_Radius + cladding_thick, fiber_lendth/2., 0.*deg, 360.*deg);
    logicCore = new G4LogicalVolume(solidCore, polystyrene, "CoreLogical");
    logicClad = new G4LogicalVolume(solidClad, PMMA, "CladLogical");

    for (size_t i=0;i<holePositions.size();++i) {
      G4ThreeVector pos = holePositions[i];
      // place cladding then core (so core is inner volume)
      physFiberClad[i] = new G4PVPlacement(nullptr, pos, logicClad, "FiberClad_" + std::to_string(i), logicBox, false, 0);
      physFiberCore[i] = new G4PVPlacement(nullptr, pos, logicCore, "FiberCore_" + std::to_string(i), logicBox, false, 0);
    }

    logicCore->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,1.0,0.0)));
    logicClad->SetVisAttributes(new G4VisAttributes(G4Colour(0.8,0.8,0.8)));


    //Place Boxes
    x0 = - ( (nbox_x-1)*box_halfx*2 ) / 2.;
    y0 = - ( (nbox_y-1)*box_halfy*2 ) / 2.;
    G4int boxCopyNo = 0;
    for(int iy=0; iy<nbox_y; iy++){
      G4double yPos = y0 + iy * box_halfy*2;
      for(int ix=0; ix<nbox_x; ix++){
        boxCopyNo++;
        G4double xPos = x0 + ix * box_halfx*2;
        G4ThreeVector position(xPos, yPos, 0);
        physFrame = new G4PVPlacement(0, position, logicBox, "physFrame", logicWorld, false, boxCopyNo);
      }
    }

    // Optical surface of fibers
    // Not needed now
    /*
    // Surface between Crystal and fiber cladding (dielectric-dielectric)
    G4OpticalSurface* surfCrystaltoClad = new G4OpticalSurface("Crystal_Clad_Surface");
    surfCrystaltoClad->SetType(dielectric_dielectric);
    //surfCrystaltoClad->SetFinish(ground); // or polished for good coupling
    surfCrystaltoClad->SetFinish(polished);
    surfCrystaltoClad->SetModel(unified);

    // Surface MPT (optional reflectivity - here assume small reflection)
    const G4int nsurf = 2;
    G4double peSurf[nsurf] = {2.0*eV, 3.5*eV};
    G4double refl[nsurf] = {0.02, 0.02};
    auto surfMPT = new G4MaterialPropertiesTable();
    surfMPT->AddProperty("REFLECTIVITY", peSurf, refl, nsurf);
    surfCrystaltoClad->SetMaterialPropertiesTable(surfMPT);

    // Attach logical border surfaces Crystal <-> Clad
    // We must attach borders between the placed volumes -- as we placed cladding in world
    // we need to create a border surface between each Crystal physical volume and its fiber cladding
    // For simplicity we create border surfaces for each fiber placement (loop)
    for (int i=0; i< (int)holePositions.size(); ++i) {
      // We need the physical volumes: Crystal phys vol is physCrystal (single), fiber placement phys volumes are in order
      // In this simplified example we didn't keep pointers to the PVs; in production keep them.
      // Here we assume names and retrieve from PV store (not shown). Instead show the intended call:
      // new G4LogicalBorderSurface("Crystal_to_FiberSurf", physCrystal, physFiberClad_i, surfCrystaltoClad);
      // For clarity: you'd attach a border surface between the physical Crystal and the physical Clad PVs.
       new G4LogicalBorderSurface("Crystal_to_FiberSurf",physCrystal, physFiberClad[i], surfCrystaltoClad );
    }
    */

    // ----------------- PMT / photon collector at one fiber end -----------------
    // Create small detector volume at +Z end of fibers to collect photons
    // (Place a small disk covering the fiber bundle area)
    auto solidPMT = new G4Tubs("PMT", 0., pmtRadius, pmtThickness/2., 0.*deg, 360.*deg);
    logicPMT = new G4LogicalVolume(solidPMT, nist->FindOrBuildMaterial("G4_AIR"), "PMTLogical");
    // Place PMT at some position slightly beyond the fiber end (z = +bz + margin)
    G4double pmtZ = (fiber_lendth/2.) + 0.01*mm;
    physPMT = new G4PVPlacement(nullptr, G4ThreeVector(0.,0., pmtZ), logicPMT, "PMTPhys", logicWorld, false, 0);
    logicPMT->SetVisAttributes(new G4VisAttributes(G4Colour(0.0,0.0,1.0)));

    return physWorld;
}

// Construct sensitive detector
void MyDetectorConstruction::ConstructSDandField()
{
    auto sdManager = G4SDManager::GetSDMpointer();
    MySensitiveDetector *sensDet_crystal = new MySensitiveDetector("CrystalModule");
    sdManager->AddNewDetector(sensDet_crystal);
    logicCrystal->SetSensitiveDetector(sensDet_crystal);

    MySensitiveDetector *sensDet_fibercore = new MySensitiveDetector("FiberCore");
    sdManager->AddNewDetector(sensDet_fibercore);
    logicCore->SetSensitiveDetector(sensDet_fibercore);

    MySensitiveDetector *sensDet_fiberclad = new MySensitiveDetector("FiberCladding");
    sdManager->AddNewDetector(sensDet_fiberclad);
    logicClad->SetSensitiveDetector(sensDet_fiberclad);    

    MySensitiveDetector *sensDet_carbonframe = new MySensitiveDetector("CarbonFrame");
    sdManager->AddNewDetector(sensDet_carbonframe);
    logicFrame->SetSensitiveDetector(sensDet_carbonframe);

}

