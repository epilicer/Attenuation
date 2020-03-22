
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

/**********************************************************************
     Version :  Geant4.10.0.p02
     Created :  25/10/2014
     Author  :  Belgin Pilicer & Ercan Pilicer
     Company :  Uludag University
**********************************************************************/

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4RunManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include "G4UIcmdWithAString.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

/////////////////////////////////////////////////////////////////////////////

DetectorConstruction::DetectorConstruction()
   : worldSolid(0), worldLogic(0), worldPhysical(0), 
     detectorSolid(0), detectorLogical(0), detectorPhysical(0),
     dummySolid(0), dummyLogical(0), dummyPhysical(0), fMaterial(0)
{
  // Detector sizes
//  detectorSizeZ 	= 1.0 *cm;

  DefineMaterials();
  SetDetectorMaterial("Silicon");

  // Messenger to change parameters of the geometry
  detectorMessenger = new DetectorMessenger(this);
}

/////////////////////////////////////////////////////////////////////////////

DetectorConstruction::~DetectorConstruction()
{ delete detectorMessenger;}

/////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* DetectorConstruction::Construct()
{ return ConstructVolumes();}

/////////////////////////////////////////////////////////////////////////////
    
void DetectorConstruction::DefineMaterials()
{
  //
  // define few Elements by hand
  //
  G4double a, z;
    
  G4Element* H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"   ,"C" , z= 6., a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium"   ,"Na", z=11., a=  22.99*g/mole);
  G4Element* Mg = new G4Element("Magnesium","Mg", z=12., a=  24.31*g/mole);
  G4Element* Si = new G4Element("Silicon"  ,"Si", z=14., a=  28.09*g/mole);
  G4Element* P  = new G4Element("Phosphor" ,"P" , z=15., a=  30.97*g/mole);
  G4Element* S  = new G4Element("Sulfur"   ,"S" , z=16., a=  32.07*g/mole);
  G4Element* Cl = new G4Element("Chlorine" ,"Cl", z=17., a=  35.45*g/mole);
  G4Element* K  = new G4Element("Potasium" ,"P" , z=19., a=  39.10*g/mole);
  G4Element* Ca = new G4Element("Calsium"  ,"Ca", z=20., a=  40.08*g/mole);
  G4Element* Fe = new G4Element("Iron"     ,"Fe", z=26., a=  55.85*g/mole);
  G4Element* Zn = new G4Element("Zinc"     ,"Zn", z=30., a=  65.39*g/mole);

  //
  // define materials
  //
  G4double density;
  G4int ncomponents, natoms;
  G4double fractionmass;  

  // water with ionisation potential 78 eV
  G4Material* H2O = 
  new G4Material("Water", density= 1.00*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

 G4Material* BLOOD =
 new G4Material("BLOOD", density= 1.069*g/cm3, ncomponents=14);
 BLOOD->AddElement(H, fractionmass=0.101866);
 BLOOD->AddElement(C,  fractionmass=0.100020);
 BLOOD->AddElement(N,  fractionmass=0.029640);
 BLOOD->AddElement(O,  fractionmass=0.759414);
 BLOOD->AddElement(Na,  fractionmass=0.001850);
 BLOOD->AddElement(Mg,  fractionmass=0.00004);
 BLOOD->AddElement(Si,  fractionmass=0.00003);
 BLOOD->AddElement(P,  fractionmass=0.00035);
 BLOOD->AddElement(S,  fractionmass=0.00185);
 BLOOD->AddElement(Cl,  fractionmass=0.00278);
 BLOOD->AddElement(K,  fractionmass=0.00163);
 BLOOD->AddElement(Ca,  fractionmass=0.00006);
 BLOOD->AddElement(Fe,  fractionmass=0.00046);
 BLOOD->AddElement(Zn,  fractionmass=0.00001);


G4Material* BONECOMP =
 new G4Material("BONECOMP", density= 1.85*g/cm3, ncomponents=8);
 BONECOMP->AddElement(H, fractionmass=0.06398400);
 BONECOMP->AddElement(C,  fractionmass=0.27800000);
 BONECOMP->AddElement(N,  fractionmass=0.02700000);
 BONECOMP->AddElement(O,  fractionmass=0.41001600);
 BONECOMP->AddElement(Mg,  fractionmass=0.00200000);
 BONECOMP->AddElement(P,  fractionmass=0.07000000);
 BONECOMP->AddElement(S,  fractionmass=0.00200000);
 BONECOMP->AddElement(Ca,  fractionmass=0.14700000);


G4Material* LUNG =
 new G4Material("LUNG", density= 1.05*g/cm3, ncomponents=13);
 LUNG->AddElement(H ,  fractionmass=0.10127800);
 LUNG->AddElement(C ,   fractionmass=0.10231000);
 LUNG->AddElement(N ,   fractionmass=0.02865000);
 LUNG->AddElement(O ,   fractionmass=0.75707200);
 LUNG->AddElement(Na,  fractionmass=0.00184000);
 LUNG->AddElement(Mg,  fractionmass=0.00073000);
 LUNG->AddElement(P ,   fractionmass=0.00080000);
 LUNG->AddElement(S ,   fractionmass=0.00225000);
 LUNG->AddElement(Cl,  fractionmass=0.00266000);
 LUNG->AddElement(K ,   fractionmass=0.00194000);
 LUNG->AddElement(Ca,  fractionmass=0.00009000);
 LUNG->AddElement(Fe,  fractionmass=0.00037000);
 LUNG->AddElement(Zn,  fractionmass=0.00001000);

G4Material* EYELENS =
 new G4Material("EYELENS", density= 1.1*g/cm3, ncomponents=4);
 EYELENS->AddElement(H, fractionmass=0.099269);
 EYELENS->AddElement(C,  fractionmass=0.193710);
 EYELENS->AddElement(N,  fractionmass=0.053270);
 EYELENS->AddElement(O,  fractionmass=0.653751);


G4Material* Adipose_ =
 new G4Material("Adipose_", density= 0.92*g/cm3, ncomponents=8);
 Adipose_->AddElement(H ,  fractionmass=0.11947700);
 Adipose_->AddElement(C ,   fractionmass=0.63724000);
 Adipose_->AddElement(N ,   fractionmass=0.00797000);
 Adipose_->AddElement(O ,   fractionmass=0.23233300);
 Adipose_->AddElement(Na,  fractionmass=0.00050000);
 Adipose_->AddElement(Mg,  fractionmass=0.00002000);
 Adipose_->AddElement(P ,   fractionmass=0.00016000);
 Adipose_->AddElement(Cl,  fractionmass=0.00119000);

 G4Material* MuscleEl =
 new G4Material("MuscleEl", density= 1.04*g/cm3, ncomponents=8);
 MuscleEl->AddElement(H, fractionmass=0.10063700);
 MuscleEl->AddElement(C,  fractionmass=0.10783000);
 MuscleEl->AddElement(N,  fractionmass=0.02768000);
 MuscleEl->AddElement(O,  fractionmass=0.75477300);
 MuscleEl->AddElement(Na,  fractionmass=0.00075000);
 MuscleEl->AddElement(Mg,  fractionmass=0.00019000);
 MuscleEl->AddElement(P,  fractionmass=0.00180000);
 MuscleEl->AddElement(S,  fractionmass=0.00241000);
 

 G4Material* TISSUEIC =
 new G4Material("TISSUEIC", density= 1.0*g/cm3, ncomponents=13);
 TISSUEIC->AddElement(H, fractionmass=0.10447200);
 TISSUEIC->AddElement(C,  fractionmass=0.23219000);
 TISSUEIC->AddElement(N,  fractionmass=0.02488000);
 TISSUEIC->AddElement(O,  fractionmass=0.63023800);
 TISSUEIC->AddElement(Na,  fractionmass=0.00113000);
 TISSUEIC->AddElement(Mg,  fractionmass=0.00013000);
 TISSUEIC->AddElement(P,  fractionmass=0.00133000);
 TISSUEIC->AddElement(S,  fractionmass=0.00199000);
 TISSUEIC->AddElement(Cl,  fractionmass=0.00134000);
 TISSUEIC->AddElement(K,  fractionmass=0.00199000);
 TISSUEIC->AddElement(Ca,  fractionmass=0.00023000);
 TISSUEIC->AddElement(Fe,  fractionmass=0.00005000);
 TISSUEIC->AddElement(Zn,  fractionmass=0.00003000);

 G4Material* BRAIN =
 new G4Material("BRAIN", density= 1.039*g/cm3, ncomponents=13);
 BRAIN->AddElement(H, fractionmass=0.11066700);
 BRAIN->AddElement(C,  fractionmass=0.12542000);
 BRAIN->AddElement(N,  fractionmass=0.01328000);
 BRAIN->AddElement(O,  fractionmass=0.73772300);
 BRAIN->AddElement(Na,  fractionmass=0.00184000);
 BRAIN->AddElement(Mg,  fractionmass=0.00015000);
 BRAIN->AddElement(P,  fractionmass=0.00354000);
 BRAIN->AddElement(S,  fractionmass=0.00177000);
 BRAIN->AddElement(Cl,  fractionmass=0.00236000);
 BRAIN->AddElement(K,  fractionmass=0.00310000);
 BRAIN->AddElement(Ca,  fractionmass=0.00009000);
 BRAIN->AddElement(Fe,  fractionmass=0.00005000);
 BRAIN->AddElement(Zn,  fractionmass=0.00001000);

 G4Material* SKIN_EL =
 new G4Material("SKIN_EL", density= 1.1*g/cm3, ncomponents=13);
 SKIN_EL->AddElement(H, fractionmass=0.10058800);
 SKIN_EL->AddElement(C,  fractionmass=0.22825000);
 SKIN_EL->AddElement(N,  fractionmass=0.04642000);
 SKIN_EL->AddElement(O,  fractionmass=0.61900200);
 SKIN_EL->AddElement(Na,  fractionmass=0.00007000);
 SKIN_EL->AddElement(Mg,  fractionmass=0.00006000);
 SKIN_EL->AddElement(P,  fractionmass=0.00033000);
 SKIN_EL->AddElement(S,  fractionmass=0.00159000);
 SKIN_EL->AddElement(Cl,  fractionmass=0.00267000);
 SKIN_EL->AddElement(K,  fractionmass=0.00085000);
 SKIN_EL->AddElement(Ca,  fractionmass=0.00015000);
 SKIN_EL->AddElement(Fe,  fractionmass=0.00001000);
 SKIN_EL->AddElement(Zn,  fractionmass=0.00001000);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

/////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // -----------------------------
  //  World volume
  //------------------------------
  
  const G4double worldX = 50.0 *cm;
  const G4double worldY = 50.0 *cm;
  const G4double worldZ = 40.0 *cm;
  
  //G4bool isotopes = false;
 
  G4Material* voidMat =  G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  worldSolid = new G4Box("World",worldX/2,worldY/2,worldZ/2);
  worldLogic = new G4LogicalVolume(worldSolid, voidMat, "World",0,0,0);
  worldPhysical = new G4PVPlacement(0,
					    G4ThreeVector(),
					    "World", 
					    worldLogic, 
					    0,false,0);
 
  //-----------
  // Dummy volume before detector
  //-----------

  const G4double dummyX = 50.0 *cm; //= worldX;
  const G4double dummyY = 50.0 *cm; //= worldY;
  const G4double dummyZ = 20.0 *cm; //= worldZ/2;

  G4Material* dummyMat =  G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  
  dummySolid = new G4Box("Dummy",dummyX/2,dummyY/2,dummyZ/2);
  dummyLogical = new G4LogicalVolume(dummySolid, dummyMat, "Dummy",0,0,0); 
  dummyPhysical = new G4PVPlacement(0,
					     G4ThreeVector(0, 0, -dummyZ/2),
					     "Dummy",
					     dummyLogical,
					     worldPhysical,
					     false,0);

  //-----------
  // Detector
  //-----------

  const G4double detectorSizeX = 50.0 *cm; //= worldX;
  const G4double detectorSizeY = 50.0 *cm; //= worldY;
  const G4double detectorSizeZ 	= 5.0 *cm;

  detectorSolid = new G4Box("Detector",detectorSizeX/2,detectorSizeY/2,detectorSizeZ/2);
  detectorLogical = new G4LogicalVolume(detectorSolid, fMaterial,"Detector",0,0,0); 
  detectorPhysical = new G4PVPlacement(0,G4ThreeVector(0, 0, detectorSizeZ/2),"Detector",detectorLogical,worldPhysical,false,0);
  
 
    G4VisAttributes* logicDetector_VisAtt = new G4VisAttributes ( G4Colour ( 1.0, 0.5, 0.5 ) );
    logicDetector_VisAtt -> SetVisibility ( true );
    logicDetector_VisAtt -> SetForceWireframe ( true );
    detectorLogical->SetVisAttributes ( logicDetector_VisAtt );

    //--------------------------------------------------------------
    // print some parameters
    //--------------------------------------------------------------
    G4cout 
    << "-----------------------------------------------------------------------"
    << G4endl
    << "Z: " << detectorSizeZ << G4endl
    << "detectorMat  : " << fMaterial->GetName() << G4endl
    << G4endl
    << "-----------------------------------------------------------------------"
    << G4endl;

   // always return the physical World
   return worldPhysical;
}

/////////////////////////////////////////////////////////////////////////////
/*
void DetectorConstruction::SetDetectorSizeZThickness(G4double value){
  detectorSizeZ=value;
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
*/
/////////////////////////////////////////////////////////////////////////////

/*
void DetectorConstruction::SetDetectorMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    detectorMat = pttoMaterial;
    if(detectorLogical) detectorLogical->SetMaterial(detectorMat);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}
*/
void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) {
   fMaterial = pttoMaterial;
   UpdateGeometry();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;  
  } 
}


/////////////////////////////////////////////////////////////////////////////

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}


