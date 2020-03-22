
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
     Author  :  Ercan Pilicer
     Company :  Uludag University
**********************************************************************/

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"

class G4Box;
class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorMessenger;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction();

public: 

  void SetDetectorSizeZThickness (G4double );
//  void SetDetectorMaterial       (const G4String& );
  void SetDetectorMaterial       (G4String);
  
  virtual G4VPhysicalVolume* Construct();

  void UpdateGeometry();
  G4Material*        GetMaterial()   {return fMaterial;};

private:

 void DefineMaterials();

  G4VPhysicalVolume* ConstructVolumes();

  G4Box*             worldSolid;
  G4LogicalVolume*   worldLogic;
  G4VPhysicalVolume* worldPhysical;

  G4Box*             detectorSolid;
  G4LogicalVolume*   detectorLogical;
  G4VPhysicalVolume* detectorPhysical;

  G4Box*             dummySolid;
  G4LogicalVolume*   dummyLogical;
  G4VPhysicalVolume* dummyPhysical;

  DetectorMessenger* detectorMessenger; 

//  G4double detectorSizeZ;

  G4Material*           fMaterial;     
};
#endif

