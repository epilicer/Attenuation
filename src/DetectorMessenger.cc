
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

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
//#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* detector)
  :Detector(detector)
{

  myDir = new G4UIdirectory("/mygeom/");
  myDir->SetGuidance("UI commands specific to my geometry.");
  
  setupDir = new G4UIdirectory("/mygeom/setup/");
  setupDir->SetGuidance("setup control.");

/*  detectorSizeZCmd = new G4UIcmdWithADoubleAndUnit("/mygeom/setup/detectorSizeZ",this);
  detectorSizeZCmd -> SetGuidance("Set detector thickness (along Z)");
  detectorSizeZCmd -> SetParameterName("Size",false);
  detectorSizeZCmd -> SetDefaultUnit("um");  
  detectorSizeZCmd -> SetUnitCandidates("um mm cm");  
  detectorSizeZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
*/

  detectorMaterCmd = new G4UIcmdWithAString("/mygeom/setup/detectorMaterial",this);
  detectorMaterCmd -> SetGuidance("Select material for the detector ( default = G4_Si )");
  detectorMaterCmd -> SetParameterName("choice",false);
  detectorMaterCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);

}

DetectorMessenger::~DetectorMessenger()
{ 
  //delete detectorSizeZCmd;
  delete detectorMaterCmd;
  delete setupDir;
  delete myDir;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  /*if( command == detectorSizeZCmd )
    { Detector -> SetDetectorSizeZThickness(detectorSizeZCmd -> GetNewDoubleValue(newValue));}*/

  if( command == detectorMaterCmd )
    { Detector -> SetDetectorMaterial(newValue);}
}

