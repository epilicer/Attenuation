
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

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4ios.hh"

#include <fstream>
#include <iomanip>

/////////////////////////////////////////////////////////////////////////////

SteppingAction::SteppingAction ( RunAction* run, EventAction* evt ) : 
runAction ( run ), eventAction ( evt )
{}

/////////////////////////////////////////////////////////////////////////////

SteppingAction::~SteppingAction()
{}

/////////////////////////////////////////////////////////////////////////////4
void SteppingAction::UserSteppingAction ( const G4Step* aStep )
{
    G4String procName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    runAction->CountProcesses ( procName );
	
//    std::ofstream File("st.out", std::ios::app | std::ios::out);
/*
	if ( aStep -> GetPostStepPoint() -> GetPhysicalVolume() != NULL ) {
File << aStep -> GetTrack() -> GetVolume() -> GetName()
 << "\t" << aStep -> GetTrack() -> GetNextVolume()-> GetName() 
 << "\t" << aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding() 
 << "\t" << aStep -> GetTrack() -> GetTrackID()
 << "\t" << aStep -> GetTrack() -> GetParentID()
 << "\t" << aStep -> GetTrack() -> GetDynamicParticle() -> GetKineticEnergy() / keV
 << "\t" << aStep -> GetPreStepPoint() -> GetKineticEnergy() / keV
 << "\t" << aStep -> GetPostStepPoint() -> GetKineticEnergy() / keV
 << "\t" << aStep -> GetPreStepPoint() -> GetPosition() / mm
 << "\t" << aStep -> GetPostStepPoint() -> GetPosition() / mm
 << "\t" << aStep -> GetPreStepPoint() -> GetMomentumDirection()
 << "\t" << aStep -> GetPostStepPoint() -> GetMomentumDirection()
 << G4endl;
	}
*/
	
	
/*File << aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding() 
 << "\t" << aStep -> GetTrack() -> GetTrackID()
 << "\t" << aStep -> GetTrack() -> GetDynamicParticle() -> GetKineticEnergy() / GeV
 << "\t" << aStep -> GetPreStepPoint() -> GetKineticEnergy() / GeV
 << "\t" << aStep -> GetPostStepPoint() -> GetKineticEnergy() / GeV
 << "\t" << aStep -> GetPostStepPoint() -> GetPosition().getX() / cm
 << "\t" << aStep -> GetPostStepPoint() -> GetPosition().getY() / cm
 << "\t" << aStep -> GetPostStepPoint() -> GetPosition().getZ() / cm
 << "\t" << aStep -> GetPostStepPoint() -> GetMomentumDirection().getX()
 << "\t" << aStep -> GetPostStepPoint() -> GetMomentumDirection().getY()
 << "\t" << aStep -> GetPostStepPoint() -> GetMomentumDirection().getZ()
 << G4endl;
*/
	
	if ( ( aStep -> GetPostStepPoint() -> GetPhysicalVolume() != NULL ) &&
         ( aStep -> GetTrack() -> GetVolume() -> GetName() == "Detector" ) &&
         ( aStep -> GetTrack() -> GetNextVolume()-> GetName() == "World") )
       {
           //counter primary particles exiting detector material
           if ( ( aStep -> GetTrack() -> GetTrackID() == 1 ) &&
			    ( aStep->GetPostStepPoint()->GetMomentumDirection().getZ() > 0.0 ) )
			    runAction -> AddFwdCountDetector();
           
           /*G4cout << aStep -> GetTrack() -> GetTrackID()
			   << "\t" << aStep -> GetTrack() -> GetParentID() 
               << "\t" << aStep -> GetTrack() -> GetDynamicParticle()->GetDefinition()->GetParticleName()
               << "\t" << aStep -> GetPreStepPoint() -> GetPosition()
               << "\t" << aStep -> GetPostStepPoint() -> GetPosition()
               << "\t" << aStep -> GetPostStepPoint() -> GetKineticEnergy() / keV
               << "\t" << aStep->GetPreStepPoint()->GetMomentumDirection()
               << "\t" << aStep->GetPostStepPoint()->GetMomentumDirection()
               << G4endl;*/
                
           eventAction-> SurfaceInfo ( aStep -> GetPostStepPoint() -> GetKineticEnergy() / keV,
                                       aStep -> GetPostStepPoint() -> GetMomentumDirection().getZ(),
                                       aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding(),
                                       aStep -> GetTrack() -> GetTrackID() );
       }


	   if ( ( aStep -> GetPostStepPoint() -> GetPhysicalVolume() != NULL ) &&
          ( aStep -> GetTrack() -> GetVolume() -> GetName() == "Detector" ) &&
          ( aStep -> GetTrack() -> GetNextVolume()-> GetName() == "Dummy") && 
		  ( aStep -> GetTrack() -> GetTrackID() == 1 ) )
		  runAction -> AddBwdCountDetector();     

       if ( ( aStep -> GetPostStepPoint() -> GetPhysicalVolume() != NULL ) && 
          ( aStep -> GetTrack() -> GetVolume() -> GetName() == "Detector" ) ) {

          G4double Ekpre = aStep -> GetPreStepPoint() -> GetKineticEnergy();
          G4double Ekpos = aStep -> GetPostStepPoint() -> GetKineticEnergy();

          G4Material* material = aStep -> GetTrack() -> GetMaterial();
          G4double density = aStep -> GetTrack() -> GetMaterial() -> GetDensity ();

		  G4ParticleDefinition* particle = aStep -> GetTrack() -> GetDefinition();
        
		  G4EmCalculator emCalculator;
		  
		  G4double amassPre = 0.0;
		  G4double amassPos = 0.0;
		  
		  G4double alenPre = 0.0;
		  G4double alenPos = 0.0;

		  if ( particle->GetParticleName() == "gamma" && procName != "Transportation" ) {

			  amassPre = emCalculator.ComputeCrossSectionPerVolume(Ekpre,particle,procName,material)/density;
              amassPos = emCalculator.ComputeCrossSectionPerVolume(Ekpos,particle,procName,material)/density;

			  alenPre = emCalculator.ComputeGammaAttenuationLength(Ekpre,material);
			  alenPos = emCalculator.ComputeGammaAttenuationLength(Ekpos,material);
/*
	          G4cout 
			  << "trkID : "      << aStep -> GetTrack() -> GetTrackID()
			  << "\tpdgID : "    << aStep -> GetTrack() -> GetParentID()
			  << "\tcurStp : "   << aStep -> GetTrack() -> GetCurrentStepNumber()
			  << "\tproc : "     << procName
			  << "\tpName : "    << particle->GetParticleName()
			  << "\tEkpre : "    << Ekpre
			  << "\tEkpos : "    << Ekpos
			  << "\tPospos : "   << aStep -> GetPostStepPoint() -> GetPosition()
			  << "\tMat : "      << material->GetName()
			  << "\tDens : "     << G4BestUnit(density,"Volumic Mass")
			  //<< "\talenPre : "  << alenPre / mm
			  << "\talenPre : "  << G4BestUnit(alenPre, "Length")
			  //<< "\talenPos : "  << alenPos / mm
			  << "\talenPos : "  << G4BestUnit(alenPos, "Length")
			  //<< "\tamassPre : " << amassPre / (mm2/g)
			  << "\tamassPre : " << G4BestUnit(amassPre, "Surface/Mass")
			  //<< "\tamassPos : " << amassPos / (mm2/g) 
			  << "\tamassPos : " << G4BestUnit(amassPos, "Surface/Mass")
			  << G4endl;
*/
			  
		  }
		
		   eventAction-> SumEvent (aStep->GetTotalEnergyDeposit() / keV, amassPre / (mm2/g), amassPos / (mm2/g)); 
  // kill event after first interaction
  //
  G4RunManager::GetRunManager()->AbortEvent();  
		
	   }

}
