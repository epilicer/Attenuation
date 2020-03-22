
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

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#include <fstream>
#include <iomanip>

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction ( RunAction* run ) : runAct ( run )
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction ( const G4Event* evt )
//void EventAction::BeginOfEventAction ( const G4Event* )
{
    G4int evtNb = evt->GetEventID();
    //printing survey
	if ( evtNb % 1 == 0 )
    //  G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
	
    Edep = amassPre = amassPos = ekin = momz = pdg = ltrk = 0;
}

void EventAction::EndOfEventAction ( const G4Event* evt )
{
    runAct->PerEvent (Edep);
    G4int evtNb = evt->GetEventID();
    /*G4cout << "---> End of Event: " << evtNb << "\t Edep " << Edep << "\t ekin " << ekin << "\t momz : " << momz 
		 << "\t pdg : " << pdg  << "\t ltrk : " << ltrk	<< "\n" << G4endl;*/

	std::ofstream File("ev.out", std::ios::app | std::ios::out);
    File << evtNb << "\t" << Edep << "\t" << ekin  << "\t" << momz << "\t" 
		<< pdg << "\t" << ltrk << "\t" << amassPre << "\t" << amassPos << G4endl;

/*		<<  G4BestUnit(amassPre, "Surface/Mass") << "\t"
		<<  G4BestUnit(amassPos, "Surface/Mass") << G4endl;*/
	
}

