
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

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include <fstream>
#include <iomanip>
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"


RunAction::RunAction(DetectorConstruction* det)
:ProcCounter(0), fDetector(det)
{ }

RunAction::~RunAction()
{ }

void RunAction::BeginOfRunAction ( const G4Run* aRun )
{
    //G4RunManager::GetRunManager()->SetRandomNumberStore (false);
    G4cout << "Run " << aRun -> GetRunID() << " starts ..." << G4endl;

    //initialize cumulative quantities
    sumedep = sum2edep = 0.;
    counterFwdDet = counterBwdDet = 0;

	if (ProcCounter) delete ProcCounter;
	ProcCounter = new ProcessesCount;
}

void RunAction::CountProcesses(G4String procName)
{
  // is the process already counted ?
  // *AS* change to std::map?!
  size_t nbProc = ProcCounter->size();
  size_t i = 0;
  while ((i<nbProc)&&((*ProcCounter)[i]->GetName()!=procName)) i++;
  if (i == nbProc) ProcCounter->push_back( new OneProcessCount(procName));
  (*ProcCounter)[i]->Count();
}

void RunAction::EndOfRunAction ( const G4Run* aRun )
{

    NbOfEvents = aRun->GetNumberOfEvent();
    if ( NbOfEvents == 0 ) return;

    //compute statistics: mean and rms
    sumedep /= NbOfEvents;
    sum2edep /= NbOfEvents;
    G4double rmsedep = sum2edep - sumedep * sumedep;
    if ( rmsedep > 0. ) rmsedep = std::sqrt ( rmsedep );
    else rmsedep = 0.;

    G4cout << " Summary of Run " << aRun -> GetRunID() << G4endl;
	
  //frequency of processes
/*  G4cout << "\n---------------------------------------------------------------------------------------";
  G4cout << "\n Process calls frequency : \n";
*/

  for (size_t i=0; i< ProcCounter->size();i++) {
     G4String procName = (*ProcCounter)[i]->GetName();
     G4int    count    = (*ProcCounter)[i]->GetCounter();
  }

  G4Material* material = fDetector->GetMaterial();
  G4double density = material->GetDensity() / (g/cm3);
  G4double thickness = 5.0;

/*	          G4cout 
			  << "\tMat : "      << material->GetName()
//			  << "\tDens : "     << G4BestUnit(density,"Volumic Mass")
			  << "\tDens : "     << density
			  << G4endl;
*/
  
   G4double  massattenu = -1. * log(counterFwdDet/NbOfEvents) / (density*thickness) ;
//   G4double  massattenu = -1. * log(counterFwdDet/NbOfEvents) / ((density)*5.0) ;
//   G4double  massattenu = -1. * log(counterFwdDet/NbOfEvents) / 8.92;

    //print
/*    G4cout << "---------------------------------------------------------------------------------------\n"
    << " Number of Events                                            :   " << NbOfEvents     << "\n"
    << " Number of primary particles exiting from the material       :   " << counterFwdDet  << "\n"
    << " Number of primary particles backscattered from the material :   " << counterBwdDet  << "\n"
    << " mass attenuation                                            :   " << massattenu  << "\n"
    << " mean deposited energy in Material                           :   " << sumedep
    << " +- "                                                              << rmsedep        << " keV\n"
    << "\n---------------------------------------------------------------------------------------\n"
    << G4endl;
*/
   
    
    std::ofstream File("rn.out", std::ios::app | std::ios::out);
    File << NbOfEvents << " " << counterFwdDet << " " << counterBwdDet << " " << massattenu << G4endl;
}

void RunAction::PerEvent ( G4double edep)
{
    //accumulate statistic
    sumedep  += edep;
    sum2edep += edep * edep;

}

void RunAction::AddFwdCountDetector()
{
    counterFwdDet += 1;
}

void RunAction::AddBwdCountDetector()
{
    counterBwdDet += 1;
}

////////////////////////////////////////////////////////////
