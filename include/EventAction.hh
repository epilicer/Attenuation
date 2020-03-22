
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

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;
class EventAction : public G4UserEventAction
{

public:
  EventAction(RunAction*);
  ~EventAction();

public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  void SumEvent(G4double m_Edep, G4double m_amassPre, G4double m_amassPos) 
       { Edep += m_Edep; amassPre += m_amassPre; amassPos += m_amassPos;}
  void SurfaceInfo(G4double m_ekin, G4double m_momz, G4int m_pdg, G4int m_ltrk) 
       {ekin=m_ekin; momz=m_momz; pdg=m_pdg; ltrk=m_ltrk ;}

private: 
  RunAction* runAct;

public:
  G4double  Edep;
  G4double  amassPre;
  G4double  amassPos;
  G4double  ekin;
  G4double  momz;
  G4int     pdg;
  G4int     ltrk;
};

#endif
