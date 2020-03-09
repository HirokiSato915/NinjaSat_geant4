//
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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// \file B3aActionInitialization.cc
/// \brief Implementation of the B3aActionInitialization class

#include "B3aActionInitialization.hh"
#include "B3PrimaryGeneratorAction.hh"
#include "B3aRunAction.hh"
#include "B3aEventAction.hh"
#include "B3SteppingAction.hh"
#include "B3DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aActionInitialization::B3aActionInitialization(B3DetectorConstruction *detConstruction)
    : G4VUserActionInitialization(),
      fDetConstruction(detConstruction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aActionInitialization::~B3aActionInitialization()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aActionInitialization::BuildForMaster() const
{
  SetUserAction(new B3aRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aActionInitialization::Build() const
{
  auto eventAction = new B3aEventAction;
  SetUserAction(new B3PrimaryGeneratorAction(eventAction));
  SetUserAction(new B3aRunAction);
  SetUserAction(eventAction);
  //auto trackingAction = new B3aTrackingAction(fDetConstruction);
  //SetUserAction(trackingAction);
  //SetUserAction(eventAction);
  SetUserAction(new B3aSteppingAction(fDetConstruction, eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
