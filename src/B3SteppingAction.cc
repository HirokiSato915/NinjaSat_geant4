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
/// \file B3SteppingAction.cc
/// \brief Implementation of the B3StackingAction class

#include "B3SteppingAction.hh"
#include "B3aEventAction.hh"
#include "B3DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aSteppingAction::B3aSteppingAction(
    const B3DetectorConstruction *detectorConstruction,
    B3aEventAction *eventAction)
    : G4UserSteppingAction(),
      fDetConstruction(detectorConstruction),
      fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aSteppingAction::~B3aSteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aSteppingAction::UserSteppingAction(const G4Step *step)
{
  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  //auto momentum = step->GetPreStepPoint()->GetMomentumDirection();

  // energy deposit
  if (volume != fDetConstruction->GetSensitivePV())
    return;

  auto momentum = step->GetTrack()->GetMomentumDirection();
  G4double theta = momentum.cosTheta(), phi = momentum.phi();
  //auto Kenergy = step->GetTrack()->GetKineticEnergy();
  std::ofstream ofss("cosT.txt", std::ios::app);
  ofss << theta << std::endl;
  std::ofstream ofsp("phi.txt", std::ios::app);
  ofsp << phi << std::endl;

  /*
  auto edep = step->GetTotalEnergyDeposit();
  fEventAction->AddSens(edep);
*/
  //fEventAction->AddSens(1);
  //fEventAction->AddSens(momentum);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......