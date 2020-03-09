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
/// \file B3aEventAction.cc
/// \brief Implementation of the B3aEventAction class

#include "B3aEventAction.hh"
#include "B3aRunAction.hh"
#include "B3Analysis.hh"
#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//#include "G4SDManager.hh"
//#include "G4HCofThisEvent.hh"
//#include "G4THitsMap.hh"

//#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::B3aEventAction()
    : G4UserEventAction(),
      fEnergySens(0.),
      fEnergy(0.)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3aEventAction::~B3aEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::BeginOfEventAction(const G4Event * /*event*/)
{
  //std::ofstream ofsp("particle.txt", std::ios::app);
  //ofsp << fEnergy << std::endl;
  // initialisation per event
  fEnergySens = 0.;
  //fEnergy = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3aEventAction::EndOfEventAction(const G4Event * /*event*/)
{
  // Accumulate statistics
  //
  //if (fEnergySens != 0)
  //{
  /*
    std::ofstream ofsp("particle.txt", std::ios::app);
    ofsp << fEnergy << std::endl;
    */

  //auto analysisManager = G4AnalysisManager::Instance();

  /*
    analysisManager->FillH1(0, fEnergy);
    analysisManager->FillH1(1, fEnergySens);
    */

  // fill ntuple
  //analysisManager->FillNtupleDColumn(0, fEnergy);
  //analysisManager->FillNtupleDColumn(1, fEnergySens);
  //analysisManager->AddNtupleRow();

  //}

  // get analysis manager
  //std::cout << initEnergy << std::endl;
  // fill histograms

  //fEnergy = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
