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
/// \file B3A.cc
/// \brief Implementation of the B3PrimaryGeneratorAction class

#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "B3Constants.hh"
#include "B3aEventAction.hh"
#include "G4RotationMatrix.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4ThreeVector.hh"
#include <iostream>
#include <string> // useful for reading and writing

#include <fstream> // ifstream, ofstream

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
B3PrimaryGeneratorAction::B3PrimaryGeneratorAction(
    B3aEventAction *eventAction)
    : G4VUserPrimaryGeneratorAction(),
      fEventAction(eventAction)
{
  particleGun = new G4GeneralParticleSource();
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *particle = particleTable->FindParticle("gamma");
  //....PARTICLE DEFINITIONS
  particleGun->SetParticleDefinition(particle);

  // DEFINE A Pow-ENERGETIC SOURCE
  G4SPSEneDistribution *eneDist = particleGun->GetCurrentSource()->GetEneDist();
  eneDist->SetEnergyDisType("Pow");
  eneDist->SetAlpha(-1.29);
  eneDist->SetEmax(200 * keV);
  eneDist->SetEmin(1 * keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
B3PrimaryGeneratorAction::~B3PrimaryGeneratorAction()
{
  delete particleGun;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void B3PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  G4double targetR = 5.4 * cm;
  G4double length = 8. * cm;

  G4ThreeVector dir1 = G4RandomDirection();
  G4ThreeVector dir2 = dir1.orthogonal();
  dir2.setR(targetR * sqrt(G4UniformRand()));
  dir2.rotate(2 * CLHEP::pi * G4UniformRand(), dir1);

  G4ThreeVector center(0, 0, 7.5 * mm);
  G4ThreeVector position = center + dir1 * length + dir2;
  G4ThreeVector mom = -dir1;

  G4SPSPosDistribution *posDist = particleGun->GetCurrentSource()->GetPosDist();
  G4SPSAngDistribution *angDist = particleGun->GetCurrentSource()->GetAngDist();

  posDist->SetPosDisType("Point");
  posDist->SetCentreCoords(position);
  angDist->SetParticleMomentumDirection(mom);

  particleGun->GeneratePrimaryVertex(anEvent);
  G4double initEnergy = particleGun->GetParticleEnergy();
  fEventAction->AddEn(initEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......