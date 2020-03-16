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
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Para.hh"
#include "G4Orb.hh"
#include "G4AssemblyVolume.hh"
#include "G4Polyhedra.hh"
#include "G4UnionSolid.hh"
#include "G4ExtrudedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "B3Constants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::B3DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fSensitivePV(0),
      fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B3DetectorConstruction::Construct()
{
    //
    // Material
    //
    G4double mmN = 14.01 * g / mole;
    G4Element *elN = new G4Element("Nitrogen", "N", 7., mmN);
    G4double mmO = 16.00 * g / mole;
    G4Element *elO = new G4Element("Oxygen", "O", 8., mmO);

    G4Material *Vacuum = new G4Material("Vacuum", universe_mean_density, 2);
    Vacuum->AddElement(elN, .7);
    Vacuum->AddElement(elO, .3);

    G4double denXe = 5.858 * mg / cm3;
    G4double mmXe = 131.29 * g / mole;
    G4Material *Xe = new G4Material("Xenon", 54., mmXe, denXe);

    G4double denC = 12.01 * g / mole;
    G4Element *elC = new G4Element("Carbon", "C", 6., denC);

    G4double denCO2 = 1.977 * mg / cm3;
    G4Material *CarbonDioxide = new G4Material("CO2", denCO2, 2);
    CarbonDioxide->AddElement(elC, 1);
    CarbonDioxide->AddElement(elO, 2);

    G4double denXe5CO2 = 6.8 * mg / cm3;
    G4Material *Xe5CO2 = new G4Material("Xe5CO2", denXe5CO2, 2);
    Xe5CO2->AddMaterial(Xe, 0.9826);
    Xe5CO2->AddMaterial(CarbonDioxide, 0.0174);

    G4double denAl = 2.70 * g / cm3;
    G4double mmAl = 26.98 * g / mole;
    G4Material *Al = new G4Material("Al", 13., mmAl, denAl);

    G4double denSn = 7.365 * g / cm3;
    G4double mmSn = 118.710 * g / mole;
    G4Material *Sn = new G4Material("Sn", 50., mmSn, denSn);

    G4double denPb = 11.34 * g / cm3;
    G4double mmPb = 207.2 * g / mole;
    G4Material *Pb = new G4Material("Pb", 82., mmPb, denPb);

    G4Element *elBe = new G4Element("Berillium", "Be", 4., 9.01 * g / mole);
    G4Material *Be = new G4Material("Be", 1.848 * g / cm3, 1);
    Be->AddElement(elBe, 1);

    G4NistManager *nist = G4NistManager::Instance();
    G4Material *St = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    G4Material *LCP = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Material *Cu = nist->FindOrBuildMaterial("G4_Cu");

    //
    // World
    //
    G4Box *solidWorld =
        new G4Box("World",                                            //its name
                  worldLength / 2, worldLength / 2, worldLength / 2); //its size

    G4LogicalVolume *logicWorld =
        new G4LogicalVolume(solidWorld, //its solid
                            Vacuum,     //its material
                            "World");   //its name

    G4VPhysicalVolume *physWorld =
        new G4PVPlacement(0,               //no rotation
                          G4ThreeVector(), //at (0,0,0)
                          logicWorld,      //its logical volume
                          "World",         //its name
                          0,               //its mother  volume
                          false,           //no boolean operation
                          0,               //copy number
                          fCheckOverlaps); // checking overlaps
                                           //

    /*    
    //                       
    // Pbシールド
    // チェンバー本体
    G4Tubs *PbTube =
        new G4Tubs("PbTube", pRMax + Alring_thickness, pRMax + Alring_thickness + Pbring_thickness, pDz, pSPhi, pDPhi);

    G4LogicalVolume *logicPbTube =
        new G4LogicalVolume(PbTube,    //its solid
                            Pb,        //its material
                            "PbTube"); //its name

    G4VPhysicalVolume *physPbTube =
        new G4PVPlacement(0,                               // no rotation
                          G4ThreeVector(0, 0, 0), // at (0,0,0)
                          logicPbTube,                     // its logical volume
                          "PbPV",                          // its name
                          logicWorld,                      // its mother  volume
                          false,                           // no boolean operations
                          0,                               // copy number
                          fCheckOverlaps);                 // checking overlaps

    //コリメータ部
    G4Tubs *PbLid1 =
        new G4Tubs("PbLid", winR + 1 * mm, winR + 1 * mm + Pbring_thickness, 0.75 * cm, pSPhi, pDPhi);
    G4LogicalVolume *logicPbLid1 =
        new G4LogicalVolume(PbLid1,   //its solid
                            Pb,       //its material
                            "PbLid"); //its name
    G4VPhysicalVolume *physPbLid1 =
        new G4PVPlacement(0,                                                // no rotation
                          G4ThreeVector(0, 0, pDz + 0.75 * cm 0), // at (0,0,0)
                          logicPbLid1,                                      // its logical volume
                          "PbLidPV1",                                       // its name
                          logicWorld,                                       // its mother  volume
                          false,                                            // no boolean operations
                          0,                                                // copy number
                          fCheckOverlaps);                                  // checking overlaps

    //下蓋
    G4Tubs *PbLid2 =
        new G4Tubs("PbLid", pRMin, pRMax + Alring_thickness + Pbring_thickness, Pbring_thickness / 2, pSPhi, pDPhi);
    G4LogicalVolume *logicPbLid2 =
        new G4LogicalVolume(PbLid2,   //its solid
                            Pb,       //its material
                            "PbLid"); //its name
    G4VPhysicalVolume *physPbLid2 =
        new G4PVPlacement(0,                                                              // no rotation
                          G4ThreeVector(0, 0, -(pDz + Pbring_thickness / 2)), // at (0,0,0)
                          logicPbLid2,                                                    // its logical volume
                          "PbLidPV2",                                                     // its name
                          logicWorld,                                                     // its mother  volume
                          false,                                                          // no boolean operations
                          0,                                                              // copy number
                          fCheckOverlaps);

    //上リング部
    G4Tubs *PbLid3 =
        new G4Tubs("PbLid", winR + 1 * mm + Pbring_thickness, pRMax + Alring_thickness + Pbring_thickness, Pbring_thickness / 2, pSPhi, pDPhi);
    G4LogicalVolume *logicPbLid3 =
        new G4LogicalVolume(PbLid3,   //its solid
                            Pb,       //its material
                            "PbLid"); //its name
    G4VPhysicalVolume *physPbLid3 =
        new G4PVPlacement(0,                                                           // no rotation
                          G4ThreeVector(0, 0, pDz + Pbring_thickness / 2), // at (0,0,0)
                          logicPbLid3,                                                 // its logical volume
                          "PbLidPV3",                                                  // its name
                          logicWorld,                                                  // its mother  volume
                          false,                                                       // no boolean operations
                          0,                                                           // copy number
                          fCheckOverlaps);                                             // checking overlaps
*/

    //
    // Al Chamber
    //
    G4Tubs *AlChamber =
        new G4Tubs("Alchamber", pRMin, pRMax + Alring_thickness, pDz, pSPhi, pDPhi);

    G4LogicalVolume *logicAlChamber =
        new G4LogicalVolume(AlChamber,    //its solid
                            Al,           //its material
                            "AlChamber"); //its name

    G4VPhysicalVolume *physAlChamber =
        new G4PVPlacement(0,                      // no rotation
                          G4ThreeVector(0, 0, 0), // at (0,0,0)
                          logicAlChamber,         // its logical volume
                          "AlPV",                 // its name
                          logicWorld,             // its mother  volume
                          false,                  // no boolean operations
                          0,                      // copy number
                          fCheckOverlaps);        // checking overlaps

    //
    // XeCO2 GasCell
    //
    G4Tubs *GasCell =
        new G4Tubs("Ring", pRMin, pRMax, 9.25 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicGasCell =
        new G4LogicalVolume(GasCell,    //its solid
                            Xe5CO2,     //its material
                            "GasCell"); //its name

    G4VPhysicalVolume *physGasCell =
        new G4PVPlacement(0,                             // no rotation
                          G4ThreeVector(0, 0, 0.5 * mm), // at (0,0,0)
                          logicGasCell,                  // its logical volume
                          "TargetPV",                    // its name
                          logicAlChamber,                // its mother  volume
                          false,                         // no boolean operations
                          0,                             // copy number
                          fCheckOverlaps);               // checking overlaps

    //
    //Sensible Volume
    //

    G4Tubs *SolidSensitive =
        new G4Tubs("Sensitive", pRMin, 33.5 * mm, 7.65 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicSensitive =
        new G4LogicalVolume(SolidSensitive, //its solid
                            Xe5CO2,         //its material
                            "Sensitive");   //its name

    fSensitivePV =
        new G4PVPlacement(0,                             // no rotation
                          G4ThreeVector(0, 0, 1.6 * mm), // at (0,0,0)
                          logicSensitive,                // its logical volume
                          "SensitivePV",                 // its name
                          logicGasCell,                  // its mother  volume
                          false,                         // no boolean operations
                          0,                             // copy number
                          fCheckOverlaps);               // checking overlaps

    //
    //GEM
    //
    G4Tubs *GEM_Cu =
        new G4Tubs("GemCu", 0, winR + 1 * mm, 0.059 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicGem_Cu =
        new G4LogicalVolume(GEM_Cu,    //its solid
                            Cu,        //its material
                            "Gem_Cu"); //its name

    G4VPhysicalVolume *physGem_Cu =
        new G4PVPlacement(0,                                         // no rotation
                          G4ThreeVector(0, 0, -(6.95 + 0.059) * mm), // at (0,0,0)
                          logicGem_Cu,                               // its logical volume
                          "Gem_CuPV",                                // its name
                          logicGasCell,                              // its mother  volume
                          false,                                     // no boolean operations
                          0,                                         // copy number
                          fCheckOverlaps);                           // checking overlaps

    G4Tubs *GEM_LCP =
        new G4Tubs("GemLcp", 0, winR + 1 * mm, 0.05 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicGem_Lcp =
        new G4LogicalVolume(GEM_LCP,    //its solid
                            LCP,        //its material
                            "Gem_Lcp"); //its name

    G4VPhysicalVolume *physGem_Lcp =
        new G4PVPlacement(0,                      // no rotation
                          G4ThreeVector(0, 0, 0), // at (0,0,0)
                          logicGem_Lcp,           // its logical volume
                          "Gem_Lcp",              // its name
                          logicGem_Cu,            // its mother  volume
                          false,                  // no boolean operations
                          0,                      // copy number
                          fCheckOverlaps);        // checking overlaps

    //
    //Vavcuum_Tube
    //
    G4Tubs *VacuumRing =
        new G4Tubs("VacuRing", pRMin, winR, 3.5 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicVacu =
        new G4LogicalVolume(VacuumRing,  //its solid
                            Vacuum,      //its material
                            "VacuRing"); //its name

    G4VPhysicalVolume *physVacu =
        new G4PVPlacement(0,                                      // no rotation
                          G4ThreeVector(0, 0, (9.75 + 3.5) * mm), // at (0,0,0)
                          logicVacu,                              // its logical volume
                          "VacuPV",                               // its name
                          logicAlChamber,                         // its mother  volume
                          false,                                  // no boolean operations
                          0,                                      // copy number
                          fCheckOverlaps);                        // checking overlaps

    //
    //Al_grid
    //
    for (int j = 1; j < 4; j++)
    {
        auto xleng = sqrt(pow(3.5, 2) - pow(1.2 * (j - 0.5) + 0.1, 2)) * cm;

        G4Box *AlGrid =
            new G4Box("Grid",                      //its name
                      xleng, 0.75 * mm, 2.5 * mm); //its size

        G4LogicalVolume *logicGrid =
            new G4LogicalVolume(AlGrid,  //its solid
                                Al,      //its material
                                "Grid"); //its name

        G4VPhysicalVolume *GridPV1 =
            new G4PVPlacement(0,                                                   //no rotation
                              G4ThreeVector(0, (1.2 * (j - 0.5)) * cm, -0.9 * mm), //at (0,0,0)
                              logicGrid,                                           //its logical volume
                              "Grid",                                              //its name
                              logicVacu,                                           //its mother  volume
                              false,                                               //no boolean operation
                              0,                                                   //copy number
                              fCheckOverlaps);                                     // checking overlaps

        G4VPhysicalVolume *GridPV2 =
            new G4PVPlacement(0,                                                      //no rotation
                              G4ThreeVector(0, -((1.2 * (j - 0.5)) * cm), -0.9 * mm), //at (0,0,0)
                              logicGrid,                                              //its logical volume
                              "Grid",                                                 //its name
                              logicVacu,                                              //its mother  volume
                              false,                                                  //no boolean operation
                              0,                                                      //copy number
                              fCheckOverlaps);                                        // checking overlaps

        G4VisAttributes *GVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
        logicGrid->SetVisAttributes(GVisAtt);

        int kk = (xleng + 0.1 * cm) / (1.0 * cm);
        for (int k = -kk; k < kk; k++)
        {
            auto xl = 5 * mm;
            auto yp = 1.2 * j * cm;
            auto xp = (1.1 * k + 0.5 + 0.05) * cm;

            if (j == 3)
            {
                xl = xl - 3.5 * mm;
                yp = yp - 3.5 * mm;
            }

            if (j == 2 && k == -3)
            {
                xl = xl - 4.5 * mm;
                yp = yp - 4.5 * mm;
            }

            if (j == 2 && k == 2)
            {
                xl = xl - 4.5 * mm;
                yp = yp - 4.5 * mm;
            }

            G4Box *AlGrid_s =
                new G4Box("Grid",                   //its name
                          0.75 * mm, xl, 2.5 * mm); //its size

            G4LogicalVolume *logicGrid_s =
                new G4LogicalVolume(AlGrid_s, //its solid
                                    Al,       //its material
                                    "Grid");  //its name

            G4VPhysicalVolume *Grid_sPV1 =
                new G4PVPlacement(0,                                //no rotation
                                  G4ThreeVector(xp, yp, -0.9 * mm), //at (0,0,0)
                                  logicGrid_s,                      //its logical volume
                                  "Grid",                           //its name
                                  logicVacu,                        //its mother  volume
                                  false,                            //no boolean operation
                                  0,                                //copy number
                                  fCheckOverlaps);                  // checking overlaps

            G4VPhysicalVolume *Grid_sPV2 =
                new G4PVPlacement(0,                                 //no rotation
                                  G4ThreeVector(xp, -yp, -0.9 * mm), //at (0,0,0)
                                  logicGrid_s,                       //its logical volume
                                  "Grid",                            //its name
                                  logicVacu,                         //its mother  volume
                                  false,                             //no boolean operation
                                  0,                                 //copy number
                                  fCheckOverlaps);                   // checking overlaps

            if (j == 1)
            {
                G4VPhysicalVolume *Grid_sPV0 =
                    new G4PVPlacement(0,                               //no rotation
                                      G4ThreeVector(xp, 0, -0.9 * mm), //at (0,0,0)
                                      logicGrid_s,                     //its logical volume
                                      "Grid",                          //its name
                                      logicVacu,                       //its mother  volume
                                      false,                           //no boolean operation
                                      0,                               //copy number
                                      fCheckOverlaps);                 // checking overlaps
            }
            logicGrid_s->SetVisAttributes(GVisAtt);
        }
    }

    //
    //Be_window
    //
    G4Tubs *Be_window =
        new G4Tubs("BeLid", 0, winR, 0.05 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicBe_window =
        new G4LogicalVolume(Be_window,    //its solid
                            Be,           //its material
                            "Be_window"); //its name

    G4VPhysicalVolume *physBe_window =
        new G4PVPlacement(0,                               // no rotation
                          G4ThreeVector(0, 0, -3.45 * mm), // at (0,0,0)
                          logicBe_window,                  // its logical volume
                          "Be_windowPV",                   // its name
                          logicVacu,                       // its mother  volume
                          false,                           // no boolean operations
                          0,                               // copy number
                          fCheckOverlaps);                 // checking overlaps

    //
    //Collimeter
    //
    G4Tubs *Collimeter =
        new G4Tubs("ColLid", 0, winR + 1 * mm, 0.75 * cm, pSPhi, pDPhi);

    G4LogicalVolume *logicCollimeter =
        new G4LogicalVolume(Collimeter,    //its solid
                            St,            //its material
                            "Collimeter"); //its name

    G4VPhysicalVolume *physCollimeter =
        new G4PVPlacement(0,                                    // no rotation
                          G4ThreeVector(0, 0, pDz + 0.75 * cm), // at (0,0,0)
                          logicCollimeter,                      // its logical volume
                          "CollimeterPV",                       // its name
                          logicWorld,                           // its mother  volume
                          false,                                // no boolean operations
                          0,                                    // copy number
                          fCheckOverlaps);                      // checking overlaps

    //
    //Vacuum_Tubes　Hexagon　
    //
    const G4double zP[] = {-0.75 * cm, 0.75 * cm};
    const G4double rI[] = {0, 0};
    auto r = 0.3 * mm;
    const G4double rO[] = {r, r};
    G4Polyhedra *Vacuum_Tube = new G4Polyhedra("VacuT", 0 * deg, 360 * deg, 6, 2, zP, rI, rO);

    G4LogicalVolume *logicVacuum_Tube =
        new G4LogicalVolume(Vacuum_Tube,    //its solids
                            Vacuum,         //its material
                            "Vacuum_Tube"); //its name

    G4RotationMatrix *rotCounter = new G4RotationMatrix;
    rotCounter->rotateZ(90. * deg);

    auto y_width = (1.8 / sqrt(3) + 0.1 * sqrt(3)) * mm;
    for (auto y_num = -28; y_num < 29; y_num++)
    {
        int h = sqrt(pow(y_num, 2));
        int x_num = sqrt(pow(35, 2) - pow(y_width * (h + 0.25), 2)) / 0.7;
        for (int i = -x_num - 1; i < x_num + 1; i++)
        {
            auto x_width = 0.7 * i * mm;
            auto hh = y_num - 0.25;

            if (i == -x_num - 1 && y_num != -29)
            {
                G4VPhysicalVolume *physVacuum_Tube_b =
                    new G4PVPlacement(rotCounter,                                          // no rotation
                                      G4ThreeVector(x_width + 0.35 * mm, y_width * hh, 0), // at (0,0,0)
                                      logicVacuum_Tube,                                    // its logical volume
                                      "Vacuum_TubePV",                                     // its name
                                      logicCollimeter,                                     // its mother  volume
                                      false,                                               // no boolean operations
                                      0,                                                   // copy number
                                      fCheckOverlaps);                                     // checking overlaps
            }

            if (i != -x_num - 1 && y_num == -29)
            {
                G4VPhysicalVolume *physVacuum_Tube_a =
                    new G4PVPlacement(rotCounter,                                          // no rotation
                                      G4ThreeVector(x_width, y_width * (y_num + 0.25), 0), // at (0,0,0)
                                      logicVacuum_Tube,                                    // its logical volume
                                      "Vacuum_TubePV",                                     // its name
                                      logicCollimeter,                                     // its mother  volume
                                      false,                                               // no boolean operations
                                      0,                                                   // copy number
                                      fCheckOverlaps);                                     // checking overlaps
            }

            if (i != -x_num - 1 && y_num != -29)
            {
                G4VPhysicalVolume *physVacuum_Tube_b =
                    new G4PVPlacement(rotCounter,                                          // no rotation
                                      G4ThreeVector(x_width + 0.35 * mm, y_width * hh, 0), // at (0,0,0)
                                      logicVacuum_Tube,                                    // its logical volume
                                      "Vacuum_TubePV",                                     // its name
                                      logicCollimeter,                                     // its mother  volume
                                      false,                                               // no boolean operations
                                      0,                                                   // copy number
                                      fCheckOverlaps);                                     // checking overlaps

                G4VPhysicalVolume *physVacuum_Tube_a =
                    new G4PVPlacement(rotCounter,                                          // no rotation
                                      G4ThreeVector(x_width, y_width * (y_num + 0.25), 0), // at (0,0,0)
                                      logicVacuum_Tube,                                    // its logical volume
                                      "Vacuum_TubePV",                                     // its name
                                      logicCollimeter,                                     // its mother  volume
                                      false,                                               // no boolean operations
                                      0,                                                   // copy number
                                      fCheckOverlaps);                                     // checking overlaps
            }
        }
    }

    //
    // Visualization attributes
    //
    G4VisAttributes *RedVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
    G4VisAttributes *OrangeVisAtt = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0));
    G4VisAttributes *MagentaVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
    G4VisAttributes *GreenVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    G4VisAttributes *BlueVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));

    logicWorld->SetVisAttributes(RedVisAtt);
    logicGasCell->SetVisAttributes(OrangeVisAtt);
    logicAlChamber->SetVisAttributes(OrangeVisAtt);
    logicSensitive->SetVisAttributes(BlueVisAtt);
    logicBe_window->SetVisAttributes(MagentaVisAtt);
    logicCollimeter->SetVisAttributes(BlueVisAtt);
    logicVacuum_Tube->SetVisAttributes(GreenVisAtt);
    logicVacu->SetVisAttributes(OrangeVisAtt);
    logicGem_Cu->SetVisAttributes(GreenVisAtt);
    logicGem_Lcp->SetVisAttributes(RedVisAtt);
    /*
    logicPbTube->SetVisAttributes(BlueVisAtt);
    logicPbLid1->SetVisAttributes(BlueVisAtt);
    logicPbLid2->SetVisAttributes(BlueVisAtt);
    logicPbLid3->SetVisAttributes(BlueVisAtt);
*/
    return physWorld;
}