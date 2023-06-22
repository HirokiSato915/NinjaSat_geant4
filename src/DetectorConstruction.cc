//
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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "Constants.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fCheckOverlaps(false)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{
    //*************************************************************************************************
    // Material
    //*************************************************************************************************
    G4NistManager *nist = G4NistManager::Instance();
    G4Material *Al = nist->FindOrBuildMaterial("G4_Al");
    G4Material *LCP = nist->FindOrBuildMaterial("G4_KAPTON");
    G4Material *Cu = nist->FindOrBuildMaterial("G4_Cu");
    G4Material *Si = nist->FindOrBuildMaterial("G4_Si");
    G4Material *Mn = nist->FindOrBuildMaterial("G4_Mn");
    G4Material *Ni = nist->FindOrBuildMaterial("G4_Ni");
    G4Material *Cr = nist->FindOrBuildMaterial("G4_Cr");
    G4Material *Fe = nist->FindOrBuildMaterial("G4_Fe");
    G4Material *Ti = nist->FindOrBuildMaterial("G4_Ti");
    G4Material *Mo = nist->FindOrBuildMaterial("G4_Mo");
    G4Material *C = nist->FindOrBuildMaterial("G4_C");
    G4Material *H = nist->FindOrBuildMaterial("G4_H");
    G4Material *O = nist->FindOrBuildMaterial("G4_O");
    G4Material *Sn = nist->FindOrBuildMaterial("G4_Sn");
    G4Material *Sb = nist->FindOrBuildMaterial("G4_Sb");

    G4double mmN = 14.01 * g / mole;
    G4Element *elN = new G4Element("Nitrogen", "N", 7., mmN);
    G4double mmO = 15.999 * g / mole;
    G4Element *elO = new G4Element("Oxygen", "O", 8., mmO);

    G4Material *Vacuum = new G4Material("Vacuum", universe_mean_density, 2);
    Vacuum->AddElement(elN, .7);
    Vacuum->AddElement(elO, .3);

    G4double mmXe = 131.29 * g / mole;
    G4Element *elXe = new G4Element("Xenon", "Xe", 54., mmXe);

    G4double mmAr = 39.948 * g / mole;
    G4Element *elAr = new G4Element("Argon", "Ar", 18., mmAr);

    G4double denC = 12.011 * g / mole;
    G4Element *elC = new G4Element("Carbon", "C", 6., denC);
    G4double denH = 1.008 * g / mole;
    G4Element *elH = new G4Element("Hydrogen", "H", 1., denH);

    G4double denXeArDME = 5.844 * mg / cm3;
    G4Material *XeArDME = new G4Material("XeArDME", denXeArDME, 5);
    XeArDME->AddElement(elXe, .9077);
    XeArDME->AddElement(elAr, .0879);
    XeArDME->AddElement(elC, .0023);
    XeArDME->AddElement(elH, .0006);
    XeArDME->AddElement(elO, .0015);

    G4Element *elBe = new G4Element("Berillium", "Be", 4., 9.01 * g / mole);
    G4Material *Be = new G4Material("Be", 1.848 * g / cm3, 1);
    Be->AddElement(elBe, 1);

    G4double denSUS304 = 7.93 * g / cm3;
    G4Material *SUS304 = new G4Material("SUS304", denSUS304, 3);
    SUS304->AddMaterial(Fe, .74);
    SUS304->AddMaterial(Cr, .18);
    SUS304->AddMaterial(Ni, .08);

    G4double denSUS310 = 7.98 * g / cm3;
    G4Material *SUS310 = new G4Material("SUS310", denSUS310, 5);
    SUS310->AddMaterial(Fe, .52);
    SUS310->AddMaterial(Cr, .23);
    SUS310->AddMaterial(Ni, .22);
    SUS310->AddMaterial(Mn, .02);
    SUS310->AddMaterial(Si, .01);

    G4Material *Peek = new G4Material("Peek", 1.31 * g / cm3, 3);
    Peek->AddMaterial(C, .76);
    Peek->AddMaterial(H, .08);
    Peek->AddMaterial(O, .16);

    G4Material *Alumina = new G4Material("Alumina", 3.8 * g / cm3, 2);
    Alumina->AddMaterial(Al, .53);
    Alumina->AddMaterial(O, .47);

    G4double denSnSb = 7.28 * g / cm3;
    G4Material *SnSb = new G4Material("SnSb", denSnSb, 2);
    SnSb->AddMaterial(Sn, .95);
    SnSb->AddMaterial(Sb, .05);

    //*************************************************************************************************
    // VisAttributes
    //*************************************************************************************************
    G4VisAttributes *BlueVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
    G4VisAttributes *OrangeVisAtt = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0));
    G4VisAttributes *WhiteVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    G4VisAttributes *GrayVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    G4VisAttributes *DeepGreenVisAtt = new G4VisAttributes(G4Colour(0., 0.31, 0.18));
    G4VisAttributes *SUS310VisAtt = new G4VisAttributes(G4Colour(0.5, 0.55, 0.53));
    G4VisAttributes *CuVisAtt = new G4VisAttributes(G4Colour(0.6, 0.4, 0.25));
    G4VisAttributes *AlVisAtt = new G4VisAttributes(G4Colour(0.76, 0.75, 0.71));
    G4VisAttributes *TiVisAtt = new G4VisAttributes(G4Colour(0.6, 0.6, 0.6));
    G4VisAttributes *BeVisAtt = new G4VisAttributes(G4Colour(0.93, 0.89, 0.83));
    G4VisAttributes *PEEKVisAtt = new G4VisAttributes(G4Colour(0.87, 0.81, 0.74));

    //*************************************************************************************************
    // World
    //*************************************************************************************************
    G4Box *solidWorld =
        new G4Box("WorldS",                                           // its name
                  worldLength / 2, worldLength / 2, worldLength / 2); // its size

    G4LogicalVolume *logicWorld =
        new G4LogicalVolume(solidWorld, // its solid
                            Vacuum,     // its material
                            "WorldLV"); // its name

    G4VPhysicalVolume *physWorld =
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(), // at (0,0,0)
                          logicWorld,      // its logical volume
                          "WorldPV",       // its name
                          0,               // its mother  volume
                          false,           // no boolean operation
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps

    //*************************************************************************************************
    // Collimeter Ear
    //*************************************************************************************************
    G4Box *solidColliEar =
        new G4Box("ColliEarS",                   // its name
                  3.5 * mm, 2.0 * mm, 7.5 * mm); // its size

    G4LogicalVolume *logicColliEar =
        new G4LogicalVolume(solidColliEar, // its solid
                            SUS304,        // its material
                            "ColliEarLV"); // its name

    G4double colli_earR = 39.5 * mm;

    for (G4int i = 0; i < 6; i++)
    {
        G4double angV = ang60 * i;
        G4RotationMatrix *rotZV = new G4RotationMatrix;
        rotZV->rotateZ(-angV);
        G4double px = colli_earR * cos(angV);
        G4double py = colli_earR * sin(angV);
        new G4PVPlacement(rotZV,
                          G4ThreeVector(px, py, chamber_half_hight + 0.75 * cm),
                          logicColliEar,   // its logical volume
                          "ColliEarPV",    // its name
                          logicWorld,      // its mother  volume
                          false,           // no boolean operation
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
    }

    //*************************************************************************************************
    // Collimeter (Hexagonal Aperture)
    //*************************************************************************************************
    G4Tubs *solidCollimeter =
        new G4Tubs("CollimeterS", 0, winR + 1 * mm, 0.75 * cm, pSPhi, pDPhi);

    G4LogicalVolume *logicCollimeter =
        new G4LogicalVolume(solidCollimeter, // its solid
                            SUS304,          // its material
                            "CollimeterLV"); // its name

    new G4PVPlacement(0,                                                   // no rotation
                      G4ThreeVector(0, 0, chamber_half_hight + 0.75 * cm), // at
                      logicCollimeter,                                     // its logical volume
                      "CollimeterPV",                                      // its name
                      logicWorld,                                          // its mother  volume
                      false,                                               // no boolean operations
                      0,                                                   // copy number
                      fCheckOverlaps);                                     // checking overlaps

    const G4double zP[] = {-0.75 * cm, 0.75 * cm};
    const G4double rI[] = {0, 0};
    auto hexr = 0.3 * mm;
    const G4double rO[] = {hexr, hexr};

    G4Polyhedra *solidHexAperture =
        new G4Polyhedra("HexApertureS", 0 * deg, 360 * deg, 6, 2, zP, rI, rO);

    G4LogicalVolume *logicHexAperture =
        new G4LogicalVolume(solidHexAperture, // its solids
                            Vacuum,           // its material
                            "HexApertureLV"); // its name

    G4RotationMatrix *rotZ90 = new G4RotationMatrix;
    rotZ90->rotateZ(ang90);
    G4double x_grid_width = 0.7;
    G4double y_grid_width = 1.8 / sqrt(3) + 0.1 * sqrt(3);
    G4int y_grid_num_a = (35 - 0.6 / sqrt(3)) / (y_grid_width);

    for (G4int i = 0; i < y_grid_num_a + 1; i++)
    {
        G4double y_grid_coordinate = y_grid_width * i;
        G4double x_grid_length = sqrt(pow(35, 2) - pow(y_grid_coordinate + 0.3 / sqrt(3), 2));
        G4int x_grid_num = x_grid_length / x_grid_width;
        for (G4int j = 0; j < x_grid_num + 1; j++)
        {
            G4double x_grid_coordinate = x_grid_width * j;
            for (G4int k = -1; k < 2; k = k + 2)
            {
                new G4PVPlacement(rotZ90,
                                  G4ThreeVector(k * x_grid_coordinate * mm,
                                                k * y_grid_coordinate * mm, 0),
                                  logicHexAperture, // its logical volume
                                  "HexAperturePV",  // its name
                                  logicCollimeter,  // its mother  volume
                                  false,            // no boolean operations
                                  0,                // copy number
                                  fCheckOverlaps);  // checking overlaps
            }

            if (i != 0 && j != 0)
            {
                for (G4int l = -1; l < 2; l = l + 2)
                {
                    new G4PVPlacement(rotZ90,
                                      G4ThreeVector(l * x_grid_coordinate * mm,
                                                    -l * y_grid_coordinate * mm, 0),
                                      logicHexAperture, // its logical volume
                                      "HexAperturePV",  // its name
                                      logicCollimeter,  // its mother  volume
                                      false,            // no boolean operations
                                      0,                // copy number
                                      fCheckOverlaps);
                }
            }
        }
    }
    G4int y_grid_num_b = (35 - (y_grid_width / 2 + 0.6 / sqrt(3))) / (y_grid_width) + 1;
    for (G4int i = 0; i < y_grid_num_b + 1; i++)
    {
        G4double y_grid_coordinate = y_grid_width * i + y_grid_width / 2;
        G4double x_grid_length = sqrt(pow(35, 2) - pow(y_grid_coordinate + 0.3 / sqrt(3), 2)) - 0.35;
        G4int x_grid_num = x_grid_length / x_grid_width;
        for (G4int j = 0; j < x_grid_num + 1; j++)
        {
            G4double x_grid_coordinate = x_grid_width * j + 0.35;
            for (G4int k = -1; k < 2; k = k + 2)
            {
                for (G4int l = -1; l < 2; l = l + 2)
                {
                    new G4PVPlacement(rotZ90,
                                      G4ThreeVector(k * x_grid_coordinate * mm,
                                                    l * y_grid_coordinate * mm, 0),
                                      logicHexAperture, // its logical volume
                                      "HexAperturePV",  // its name
                                      logicCollimeter,  // its mother  volume
                                      false,            // no boolean operations
                                      0,                // copy number
                                      fCheckOverlaps);
                }
            }
        }
    }

    //*************************************************************************************************
    // Al Chamber
    //*************************************************************************************************
    G4Tubs *solidAlChamber =
        new G4Tubs("AlchamberS", pRMin, pRMax + chamber_side_thickness, chamber_half_hight, pSPhi, pDPhi);

    G4LogicalVolume *logicAlChamber =
        new G4LogicalVolume(solidAlChamber, // its solid
                            Al,             // its material
                            "AlChamberLV"); // its name

    new G4PVPlacement(0,                      // no rotation
                      G4ThreeVector(0, 0, 0), // at (0,0,0)
                      logicAlChamber,         // its logical volume
                      "AlChamberPV",          // its name
                      logicWorld,             // its mother  volume
                      false,                  // no boolean operations
                      0,                      // copy number
                      fCheckOverlaps);        // checking overlaps

    //*************************************************************************************************
    // SUS310 Plate (SUS310 & Al)
    //*************************************************************************************************
    G4Tubs *solidSUS310Plate =
        new G4Tubs("SUS310PlateS", 37.5 * mm, Be_window_R, 3.25 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicSUS310Plate =
        new G4LogicalVolume(solidSUS310Plate, // its solid
                            SUS310,           // its material
                            "SUS310PlateLV"); // its name

    new G4PVPlacement(0,                              // no rotation
                      G4ThreeVector(0, 0, 13.5 * mm), // at
                      logicSUS310Plate,               // its logical volume
                      "SUS310PlatePV",                // its name
                      logicAlChamber,                 // its mother  volume
                      false,                          // no boolean operations
                      0,                              // copy number
                      fCheckOverlaps);                // checking overlaps

    G4Tubs *solidSUS310Plate_Al =
        new G4Tubs("SUS310Plate_AlS", 37.5 * mm, 42.5 * mm, 0.75 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicSUS310Plate_Al =
        new G4LogicalVolume(solidSUS310Plate_Al, // its solid
                            Al,                  // its material
                            "SUS310Plate_AlLV"); // its name

    new G4PVPlacement(0,                             // no rotation
                      G4ThreeVector(0, 0, 2.5 * mm), // at
                      logicSUS310Plate_Al,           // its logical volume
                      "SUS310Plate_AlPV",            // its name
                      logicSUS310Plate,              // its mother  volume
                      false,                         // no boolean operations
                      0,                             // copy number
                      fCheckOverlaps);               // checking overlaps

    //*************************************************************************************************
    // SUS304 Transition　(Upper & Lower)
    //*************************************************************************************************
    G4Tubs *solidSUS304TraU =
        new G4Tubs("SUS304TraUS", Be_window_R, pRMax + chamber_side_thickness, chamber_top_hight / 2, pSPhi, pDPhi);

    G4LogicalVolume *logicSUS304TraU =
        new G4LogicalVolume(solidSUS304TraU, // its solid
                            SUS304,          // its material
                            "SUS304TraULV"); // its name

    new G4PVPlacement(0,                              // no rotation
                      G4ThreeVector(0, 0, 13.4 * mm), // at
                      logicSUS304TraU,                // its logical volume
                      "SUS304TraUPV",                 // its name
                      logicAlChamber,                 // its mother  volume
                      false,                          // no boolean operations
                      0,                              // copy number
                      fCheckOverlaps);                // checking overlaps

    G4Tubs *solidSUS304TraL =
        new G4Tubs("SUS304TraLS", pRMax, pRMax + chamber_side_thickness, 1.65 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicSUS304TraL =
        new G4LogicalVolume(solidSUS304TraL, // its solid
                            SUS304,          // its material
                            "SUS304TraLLV"); // its name

    new G4PVPlacement(0,                             // no rotation
                      G4ThreeVector(0, 0, 8.4 * mm), // at
                      logicSUS304TraL,               // its logical volume
                      "SUS304TraLPV",                // its name
                      logicAlChamber,                // its mother  volume
                      false,                         // no boolean operations
                      0,                             // copy number
                      fCheckOverlaps);               // checking overlaps

    //*************************************************************************************************
    // Ti Transition (1mm)
    //*************************************************************************************************
    G4Tubs *solidTiTra =
        new G4Tubs("TiTraS", pRMax, pRMax + chamber_side_thickness, 0.5 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicTiTra =
        new G4LogicalVolume(solidTiTra, // its solid
                            Ti,         // its material
                            "TiTraLV"); // its name

    new G4PVPlacement(0,                              // no rotation
                      G4ThreeVector(0, 0, 6.25 * mm), // at
                      logicTiTra,                     // its logical volume
                      "TiTraPV",                      // its name
                      logicAlChamber,                 // its mother  volume
                      false,                          // no boolean operations
                      0,                              // copy number
                      fCheckOverlaps);                // checking overlaps

    //*************************************************************************************************
    // Al Lattice Window
    //*************************************************************************************************
    auto LWRH = 1.5 * mm;
    G4Tubs *solidLWApertureRing =
        new G4Tubs("LWApertureRingS", pRMin, winR, LWRH / 2, pSPhi, pDPhi);

    G4LogicalVolume *logicLWApertureRing =
        new G4LogicalVolume(solidLWApertureRing, // its solid
                            Vacuum,              // its material
                            "LWApertureRingLV"); // its name

    new G4PVPlacement(0,                                                  // no rotation
                      G4ThreeVector(0, 0, chamber_half_hight - LWRH / 2), // at
                      logicLWApertureRing,                                // its logical volume
                      "LWApertureRingPV",                                 // its name
                      logicAlChamber,                                     // its mother  volume
                      false,                                              // no boolean operations
                      0,                                                  // copy number
                      fCheckOverlaps);                                    // checking overlaps

    auto LWH = chamber_top_hight - LWRH - Be_window_thickness;
    G4Tubs *solidLWAperture =
        new G4Tubs("LWApertureS", pRMin, winR + 0.7 * mm, LWH / 2, pSPhi, pDPhi);

    G4LogicalVolume *logicLWAperture =
        new G4LogicalVolume(solidLWAperture, // its solid
                            Vacuum,          // its material
                            "LWApertureLV"); // its name

    new G4PVPlacement(0,                                                        // no rotation
                      G4ThreeVector(0, 0, chamber_half_hight - LWRH - LWH / 2), // at
                      logicLWAperture,                                          // its logical volume
                      "LWAperturePV",                                           // its name
                      logicAlChamber,                                           // its mother  volume
                      false,                                                    // no boolean operations
                      0,                                                        // copy number
                      fCheckOverlaps);                                          // checking overlaps

    for (G4int j = 1; j < 4; j++)
    {
        auto xlength = sqrt(pow(35 + 0.7, 2) - pow(11.5 * (j - 1) + 5.75 + 0.75, 2)) * mm;
        G4Box *solidAlLatticeX =
            new G4Box("AlLatticeXS",                // its name
                      xlength, 0.75 * mm, LWH / 2); // its size
        G4LogicalVolume *logicAlLatticeX =
            new G4LogicalVolume(solidAlLatticeX, // its solid
                                Al,              // its material
                                "AlLatticeXLV"); // its name

        for (G4int k = -1; k < 2; k = k + 2)
        {
            new G4PVPlacement(0,                                                     // no rotation
                              G4ThreeVector(0, k * (11.5 * (j - 1) + 5.75) * mm, 0), // at
                              logicAlLatticeX,                                       // its logical volume
                              "AlLatticePV",                                         // its name
                              logicLWAperture,                                       // its mother  volume
                              false,                                                 // no boolean operation
                              0,                                                     // copy number
                              fCheckOverlaps);                                       // checking overlaps
        }

        logicAlLatticeX->SetVisAttributes(GrayVisAtt);

        G4int y_grid_num = (xlength + 0.1 * cm) / (1.0 * cm);
        for (G4int k = -y_grid_num; k < y_grid_num; k++)
        {
            auto xl = 5. * mm;
            auto yp = 11.5 * j * mm;
            auto xp = (11.5 * k + 5.75) * mm;

            if (j == 3)
            {
                if (k == -1 || k == 0)
                {
                    xl = xl - 2.20 * mm;
                    yp = yp - 2.20 * mm;
                }

                if (k == -2 || k == 1)
                {
                    xl = xl - 4.34 * mm;
                    yp = yp - 4.34 * mm;
                }
            }

            if (j == 2 && (k == -3 || k == 2))
            {
                xl = xl - 3.95 * mm;
                yp = yp - 3.95 * mm;
            }

            G4Box *solidAlLatticeY =
                new G4Box("AlLatticeYS",           // its name
                          0.75 * mm, xl, LWH / 2); // its size

            G4LogicalVolume *logicAlLatticeY =
                new G4LogicalVolume(solidAlLatticeY, // its solid
                                    Al,              // its material
                                    "AlLatticeYLV"); // its name

            for (G4int l = -1; l < 2; l = l + 2)
            {
                new G4PVPlacement(0,                            // no rotation
                                  G4ThreeVector(xp, l * yp, 0), // at
                                  logicAlLatticeY,              // its logical volume
                                  "AlLatticeYPV",               // its name
                                  logicLWAperture,              // its mother  volume
                                  false,                        // no boolean operation
                                  0,                            // copy number
                                  fCheckOverlaps);              // checking overlaps

                if (j == 1 && l == 1)
                {
                    new G4PVPlacement(0,                       // no rotation
                                      G4ThreeVector(xp, 0, 0), // at
                                      logicAlLatticeY,         // its logical volume
                                      "AlLatticeY0PV",         // its name
                                      logicLWAperture,         // its mother  volume
                                      false,                   // no boolean operation
                                      0,                       // copy number
                                      fCheckOverlaps);         // checking overlaps
                }
            }
            logicAlLatticeY->SetVisAttributes(GrayVisAtt);
        }
    }

    //*************************************************************************************************
    // Be Window
    //*************************************************************************************************
    G4Tubs *solidBeWindow =
        new G4Tubs("BeWindowS", 0, Be_window_R, Be_window_thickness / 2, pSPhi, pDPhi);

    G4LogicalVolume *logicBeWindow =
        new G4LogicalVolume(solidBeWindow, // its solid
                            Be,            // its material
                            "BeWindowLV"); // its name

    new G4PVPlacement(0,                                                                                     // no rotation
                      G4ThreeVector(0, 0, chamber_half_hight - chamber_top_hight + Be_window_thickness / 2), // at
                      logicBeWindow,                                                                         // its logical volume
                      "BeWindowPV",                                                                          // its name
                      logicAlChamber,                                                                        // its mother  volume
                      false,                                                                                 // no boolean operations
                      0,                                                                                     // copy number
                      fCheckOverlaps);                                                                       // checking overlaps

    //*****************************************************************************************************
    // XeArDME GasCell
    //*****************************************************************************************************
    auto gas_cell_hight = 18.8 * mm;
    G4Tubs *solidGasCell =
        new G4Tubs("GasCellS", pRMin, pRMax, gas_cell_hight / 2, pSPhi, pDPhi);

    G4LogicalVolume *logicGasCell =
        new G4LogicalVolume(solidGasCell, // its solid
                            XeArDME,      // its material
                            "GasCellLV"); // its name

    fGasCellPV = new G4PVPlacement(0,                                                                                // no rotation
                                   G4ThreeVector(0, 0, chamber_half_hight - chamber_top_hight - gas_cell_hight / 2), // at
                                   logicGasCell,                                                                     // its logical volume
                                   "GasCellPV",                                                                      // its name
                                   logicAlChamber,                                                                   // its mother  volume
                                   false,                                                                            // no boolean operations
                                   0,                                                                                // copy number
                                   fCheckOverlaps);                                                                  // checking overlaps

    G4Region *regionA = new G4Region("regionA");
    // Attach a logical volume to the region
    regionA->AddRootLogicalVolume(logicGasCell);

    //*****************************************************************************************************
    //　Sensible Volume (Inner & Outer)
    //*****************************************************************************************************
    G4Tubs *solidSensitiveInner =
        new G4Tubs("SensitiveInnerS", pRMin, InnerR, sensitive_layer_thickness / 2, pSPhi, pDPhi);

    G4LogicalVolume *logicSensitiveInner =
        new G4LogicalVolume(solidSensitiveInner, // its solid
                            XeArDME,             // its material
                            "SensitiveInnerLV"); // its name

    fSensitiveInnerPV = new G4PVPlacement(0,                                                                     // no rotation
                                          G4ThreeVector(0, 0, (gas_cell_hight - sensitive_layer_thickness) / 2), // at
                                          logicSensitiveInner,                                                   // its logical volume
                                          "SensitiveInnerPV",                                                    // its name
                                          logicGasCell,                                                          // its mother  volume
                                          false,                                                                 // no boolean operations
                                          0,                                                                     // copy number
                                          fCheckOverlaps);                                                       // checking overlaps

    G4Tubs *solidSensitiveOuter =
        new G4Tubs("SensitiveOuterS", InnerR, OuterR, sensitive_layer_thickness / 2, pSPhi, pDPhi);

    G4LogicalVolume *logicSensitiveOuter =
        new G4LogicalVolume(solidSensitiveOuter, // its solid
                            XeArDME,             // its material
                            "SensitiveOuterLV"); // its name

    fSensitiveOuterPV = new G4PVPlacement(0,                                                                     // no rotation
                                          G4ThreeVector(0, 0, (gas_cell_hight - sensitive_layer_thickness) / 2), // at
                                          logicSensitiveOuter,                                                   // its logical volume
                                          "SensitiveOuterPV",                                                    // its name
                                          logicGasCell,                                                          // its mother  volume
                                          false,                                                                 // no boolean operations
                                          0,                                                                     // copy number
                                          fCheckOverlaps);                                                       // checking overlaps

    //*****************************************************************************************************
    //　GEM (Cu & LCP)
    //*****************************************************************************************************
    G4Tubs *solidGEMCu =
        new G4Tubs("GEMCuS", 0, winR + 1 * mm, 0.059 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicGEMCu =
        new G4LogicalVolume(solidGEMCu, // its solid
                            Cu,         // its material
                            "GemCuLV"); // its name

    new G4PVPlacement(0,                                         // no rotation
                      G4ThreeVector(0, 0, -(6.05 + 0.059) * mm), // at
                      logicGEMCu,                                // its logical volume
                      "GEMCuPV",                                 // its name
                      logicGasCell,                              // its mother  volume
                      false,                                     // no boolean operations
                      0,                                         // copy number
                      fCheckOverlaps);                           // checking overlaps

    G4Tubs *solidGEMLCP =
        new G4Tubs("GEMLCPS", 0, winR + 1 * mm, 0.05 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicGEMLCP =
        new G4LogicalVolume(solidGEMLCP, // its solid
                            LCP,         // its material
                            "GEMLCPLV"); // its name

    new G4PVPlacement(0,                      // no rotation
                      G4ThreeVector(0, 0, 0), // at (0,0,0)
                      logicGEMLCP,            // its logical volume
                      "GEMLCPPV",             // its name
                      logicGEMCu,             // its mother  volume
                      false,                  // no boolean operations
                      0,                      // copy number
                      fCheckOverlaps);        // checking overlaps

    //*****************************************************************************************************
    //　Readout Pad
    //*****************************************************************************************************
    G4Tubs *solidReadoutPad =
        new G4Tubs("ReadoutPadS", 0, 41. * mm, 0.5 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicReadoutPad =
        new G4LogicalVolume(solidReadoutPad, // its solid
                            Alumina,         // its material
                            "ReadoutPadLV"); // its name

    new G4PVPlacement(0,                               // no rotation
                      G4ThreeVector(0, 0, -7.98 * mm), // at (0,0,0)
                      logicReadoutPad,                 // its logical volume
                      "ReadoutPadPV",                  // its name
                      logicGasCell,                    // its mother  volume
                      false,                           // no boolean operations
                      0,                               // copy number
                      fCheckOverlaps);                 // checking overlaps

    //*****************************************************************************************************
    // Al-Cu Transisiton
    //*****************************************************************************************************
    G4Tubs *solidCuTra =
        new G4Tubs("CuTraS", 3.75 * mm, 13.5 * mm, 2.5 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicCuTra =
        new G4LogicalVolume(solidCuTra, // its solid
                            Cu,         // its material
                            "CuTraLV"); // its name

    new G4PVPlacement(0,                                                     // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + 2.5 * mm)), // at
                      logicCuTra,                                            // its logical volume
                      "CuTraPV",                                             // its name
                      logicWorld,                                            // its mother  volume
                      false,                                                 // no boolean operations
                      0,                                                     // copy number
                      fCheckOverlaps);                                       // checking overlaps

    G4Tubs *solidAlTra =
        new G4Tubs("AlTraS", 3.75 * mm, 13.5 * mm, 0.1 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicAlTra =
        new G4LogicalVolume(solidAlTra, // its solid
                            Al,         // its material
                            "AlTraLV"); // its name

    new G4PVPlacement(0,                                     // no rotation
                      G4ThreeVector(0, 0, (2.5 - 0.1) * mm), // at
                      logicAlTra,                            // its logical volume
                      "AlTraPV",                             // its name
                      logicCuTra,                            // its mother  volume
                      false,                                 // no boolean operations
                      0,                                     // copy number
                      fCheckOverlaps);                       // checking overlaps

    //*****************************************************************************************************
    // Cu Pipe
    //*****************************************************************************************************
    G4Tubs *solidCuPipe =
        new G4Tubs("CuPipeS", 3.75 * mm, 4.765 * mm, 14. * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicCuPipe =
        new G4LogicalVolume(solidCuPipe, // its solid
                            Cu,          // its material
                            "CuPipeLV"); // its name

    new G4PVPlacement(0,                                                     // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + 19. * mm)), // at
                      logicCuPipe,                                           // its logical volume
                      "CuPipePV",                                            // its name
                      logicWorld,                                            // its mother  volume
                      false,                                                 // no boolean operations
                      0,                                                     // copy number
                      fCheckOverlaps);                                       // checking overlaps

    G4Cons *solidCuCon =
        new G4Cons("CuConS", pRMin, 2.5 * mm, 0. * mm, 4.765 * mm, 0.75 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicCuCon =
        new G4LogicalVolume(solidCuCon, // its solid
                            Cu,         // its material
                            "CuConLV"); // its name

    new G4PVPlacement(0,                                                              // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + (33. + 0.75) * mm)), // at
                      logicCuCon,                                                     // its logical volume
                      "CuConPV",                                                      // its name
                      logicWorld,                                                     // its mother  volume
                      false,                                                          // no boolean operations
                      0,                                                              // copy number
                      fCheckOverlaps);                                                // checking overlaps

    //*****************************************************************************************************
    // XeArDME Gas in Cu Pipe
    //*****************************************************************************************************
    G4Tubs *solidGasTube1 =
        new G4Tubs("GasTube1S", pRMin, 3.75 * mm, 16.5 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicGasTube1 =
        new G4LogicalVolume(solidGasTube1, // its solid
                            XeArDME,       // its material
                            "GasTube1LV"); // its name

    new G4PVPlacement(0,                                                      // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + 16.5 * mm)), // at
                      logicGasTube1,                                          // its logical volume
                      "GasTube1PV",                                           // its name
                      logicWorld,                                             // its mother  volume
                      false,                                                  // no boolean operations
                      0,                                                      // copy number
                      fCheckOverlaps);                                        // checking overlaps

    G4Tubs *solidGasTube2 =
        new G4Tubs("GasTube2S", pRMin, 3.75 * mm, 4.0 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicGasTube2 =
        new G4LogicalVolume(solidGasTube2, // its solid
                            XeArDME,       // its material
                            "GasTube2LV"); // its name

    new G4PVPlacement(0,                                                     // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight - 4.0 * mm)), // at
                      logicGasTube2,                                         // its logical volume
                      "GasTube2PV",                                          // its name
                      logicAlChamber,                                        // its mother  volume
                      false,                                                 // no boolean operations
                      0,                                                     // copy number
                      fCheckOverlaps);                                       // checking overlaps

    //*****************************************************************************************************
    // Cu Electrode (Upper & Lower)
    //*****************************************************************************************************
    G4Tubs *solidCuElectrode1 =
        new G4Tubs("CuElectrode1S", pRMin, 1.6 * mm, 3. * mm, pSPhi, pDPhi);
    G4Tubs *solidCuElectrode2 =
        new G4Tubs("CuElectrode2S", pRMin, 3. * mm, 1. * mm, pSPhi, pDPhi);
    G4Tubs *solidCuElectrode3 =
        new G4Tubs("CuElectrode3S", pRMin, 3. * mm, 1.5 * mm, pSPhi, pDPhi);
    G4Box *solidCuElectrode4 =
        new G4Box("CuElectrode4S", 2.25 * mm, 1. * mm, 5.5 * mm);
    G4VSolid *solidCuElectrodeU =
        new G4UnionSolid("CuElectrodeUS", solidCuElectrode2, solidCuElectrode1,
                         0, G4ThreeVector(0, 0, 4. * mm));
    G4VSolid *solidCuElectrodeL =
        new G4UnionSolid("CuElectrodeLS", solidCuElectrode3, solidCuElectrode4,
                         0, G4ThreeVector(0, 0, -7. * mm));
    G4LogicalVolume *logicCuElectrodeU =
        new G4LogicalVolume(solidCuElectrodeU, // its solid
                            Cu,                // its material
                            "CuElectrodeULV"); // its name
    G4LogicalVolume *logicCuElectrodeL =
        new G4LogicalVolume(solidCuElectrodeL, // its solid
                            Cu,                // its material
                            "CuElectrodeLLV"); // its name

    //*****************************************************************************************************
    // Feedthrough　Peek
    //*****************************************************************************************************
    G4Tubs *solidPeekParts1 =
        new G4Tubs("PeekParts1S", 1.6 * mm, 2.75 * mm, 0.75 * mm, pSPhi, pDPhi);
    G4LogicalVolume *logicPeekParts1 =
        new G4LogicalVolume(solidPeekParts1, // its solid
                            Peek,            // its material
                            "PeekParts1LV"); // its name

    G4Tubs *solidPeekParts2 =
        new G4Tubs("PeekParts2S", 1.6 * mm, 4.0 * mm, 1.25 * mm, pSPhi, pDPhi);
    G4LogicalVolume *logicPeekParts2 =
        new G4LogicalVolume(solidPeekParts2, // its solid
                            Peek,            // its material
                            "PeekParts2LV"); // its name

    G4Box *solidPeekParts3Base =
        new G4Box("PeekParts3BaseS",               // its name
                  7.75 * mm, 7.75 * mm, 2.0 * mm); // its size
    G4Tubs *PeekParts3Frame =
        new G4Tubs("PeekParts3FrameS", pRMin, 9.5 * mm, 2.0 * mm, pSPhi, pDPhi);
    G4VSolid *solidPeekParts3Fit =
        new G4IntersectionSolid("PeekParts3FitS", solidPeekParts3Base, PeekParts3Frame);
    G4VSolid *solidPeekParts3 =
        new G4SubtractionSolid("PeekParts3S", solidPeekParts3Fit, solidCuElectrodeU,
                               0, G4ThreeVector(0, 0, -1. * mm));
    G4LogicalVolume *logicPeekParts3 =
        new G4LogicalVolume(solidPeekParts3, // its solid
                            Peek,            // its material
                            "PeekParts3LV"); // its name

    G4Box *solidPeekParts4Base =
        new G4Box("PeekParts4BaseS",               // its name
                  10.0 * mm, 10.0 * mm, 3.0 * mm); // its size
    G4VSolid *PeekParts4 =
        new G4SubtractionSolid("PeekParts4S", solidPeekParts4Base, solidCuElectrodeL,
                               0, G4ThreeVector(0, 0, 1.5 * mm));
    G4LogicalVolume *logicPeekParts4 =
        new G4LogicalVolume(PeekParts4,      // its solid
                            Peek,            // its material
                            "PeekParts4LV"); // its name

    //*****************************************************************************************************
    // Feedthrough Screw (Screw & Head)
    //*****************************************************************************************************
    G4Tubs *solidScrewScrew =
        new G4Tubs("ScrewScrewS", pRMin, 1.0 * mm, 2.5 * mm, pSPhi, pDPhi);
    G4Tubs *solidScrewHead =
        new G4Tubs("ScrewHeadS", pRMin, 1.75 * mm, 0.5 * mm, pSPhi, pDPhi);
    G4LogicalVolume *logicScrewScrew =
        new G4LogicalVolume(solidScrewScrew, // its solid
                            SUS304,          // its material
                            "ScrewScrewLV"); // its name
    G4LogicalVolume *logicScrewHead =
        new G4LogicalVolume(solidScrewHead,  // its solid
                            SUS304,          // its material
                            "Screw_HeadLV"); // its name

    //*****************************************************************************************************
    // Chamber Ear
    //*****************************************************************************************************
    G4Box *solidChamberEar =
        new G4Box("ChamberEarS",                   // its name
                  3.75 * mm, 10.0 * mm, 2.0 * mm); // its size
    G4LogicalVolume *logicChamberEar =
        new G4LogicalVolume(solidChamberEar, // its solid
                            Al,              // its material
                            "ChamberEarLV"); // its name

    //*****************************************************************************************************
    // Helisert
    //*****************************************************************************************************
    G4Box *solidHelisert =
        new G4Box("HelisertS",                     // its name
                  16.15 * mm, 4.0 * mm, 0.5 * mm); // its size
    G4LogicalVolume *logicHelisert =
        new G4LogicalVolume(solidHelisert, // its solid
                            XeArDME,       // its material
                            "HelisertLV"); // its name

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // placement
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    G4double peekR = 38. * mm;
    G4double screwR = 7.0 * sqrt(2) * mm;
    G4double earR = pRMax + chamber_side_thickness + 3.75 * mm;
    G4double helisertR = 26.65 * mm;
    for (G4int p = 0; p < 4; p++)
    {
        G4double pT = ang45 + ang90 * p;
        G4double px = peekR * cos(pT);
        G4double py = peekR * sin(pT);
        G4double ex = earR * cos(pT);
        G4double ey = earR * sin(pT);
        G4double hT = ang90 * p;
        G4double hx = helisertR * cos(hT);
        G4double hy = helisertR * sin(hT);

        G4RotationMatrix *rotZV45 = new G4RotationMatrix;
        rotZV45->rotateZ(pT + ang90);
        G4RotationMatrix *rotZV90 = new G4RotationMatrix;
        rotZV90->rotateZ(hT);

        //
        // Cu Electrode
        //
        new G4PVPlacement(0,                                                    // no rotation
                          G4ThreeVector(px, py, -chamber_half_hight + 1. * mm), // at
                          logicCuElectrodeU,                                    // its logical volume
                          "CuElectrodeUPV",                                     // its name
                          logicAlChamber,                                       // its mother  volume
                          false,                                                // no boolean operations
                          0,                                                    // copy number
                          fCheckOverlaps);                                      // checking overlaps

        new G4PVPlacement(rotZV45,                                               // rotation
                          G4ThreeVector(px, py, -chamber_half_hight - 1.5 * mm), // at
                          logicCuElectrodeL,                                     // its logical volume
                          "CuElectrodeLPV",                                      // its name
                          logicWorld,                                            // its mother  volume
                          false,                                                 // no boolean operations
                          0,                                                     // copy number
                          fCheckOverlaps);                                       // checking overlaps

        //
        // Feedthrough　Peek
        //
        new G4PVPlacement(0,                                                                // no rotation
                          G4ThreeVector(px, py, -(chamber_half_hight - (6.5 + 0.75) * mm)), // at
                          logicPeekParts1,                                                  // its logical volume
                          "PeekParts1PV",                                                   // its name
                          logicAlChamber,                                                   // its mother  volume
                          false,                                                            // no boolean operations
                          0,                                                                // copy number
                          fCheckOverlaps);                                                  // checking overlaps

        new G4PVPlacement(0,                                                                // no rotation
                          G4ThreeVector(px, py, -(chamber_half_hight - (4.0 + 1.25) * mm)), // at
                          logicPeekParts2,                                                  // its logical volume
                          "PeekParts2PV",                                                   // its name
                          logicAlChamber,                                                   // its mother  volume
                          false,                                                            // no boolean operations
                          0,                                                                // copy number
                          fCheckOverlaps);                                                  // checking overlaps

        new G4PVPlacement(0,                                                       // no rotation
                          G4ThreeVector(px, py, -(chamber_half_hight - 2.0 * mm)), // at
                          logicPeekParts3,                                         // its logical volume
                          "PeekParts3PV",                                          // its name
                          logicAlChamber,                                          // its mother  volume
                          false,                                                   // no boolean operations
                          0,                                                       // copy number
                          fCheckOverlaps);                                         // checking overlap

        new G4PVPlacement(rotZV45,                                                 // rotation
                          G4ThreeVector(px, py, -(chamber_half_hight + 3.0 * mm)), // at
                          logicPeekParts4,                                         // its logical volume
                          "PeekParts4PV",                                          // its name
                          logicWorld,                                              // its mother  volume
                          false,                                                   // no boolean operations
                          0,                                                       // copy number
                          fCheckOverlaps);                                         // checking overlaps

        //
        // Chamber Ear
        //
        new G4PVPlacement(rotZV45,                            // rotation
                          G4ThreeVector(ex, ey, -14.75 * mm), // at
                          logicChamberEar,                    // its logical volume
                          "ChamberEarPV",                     // its name
                          logicWorld,                         // its mother  volume
                          false,                              // no boolean operations
                          0,                                  // copy number
                          fCheckOverlaps);                    // checking overlaps

        //
        // Helisert
        //
        new G4PVPlacement(rotZV90,                                                 // rotation
                          G4ThreeVector(hx, hy, -(chamber_half_hight - 7.5 * mm)), // at
                          logicHelisert,                                           // its logical volume
                          "HelisertPV",                                            // its name
                          logicAlChamber,                                          // its mother  volume
                          false,                                                   // no boolean operations
                          0,                                                       // copy number
                          fCheckOverlaps);                                         // checking overlaps
    }

    //
    // Feedthrough Screw
    //
    for (G4int i = 0; i < 4; i++)
    {
        G4double angV = ang45 + ang90 * i;
        G4double sx = screwR * cos(angV);
        G4double sy = screwR * sin(angV);

        new G4PVPlacement(0,                               // no rotation
                          G4ThreeVector(sx, sy, 0.5 * mm), // at
                          logicScrewScrew,                 // its logical volume
                          "ScrewScrewPV",                  // its name
                          logicPeekParts4,                 // its mother  volume
                          false,                           // no boolean operations
                          0,                               // copy number
                          fCheckOverlaps);

        new G4PVPlacement(0,                                // no rotation
                          G4ThreeVector(sx, sy, -2.5 * mm), // at
                          logicScrewHead,                   // its logical volume
                          "ScrewHeadPV",                    // its name
                          logicPeekParts4,                  // its mother  volume
                          false,                            // no boolean operations
                          0,                                // copy number
                          fCheckOverlaps);
    }

    //*************************************************************************************************
    // Mo & Sn Shield @ Chamber Top
    //*************************************************************************************************
    G4Tubs *solidMoShieldChamTopBase =
        new G4Tubs("MoShieldChamTopBaseS", winR + 1 * mm,
                   pRMax + chamber_side_thickness, Mo_thickness / 2, pSPhi, pDPhi);

    G4VSolid *solidMoShieldChamTop =
        new G4SubtractionSolid("MoShieldChamTopS", solidMoShieldChamTopBase, solidColliEar, 0,
                               G4ThreeVector(colli_earR, 0, 0));

    G4Tubs *solidSnShieldChamTopBase =
        new G4Tubs("SnShieldChamTopBaseS", winR + 1 * mm,
                   pRMax + chamber_side_thickness + Mo_thickness + chamber_side_Sn_thickness,
                   chamber_top_Sn_thickness / 2, pSPhi, pDPhi);

    G4VSolid *solidSnShieldChamTop =
        new G4SubtractionSolid("SnShieldChamTopS", solidSnShieldChamTopBase, solidColliEar, 0,
                               G4ThreeVector(colli_earR, 0, 0));

    for (G4int i = 1; i < 6; i++)
    {
        G4RotationMatrix *rotZV = new G4RotationMatrix;
        G4double angV = ang60 * i;
        rotZV->rotateZ(-angV);

        G4VSolid *solidMoShieldChamTopHollow =
            new G4SubtractionSolid("MoShieldChamTopHollowS", solidMoShieldChamTop, solidColliEar, rotZV,
                                   G4ThreeVector(colli_earR * cos(angV), colli_earR * sin(angV), 0));

        solidMoShieldChamTop = solidMoShieldChamTopHollow;

        G4VSolid *solidSnShieldChamTopHollow =
            new G4SubtractionSolid("SnShieldChamTopHollowS", solidSnShieldChamTop, solidColliEar, rotZV,
                                   G4ThreeVector(colli_earR * cos(angV),
                                                 colli_earR * sin(angV), 0));

        solidSnShieldChamTop = solidSnShieldChamTopHollow;
    }

    G4LogicalVolume *logicMoShieldChamTop = new G4LogicalVolume(solidMoShieldChamTop, // its solid
                                                                Mo,                   // its material
                                                                "MoShieldChamTopLV"); // its name

    new G4PVPlacement(0,                                                          // no rotation
                      G4ThreeVector(0, 0, chamber_half_hight + Mo_thickness / 2), // at
                      logicMoShieldChamTop,                                       // its logical volume
                      "MoShieldChamTopPV",                                        // its name
                      logicWorld,                                                 // its mother  volume
                      false,                                                      // no boolean operations
                      0,                                                          // copy number
                      fCheckOverlaps);                                            // checking overlaps

    G4LogicalVolume *logicSnShieldChamTop = new G4LogicalVolume(solidSnShieldChamTop, // its solid
                                                                SnSb,                 // its material
                                                                "SnShieldChamTopLV"); // its name

    new G4PVPlacement(0, // no rotation
                      G4ThreeVector(0, 0,
                                    chamber_half_hight + Mo_thickness + chamber_top_Sn_thickness / 2),
                      logicSnShieldChamTop, // its logical volume
                      "SnShieldChamTopPV",  // its name
                      logicWorld,           // its mother  volume
                      false,                // no boolean operations
                      0,                    // copy number
                      fCheckOverlaps);      // checking overlaps

    //*************************************************************************************************
    // Mo Shield @ Chamber Side
    //*************************************************************************************************
    G4Tubs *solidMoShieldChamSide =
        new G4Tubs("MoShieldChamSideS", pRMax + chamber_side_thickness,
                   pRMax + chamber_side_thickness + Mo_thickness, 9.75 * mm, pSPhi, pDPhi);

    G4LogicalVolume *logicMoShieldChamSide =
        new G4LogicalVolume(solidMoShieldChamSide, // its solid
                            Mo,                    // its material
                            "MoShieldChamSideLV"); // its name

    new G4PVPlacement(0,                              // no rotation
                      G4ThreeVector(0, 0, -3.0 * mm), // at
                      logicMoShieldChamSide,          // its logical volume
                      "MoShieldChamSidePV",           // its name
                      logicWorld,                     // its mother  volume
                      false,                          // no boolean operations
                      0,                              // copy number
                      fCheckOverlaps);                // checking overlaps

    //*************************************************************************************************
    // Sn Shield @ Chamber Side
    //*************************************************************************************************
    G4Tubs *solidSnShieldChamSideBase =
        new G4Tubs("SnShieldChamSideBaseS",
                   pRMax + chamber_side_thickness + Mo_thickness,
                   pRMax + chamber_side_thickness + Mo_thickness + chamber_side_Sn_thickness,
                   chamber_half_hight, pSPhi, pDPhi);

    G4RotationMatrix *rotZ135 = new G4RotationMatrix;
    rotZ135->rotateZ(-ang45);
    G4double cham_ear_rotR = earR - 0.25 * mm;
    G4double cham_ear_z = -chamber_half_hight + 2.0 * mm;

    G4VSolid *solidSnShiledChamSide =
        new G4SubtractionSolid("SnShiledChamSideS", solidSnShieldChamSideBase, solidChamberEar, rotZ135,
                               G4ThreeVector(cham_ear_rotR * cos(ang45),
                                             cham_ear_rotR * sin(ang45), cham_ear_z));

    for (G4int i = 1; i < 4; i++)
    {
        G4RotationMatrix *rotZV = new G4RotationMatrix;
        rotZV->rotateZ(ang45 * pow(-1, (i + 1)));
        G4double angV = ang45 + ang90 * i;
        G4VSolid *solidSnShiledChamSideRe =
            new G4SubtractionSolid("SnShiledChamSideReS", solidSnShiledChamSide, solidChamberEar, rotZV,
                                   G4ThreeVector(cham_ear_rotR * cos(angV), cham_ear_rotR * sin(angV), cham_ear_z));
        solidSnShiledChamSide = solidSnShiledChamSideRe;
    }

    G4LogicalVolume *logicSnShieldChamSide =
        new G4LogicalVolume(solidSnShiledChamSide, // its solid
                            SnSb,                  // its material
                            "SnShiledChamSideLV"); // its name

    new G4PVPlacement(0,                      // no rotation
                      G4ThreeVector(0, 0, 0), // at (0,0,0)
                      logicSnShieldChamSide,  // its logical volume
                      "SnShiledChamSideSnPV", // its name
                      logicWorld,             // its mother  volume
                      false,                  // no boolean operations
                      0,                      // copy number
                      fCheckOverlaps);        // checking overlaps

    //*************************************************************************************************
    // Mo & Sn Shield @ Chamber Bottom
    //*************************************************************************************************
    G4Box *solidPeekParts4Model =
        new G4Box("PeekParts4ModelS",                                          // its name
                  10.0 * mm, 10.0 * mm, chamber_bottom_Sn_thickness / 2 * mm); // its size
    G4RotationMatrix *rotC45 = new G4RotationMatrix;
    rotC45->rotateZ(ang45);

    G4Tubs *solidMoShieldChamBottomBase =
        new G4Tubs("MoShieldChamBottomBaseS", 13.5 * mm, 42.5 * mm, Mo_thickness / 2, pSPhi, pDPhi);
    G4VSolid *solidMoShieldChamBottom =
        new G4SubtractionSolid("MoShieldChamBottomS", solidMoShieldChamBottomBase,
                               solidPeekParts4Model, rotC45,
                               G4ThreeVector(peekR * cos(ang45), peekR * sin(ang45), 0));

    G4Tubs *solidSnShieldChamBottomBase =
        new G4Tubs("SnShieldChamBottomBaseS", 13.5 * mm, 42.5 * mm, chamber_bottom_Sn_thickness / 2, pSPhi, pDPhi);
    G4VSolid *solidSnShieldChamBottom =
        new G4SubtractionSolid("SnShieldChamBottomS", solidSnShieldChamBottomBase, solidPeekParts4Model, rotC45,
                               G4ThreeVector(peekR * cos(ang45), peekR * sin(ang45), 0));

    for (G4int i = 1; i < 4; i++)
    {
        G4double angV = ang45 + ang90 * i;
        G4VSolid *solidMoShieldChamBottomHollow =
            new G4SubtractionSolid("MoShieldChamBottomHollowS",
                                   solidMoShieldChamBottom, solidPeekParts4Model, rotC45,
                                   G4ThreeVector(peekR * cos(angV), peekR * sin(angV), 0));

        solidMoShieldChamBottom = solidMoShieldChamBottomHollow;

        G4VSolid *solidSnShieldChamBottomHollow =
            new G4SubtractionSolid("SnShieldChamBottomHollowS", solidSnShieldChamBottom,
                                   solidPeekParts4Model, rotC45,
                                   G4ThreeVector(peekR * cos(angV), peekR * sin(angV), 0));

        solidSnShieldChamBottom = solidSnShieldChamBottomHollow;
    }

    G4LogicalVolume *logicMoShieldChamBottom = new G4LogicalVolume(solidMoShieldChamBottom, // its solid
                                                                   Mo,                      // its material
                                                                   "MoShieldChamBottomLV"); // its name

    new G4PVPlacement(0,                                                             // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + Mo_thickness / 2)), // at
                      logicMoShieldChamBottom,                                       // its logical volume
                      "MoShieldChamBottomPV",                                        // its name
                      logicWorld,                                                    // its mother  volume
                      false,                                                         // no boolean operations
                      0,                                                             // copy number
                      fCheckOverlaps);                                               // checking overlaps

    G4LogicalVolume *logicSnShieldChamBottom = new G4LogicalVolume(solidSnShieldChamBottom, // its solid
                                                                   Sn,                      // its material
                                                                   "SnShieldChamBottomLV"); // its name

    new G4PVPlacement(0, // no rotation
                      G4ThreeVector(0, 0,
                                    -(chamber_half_hight + Mo_thickness + chamber_bottom_Sn_thickness / 2)),
                      logicSnShieldChamBottom, // its logical volume
                      "SnShieldChamBottomPV",  // its name
                      logicWorld,              // its mother  volume
                      false,                   // no boolean operations
                      0,                       // copy number
                      fCheckOverlaps);         // checking overlaps

    //*************************************************************************************************
    // Base Plate
    //*************************************************************************************************
    const G4double baseplate_size = 91.0 * mm;
    G4Box *solidBasePlateBase =
        new G4Box("BasePlateBaseS",                                  // its name
                  baseplate_size / 2, baseplate_size / 2, 3.5 * mm); // its size
    G4Tubs *solidBasePlateCircle =
        new G4Tubs("BasePlateCircleS", pRMin, 42.5 * mm, 3.5 * mm, pSPhi, pDPhi);
    G4Box *solidBasePlateBox =
        new G4Box("BasePlateBaseBoxS",                                 // its name
                  10.0 * mm + 0.5 * mm, 10.0 * mm + 0.5 * mm, 5 * mm); // its size

    G4VSolid *solidBasePlate =
        new G4SubtractionSolid("BasePlateS", solidBasePlateBase, solidBasePlateCircle,
                               0, G4ThreeVector(0, 0, 0));

    for (G4int i = 0; i < 4; i++)
    {
        G4double angV = ang45 + ang90 * i;
        G4VSolid *solidBasePlateHollow =
            new G4SubtractionSolid("BasePlateHollowS", solidBasePlate, solidBasePlateBox, rotC45,
                                   G4ThreeVector(peekR * cos(angV),
                                                 peekR * sin(angV), 0));

        solidBasePlate = solidBasePlateHollow;
    }

    G4LogicalVolume *logicBasePlate = new G4LogicalVolume(solidBasePlate, // its solid
                                                          Al,             // its material
                                                          "BasePlateLV"); // its name

    new G4PVPlacement(0,                                                     // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + 3.5 * mm)), // at
                      logicBasePlate,                                        // its logical volume
                      "BasePlatePV",                                         // its name
                      logicWorld,                                            // its mother  volume
                      false,                                                 // no boolean operations
                      0,                                                     // copy number
                      fCheckOverlaps);                                       // checking overlaps

    //*************************************************************************************************
    // FEC
    //*************************************************************************************************
    G4Box *solidFECBase =
        new G4Box("FECBaeS",                       // its name
                  42.0 * mm, 44.5 * mm, 0.8 * mm); // its size
    G4Tubs *solidFECHole1 =
        new G4Tubs("FECHole1S", pRMin, 10.0 * mm, 1.0 * mm, pSPhi, pDPhi);
    G4Tubs *solidFECHole2 =
        new G4Tubs("FECHole2S", pRMin, 5.0 * mm, 1.0 * mm, pSPhi, pDPhi);
    G4VSolid *solidFEC =
        new G4SubtractionSolid("FECS", solidFECBase, solidFECHole1, 0, G4ThreeVector(0, 0, 0));
    for (G4int i = 0; i < 4; i++)
    {
        G4double angV = ang45 + ang90 * i;
        G4VSolid *solidFECHollow =
            new G4SubtractionSolid("FECHollowS", solidFEC, solidFECHole2, 0,
                                   G4ThreeVector(peekR * cos(angV), peekR * sin(angV), 0));

        solidFEC = solidFECHollow;
    }

    G4LogicalVolume *logicFEC = new G4LogicalVolume(solidFEC, // its solid
                                                    Si,       // its material
                                                    "FECLV"); // its name

    new G4PVPlacement(0,                                                             // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + (7.0 + 0.8) * mm)), // at
                      logicFEC,                                                      // its logical volume
                      "FECPV",                                                       // its name
                      logicWorld,                                                    // its mother  volume
                      false,                                                         // no boolean operations
                      0,                                                             // copy number
                      fCheckOverlaps);                                               // checking overlaps

    //*************************************************************************************************
    // EMI Shield
    //*************************************************************************************************
    G4Box *solidEMIBase =
        new G4Box("EMIBoxS",                                         // its name
                  baseplate_size / 2, baseplate_size / 2, 9.5 * mm); // its size
    G4Box *solidEMIBoxHollow =
        new G4Box("EMIBoxHollowS",                 // its name
                  42.0 * mm, 44.5 * mm, 9.0 * mm); // its size

    G4Box *solidEMIBltomHole_square =
        new G4Box("EMIBottomHole_squareS",          // its name
                  4.25 * mm, 10.0 * mm, 0.51 * mm); // its size

    G4Tubs *solidEMIBottomHole =
        new G4Tubs("EMIBottomHoleS", pRMin, 10.0 * mm, 10.0 * mm, pSPhi, pDPhi);
    G4VSolid *solidEMIHollow1 =
        new G4SubtractionSolid("EMIHollow1S", solidEMIBase, solidEMIBoxHollow,
                               0, G4ThreeVector(0, 0, 0.5 * mm));
    G4VSolid *solidEMIHollow2 =
        new G4SubtractionSolid("EMIHollow2S", solidEMIHollow1, solidEMIBottomHole,
                               0, G4ThreeVector(0, 0, 0));
    G4VSolid *solidEMIHollow3 =
        new G4SubtractionSolid("EMIHollow3S", solidEMIHollow2, solidEMIBltomHole_square,
                               0, G4ThreeVector(-41.25 * mm, 22.0 * mm, -9.0 * mm));

    G4VSolid *solidEMI =
        new G4SubtractionSolid("EMIS", solidEMIHollow3, solidEMIBltomHole_square,
                               0, G4ThreeVector(-41.25 * mm, -22.0 * mm, -9.0 * mm));

    G4LogicalVolume *logicEMI = new G4LogicalVolume(solidEMI, // its solid
                                                    Al,       // its material
                                                    "EMILV"); // its name

    new G4PVPlacement(0,                                                             // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + (7.0 + 9.5) * mm)), // at
                      logicEMI,                                                      // its logical volume
                      "EMIPV",                                                       // its name
                      logicWorld,                                                    // its mother  volume
                      false,                                                         // no boolean operations
                      0,                                                             // copy number
                      fCheckOverlaps);                                               // checking overlaps

    //*************************************************************************************************
    // Sn Shield @ Base Plate Side
    //*************************************************************************************************
    G4Box *solidSnShieldBPSide =
        new G4Box("SnShieldBPSideS", // its name
                  35.0 * mm,
                  baseplate_side_Sn_thickness / 2, 3.5 * mm);
    G4LogicalVolume *logicSnShieldBPSide = new G4LogicalVolume(solidSnShieldBPSide, // its solid
                                                               SnSb,                // its material
                                                               "SnShieldBPSideLV"); // its name
    G4double bp_shieldR = baseplate_size / 2 + baseplate_side_Sn_thickness / 2;
    for (G4int i = 0; i < 4; i++)
    {
        G4double angV = ang90 * i;
        G4RotationMatrix *rotZV = new G4RotationMatrix;
        rotZV->rotateZ(angV + ang90);

        new G4PVPlacement(rotZV, // rotation
                          G4ThreeVector(bp_shieldR * cos(angV), bp_shieldR * sin(angV), -(chamber_half_hight + 3.5 * mm)),
                          logicSnShieldBPSide, // its logical volume
                          "SnShieldBPSidePV",  // its name
                          logicWorld,          // its mother  volume
                          false,               // no boolean operations
                          0,                   // copy number
                          fCheckOverlaps);     // checking overlaps
    }

    //*************************************************************************************************
    // Sn Shield @ EMI Shield
    //*************************************************************************************************
    G4Box *solidSnShieldEMISide =
        new G4Box("SnShieldEMISideS", // its name
                  41.9 * mm,
                  EMI_side_Sn_thickness / 2, 8.7 * mm);
    G4LogicalVolume *logicSnShieldEMISide1 = new G4LogicalVolume(solidSnShieldEMISide, // its solid
                                                                 SnSb,                 // its material
                                                                 "SnShieldEMISideLV"); // its name
    G4LogicalVolume *logicSnShieldEMISide2 = new G4LogicalVolume(solidSnShieldEMISide, // its solid
                                                                 SnSb,                 // its material
                                                                 "SnShieldEMISideLV"); // its name
    G4Box *solidSnShieldEMISideHole =
        new G4Box("SnShieldEMISideHoleS", // its name
                  3.0 * mm,
                  EMI_side_Sn_thickness / 2, 4.3 * mm);
    G4LogicalVolume *logicSnShieldEMISideHole = new G4LogicalVolume(solidSnShieldEMISideHole, // its solid
                                                                    Al,                       // its material
                                                                    "SnShieldEMISideHoleLV"); // its name

    for (G4int i = 1; i < 3; i++)
    {
        new G4PVPlacement(rotZ90, // rotation
                          G4ThreeVector(pow(-1, i) * bp_shieldR, 0, -(chamber_half_hight + (7.0 + 8.7) * mm)),
                          logicSnShieldEMISide1, // its logical volume
                          "SnShieldEMISide1PV",  // its name
                          logicWorld,            // its mother  volume
                          false,                 // no boolean operations
                          0,                     // copy number
                          fCheckOverlaps);       // checking overlaps
        new G4PVPlacement(0,                     // rotation
                          G4ThreeVector(0, pow(-1, i) * bp_shieldR, -(chamber_half_hight + (7.0 + 8.7) * mm)),
                          logicSnShieldEMISide2, // its logical volume
                          "SnShieldEMISide2PV",  // its name
                          logicWorld,            // its mother  volume
                          false,                 // no boolean operations
                          0,                     // copy number
                          fCheckOverlaps);       // checking overlaps
    }
    new G4PVPlacement(0, // rotation
                      G4ThreeVector(0, 0, -4.4 * mm),
                      logicSnShieldEMISideHole, // its logical volume
                      "SnShieldEMISideHolePV",  // its name
                      logicSnShieldEMISide1,    // its mother  volume
                      false,                    // no boolean operations
                      0,                        // copy number
                      fCheckOverlaps);          // checking overlaps

    G4Box *solidSnShieldEMIBottomBase =
        new G4Box("SnShieldEMIBottomBaseS",                           // its name
                  45.0 * mm, 45.0 * mm, EMI_bottom_Sn_thickness / 2); // its size
    G4VSolid *solidSnShieldEMIBottomHollow1 =
        new G4SubtractionSolid("SnShieldEMIBottomHollow1S", solidSnShieldEMIBottomBase, solidEMIBottomHole,
                               0, G4ThreeVector(0, 0, 0));
    G4VSolid *solidSnShieldEMIBottomHollow2 =
        new G4SubtractionSolid("SnShieldEMIBottomHollow2S", solidSnShieldEMIBottomHollow1, solidEMIBltomHole_square,
                               0, G4ThreeVector(-41.25 * mm, 22.0 * mm, 0));
    G4VSolid *solidSnShieldEMIBottom =
        new G4SubtractionSolid("solidSnShieldEMIBottom", solidSnShieldEMIBottomHollow2, solidEMIBltomHole_square,
                               0, G4ThreeVector(-41.25 * mm, -22.0 * mm, 0));

    G4LogicalVolume *logicSnShieldEMIBottom = new G4LogicalVolume(solidSnShieldEMIBottom, // its solid
                                                                  SnSb,                   // its material
                                                                  "SnShieldEMIBottomLV"); // its name

    new G4PVPlacement(0,                                                                                            // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + (7.0 + 19.0 + EMI_bottom_Sn_thickness / 2) * mm)), // at
                      logicSnShieldEMIBottom,                                                                       // its logical volume
                      "SnShieldEMIBottomPV",                                                                        // its name
                      logicWorld,                                                                                   // its mother  volume
                      false,                                                                                        // no boolean operations
                      0,                                                                                            // copy number
                      fCheckOverlaps);                                                                              // checking overlaps

    G4Box *solidSnShieldEMIBottomHole1 =
        new G4Box("SnShieldEMIBottomHole1S",                                                        // its name
                  2.0 * mm, 2.0 * mm, EMI_bottom_Sn_thickness / 2);                                 // its size
    G4LogicalVolume *logicSnShieldEMIBottomHole1 = new G4LogicalVolume(solidSnShieldEMIBottomHole1, // its solid
                                                                       Vacuum,                      // its material
                                                                       "SnShieldEMIBottomHole1LV"); // its name
    for (G4int i = 0; i < 4; i++)
    {
        G4double EMI_holeR = (sqrt(2) * 45 - (sqrt(2) * 4) * 1 / 2) * mm;
        G4double angV = ang90 * i + ang45;
        new G4PVPlacement(0, // rotation
                          G4ThreeVector(EMI_holeR * cos(angV), EMI_holeR * sin(angV), 0),
                          logicSnShieldEMIBottomHole1, // its logical volume
                          "SnShieldEMIBottomHolePV",   // its name
                          logicSnShieldEMIBottom,      // its mother  volume
                          false,                       // no boolean operations
                          0,                           // copy number
                          fCheckOverlaps);             // checking overlaps
    }
    G4Box *solidSnShieldEMIBottomHole2 =
        new G4Box("SnShieldEMIBottomHole2S",                                                        // its name
                  2.0 * mm, 3.5 * mm, EMI_bottom_Sn_thickness / 2);                                 // its size
    G4LogicalVolume *logicSnShieldEMIBottomHole2 = new G4LogicalVolume(solidSnShieldEMIBottomHole2, // its solid
                                                                       Vacuum,                      // its material
                                                                       "SnShieldEMIBottomHole2LV"); // its name
    for (G4int i = 1; i < 3; i++)
    {
        new G4PVPlacement(0, // rotation
                          G4ThreeVector(pow(-1, i) * 43.0 * mm, 0, 0),
                          logicSnShieldEMIBottomHole2, // its logical volume
                          "SnShieldEMISide2PV",        // its name
                          logicSnShieldEMIBottom,      // its mother  volume
                          false,                       // no boolean operations
                          0,                           // copy number
                          fCheckOverlaps);             // checking overlaps
    }

    //*************************************************************************************************
    // DAQ
    //*************************************************************************************************
    G4Box *solidDAQBase =
        new G4Box("DAQBaseS",                                     // its name
                  45.0 * mm, 45.0 * mm, 0.8 * mm);                // its size
    G4LogicalVolume *logicDAQ = new G4LogicalVolume(solidDAQBase, // its solid
                                                    Si,           // its material
                                                    "DAQLV");     // its name

    new G4PVPlacement(0,                                                      // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + 38.7 * mm)), // at
                      logicDAQ,                                               // its logical volume
                      "DAQPV",                                                // its name
                      logicWorld,                                             // its mother  volume
                      false,                                                  // no boolean operations
                      0,                                                      // copy number
                      fCheckOverlaps);                                        // checking overlaps

    //*************************************************************************************************
    // DAQ housing
    //*************************************************************************************************
    G4Box *solidDAQHousingBase =
        new G4Box("DAQHousingBaseS",                // its name
                  47.0 * mm, 47.0 * mm, 9.25 * mm); // its size
    G4Box *solidDAQHousingBoxHollow =
        new G4Box("DAQHousingBoxHollowS",           // its name
                  46.0 * mm, 46.0 * mm, 8.75 * mm); // its size
    G4Box *solidDAQHousingHole =
        new G4Box("DAQHousingHoleS",              // its name
                  37.0 * mm, 5.0 * mm, 2.5 * mm); // its size

    G4VSolid *solidDAQHousinHollow1 =
        new G4SubtractionSolid("DAQHousinHollow1S", solidDAQHousingBase, solidDAQHousingBoxHollow,
                               0, G4ThreeVector(0, 0, 0.5 * mm));
    G4VSolid *solidDAQHousinHollow2 =
        new G4SubtractionSolid("DAQHousinHollow2S", solidDAQHousinHollow1, solidDAQHousingHole,
                               0, G4ThreeVector(0, -42.0 * mm, -6.75 * mm));

    G4LogicalVolume *logicDAQHousing = new G4LogicalVolume(solidDAQHousinHollow2, // its solid
                                                           Al,                    // its material
                                                           "DAQHousingLV");       // its name

    new G4PVPlacement(0,                                                       // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + 35.25 * mm)), // at
                      logicDAQHousing,                                         // its logical volume
                      "DAQHousingPV",                                          // its name
                      logicWorld,                                              // its mother  volume
                      false,                                                   // no boolean operations
                      0,                                                       // copy number
                      fCheckOverlaps);                                         // checking overlaps

    G4Box *solidDAQHousingBottom =
        new G4Box("DAQHousingBottomS",                                                  // its name
                  20.0 * mm, 14.0 * mm, 2.0 * mm);                                      // its size
    G4LogicalVolume *logicDAQHousingBottom = new G4LogicalVolume(solidDAQHousingBottom, // its solid
                                                                 Al,                    // its material
                                                                 "DAQHousingBottomLV"); // its name
    new G4PVPlacement(0,                                                                // no rotation
                      G4ThreeVector(0, 0, -(chamber_half_hight + 46.5 * mm)),           // at
                      logicDAQHousingBottom,                                            // its logical volume
                      "DAQHousingBottomPV",                                             // its name
                      logicWorld,                                                       // its mother  volume
                      false,                                                            // no boolean operations
                      0,                                                                // copy number
                      fCheckOverlaps);                                                  // checking overlaps

    //*************************************************************************************************
    // Visualization attributes
    //*************************************************************************************************
    logicWorld->SetVisAttributes(WhiteVisAtt);
    logicCollimeter->SetVisAttributes(GrayVisAtt);
    logicHexAperture->SetVisAttributes(GrayVisAtt);
    logicColliEar->SetVisAttributes(GrayVisAtt);
    logicGasCell->SetVisAttributes(BlueVisAtt);
    logicAlChamber->SetVisAttributes(AlVisAtt);
    logicSensitiveInner->SetVisAttributes(BlueVisAtt);
    logicSensitiveOuter->SetVisAttributes(BlueVisAtt);
    logicBeWindow->SetVisAttributes(BeVisAtt);
    logicLWApertureRing->SetVisAttributes(OrangeVisAtt);
    logicLWAperture->SetVisAttributes(OrangeVisAtt);
    logicGEMCu->SetVisAttributes(CuVisAtt);
    logicGEMLCP->SetVisAttributes(PEEKVisAtt);
    logicReadoutPad->SetVisAttributes(OrangeVisAtt);
    logicTiTra->SetVisAttributes(TiVisAtt);
    logicSUS304TraU->SetVisAttributes(GrayVisAtt);
    logicSUS304TraL->SetVisAttributes(GrayVisAtt);
    logicSUS310Plate_Al->SetVisAttributes(AlVisAtt);
    logicSUS310Plate->SetVisAttributes(SUS310VisAtt);
    logicAlTra->SetVisAttributes(AlVisAtt);
    logicCuPipe->SetVisAttributes(CuVisAtt);
    logicCuTra->SetVisAttributes(CuVisAtt);
    logicGasTube1->SetVisAttributes(BlueVisAtt);
    logicGasTube2->SetVisAttributes(BlueVisAtt);
    logicHelisert->SetVisAttributes(BlueVisAtt);
    logicPeekParts1->SetVisAttributes(PEEKVisAtt);
    logicPeekParts2->SetVisAttributes(PEEKVisAtt);
    logicPeekParts3->SetVisAttributes(PEEKVisAtt);
    logicPeekParts4->SetVisAttributes(PEEKVisAtt);
    logicCuCon->SetVisAttributes(CuVisAtt);
    logicBasePlate->SetVisAttributes(AlVisAtt);
    logicFEC->SetVisAttributes(DeepGreenVisAtt);
    logicSnShieldEMISideHole->SetVisAttributes(AlVisAtt);
    logicDAQ->SetVisAttributes(DeepGreenVisAtt);
    logicDAQHousing->SetVisAttributes(AlVisAtt);
    logicDAQHousingBottom->SetVisAttributes(AlVisAtt);
    logicEMI->SetVisAttributes(AlVisAtt);
    logicSnShieldEMIBottomHole1->SetVisAttributes(AlVisAtt);
    logicSnShieldEMIBottomHole2->SetVisAttributes(AlVisAtt);
    logicCuElectrodeU->SetVisAttributes(CuVisAtt);
    logicCuElectrodeL->SetVisAttributes(CuVisAtt);
    logicScrewScrew->SetVisAttributes(GrayVisAtt);
    logicScrewHead->SetVisAttributes(GrayVisAtt);
    logicChamberEar->SetVisAttributes(AlVisAtt);

    return physWorld;
}
