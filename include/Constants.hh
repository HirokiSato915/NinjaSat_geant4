#ifndef Constants_h
#define Constants_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
//
const G4double worldLength = 0.3 * m;
const G4double pRMin = 0. * cm;
const G4double pRMax = 4.2 * cm;
const G4double chamber_half_hight = 16.75 * mm;
const G4double pSPhi = 0. * deg;
const G4double pDPhi = 360 * deg;
const G4double winR = 3.5 * cm;
const G4double Be_window_thickness = 0.15 * mm;
const G4double Be_window_R = 4.5 * cm;
const G4double chamber_side_thickness = 5.5 * mm;
const G4double chamber_top_hight = 6.7 * mm; // the length from chamber top to Be window bottom
const G4double sensitive_layer_thickness = 15.45 * mm;
const G4double OuterR_thickness = 8.5 * mm;
const G4double OuterR = 33.5 * mm;
const G4double InnerR = OuterR - OuterR_thickness;
const G4double Mo_thickness = 0.05 * mm;
const G4double chamber_top_Sn_thickness = 0.5 * mm;
const G4double chamber_side_Sn_thickness = 0.8 * mm;
const G4double chamber_bottom_Sn_thickness = 0.5 * mm;
const G4double baseplate_side_Sn_thickness = 0.5 * mm;
const G4double EMI_side_Sn_thickness = 0.5 * mm;
const G4double EMI_bottom_Sn_thickness = 0.5 * mm;
const G4double ang45 = 45. * deg;
const G4double ang60 = 60. * deg;
const G4double ang90 = 90. * deg;
//

#endif