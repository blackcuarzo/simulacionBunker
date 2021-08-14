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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4EllipticalTube.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4DisplacedSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),  
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  //
  G4bool checkOverlaps = true;
  G4bool conPared = true;
  G4bool conFantoma1 = true;



// Elementos

  G4double density,      // density
    a,                   // atomic mass
    z;                   // atomic number
  G4String name,         // name
    symbol;              // symbol
  G4int ncomponents;     // n components

  G4Element* Si = new G4Element
    (name="Silicon",symbol="Si" , z= 14., a=28.09*g/mole);
  G4Element* O  = new G4Element
    (name="Oxygen"  ,symbol="O" , z= 8., a=16.00*g/mole);
  G4Element* H = new G4Element
    (name="Hydrogen",symbol="H" , z= 1., a=1.00794*g/mole);
  G4Element* Ca = new G4Element
    (name="Calcium",symbol="Ca" , z= 20., a=40.078*g/mole);
  G4Element* Al = new G4Element
    (name="Aluminium"  ,symbol="Al" , z= 13., a=26.98*g/mole);  
  G4Element* Fe = new G4Element
    (name="Iron"  ,symbol="Fe" , z= 26., a=55.85*g/mole);  
 

  G4Element* elO = new G4Element(name="Oxygen",symbol="O2",z=8.,a=16.*g/mole);
  G4Element* elSi=new G4Element(name="Silicon",symbol="Si",z=14.,a=28.085*g/mole);
  G4Element* elCa=new G4Element(name="Calcium",symbol="Ca",z=20.,a=40.08*g/mole);
  G4Element* elNa=new G4Element(name="Sodium",symbol="Na",z=11.,a=22.99*g/mole);
  G4Element* elFe=new G4Element(name="Iron",symbol="Fe",z=26.,a=55.850*g/mole);
  G4Element* elAl=new G4Element(name="Aluminium",symbol="Al",z=13.,a=26.98*g/mole);

  G4Element* elH=new G4Element(name="Hydrogen",symbol="H2",z=1.,a = 1.01*g/mole);
  G4Element* elC=new G4Element(name="Carbon",symbol="C",z=6.,a = 12.01*g/mole);
  G4Element* elN=new G4Element(name="Nitrogen",symbol="N2",z=7.,a = 14.01*g/mole);

// Materiales 
  G4Material* concreto1 = new G4Material
    (name="Concrete", density=2.3*g/cm3, ncomponents=6);

  concreto1->AddElement(Si, 0.227915);
  concreto1->AddElement(O, 0.60541);
  concreto1->AddElement(H, 0.09972);
  concreto1->AddElement(Ca, 0.04986);
  concreto1->AddElement(Al, 0.014245);
  concreto1->AddElement(Fe, 0.00285);


  G4Material* water = new G4Material(name="water", density=1.00*g/cm3, ncomponents=2);
  water->AddElement(H , 2);
  water->AddElement(O , 1);

  G4double fractionmass;
  G4Material* concreto2 = new G4Material
    (name="concreto2",density = 2.35*g/cm3,ncomponents=6); // Se escoje una densidad de 2.35
  concreto2->AddElement(elO, fractionmass = 0.52);
  concreto2->AddElement(elSi, fractionmass = 0.325);
  concreto2->AddElement(elCa, fractionmass = 0.06);
  concreto2->AddElement(elNa, fractionmass = 0.015);
  concreto2->AddElement(elFe, fractionmass = 0.04);
  concreto2->AddElement(elAl, fractionmass = 0.04);
  
  G4Material* tejidoICRU = new G4Material (name="tejidoICRU", density = 1.00*g/cm3, ncomponents = 4);

  tejidoICRU -> AddElement (elH, 0.101174);
  tejidoICRU -> AddElement (elC, 0.111000);
  tejidoICRU -> AddElement (elN, 0.026000);
  tejidoICRU -> AddElement (elO, 0.761826);



// World
  //
  G4double world_size = 10.*m;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_size, world_size, world_size);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
// Concreto

// Pared Concreto

  G4double halfZPared = (157./2.0)*cm; //  1 TVL concreto 37cm, TVLe 33cm
  G4double halfXPared = (3./2.0)*m; 
  G4double halfYPared = (3./2.0)*m; // La mitad de la verdadera altura
  G4Box* pared = new G4Box("ParedConcreto", halfXPared, halfYPared, halfZPared);
  
  G4Material* materialPared = concreto2; 
  G4ThreeVector posicionZPared = G4ThreeVector(0.*cm, 0.*cm, - halfZPared + 5.*m);
  G4RotationMatrix* rotPared = new G4RotationMatrix();
  //rotPared->rotateX(M_PI/2.*rad);
  G4Transform3D trPared = G4Translate3D(posicionZPared)*G4Rotate3D(*rotPared);
  
  /*
  if (conPared == true){
  
    G4LogicalVolume* logicParedConcreto =                         
      new G4LogicalVolume(pared,         //its solid
                        materialPared,          //its material
                        "BarreraConcreto");    // Este nombre sirve para vis.mac
          
    new G4PVPlacement (trPared,                      //Traslación 
                    logicParedConcreto,              //its logical volume
                    "Barrera de Concreto",           //its name
                    logicWorld,                      //its mother  volume
                    false,                           //no boolean operation
                    0,                               //copy number
                    checkOverlaps);                  //overlaps checking
  }
  */

// Bunker

 


  G4double halfXCuboG = 5.*m,  halfYCuboG = (3./2.0 + 0.5)*m, halfZCuboG = 5*m; 
  G4Box* cuboGrande = new G4Box("cuboGrande", halfXCuboG, halfYCuboG, halfZCuboG);
  G4Box* cuboPeq = new G4Box("cuboPeq", halfXCuboG - 50*cm, halfYCuboG - 50*cm, halfZCuboG - 50*cm);

  G4SubtractionSolid* subtraction = new G4SubtractionSolid("subtracted_boxes",cuboGrande,cuboPeq);

  G4VSolid* paredLocalizada = new G4DisplacedSolid("ParedPrimaria",pared,trPared);

  
  G4UnionSolid* unionBunker = new G4UnionSolid("Box+Cylinder", subtraction, paredLocalizada);

  G4LogicalVolume* logicBunker = new G4LogicalVolume(unionBunker,materialPared, "Bunker");
  G4ThreeVector posBunker = G4ThreeVector(0.*cm, -20.*cm, 0.*cm); // se pone el menos para que la cara de con los 8m
  G4Transform3D trBunker = G4Rotate3D() * G4Translate3D(posBunker);
  if (conPared == true){
  
    new G4PVPlacement (trBunker,                      //Traslación 
                    logicBunker,              //its logical volume
                    "Bunker de Concreto",           //its name
                    logicWorld,                      //its mother  volume
                    false,                           //no boolean operation
                    0,                               //copy number
                    checkOverlaps);                  //overlaps checking
  }



  //G4RotationMatrix* rotParedSecundaria = G4RotationMatrix();
  //rotParedSuperior->rotateX(M_PI/2.*rad);
  //G4Transform3D trParedSuperior = G4Translate3D(posFantoma2)*G4Rotate3D(*rotFantoma2); // Vector de transformación


// Detectores
// Fantoma ICRU


//Fantoma 1: Irradiado

  G4double halfXFantom = (40./2.0)*cm; //Ancho del fantoma
  G4double halfYFantom = (20./2.0)*cm; //Profundidad fantoma
  G4double halfZFantom = (70./2.0)*cm; //Altura fantoma

 
  G4EllipticalTube* geomFantoma1 = new G4EllipticalTube("Fantoma1", halfXFantom, halfYFantom, halfZFantom);
 
  G4ThreeVector posFantoma1 = G4ThreeVector(0.*cm,0.*cm, 0.*cm);
  G4RotationMatrix* rotFantoma1 =  new G4RotationMatrix();
  rotFantoma1->rotateY(M_PI/2.*rad);
  //G4RotationMatrix rotFantoma1 =  G4RotationMatrix(); // También funciona
  //rotFantoma1.rotateY(M_PI/2.*rad);



  G4Transform3D trFantoma1 = G4Rotate3D(*rotFantoma1)* G4Translate3D(posFantoma1); // Vector de transformación
  //G4Transform3D trFantoma1* = G4Rotate3D(G4RotationMatrix()) * G4Translate3D (G4ThreeVector()) // Otra forma de obtener una transf 0
  //G4Transform3D trFantoma1* = new G4Transform3D() // Otra forma de obtener una transf 0

  G4Material* materialFantoma1 = tejidoICRU;

  G4LogicalVolume* logicFantoma1 = new G4LogicalVolume(geomFantoma1, materialFantoma1,"Fantoma1"); 
  if (conFantoma1 == true)
  {
    new G4PVPlacement (trFantoma1,                // vector transformación
                    logicFantoma1,              //its logical volume
                    "Fantoma fisico",           //its name
                    logicWorld,                      //its mother  volume
                    false,                           //no boolean operation
                    0,                               //copy number
                    checkOverlaps);                  //overlaps checking
  }

  
  
// Fantoma 2
 
  //G4EllipticalTube* geomFantoma2 = new G4EllipticalTube("Fantoma2", halfXFantom, halfYFantom, halfZFantom);
 
  G4ThreeVector posFantoma2 = G4ThreeVector(0.*m,0.*m, 5.3*m);
  G4RotationMatrix* rotFantoma2 =  new G4RotationMatrix();
  rotFantoma2->rotateX(M_PI/2.*rad);
  //G4RotationMatrix rotFantoma2 =  G4RotationMatrix(); // También funciona
  //rotFantoma2.rotateY(M_PI/2.*rad);


  //En este caso primero roto y después translado (primer caso que me toca pensar en esto)
  G4Transform3D trFantoma2 = G4Translate3D(posFantoma2)*G4Rotate3D(*rotFantoma2); // Vector de transformación

  G4Material* materialFantoma2 = tejidoICRU;

  G4LogicalVolume* logicFantoma2 = new G4LogicalVolume(geomFantoma1, materialFantoma2,"Fantoma2"); 

  new G4PVPlacement (trFantoma2,                // vector transformación
                    logicFantoma2,              //its logical volume
                    "Fantoma fisico 2",           //its name
                    logicWorld,                      //its mother  volume
                    false,                           //no boolean operation
                    0,                               //copy number
                    checkOverlaps);                  //overlaps checking

// Detector ICRU dentro del fantoma 1 

  // Caja de detector/
  G4double halfDetector = (10./2.0)*cm; //Lado 10 cm
  G4Box* cajaDetector = new G4Box("BoxA", halfDetector, halfDetector, halfDetector);
 
  //Posición detector en el centro del fantoma1
  G4ThreeVector posCajaDetector = G4ThreeVector(0.*cm,0.*cm, 0.*cm);

  G4Material* materialDetector = tejidoICRU;

  G4LogicalVolume* logicCajaDetector = new G4LogicalVolume(cajaDetector, materialDetector,"CajaDetector"); 
  
  if (conFantoma1 == true)
  {
    new G4PVPlacement(trFantoma1, logicCajaDetector, "CajaDetector", logicFantoma1, false, 0,  checkOverlaps);  
  }


  // oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

// Bordes, en caso de querer ponerlos
/*
  G4Material* matBordes =  nist->FindOrBuildMaterial("G4_Pb");

  G4double halfZBordes = (4./2.0)*cm; 
  G4double halfGrosorBorde = (1./2.0)*cm;

  G4double halfYBordes = halfDetector + halfGrosorBorde*2.0;
  G4double halfXBordes = halfGrosorBorde; 
  G4Box* borde1 = new G4Box("BordePlomo", halfDetector, halfXBordes, halfZBordes);
  G4Box* borde2 = new G4Box("BordePlomo", halfYBordes, halfXBordes, halfZBordes);
  G4Box* borde3 = new G4Box("BordePlomo", halfDetector, halfXBordes, halfZBordes);
  G4Box* borde4 = new G4Box("BordePlomo", halfYBordes, halfXBordes, halfZBordes);

  G4ThreeVector posBorde1 = G4ThreeVector( 0.*cm, halfDetector + halfGrosorBorde, 0.*cm);  
  G4ThreeVector posBorde2 = G4ThreeVector( - halfDetector - halfGrosorBorde, 0.*cm , 0*cm);  
  G4ThreeVector posBorde3 = G4ThreeVector( 0.*cm, - halfDetector - halfGrosorBorde , 0.*cm);  
  G4ThreeVector posBorde4 = G4ThreeVector( halfDetector + halfGrosorBorde, 0.*cm , 0.*cm);  

  G4RotationMatrix * zRot = new G4RotationMatrix;  // Rotates XandY axes only
  //zRot->rotateZ(M_PI/2.*rad); //rot 90 degrees
  G4Transform3D trBorde1 = G4Translate3D(posBorde1)*trDetector*G4Rotate3D (*zRot); //Se mueve con respecto al detector
  zRot->rotateZ(M_PI/2.*rad); //rot 90 degrees
  G4Transform3D trBorde2 = G4Translate3D(posBorde2)*trDetector*G4Rotate3D (*zRot); //Se mueve con respecto al detector
  zRot->rotateZ(M_PI/2.*rad); //rot 90 degrees
  G4Transform3D trBorde3 = G4Translate3D(posBorde3)*trDetector*G4Rotate3D (*zRot); //Se mueve con respecto al detector
  zRot->rotateZ(M_PI/2.*rad); //rot 90 degrees
  G4Transform3D trBorde4 = G4Translate3D(posBorde4)*trDetector*G4Rotate3D (*zRot); //Se mueve con respecto al detector

  /*
  G4LogicalVolume* logicBorde1 =  new G4LogicalVolume(borde1,matBordes,"Borde1");    // Este nombre sirve para vis.mac
  G4LogicalVolume* logicBorde2 =  new G4LogicalVolume(borde2,matBordes,"Borde1");    // Este nombre sirve para vis.mac
  G4LogicalVolume* logicBorde3 =  new G4LogicalVolume(borde3,matBordes,"Borde1");    // Este nombre sirve para vis.mac
  G4LogicalVolume* logicBorde4 =  new G4LogicalVolume(borde4,matBordes,"Borde1");    // Este nombre sirve para vis.mac
          
  new G4PVPlacement(trBorde1, logicBorde1, "Borde de Plomo", logicWorld, false, 0, checkOverlaps);
  new G4PVPlacement(trBorde2, logicBorde2, "Borde de Plomo", logicWorld, false, 0, checkOverlaps);
  new G4PVPlacement(trBorde3, logicBorde3, "Borde de Plomo", logicWorld, false, 0, checkOverlaps);
  new G4PVPlacement(trBorde4, logicBorde4, "Borde de Plomo", logicWorld, false, 0, checkOverlaps);

  */
  
//Scoring volumes
  fScoringVolume = logicFantoma2;


  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
