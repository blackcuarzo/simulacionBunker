#include "B1PrimaryGeneratorAction.hh"

#include "G4ParticleTable.hh"
//#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4Event.hh"

#include "G4GeneralParticleSource.hh" // Necesario para usar el particleSource
//#include "globals.hh"
#include "B1DetectorConstruction.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
//: G4VUserPrimaryGeneratorAction(),
  //particleSource(0) 
{
    // Dan: se crea el objeto GPS
    particleSource = new G4GeneralParticleSource();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete particleSource;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    //Dan: Esto está de acuerdo al Gammaknife
    particleSource->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

