/// Date created: Jan 28, 2021
/// Authors: hexc, Ernesto, Arfa, Jarred Lockie
///   The first implementation of ActionInitialization for the ECRS simulation

#ifndef ECRSActionInitialization_h
#define ECRSActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "ECRSDetectorConstruction.hh"

/// Action initialization class.
///

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ECRSActionInitialization : public G4VUserActionInitialization
{
public:
  ECRSActionInitialization(ECRSDetectorConstruction*);
  virtual ~ECRSActionInitialization();
  
  virtual void BuildForMaster() const;
  virtual void Build() const;

private:
  ECRSDetectorConstruction *detector;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
