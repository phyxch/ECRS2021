/// Date created: Jan 28, 2021
/// Authors: hexc, Ernesto, Arfa, Jarred Lockie
///   The first implementation of ActionInitialization for the ECRS simulation

#include "ECRSActionInitialization.hh"
#include "ECRSRunAction.hh"
#include "ECRSEventAction.hh"
#include "ECRSSteppingAction.hh"
#include "ECRSTrackingAction.hh"
#include "ECRSPrimaryGeneratorAction.hh"
#include "ECRSStackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ECRSActionInitialization::ECRSActionInitialization(ECRSDetectorConstruction *detectorIn)
 : G4VUserActionInitialization()
{
  detector = detectorIn;    // set the detector object
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ECRSActionInitialization::~ECRSActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ECRSActionInitialization::BuildForMaster() const
{
  SetUserAction(new ECRSRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ECRSActionInitialization::Build() const
{
  ECRSRunAction* runAction = new ECRSRunAction();
  SetUserAction(runAction);

  ECRSEventAction* evtAction = new ECRSEventAction();
  SetUserAction(evtAction);

  ECRSPrimaryGeneratorAction* primaryGenAction = new ECRSPrimaryGeneratorAction(detector, runAction, evtAction);
  SetUserAction(primaryGenAction);
  
  SetUserAction(new ECRSStackingAction);
  SetUserAction(new ECRSSteppingAction(detector, runAction, evtAction));
  SetUserAction(new ECRSTrackingAction(runAction, evtAction));
  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
