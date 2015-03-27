//
// $Id: G4BogaShamRK23.hh 2015-03-27 09:02:32Z gcosmo $
//
// The Bogacki-Shampine 2/3 method is an embedded second
// order method (giving third-order accuracy) for the solution of an ODE.
// Two different second order estimates are calculated; their difference
// gives an error estimate. [ref. Numerical Recipes in C, 2nd Edition]
// It is used to integrate the equations of the motion of a particle 
// in a magnetic field.
//
//  [ref. Numerical Recipes in C, 2nd Edition]
//
//	Author: Alfred Ajay Aureate .R
//
// -------------------------------------------------------------------

#ifndef G4BogaSham_RK23
#define G4BogaSham_RK23

#include "G4MagIntegratorStepper.hh"

class G4BogaShamRK23 : public G4MagIntegratorStepper
{

  public:  // with description

    G4BogaShamRK23( G4EquationOfMotion *EqRhs,
                     G4int numberOfVariables = 6,
                     G4bool primary= true ) ;
   ~G4BogaShamRK23() ;

    void Stepper( const G4double y[],
                  const G4double dydx[],
                        G4double h,
                        G4double yout[],
                        G4double yerr[] ) ;

  public:  // without description

    G4double  DistChord()   const; 
    G4int IntegratorOrder() const { return 2; }

  private:

    void StepWithEst( const G4double yIn[],
                      const G4double dydx[],
                            G4double Step,
                            G4double yOut[],
                            G4double& alpha2,
                            G4double& beta2,
                      const G4double B1[],
                            G4double B2[]    );  
      // No longer used. Obsolete.

    G4BogaShamRK23(const G4BogaShamRK23&);
    G4BogaShamRK23& operator=(const G4BogaShamRK23&);
      // Private copy constructor and assignment operator.

  private:

    G4double *ak2, *ak3, *ak4, *ak5, *yTemp, *yIn;		// modified by orgate
      // scratch space

    G4double fLastStepLength;
    G4double *fLastInitialVector, *fLastFinalVector,
             *fLastDyDx, *fMidVector, *fMidError;
      // for DistChord calculations

    G4BogaShamRK23* fAuxStepper; 

};

#endif /* G4BogaSham_RK23 */
# GSoC_CERN_SFT
