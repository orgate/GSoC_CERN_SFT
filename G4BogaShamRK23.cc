//
// $Id: G4BogaShamRK23.cc 2015-03-27 09:02:32Z gcosmo $
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

#include "G4BogaShamRK23.hh"
#include "G4LineSection.hh"

/////////////////////////////////////////////////////////////////////
//
// Constructor

G4BogaShamRK23::G4BogaShamRK23(G4EquationOfMotion *EqRhs, 
				 G4int noIntegrationVariables, 
				 G4bool primary)
  : G4MagIntegratorStepper(EqRhs, noIntegrationVariables),
    fLastStepLength(0.), fAuxStepper(0)
{
  const G4int numberOfVariables = noIntegrationVariables;

  ak2 = new G4double[numberOfVariables] ;  
  ak3 = new G4double[numberOfVariables] ; 
  ak4 = new G4double[numberOfVariables] ; 
//  ak5 = new G4double[numberOfVariables] ; 	// commented by orgate
  ak5 = 0; 					// added by orgate
//  ak6 = new G4double[numberOfVariables] ; 	// commented by orgate
//  ak7 = 0;					// commented by orgate
  yTemp = new G4double[numberOfVariables] ; 
  yIn = new G4double[numberOfVariables] ;

  fLastInitialVector = new G4double[numberOfVariables] ;
  fLastFinalVector = new G4double[numberOfVariables] ;
  fLastDyDx = new G4double[numberOfVariables];

  fMidVector = new G4double[numberOfVariables];
  fMidError =  new G4double[numberOfVariables];
  if( primary )
  { 
    fAuxStepper = new G4BogaShamRK23(EqRhs, numberOfVariables, !primary);
  }
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4BogaShamRK23::~G4BogaShamRK23()
{
  delete[] ak2;
  delete[] ak3;
  delete[] ak4;
  delete[] ak5;
//  delete[] ak6;				// commented by orgate
  // delete[] ak7;
  delete[] yTemp;
  delete[] yIn;

  delete[] fLastInitialVector;
  delete[] fLastFinalVector;
  delete[] fLastDyDx;
  delete[] fMidVector;
  delete[] fMidError; 

  delete fAuxStepper;
}

//////////////////////////////////////////////////////////////////////
//
// Given values for n = 6 variables yIn[0,...,n-1] 
// known  at x, use the third-order Bogacki-Shampine-2-3 method
// to advance the solution over an interval
// Step and return the incremented variables as yOut[0,...,n-1]. Also
// return an estimate of the local truncation error yErr[] using the
// embedded 2nd-order method. The user supplies routine
// RightHandSide(y,dydx), which returns derivatives dydx for y .

void
G4BogaShamRK23::Stepper(const G4double yInput[],
                         const G4double dydx[],
                               G4double Step,
                               G4double yOut[],
                               G4double yErr[])
{
  // const G4int nvar = 4 ;
  // const G4double a2 = 0.5 , a3 = 0.75 , a4 = 1;
 G4int i;

 const G4double  b21 = 0.5 ,
                 b31 = 0.0 , b32 = 0.75 ,
                 b41 = 2.0/9.0 , b42 = 1.0/3.0 , b43 = 4.0/9.0 ,

//                 b51 = -11.0/54.0 , b52 = 2.5 , b53 = -70.0/27.0 , // commented by orgate
//                 b54 = 35.0/27.0 ,				     // commented by orgate

//                 b61 = 1631.0/55296.0 , b62 =   175.0/512.0 ,	     // commented by orgate
//                 b63 =  575.0/13824.0 , b64 = 44275.0/110592.0 ,   // commented by orgate
//                 b65 =  253.0/4096.0 ,

                 c1 = 2.0/9.0 , c2 = 1.0/3.0 , c3 = 4.0/9.0 ,
                 c4 = 0.0 ;
//                                          dc5 = -277.0/14336.0 ;   // commented by orgate

 const G4double dc1 = c1 - 7.0/24.0 ,  dc2 = c2 - 1.0/4.0 ,
    dc3 = c3 - 1.0/3.0 , dc4 = c4 - 1.0/8.0 ;

 // Initialise time to t0, needed when it is not updated by the integration.
 //        [ Note: Only for time dependent fields (usually electric) 
 //                  is it neccessary to integrate the time.] 
 yOut[7] = yTemp[7]   = yIn[7]; 

 const G4int numberOfVariables= this->GetNumberOfVariables(); 
 // The number of variables to be integrated over

   //  Saving yInput because yInput and yOut can be aliases for same array

   for(i=0;i<numberOfVariables;i++) 
   {
     yIn[i]=yInput[i];
   }
 // RightHandSide(yIn, dydx) ;              // 1st Step

 for(i=0;i<numberOfVariables;i++) 
 {
   yTemp[i] = yIn[i] + b21*Step*dydx[i] ;
 }
 RightHandSide(yTemp, ak2) ;              // 2nd Step

 for(i=0;i<numberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b31*dydx[i] + b32*ak2[i]) ;
 }
 RightHandSide(yTemp, ak3) ;              // 3rd Step

 for(i=0;i<numberOfVariables;i++)
 {
    yOut[i] = yIn[i] + Step*(c1*dydx[i] + c2*ak2[i] + c3*ak3[i]) ;
 }										//O(h^3) output step

 RightHandSide(yOut, ak4) ;              // 4th Step	// This ak4 should be stored to separately to exploit FSAL

// commented by orgate
/*
 for(i=0;i<numberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] +
                      b54*ak4[i]) ;
 }
 RightHandSide(yTemp, ak5) ;              // 5th Step

 for(i=0;i<numberOfVariables;i++)
 {
    yTemp[i] = yIn[i] + Step*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] +
                      b64*ak4[i] + b65*ak5[i]) ;
 }
 RightHandSide(yTemp, ak6) ;              // 6th Step
*/

 for(i=0;i<numberOfVariables;i++)
 {
    // Accumulate increments with proper weights

//    yOut[i] = yIn[i] + Step*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]) ; // commented by orgate

    // Estimate error as difference between 4th and
    // 5th order methods

    yErr[i] = Step*(dc1*dydx[i] + dc2*ak2[i] + dc3*ak3[i] + dc4*ak4[i]) ;

    // Store Input and Final values, for possible use in calculating chord
    fLastInitialVector[i] = yIn[i] ;
    fLastFinalVector[i]   = yOut[i];
    fLastDyDx[i]          = dydx[i];
 }
 // NormaliseTangentVector( yOut ); // Not wanted

 fLastStepLength =Step;

 return ;
} 

///////////////////////////////////////////////////////////////////////////////

void
G4BogaShamRK23::StepWithEst( const G4double*,
                              const G4double*,
                                    G4double,
                                    G4double*,
                                    G4double&,
                                    G4double&,
                              const G4double*,
                                    G4double*  )    
{
  G4Exception("G4BogaShamRK23::StepWithEst()", "GeomField0001",
              FatalException, "Method no longer used.");
  return ;
}

/////////////////////////////////////////////////////////////////

G4double  G4BogaShamRK23::DistChord() const
{
  G4double distLine, distChord; 
  G4ThreeVector initialPoint, finalPoint, midPoint;

  // Store last initial and final points (they will be overwritten in self-Stepper call!)
  initialPoint = G4ThreeVector( fLastInitialVector[0], 
                                fLastInitialVector[1], fLastInitialVector[2]); 
  finalPoint   = G4ThreeVector( fLastFinalVector[0],  
                                fLastFinalVector[1],  fLastFinalVector[2]); 

  // Do half a step using StepNoErr

  fAuxStepper->Stepper( fLastInitialVector, fLastDyDx, 0.5 * fLastStepLength, 
           fMidVector,   fMidError );

  midPoint = G4ThreeVector( fMidVector[0], fMidVector[1], fMidVector[2]);       

  // Use stored values of Initial and Endpoint + new Midpoint to evaluate
  //  distance of Chord


  if (initialPoint != finalPoint) 
  {
     distLine  = G4LineSection::Distline( midPoint, initialPoint, finalPoint );
     distChord = distLine;
  }
  else
  {
     distChord = (midPoint-initialPoint).mag();
  }
  return distChord;
}


# GSoC_CERN_SFT
