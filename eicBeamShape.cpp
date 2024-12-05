// Modification of Pythia Class BeamShape for EIC Beam Conditions

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"
#include "eicBeamShape.hh"

//#include "TMath.h"

using namespace Pythia8;

eicBeamShape::eicBeamShape(int config, double ion, double lepton, double xAngle) 
  {
    mDivAcc = config;
    mIonBeamEnergy = ion;
    mLeptonBeamEnergy = lepton;
    mXAngle = xAngle;
    mXAngleY = 0.000100; // Fixed for IP6
    mKill == 0;

    // Ensure Beam Energies Correspond to Those Presented in CDR
    if(ion != 275.0 && ion != 110.0 && ion != 100.0 && ion != 41.0) 
      {
	cout << ion << " is not a valid Hadron Beam Energy!!" << endl;
	cout << "Valid Energies are 275.0, 110.0, 100.0, and 41.0" << endl;
	cout << "Turning off all beam effects" << endl;
	mKill = 1;
      }
    if(lepton != 18.0 && lepton != 10.0 && lepton != 5.0) 
      {
	cout << lepton << " is not a valid Lepton Beam Energy!!" << endl;
	cout << "Valid Energies are 18.0, 10.0, and 5.0" << endl;
	cout << "Turning off all beam effects" << endl;
	mKill = 1;
      }

    // Ensure Beam Energy Combinations Correspond to Those Presented in CDR
    if(ion == 275.0 && (lepton != 18.0 && lepton != 10.0))
      {
	cout << lepton << "x" << ion << " is not a valid energy combination!!" << endl;
	cout << "Valid (ep) Combinations are 18x275, 10x275, 10x100, 5x100, and 5x41" << endl;
	cout << "Valid (eA) Combinations are 18x110, 10x110, 5x110, and 5x41" << endl;
	cout << "Turning off all beam effects" << endl;
	mKill = 1;
      }
    if(ion == 110.0 && (lepton != 18.0 && lepton != 10.0 && lepton != 5.0))
      {
	cout << lepton << "x" << ion << " is not a valid energy combination!!" << endl;
	cout << "Valid (ep) Combinations are 18x275, 10x275, 10x100, 5x100, and 5x41" << endl;
	cout << "Valid (eA) Combinations are 18x110, 10x110, 5x110, and 5x41" << endl;
	cout << "Turning off all beam effects" << endl;
	mKill = 1;
      }
    if(ion == 100.0 && (lepton != 10.0 && lepton != 5.0))
      {
	cout << lepton << "x" << ion << " is not a valid energy combination!!" << endl;
	cout << "Valid (ep) Combinations are 18x275, 10x275, 10x100, 5x100, and 5x41" << endl;
	cout << "Valid (eA) Combinations are 18x110, 10x110, 5x110, and 5x41" << endl;
	cout << "Turning off all beam effects" << endl;
	mKill = 1;
      }
    if(ion == 41.0 && lepton != 5.0)
      {
	cout << lepton << "x" << ion << " is not a valid energy combination!!" << endl;
	cout << "Valid (ep) Combinations are 18x275, 10x275, 10x100, 5x100, and 5x41" << endl;
	cout << "Valid (eA) Combinations are 18x110, 10x110, 5x110, and 5x41" << endl;
	cout << "Turning off all beam effects" << endl;
	mKill = 1;
      }
  }

void eicBeamShape::pick() {

  //=============================================
  // (Re)set Values
  //=============================================

  deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB
    = vertexX = vertexY = vertexZ = vertexT = 0.;


  //=============================================
  // Set Collision Vertex and Time
  //=============================================

  // Collision vertex and time are determined by randomly sampling the z-positions of the interacting lepton and hadron within their respective bunches. Assume z=0 is middle of bunch, hadron bunch moves left to right and lepton bunch moves right to left, positive z within hadron bunch is in direction of travel and positive z within lepton bunch is opposite direction of travel

  //double xIPOffset = 0.;
  double hadronPartPos = 0.;
  double leptonPartPos = 0.;

  if(allowVertexSpread)
    {
      // RMS Bunch Length [mm] (CDR Table 3.3 & 3.4)
      double hadronBL = 0.;
      double leptonBL = 0.;
      
      // Set different values depending on energy
      if(mIonBeamEnergy == 275.0) hadronBL = 60.0;
      if(mIonBeamEnergy == 110.0) hadronBL = 70.0;
      if(mIonBeamEnergy == 100.0) hadronBL = 70.0;
      if(mIonBeamEnergy == 41.0 && mDivAcc != 3) hadronBL = 75.0;
      if(mIonBeamEnergy == 41.0 && mDivAcc == 3) hadronBL = 116.0; // for eA
      if(mLeptonBeamEnergy == 18.0) leptonBL = 9.0;
      if(mLeptonBeamEnergy == 10.0) leptonBL = 7.0;
      if(mLeptonBeamEnergy == 5.0)  leptonBL = 7.0;
      
      // Set particle positions
      hadronPartPos = hadronBL*rndmPtr->gauss();
      leptonPartPos = leptonBL*rndmPtr->gauss();

      // Find Collision time, z position of collision, z position of center of hadron bunch at collision time and x position of collision
      // Quantities are found via a system of parametric equations
      // Z position of colliding hadron as a function of time: Z_h = Cos(0.5*theta_c)*t + Z_h_offset (Z_h_offset = distance from center of bunch)
      // Z position of colliding lepton as a function of time: Z_l = -Cos(0.5*theta_c)*t + Z_l_offset
      // Collision time when Z_h = Z_l -> t_int = (Z_l_offset - Z_h_offset)/(2*Cos(0.5*theta_c))
      // Z position of collision: Z_int = (Z_l_offset + Z_h_offset)/2
      // Z position of center of hadron bunch a collision time - this gives x position of collision via relation x = z*Tan(0.5*theta_c)
      // Z_bunch_int = Cos(0.5*theta_c)*t_int
      // X_int = Z_bunch_int * Tan(0.5*theta_c)

      double c_c = std::cos(mXAngle/2.0);
      double s_c = std::sin(mXAngle/2.0);
      double t_c = std::tan(mXAngle/2.0);

      double t_int = (leptonPartPos - hadronPartPos)/(2.0*c_c);
      double z_int = (leptonPartPos + hadronPartPos)/2.0;
      double z_bunch_int = c_c*t_int;
      double x_int = z_bunch_int*t_c;

      // x_int is the x position of collision assuming the colliding particles are in the center of the bunch. Sample random x position according to x-width of bunches. x_int is then an offset to this. Get y position as well.

      double y_int = 0.;
      if(sigmaVertexX > 0.)
	{
	  x_int += sigmaVertexX * rndmPtr->gauss();
	}

      if(sigmaVertexY > 0.)
	{
	  y_int += sigmaVertexY * rndmPtr->gauss();
	}

      // We now have the x-y-z position of the collision in the accelerator frame, but we want it in the detector frame. Rotate by 0.5*theta_c to get to accelerator frame

      double tmpVtxX, tmpVtxY, tmpVtxZ;
      tmpVtxX = tmpVtxY = tmpVtxZ = 0.;

      RotY(mXAngle/2.0,x_int,y_int,z_int,&tmpVtxX,&tmpVtxY,&tmpVtxZ);

      vertexT = t_int;
      vertexX = tmpVtxX;
      vertexY = tmpVtxY;
      vertexZ = tmpVtxZ;
    }


  //=============================================
  // Smear Beam Energy
  //=============================================

  // Change in Beam Pz
  double tmpPzA, tmpPzB;
  tmpPzA = tmpPzB = 0.;

  if(allowMomentumSpread)
    {
      double gaussZA, gaussZB;
      gaussZA = gaussZB = 0.;
      if(sigmaPzA > 0.)
	{
	  gaussZA = rndmPtr->gauss();
	  tmpPzA = mIonBeamEnergy * (sigmaPzA * gaussZA);
	}
      if(sigmaPzB > 0.)
	{
	  gaussZB = rndmPtr->gauss();
	  tmpPzB = mLeptonBeamEnergy * (sigmaPzB * gaussZB);
	}
    }


  //=============================================
  // Crossing Angle Effects
  //=============================================

  if(allowMomentumSpread)
    {
      // Modify Px Due to Crossing Angle
      deltaPxA += (mIonBeamEnergy + tmpPzA)*std::sin(mXAngle);

      // Modify Py Due to Crossing Angle
      deltaPyA += (mIonBeamEnergy + tmpPzA)*std::sin(mXAngleY);

      // Modify Pz Due to Crossing Angle
      deltaPzA += (mIonBeamEnergy + tmpPzA)*std::cos(mXAngle) - mIonBeamEnergy; // Ignore Pz modification due to crossing angle in Y
      //deltaPzA += tmpPzA;

      // Electron Beam Has 0 Crossing Angle
      deltaPzB += tmpPzB;
    }


  //=============================================
  // Divergence Effects
  //=============================================

  if(allowMomentumSpread) // allowMomentumSpread
    {
      double gaussXA, gaussYA, gaussXB, gaussYB;
      double pxLocal = 0.;
      double pyLocal = 0.;

      if(sigmaPxA > 0.)
	{
	  // Ion Beam X Divergence
	  gaussXA = rndmPtr->gauss();
	  double div = sigmaPxA * gaussXA;
	  pxLocal = (mIonBeamEnergy + tmpPzA)*std::sin(div); // Dispersion in Beam Frame
	}
      if(sigmaPyA > 0.)
	{
	  // Ion Beam Y Divergence
	  gaussYA = rndmPtr->gauss();
	  double div = sigmaPyA * gaussYA;
	  pyLocal = (mIonBeamEnergy + tmpPzA)*std::sin(div);
	}

      // Rotate into Detector Frame
      double divPxA, divPyA, divPzA;
      divPxA = divPyA = divPzA = 0.;

      RotXY(mXAngle,-1.0*mXAngleY,pxLocal,pyLocal,0.,&divPxA,&divPyA,&divPzA);

      deltaPxA += divPxA;
      deltaPyA += divPyA;
      deltaPzA += divPzA;

      if(sigmaPxB > 0.)
	{
	  // Lepton Beam X Divergence
	  gaussXB = rndmPtr->gauss();
	  double div = sigmaPxB * gaussXB;
	  deltaPxB += (mLeptonBeamEnergy + tmpPzB)*std::sin(div); // Electron Beam along Z Axis
	}

      if(sigmaPyB > 0.)
	{
	  // Lepton Beam Y Divergence
	  gaussYB = rndmPtr->gauss();
	  double div = sigmaPyB * gaussYB;
	  deltaPyB += (mLeptonBeamEnergy + tmpPzB)*std::sin(div);
	}
    }


  //=============================================
  // Crabbing Momentum Kicks
  //=============================================

  if(allowMomentumSpread)
    {
      // Assign Crab and IP Beta Functions for Different Beam Energies [mm] (CDR Tab 3.3 & 3.4 and Elke)
      double betaCrabHad = 0.;
      double betaStarHad = 0.; // In horizontal direction

      if(mIonBeamEnergy == 275.0) betaCrabHad = 1300000.0;
      if(mIonBeamEnergy == 110.0) betaCrabHad = 500000.0; // Estimate
      if(mIonBeamEnergy == 100.0) betaCrabHad = 500000.0;
      if(mIonBeamEnergy == 41.0)  betaCrabHad = 200000.0;
      if(mDivAcc == 1) // High Divergence Config - CDR Table 3.3
	{
	  if(mIonBeamEnergy == 275.0) betaStarHad = 800.0;
	  if(mIonBeamEnergy == 100.0 && mLeptonBeamEnergy == 10.0) betaStarHad = 630.0; // For root[s] = 63.2
	  if(mIonBeamEnergy == 100.0 && mLeptonBeamEnergy == 5.0)  betaStarHad = 610.0; // For root[s] = 44.7
	  if(mIonBeamEnergy == 41.0)  betaStarHad = 900.0;
	}
      if(mDivAcc == 2) // High Acceptance Config - CDR Table 3.4
	{
	  if(mIonBeamEnergy == 275.0 && mLeptonBeamEnergy == 18.0) betaStarHad = 4170.0; // For root[s] = 140.7
	  if(mIonBeamEnergy == 275.0 && mLeptonBeamEnergy == 10.0) betaStarHad = 2650.0; // For root[s] = 104.9
	  if(mIonBeamEnergy == 100.0 && mLeptonBeamEnergy == 10.0) betaStarHad = 940.0; // For root[s] = 63.2
	  if(mIonBeamEnergy == 100.0 && mLeptonBeamEnergy == 5.0)  betaStarHad = 800.0; // For root[s] = 44.7
	  if(mIonBeamEnergy == 41.0 && mLeptonBeamEnergy == 5.0)   betaStarHad = 900.0; // For root[s] = 28.6
	}
      if(mDivAcc == 3) // eA Config - CDR Table 3.5
	{
	  if(mIonBeamEnergy == 110.0) betaStarHad = 910.0; 
	  if(mIonBeamEnergy == 41.0)  betaStarHad = 900.0;
	}

      double betaCrabLep = 0.;
      double betaStarLep = 0.;
      if(mLeptonBeamEnergy == 18.0) betaCrabLep = 150000.0;
      if(mLeptonBeamEnergy == 10.0) betaCrabLep = 150000.0;
      if(mLeptonBeamEnergy == 5.0)  betaCrabLep = 150000.0;
      if(mDivAcc == 1) // High Divergence Config - CDR Table 3.3
	{
	  if(mLeptonBeamEnergy == 18.0) betaStarLep = 590.0;
	  if(mLeptonBeamEnergy == 10.0 && mIonBeamEnergy == 275.0) betaStarLep = 450.0; // For root[s] = 104.9
	  if(mLeptonBeamEnergy == 10.0 && mIonBeamEnergy == 100.0) betaStarLep = 960.0; // For root[s] = 63.2
	  if(mLeptonBeamEnergy == 5.0 && mIonBeamEnergy == 100.0)  betaStarLep = 780.0; // For root[s] = 44.7
	  if(mLeptonBeamEnergy == 5.0 && mIonBeamEnergy == 41.0)   betaStarLep = 1960.0; // For root[s] = 28.6
	}
      if(mDivAcc == 2) // High Acceptance Config - CDR Table 3.4
	{
	  if(mLeptonBeamEnergy == 18.0) betaStarLep = 3060.0;
	  if(mLeptonBeamEnergy == 10.0 && mIonBeamEnergy == 275.0) betaStarLep = 1490.0; // For root[s] = 104.9
	  if(mLeptonBeamEnergy == 10.0 && mIonBeamEnergy == 100.0) betaStarLep = 1430.0; // For root[s] = 63.2
	  if(mLeptonBeamEnergy == 5.0 && mIonBeamEnergy == 100.0)  betaStarLep = 1030.0; // For root[s] = 44.7
	  if(mLeptonBeamEnergy == 5.0 && mIonBeamEnergy == 41.0)   betaStarLep = 1960.0; // For root[s] = 28.6 
	}
      if(mDivAcc == 3)
	{
	  if(mLeptonBeamEnergy == 18.0) betaStarLep = 1960.0;
	  if(mLeptonBeamEnergy == 10.0) betaStarLep = 1930.0;
	  if(mLeptonBeamEnergy == 5.0 && mIonBeamEnergy == 110.0) betaStarLep = 1930.0;
	  if(mLeptonBeamEnergy == 5.0 && mIonBeamEnergy == 41.0)  betaStarLep = 3070.0;
	}

      //double betaCrabHad = 1300000.0; // Beta Crab in mm for 275 hadron beam
      //double betaStarHad = 800.0; // Beta Star in mm for 275 hadron beam

      //double betaCrabLep = 150000.0; // Beta Crab in mm for 18 lepton beam
      //double betaStarLep = 590.0; // Beta Star in mm for 18 lepton beam

      // Calculate angular deflection
      double crabAngHad = ((mXAngle/2.0)*hadronPartPos)/(std::sqrt(betaCrabHad*betaStarHad));
      double crabAngLep = ((mXAngle/2.0)*leptonPartPos)/(std::sqrt(betaCrabLep*betaStarLep));

      // Calculate Magnitude of Px kick
      double crabKickHad = (mIonBeamEnergy + tmpPzA)*std::sin(crabAngHad);
      double crabKickLep = (-1.0)*(mLeptonBeamEnergy + tmpPzB)*std::sin(crabAngLep); // 

      // Rotate Momentum Kick into Detector Frame
      double tmpVertPxA, tmpVertPyA, tmpVertPzA;
      double tmpVertPxB, tmpVertPyB, tmpVertPzB;
      tmpVertPxA = tmpVertPyA = tmpVertPzA = tmpVertPxB = tmpVertPyB = tmpVertPzB = 0.;

      RotY(mXAngle/2.0,crabKickHad,0.,0.,&tmpVertPxA,&tmpVertPyA,&tmpVertPzA);
      RotY(mXAngle/2.0,crabKickLep,0.,0.,&tmpVertPxB,&tmpVertPyB,&tmpVertPzB);

      deltaPxA += tmpVertPxA;
      deltaPzA += tmpVertPzA;
      deltaPxB += tmpVertPxB;
      deltaPzB += tmpVertPzB;
    }


  //=============================================
  // Reset for Input to Steering File Mismatch
  //=============================================
  int localKill = 0;
  if(mIonBeamEnergy == 275.0 && mLeptonBeamEnergy == 18.0)
    {
      if(mDivAcc == 1)
	{
	  if(sigmaPxA != 0.000150 || sigmaPyA != 0.000150 || sigmaPzA != 0.00068) localKill = 1;
	  if(sigmaPxB != 0.000202 || sigmaPyB != 0.000187 || sigmaPzB != 0.00109) localKill = 1;
	  if(sigmaVertexX != 0.084 || sigmaVertexY != 0.008) localKill = 1;
	}
      if(mDivAcc == 2)
	{
	  if(sigmaPxA != 0.000065 || sigmaPyA != 0.000065 || sigmaPzA != 0.00068) localKill = 1;
	  if(sigmaPxB != 0.000089 || sigmaPyB != 0.000082 || sigmaPzB != 0.00109) localKill = 1;
	  if(sigmaVertexX != 0.192 || sigmaVertexY != 0.017) localKill = 1;
	}
      if(mDivAcc == 3) localKill = 1;
    }
  if(mIonBeamEnergy == 275.0 && mLeptonBeamEnergy == 10.0)
    {
      if(mDivAcc == 1)
	{
	  if(sigmaPxA != 0.000119 || sigmaPyA != 0.000119 || sigmaPzA != 0.00068) localKill = 1;
	  if(sigmaPxB != 0.000211 || sigmaPyB != 0.000152 || sigmaPzB != 0.00058) localKill = 1;
	  if(sigmaVertexX != 0.067 || sigmaVertexY != 0.006) localKill = 1;
	}
      if(mDivAcc == 2)
	{
	  if(sigmaPxA != 0.000065 || sigmaPyA != 0.000065 || sigmaPzA != 0.00068) localKill = 1;
	  if(sigmaPxB != 0.000116 || sigmaPyB != 0.000084 || sigmaPzB != 0.00058) localKill = 1;
	  if(sigmaVertexX != 0.122 || sigmaVertexY != 0.011) localKill = 1;
	}
      if(mDivAcc == 3) localKill = 1;
    }
  if(mIonBeamEnergy == 110.0 && mLeptonBeamEnergy == 18.0)
    {
      if(mDivAcc == 1) localKill = 1;
      if(mDivAcc == 2) localKill = 1;
      if(mDivAcc == 3) 
	{
	  if(sigmaPxA != 0.000218 || sigmaPyA != 0.000379 || sigmaPzA != 0.00062) localKill = 1;
	  if(sigmaPxB != 0.000101 || sigmaPyB != 0.000037 || sigmaPzB != 0.00109) localKill = 1;
	  if(sigmaVertexX != 0.140 || sigmaVertexY != 0.011) localKill = 1;
	}
    }
  if(mIonBeamEnergy == 110.0 && mLeptonBeamEnergy == 10.0)
    {
      if(mDivAcc == 1) localKill = 1;
      if(mDivAcc == 2) localKill = 1;
      if(mDivAcc == 3) 
	{
	  if(sigmaPxA != 0.000216 || sigmaPyA != 0.000274 || sigmaPzA != 0.00062) localKill = 1;
	  if(sigmaPxB != 0.000102 || sigmaPyB != 0.000092 || sigmaPzB != 0.00058) localKill = 1;
	  if(sigmaVertexX != 0.139 || sigmaVertexY != 0.008) localKill = 1;
	}
    }
  if(mIonBeamEnergy == 110.0 && mLeptonBeamEnergy == 5.0)
    {
      if(mDivAcc == 1) localKill = 1;
      if(mDivAcc == 2) localKill = 1;
      if(mDivAcc == 3) 
	{
	  if(sigmaPxA != 0.000215 || sigmaPyA != 0.000275 || sigmaPzA != 0.00062) localKill = 1;
	  if(sigmaPxB != 0.000102 || sigmaPyB != 0.000185 || sigmaPzB != 0.00068) localKill = 1;
	  if(sigmaVertexX != 0.139 || sigmaVertexY != 0.008) localKill = 1;
	}
    }
  if(mIonBeamEnergy == 41.0 && mLeptonBeamEnergy == 5.0 && mDivAcc == 3)
    {
      if(mDivAcc == 1) localKill = 1;
      if(mDivAcc == 2) localKill = 1;
      if(mDivAcc == 3) 
	{
	  if(sigmaPxA != 0.000275 || sigmaPyA != 0.000377 || sigmaPzA != 0.00100) localKill = 1;
	  if(sigmaPxB != 0.000081 || sigmaPyB != 0.000136 || sigmaPzB != 0.00068) localKill = 1;
	  if(sigmaVertexX != 0.175 || sigmaVertexY != 0.011) localKill = 1;
	}
    }
  if(mIonBeamEnergy == 100.0 && mLeptonBeamEnergy == 10.0)
    {
      if(mDivAcc == 1)
	{
	  if(sigmaPxA != 0.000220 || sigmaPyA != 0.000220 || sigmaPzA != 0.00097) localKill = 1;
	  if(sigmaPxB != 0.000145 || sigmaPyB != 0.000105 || sigmaPzB != 0.00058) localKill = 1;
	  if(sigmaVertexX != 0.098 || sigmaVertexY != 0.008) localKill = 1;
	}
      if(mDivAcc == 2)
	{
	  if(sigmaPxA != 0.000180 || sigmaPyA != 0.000180 || sigmaPzA != 0.00097) localKill = 1;
	  if(sigmaPxB != 0.000118 || sigmaPyB != 0.000086 || sigmaPzB != 0.00058) localKill = 1;
	  if(sigmaVertexX != 0.120 || sigmaVertexY != 0.011) localKill = 1;
	}
      if(mDivAcc == 3) localKill = 1;
    }
  if(mIonBeamEnergy == 100.0 && mLeptonBeamEnergy == 5.0)
    {
      if(mDivAcc == 1)
	{
	  if(sigmaPxA != 0.000206 || sigmaPyA != 0.000206 || sigmaPzA != 0.00097) localKill = 1;
	  if(sigmaPxB != 0.000160 || sigmaPyB != 0.000160 || sigmaPzB != 0.00068) localKill = 1;
	  if(sigmaVertexX != 0.088 || sigmaVertexY != 0.008) localKill = 1;
	}
      if(mDivAcc == 2)
	{
	  if(sigmaPxA != 0.000180 || sigmaPyA != 0.000180 || sigmaPzA != 0.00097) localKill = 1;
	  if(sigmaPxB != 0.000140 || sigmaPyB != 0.000140 || sigmaPzB != 0.00068) localKill = 1;
	  if(sigmaVertexX != 0.101 || sigmaVertexY != 0.009) localKill = 1;
	}
      if(mDivAcc == 3) localKill = 1;
    }
  if(mIonBeamEnergy == 41.0 && mLeptonBeamEnergy == 5.0 && mDivAcc != 3)
    {
      if(mDivAcc == 1)
	{
	  if(sigmaPxA != 0.000220 || sigmaPyA != 0.000380 || sigmaPzA != 0.00103) localKill = 1;
	  if(sigmaPxB != 0.000101 || sigmaPyB != 0.000129 || sigmaPzB != 0.00068) localKill = 1;
	  if(sigmaVertexX != 0.140 || sigmaVertexY != 0.019) localKill = 1;
	}
      if(mDivAcc == 2)
	{
	  if(sigmaPxA != 0.000220 || sigmaPyA != 0.000380 || sigmaPzA != 0.00103) localKill = 1;
	  if(sigmaPxB != 0.000101 || sigmaPyB != 0.000129 || sigmaPzB != 0.00068) localKill = 1;
	  if(sigmaVertexX != 0.140 || sigmaVertexY != 0.019) localKill = 1;
	}
      if(mDivAcc == 3) localKill = 1;
    }

  if(localKill == 1)
    {
      cout << "Steering File and Input Energy Mismatch" << endl;
      cout << "Turning off all beam effects" << endl;
    }


  //=============================================
  // Reset for Bad Energy Combinations
  //=============================================
  if(mKill == 1 || localKill == 1)
    {
      deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB = vertexX = vertexY = vertexZ = vertexT = 0.;
    }
}


void eicBeamShape::RotY(double theta, double xin, double yin, double zin, double *xout, double *yout, double *zout) {

  *xout = xin*std::cos(theta) + zin*std::sin(theta);
  *yout = yin;
  *zout = zin*std::cos(theta) - xin*std::sin(theta);
}


void eicBeamShape::RotXY(double theta, double phi, double xin, double yin, double zin, double *xout, double *yout, double *zout) {

  *xout = xin*std::cos(theta) + zin*std::sin(theta);
  *yout = xin*std::sin(phi)*std::sin(theta) + yin*std::cos(phi) - zin*std::sin(phi)*std::cos(theta);
  *zout = yin*std::sin(phi) - xin*std::cos(phi)*std::sin(theta) + zin*std::cos(phi)*std::cos(theta);
}
