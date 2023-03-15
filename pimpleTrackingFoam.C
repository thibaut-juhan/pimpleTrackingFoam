/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
// #include "Random.H"
#include "basicKinematicCollidingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
bool isConverged(Foam::volScalarField maxDiff,double epsilon)
{
   bool res=false;
   if(max(maxDiff).value() < epsilon)
   {
      res = true;
      // If it's true => add particle inside the domain, we have to change the SOI => injectionSteadyState.H	   
   }
   else
   {
      Info << " maximum of temporal derivative : " << max(maxDiff).value() << endl;
      res = false;  
   }
   return res;
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );
    #include "postProcess.H" 	
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    bool convergeCheck = false;
    int  countTimeStep(0);     // timestep counter
    int  secCountTimeStep(1); //  second timeStep counter 
   
    volVectorField U0("U0",U);

    label cellCount = 0;
    List<label>  listCellCount;     // label list of cell
    List<vector> listVectorPosCell; // vector list of cell

    #include "injectionSteadyState.H"

    // find the label of cell which corresponds to the if-statement:
    forAll(mesh.C(),cellI)
    {
	    scalar xC      = mesh.C()[cellI].component(0);
	    scalar yC      = mesh.C()[cellI].component(1);
	    scalar zC      = mesh.C()[cellI].component(2);	  
            scalar R1      = (xC-xC1.value())*(xC-xC1.value()) + (yC-yC1.value())*(yC-yC1.value());
            scalar R2      = (xC-xC2.value())*(xC-xC2.value()) + (yC-yC2.value())*(yC-yC2.value());	
	    cellCount++;
	    // compute the radial coordinate of cell	
	    if( ( R1 < R.value()*R.value() || R2 < R.value()*R.value() ) && (zC < L.value()/10.) )
   	    {
    		listCellCount.append(cellCount); // append the list of cell
	    }
    }
    // fill the vector list => (xC,yC,zC)
    forAll(listCellCount,i)
    {
	scalar xC      = mesh.C()[i].component(0);
	scalar yC      = mesh.C()[i].component(1);
	scalar zC      = mesh.C()[i].component(2);	  	
	listVectorPosCell.append( vector (xC,yC,zC) );

    }

    // create the dictionnary will be filled during the reinjection process : 
    scalar numbProc = UPstream::myProcNo();
    Random randObj; 
    /*
    uMaxCell.set("",listVectorPosCell);
    uMaxCell.regIOobject::write();
    */

    scalar massInside;
    scalar massInsideTemp;

    label nParcelsInsideTemp;
    label nParcelsInside;

    Info<< "\nStarting time loop\n" << endl;
   
    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
	
	countTimeStep++;
        ++runTime;
        
        Info<< "Time = "   << runTime.timeName() << nl << endl;
	Info<< "counter = "<< countTimeStep << endl; 
	
	if (convergeCheck == false)
        {
           // --- Pressure-velocity PIMPLE corrector loop
           while (pimple.loop())
           {   
               if (pimple.firstIter() || moveMeshOuterCorrectors)
               {
                   // Do any mesh changes
                   mesh.controlledUpdate();

                   if (mesh.changing())
                   {
                       MRF.update();

                       if (correctPhi)
                       {
                           // Calculate absolute flux
                           // from the mapped surface velocity
                           phi = mesh.Sf() & Uf();
                           #include "correctPhi.H"
                           // Make the flux relative to the mesh motion
                           fvc::makeRelative(phi, U);
                        }
                        if (checkMeshCourantNo)
                        {
                           #include "meshCourantNo.H"
                        }
                   }
               }

               #include "UEqn.H"

               // --- Pressure corrector loop
               while (pimple.correct())
               {
                  #include "pEqn.H"
               }
               if (pimple.turbCorr())
               {
                   laminarTransport.correct();
                   turbulence->correct();
               }
           }

           // Stored the time every 100 timeStep:
	       if( countTimeStep%100 == 0) 
	       {  
	            forAll(U,cellI) // Loop over cell => compute dC/dt with the backward (schemes) 
	            {
      		        U0[cellI] = U[cellI];  
	            }  
	       }

           // Compare U0 with U 150 dt later		
           if ( countTimeStep == secCountTimeStep*100 + 150 )
           {
                forAll(U,cellI) // loop over cell
                {
                    maxDifference[cellI] = mag ( U[cellI] - U0[cellI] );      
                }

                Info<< "second counter = " << secCountTimeStep << endl;
                convergeCheck = isConverged( maxDifference,1.e-2); 
                // check convergence and inject particles inside the domain in 10 dt later	   
                secCountTimeStep++;	   
           }  		   
        }
        else // steady state condition => add particle inside 
        {

	    Random randObj;
	    // kinematicCloud state before evolving
	    nParcelsInsideTemp     = kinematicCloud.nTotParcels();     // number of parcels inside the domain      
            massInsideTemp         = kinematicCloud.massInSystem(); // mass inside the system before evolve    
	    Info << "parcels inside (before evolving) : " << nParcelsInsideTemp << endl;
	    // kinematicCloud state after evolving
	    kinematicCloud.evolve(); // compute the number of particle outside		
	    // mass + nParcels after evolving :
	    massInside             = kinematicCloud.massInSystem();
	    nParcelsInside         = kinematicCloud.nTotParcels();
            label nParcelsOutside  = mag( nParcelsInsideTemp - nParcelsInside ); // parcel number outside 
	    Info << "parcels inside : " << nParcelsInside << endl;
	    if(nParcelsOutside != 0)
            {
                // loop over parcel outside the domain
                for (int nI ; nI < nParcelsOutside; nI++)
                {
		    if(listCellCount.size() !=0) 
		    {
			    label randList = randObj.position( 0,listVectorPosCell.size() );
	    	    Info << "rand Cell   "<< listCellCount[randList] << endl;		
			    Info << "rand Vector "<< listVectorPosCell[randList] << endl;

			// add particles inside the cloud => addParticle(parceltype *ptr)
		    }

		    // Call constructor of particle => from vector position + cell label
   		    // Need to add the good template => CollidingParcel > KinematicParcel > particle  

	   	    // Info << listVectorPosCell [randPosListVector]  << endl;
	            // Info << listCellCount     [randPosListVector]  << endl;
                }
            }

	    }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
