/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Ibai Mugica
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    deCoupledAmrFoam

Description
     Solver for the thermodynamic cycle of an Active Magnetic Regenerator. 
     Fluid velocity and magnetic fields are precalculated and read from a file.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "regionProperties.H"
#include "fixedGradientFvPatchFields.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "pimpleMultiRegionControl.H"
#include "magneticProp.H"
#include "therCycleControl.H"
#include <iostream>
#include <iomanip>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H

    #include "postProcess.H"
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    //Create pimple controls for each region
    pimpleMultiRegionControl pimples(fluidRegions, solidRegions);

    //Create linear temperature gradient in each region
    #include "initTGrad.H"

    bool cyclesTolReached=false;

    while (!cyclesTolReached && pimples.run(runTime))
    {
        runTime++;
        #include "updateTimes.H"
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update controls
   	#include "updateTimeControls.H"
    	#include "readCycleControl.H"
	pimples.read();

        // Update preCalculated input Fields
        #include "readInputFields.H"

        // Coupled solution loop
	while (pimples.loop())
	{
            	    forAll(fluidRegions, i)
            	    {
            	    	Info<< "Solving the fluid region: " << fluidRegions[i].name() << endl;
                	#include "setRegionFluidFields.H"
                	#include "solveFluid.H"
            	    }

        	    forAll(solidRegions, i)
        	    {
        	        Info<< "Solving the solid region: "<< solidRegions[i].name() << endl;
			#include "updateSolidVars.H"
        	        #include "setRegionSolidFields.H"
        	        #include "solveSolid.H"
        	    }
	  }

	#include "updateTHEx.H"
	#include "calcWmag.H"

	#include "checkTolCycles.H"

	#include "writeIntegralVal.H"
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
