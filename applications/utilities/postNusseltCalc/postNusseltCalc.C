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
    postNusseltCalc

Description
    Calculates and writes the average Nusselt number of the boundaries.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        #include "createFields.H"

        surfaceScalarField heatFlux
        (
            fvc::snGrad(T)*Lref/(T0-Tinf)
        );

        const surfaceScalarField::Boundary& patchHeatFlux =
            heatFlux.boundaryField();

        const surfaceScalarField::Boundary& magSf =
            mesh.magSf().boundaryField();

        Info<< "\n Average Nu [-]:" << endl;
        forAll(patchHeatFlux, patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                scalar convFlux = gSum(patchHeatFlux[patchi]/magSf[patchi]);

                Info<< mesh.boundary()[patchi].name() << endl
                    <<  convFlux << endl;
            }
        }
        Info<< endl;


        volScalarField Nu
        (
            IOobject
            (
                "Nu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Nu", heatFlux.dimensions(), 0.0)
        );

        volScalarField::Boundary& NuBf =
            Nu.boundaryFieldRef();

        forAll(NuBf, patchi)
        {
            NuBf[patchi] = patchHeatFlux[patchi];
        }

        Nu.write();

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
