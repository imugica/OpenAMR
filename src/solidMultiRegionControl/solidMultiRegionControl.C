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
 
 \*---------------------------------------------------------------------------*/

#include "solidMultiRegionControl.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidMultiRegionControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::Time& Foam::solidMultiRegionControl::time
(
    const PtrList<fvMesh>& solidMeshes
)
{
    if (solidMeshes.empty())
    {
        FatalErrorInFunction
            << "There needs to be at least one region"
            << exit(FatalError);
    }

    return solidMeshes[0].time();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidMultiRegionControl::solidMultiRegionControl
(
    PtrList<fvMesh>& solidMeshes,
    const word& algorithmName
)
:
    multiRegionSolutionControl(time(solidMeshes), algorithmName),
    pimpleLoop(static_cast<solutionControl&>(*this)),
    convergenceControl(static_cast<solutionControl&>(*this)),
    correctorConvergenceControl
    (
        static_cast<solutionControl&>(*this),
        "outerCorrector"
    ),
    solidControls_()
{
    bool allSteady = true, allTransient = true;

    forAll(solidMeshes, i)
    {
        solidControls_.append
        (
            new solidNoLoopControl(solidMeshes[i], algorithmName, *this)
        );

        allSteady = allSteady && solidMeshes[i].steady();
        allTransient = allTransient && solidMeshes[i].transient();
    }

    read();

    forAll(solidMeshes, i)
    {
        Info<< nl << algorithmName << ": Region " << solidMeshes[i].name();
        solidControls_[i].printResidualControls();

        if (nCorrPimple_ > 1)
        {
            Info<< nl << algorithmName << ": Region " << solidMeshes[i].name();
            solidControls_[i].printCorrResidualControls(nCorrPimple_);
        }
    }

    Info<< nl << algorithmName << ": Operating solver in "
        << (allSteady ? "steady-state" : allTransient ? "transient" :
            "mixed steady-state/transient") << " mode with " << nCorrPimple_
        << " outer corrector" << (nCorrPimple_ == 1 ? "" : "s") << nl;

    if ((allSteady || allTransient) && nCorrPimple_ == 1)
    {
        Info<< algorithmName << ": Operating solver in "
            << (allSteady ? "SIMPLE" : "PISO") << " mode" << nl;
    }

    Info<< nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidMultiRegionControl::~solidMultiRegionControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::solidMultiRegionControl::read()
{
    forAll(solidControls_, i)
    {
        if (!solidControls_[i].read())
        {
            return false;
        }
    }

    const dictionary& solutionDict = dict();

    nCorrPimple_ = solutionDict.lookupOrDefault<label>("nOuterCorrectors", 1);

    return true;
}


bool Foam::solidMultiRegionControl::hasResidualControls() const
{
    bool result = true;

    forAll(solidControls_, i)
    {
        result = result && solidControls_[i].hasResidualControls();
    }

    return result;
}


bool Foam::solidMultiRegionControl::hasCorrResidualControls() const
{
    bool result = true;

    forAll(solidControls_, i)
    {
        result = result && solidControls_[i].hasCorrResidualControls();
    }

    return result;
}


bool Foam::solidMultiRegionControl::criteriaSatisfied() const
{
    bool result = true;

    forAll(solidControls_, i)
    {
        result = solidControls_[i].criteriaSatisfied() && result;
    }

    return result;
}


bool Foam::solidMultiRegionControl::corrCriteriaSatisfied() const
{
    bool result = true;

    forAll(solidControls_, i)
    {
        result = solidControls_[i].corrCriteriaSatisfied() && result;
    }

    return result;
}


void Foam::solidMultiRegionControl::resetCorrSolveIndex()
{
 
    forAll(solidControls_, i)
    {
        solidControls_[i].resetCorrSolveIndex();
    }
}


void Foam::solidMultiRegionControl::updateCorrSolveIndex()
{
 
    forAll(solidControls_, i)
    {
        solidControls_[i].updateCorrSolveIndex();
    }
}


bool Foam::solidMultiRegionControl::loop()
{
    read();

    if (!pimpleLoop::loop(*this))
    {
        forAll(solidControls_, i)
        {
            solidControls_[i].updateFinal();
        }

        return false;
    }
    forAll(solidControls_, i)
    {
        solidControls_[i].storePrevIterFields();
    }

    forAll(solidControls_, i)
    {
        solidControls_[i].updateFinal();
    }

    return true;
}


bool Foam::solidMultiRegionControl::run(Time& time)
{
    read();

    if (!endIfConverged(time))
    {
        forAll(solidControls_, i)
        {
            solidControls_[i].storePrevIterFields();
        }
    }

    return time.run();
}


bool Foam::solidMultiRegionControl::loop(Time& time)
{
    read();

    if (!endIfConverged(time))
    {

        forAll(solidControls_, i)
        {
            solidControls_[i].storePrevIterFields();
        }
    }

    return time.loop();
}


// ************************************************************************* //
