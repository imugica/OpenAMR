/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

scalar deltaT =
    runTime.controlDict().lookupOrDefault<scalar>("deltaT", GREAT);

runTime.setDeltaT (deltaT);

// ************************************************************************* //
