/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

scalar CycleSsTol = runTime.controlDict().lookupOrDefault<scalar>("CycleSsTolerance", 1E-3);

forAll(fluidRegions, i)
{
      caloricCycleControl[i].setCycleTol( CycleSsTol );
}



// ************************************************************************* //
