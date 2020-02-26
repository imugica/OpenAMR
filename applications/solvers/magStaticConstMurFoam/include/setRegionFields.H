    fvMesh& mesh = magneticRegions[i];

    volScalarField& psi = psiFields[i];
    dimensionedScalar& mur = murValues[i];
    volVectorField& H = HFields[i];
    solidNoLoopControl& pimple = pimples.solid(i);

