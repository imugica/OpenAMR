phiFluid.clear();
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
forAll(fluidRegions, i)
{
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi"+caloricCycleControl[i].fileN(),
		caloricCycleControl[i].firstTime(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i]
        )
    );

}

///////////////////////////////////////////////////////////


HSolid.clear();
PtrList<volVectorField> HSolid(solidRegions.size());
forAll(solidRegions, i)
{
    HSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "H"+magneticVars[i].fileN(),
		magneticVars[i].firstTime(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i]
        )
    );

    magneticVars[i].calcHmag(HSolid[i]);
}
       


