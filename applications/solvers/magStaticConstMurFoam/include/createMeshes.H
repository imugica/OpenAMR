    regionProperties rp(runTime);

    const wordList magneticNames(rp["magnetic"]);

    PtrList<fvMesh> magneticRegions(magneticNames.size());

    forAll(magneticNames, i)
    {
        Info<< "Create magnetic mesh for region " << magneticNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        magneticRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    magneticNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

    }
