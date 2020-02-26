#include "createMeshes.H"

if (!magneticRegions.size())
{
    FatalErrorIn(args.executable())
        << "No magneticRegions meshes present" << exit(FatalError);
}

fvMesh& mesh = magneticRegions[0];
