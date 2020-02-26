volScalarField& T = Tsolid[i];

volScalarField& CH = magneticVars[i].CH;
volScalarField& CT = magneticVars[i].CT;
volScalarField& Hmag = magneticVars[i].Hmag;

solidNoLoopControl& pimple = pimples.solid(i);

fv::options& fvOptions = solidHeatSources[i];
