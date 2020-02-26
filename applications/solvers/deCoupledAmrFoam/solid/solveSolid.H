while(pimple.correctNonOrthogonal())
{
    fvScalarMatrix TEqn
    (
        CH*fvm::ddt(T) + CT*fvc::ddt(Hmag)
      
      - fvm::laplacian(lambdaSolid[i]/rhoSolid[i], T)
     == 
        CH*fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);
}

Info<< "Min/max T:" << min(T()).value() << ' ' << max(T()).value() << endl;
Info <<endl;


