while(pimple.correctNonOrthogonal())
{

    fvScalarMatrix EEqn
    (
        fvm::ddt(T) + fvm::div(phi, T)
      - fvm::laplacian(lambdaFluid[i]/(rhoFluid[i]*cFluid[i]), T)
     == 
      	fvOptions(T)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(T);
}
    Info<< "Min/max T:" << min(T()).value() << ' ' << max(T()).value() << endl;
    Info << endl;
