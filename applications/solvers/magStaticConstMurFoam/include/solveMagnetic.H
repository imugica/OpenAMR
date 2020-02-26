
while (pimple.correctNonOrthogonal())
{
     //Solve psi
    fvScalarMatrix psiEqn
    (
	fvm::laplacian(mur, psi)   
    );

    psiEqn.relax();

    psiEqn.solve();

     //reconstruct H
     H=constant::electromagnetic::mu0*fvc::reconstruct(fvc::snGrad(psi)*mesh.magSf());
}

Info<< "Min/max psi:" << min(psi()).value() << ' '  << max(psi()).value() << endl;

