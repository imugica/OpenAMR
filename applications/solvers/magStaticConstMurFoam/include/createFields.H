// Initialise fluid field pointer lists

PtrList<volScalarField> psiFields(magneticRegions.size());
PtrList<volVectorField> HFields(magneticRegions.size());

PtrList<IOdictionary> magneticProperties(magneticRegions.size());
PtrList<dimensionedScalar> murValues(magneticRegions.size());

//Magnet dictionary
PtrList<IOdictionary> magnetDict(magneticRegions.size());

////////////////////////////////////////////
forAll(magneticRegions, i)
{

    Info<< "    Creating psi field for region: "<< magneticNames[i] <<"\n" << endl;

    psiFields.set
    (
        i,
        new volScalarField
    	(
        	IOobject
        	(
        	    "psi",
        	    runTime.timeName(),
        	    magneticRegions[i],
        	    IOobject::MUST_READ,
        	    IOobject::AUTO_WRITE
        	),
        	magneticRegions[i]
    	)
    );

    Info<< "    Creating H fieldfor region: "<< magneticNames[i] <<"\n" << endl;

    HFields.set
    (
        i,
        new volVectorField
    	(
        	IOobject
        	(
        	    "H",
        	    runTime.timeName(),
        	    magneticRegions[i],
        	    IOobject::NO_READ,
        	    IOobject::AUTO_WRITE
        	),
		constant::electromagnetic::mu0*fvc::reconstruct(fvc::snGrad(psiFields[i])*magneticRegions[i].magSf())
    	)
    );


    magneticProperties.set
    (
	i,
	new IOdictionary
	(
		IOobject
		(
         	"magneticProperties",
         	runTime.constant(),
         	magneticRegions[i],
         	IOobject::MUST_READ_IF_MODIFIED,
         	IOobject::NO_WRITE
        	)
        )
    );

  

    Info<< "    Reading mur for region: "<< magneticNames[i] <<"\n" << endl;
	murValues.set(
			i,
		new dimensionedScalar (								
					magneticProperties[i].lookup("mur")
				     )
		);
			


    magnetDict.set
    (
	i,
	new IOdictionary
	(
		IOobject
		(
         	"magnet",
         	runTime.constant(),
         	magneticRegions[i],
         	IOobject::READ_IF_PRESENT,
         	IOobject::NO_WRITE
        	)
        )
    );

}

