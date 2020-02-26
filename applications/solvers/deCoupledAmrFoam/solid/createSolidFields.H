// Initialise solid field pointer lists
PtrList<volScalarField> Tsolid(solidRegions.size());

//Magnetic Class
PtrList<magneticProp> magneticVars(solidRegions.size());
PtrList<IOdictionary> HDictTimeLine(solidRegions.size());
PtrList<volVectorField> HSolid(solidRegions.size());
/////////////////////////////////////////////

///Properties
PtrList<IOdictionary> transportPropertiesSolid(solidRegions.size());
PtrList<dimensionedScalar> rhoSolid(solidRegions.size());
PtrList<dimensionedScalar> lambdaSolid(solidRegions.size());

//fvOptions
PtrList<fv::options> solidHeatSources(solidRegions.size());


// Populate solid field pointer lists
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

forAll(solidRegions, i)
{
    Info<< "    Creating solid fields for region: "<< solidNames[i] << endl;


    HDictTimeLine.set
    (
	i,
	new IOdictionary
	(
		IOobject
		(
         	"HTimeLineDict",
         	runTime.constant(),
         	solidRegions[i],
         	IOobject::MUST_READ_IF_MODIFIED,
         	IOobject::NO_WRITE
        	)
        )
    );


/////////////////////////////////////////////

 magneticVars.set( i, new magneticProp( 
                                 solidRegions[i], 
                                 HDictTimeLine[i], 
                                 runTime.timeName(),
                                 rhoSolid[i].value(),
                                 runTime.deltaT().value()
                                )
                  );

//////////////////////////////////////////////////////////////////////////////////

    Tsolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

/////////////////////////////////////////////

    Info<< "    Reading magnetic field for region: "<< solidNames[i] << endl;

    HSolid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "H"+magneticVars[i].fileN(),
		runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i]
        )
    );


//////////////////////////////////////////////////////////
    Info<< "    Reading solid properties for region: "<< solidNames[i] << endl;

    transportPropertiesSolid.set
    (
	i,
	new IOdictionary
	(
		IOobject
		(
         	"transportProperties",
         	runTime.constant(),
         	solidRegions[i],
         	IOobject::MUST_READ_IF_MODIFIED,
         	IOobject::NO_WRITE
        	)
        )
    );


	rhoSolid.set(
			i,
		new dimensionedScalar (								
					transportPropertiesSolid[i].lookup("rho")
				     )
		);	

	lambdaSolid.set(
			i,
		new dimensionedScalar (								
					transportPropertiesSolid[i].lookup("lambda")
				     )
		);


//////////////////////////////////////////////////////////

        Info<< "    Adding Solid fvOptions\n" << endl;
        solidHeatSources.set
        (
            i,
            new fv::options(solidRegions[i])
        );

			        
}
