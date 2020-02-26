// Initialise fluid field pointer lists
PtrList<volScalarField> Tfluid(fluidRegions.size());

//Flow Class
PtrList<IOdictionary> phiDictTimeLine(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());

//HeatExchanger Class
PtrList<IOdictionary> cycleDict(fluidRegions.size());
PtrList<therCycleControl> caloricCycleControl(fluidRegions.size());

///Properties
PtrList<IOdictionary> transportPropertiesFluid(fluidRegions.size());
PtrList<dimensionedScalar> rhoFluid(fluidRegions.size());
PtrList<dimensionedScalar> cFluid(fluidRegions.size());
PtrList<dimensionedScalar> lambdaFluid(fluidRegions.size());

//fvOptions
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

// Populate fluid field pointer lists
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


forAll(fluidRegions, i)
{
    Info<< "    Creating fluid fields for region: "<< fluidNames[i] << endl;

     phiDictTimeLine.set( i, new IOdictionary
	(
		IOobject
		(
         	"phiTimeLineDict",
         	runTime.constant(),
         	fluidRegions[i],
         	IOobject::MUST_READ_IF_MODIFIED,
         	IOobject::NO_WRITE
        	)
        )
    );


    Info<< "    Adding to Tfluid\n" << endl;
    Tfluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Info<< "    Reading fluid transport properties for region: "<< fluidNames[i] << endl;

    transportPropertiesFluid.set
    (
	i,
	new IOdictionary
	(
		IOobject
		(
         	"transportProperties",
         	runTime.constant(),
         	fluidRegions[i],
         	IOobject::MUST_READ_IF_MODIFIED,
         	IOobject::NO_WRITE
        	)
        )
    );

  
	rhoFluid.set(
			i,
		new dimensionedScalar (								
					transportPropertiesFluid[i].lookup("rho")
				     )
		);
			

	cFluid.set(
			i,
		new dimensionedScalar (								
					transportPropertiesFluid[i].lookup("c")
				     )
		);
			

	lambdaFluid.set(
			i,
		new dimensionedScalar (								
					transportPropertiesFluid[i].lookup("lambda")
				     )
		);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Info<< "    Reading thermodynamic cycle dictionary for region: "<< fluidNames[i]  << endl;

   cycleDict.set
    (
	i,
	new IOdictionary
	(
		IOobject
		(
         	"cycleDict",
         	runTime.constant(),
         	fluidRegions[i],
         	IOobject::MUST_READ_IF_MODIFIED,
         	IOobject::NO_WRITE
        	)
        )
    );

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    caloricCycleControl.set( i, new therCycleControl(  IOobject(
                                                                   "caloricCycleControl",
                                                                   runTime.timeName(),
                                                                   fluidRegions[i],
                                                                   IOobject::NO_READ,
                                                                   IOobject::NO_WRITE
                                                               ),
				                      cycleDict[i], 
				                      transportPropertiesFluid[i], 
                                                      phiDictTimeLine[i].lookup("TimeFrames"), 
                                                      runTime.timeName(),
                                                      runTime.deltaT().value(),
						      fluidNames[i]
                                                    )
                           );			
   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Info<< "    Reading fluid flow for region: "<< fluidNames[i] << endl;

    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi"+caloricCycleControl[i].fileN(),
		runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i]
        )
    );

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


 Info<< "    Adding fvOptions for region: "<< fluidNames[i]  << endl;

    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );



}



