
volScalarField T
(
	IOobject
	(
		"T",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	),
	mesh
);

IOdictionary reference
	(
		IOobject
		(
         	"NuRefDict",
         	runTime.constant(),
         	mesh,
         	IOobject::MUST_READ_IF_MODIFIED,
         	IOobject::NO_WRITE
        	)
	);


dimensionedScalar T0=reference.lookup("T0") ;
dimensionedScalar Tinf=reference.lookup("Tinf") ;
dimensionedScalar Lref=reference.lookup("Lref") ;


