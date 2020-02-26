 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     | 
     \\  /    A nd           | Copyright (C) 2020 Ibai Mugica
      \\/     M anipulation  |
 -------------------------------------------------------------------------------
 License
     This file is derivative work of OpenFOAM.
 
     OpenFOAM is free software: you can redistribute it and/or modify it
     under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
 
     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
     for more details.
 
     You should have received a copy of the GNU General Public License
     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 
 \*---------------------------------------------------------------------------*/
 
 #include "magneticProp.H"

 // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
 
 namespace Foam
 {
     defineTypeNameAndDebug(magneticProp, 0);
 }
 
 
 // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
Foam::magneticProp::magneticProp( const fvMesh& solidMesh, const IOdictionary& Hdict, const word& firstTime, const double rhoSolid, const double dt):

   cycleIDvars(Hdict.lookup("TimeFrames"), firstTime, dt, Hdict.lookup("symmetric")),
   symH_(Hdict.lookup("symmetric")),
   CHtable_("constant/"+solidMesh.name()+"/CHtable"),
   dMdTtable_("constant/"+solidMesh.name()+"/dMdTtable"),
   Mtable_("constant/"+solidMesh.name()+"/Mtable"),
   firstTime_(firstTime),
   solidMesh_(solidMesh),
   rho_(rhoSolid),
   Hmag
        (
            IOobject
            (
                "Hmag",
		firstTime_,
                solidMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh_,
	    dimensionedScalar("Hmag",dimensionSet(1,0,-2,0,0,-1,0),0.0)
        ),
    CH
        (
            IOobject
            (
                "CH",
		firstTime_,
                solidMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh_,
	    dimensionedScalar("CH",dimensionSet(0,2,-2,-1,0,0,0),0.0)
        ),

    CT
        (
            IOobject
            (
                "CT",
                firstTime_,
                solidMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh_,
	    dimensionedScalar("CT",dimensionSet(-1,2,0,0,0,1,0),0.0)
        ),
   Wmag
        (
            IOobject
            (
                "Wmag",
                firstTime_,
                solidMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh_,
             0.0
        )
{ 
}
 
  
Foam::magneticProp::magneticProp( const fvMesh& solidMesh, const word& firstTime, const double rhoSolid, const double dt):
   cycleIDvars(),
   symH_(0),
   CHtable_("constant/"+solidMesh.name()+"/CHtable"),
   dMdTtable_("constant/"+solidMesh.name()+"/dMdTtable"),
   Mtable_("constant/"+solidMesh.name()+"/Mtable"),
   firstTime_(firstTime),
   solidMesh_(solidMesh),
   rho_(rhoSolid),
   Hmag
        (
            IOobject
            (
                "Hmag",
		firstTime_,
                solidMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh_,
	    dimensionedScalar("Hmag",dimensionSet(1,0,-2,0,0,-1,0),0.0)
        ),
    CH
        (
            IOobject
            (
                "CH",
		firstTime_,
                solidMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh_,
	    dimensionedScalar("CH",dimensionSet(0,2,-2,-1,0,0,0),0.0)
        ),

    CT
        (
            IOobject
            (
                "CT",
                firstTime_,
                solidMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh_,
	    dimensionedScalar("CT",dimensionSet(-1,2,0,0,0,1,0),0.0)
        ),
   Wmag
        (
            IOobject
            (
                "Wmag",
                firstTime_,
                solidMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidMesh_,
             0.0
        )

{ 
}
 
 // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
 
 Foam::magneticProp::~magneticProp()
 {}
 
 
 // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


 void Foam::magneticProp::calcHmag(const volVectorField& Hvector)
{

      forAll(Hvector.internalField(),cellI)
      {
	 scalar& localHmag= Hmag.ref()[cellI];
	 localHmag = mag(Hvector.internalField()[cellI]);
      }

      forAll(Hvector.boundaryField(),patchID)
      {
	      forAll(Hvector.boundaryField()[patchID],cellI)
	      {
		    scalar& localHmag=  Hmag.boundaryFieldRef()[patchID][cellI];
		    localHmag = mag(Hvector.boundaryField()[patchID][cellI]);
	      }
      }
}



void Foam::magneticProp::updateCH(const volScalarField& T)
{
      forAll(CH.internalField(),cellI)
      {
	  scalar& localValue= CH.ref()[cellI];
	  localValue = CHtable_(Hmag.internalField()[cellI],T.internalField()[cellI]);
      }

      forAll(CH.boundaryField(),patchID)
      {
	    forAll(CH.boundaryField()[patchID],cellI)
	    {
		scalar& localValue= CH.boundaryFieldRef()[patchID][cellI];
		localValue = CHtable_(Hmag.boundaryField()[patchID][cellI],T.boundaryField()[patchID][cellI]);
            }
      }
}



void Foam::magneticProp::updateCT(const volScalarField& T)
{
      forAll(CT.internalField(),cellI)
      {
	  scalar& localValue= CT.ref()[cellI];
	  localValue = dMdTtable_(Hmag.internalField()[cellI],T.internalField()[cellI])*T.internalField()[cellI];
      }

      forAll(CT.boundaryField(),patchID)
      {
	    forAll(CH.boundaryField()[patchID],cellI)
	    {
		scalar& localValue= CT.boundaryFieldRef()[patchID][cellI];
		localValue = dMdTtable_(Hmag.boundaryField()[patchID][cellI],T.boundaryField()[patchID][cellI])*T.boundaryField()[patchID][cellI];
            }
      }
}



void Foam::magneticProp::updateWmag(const volScalarField& T)
{
      forAll(Wmag.internalField(),cellI)
      {
	  scalar& localValue= Wmag.ref()[cellI];
          const double M1=Mtable_(Hmag.oldTime().internalField()[cellI],T.oldTime().internalField()[cellI]);
          const double M2=Mtable_(Hmag.internalField()[cellI],T.internalField()[cellI]);
	  localValue = solidMesh_.V()[cellI]*(Hmag.internalField()[cellI]+Hmag.oldTime().internalField()[cellI])*(M2-M1)*rho_/2;

      }

      forAll(Wmag.boundaryField(),patchID)
      {
	    forAll(Wmag.boundaryField()[patchID],cellI)
	    {
		scalar& localValue= Wmag.boundaryFieldRef()[patchID][cellI];
                const double M1=Mtable_(Hmag.oldTime().boundaryField()[patchID][cellI],T.oldTime().boundaryField()[patchID][cellI]);
                const double M2=Mtable_(Hmag.boundaryField()[patchID][cellI],T.boundaryField()[patchID][cellI]);
		localValue = solidMesh_.V()[cellI]*(Hmag.boundaryField()[patchID][cellI]+Hmag.oldTime().boundaryField()[patchID][cellI])*(M2-M1)*rho_/2;


            }
      }
}
 
 // ************************************************************************* //
