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
 
 #include "therCycleControl.H"

 // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
 
 namespace Foam
 {
     defineTypeNameAndDebug(therCycleControl, 0);
 }
 
 
 // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
Foam::therCycleControl::therCycleControl(const IOobject& io, const IOdictionary& cycleDict, const IOdictionary& transportProp, const scalarField& flowTimeLine, const word& firstTime, const double dt, const word& fluidName):
   cycleIDvars(flowTimeLine, firstTime, dt, 0),
   regIOobject(io),
   gradDict_(cycleDict.subDict("initialGrad")),
   hexDict_(cycleDict.subDict("heatExchangers")),
   rhof_(transportProp.lookup("rho")),
   cf_(transportProp.lookup("c")),
   fluidName_(fluidName),
   makeGrad(cycleDict.lookup("makeGrad")),
   coldName (hexDict_.lookup("coldHExName")),
   hotName (hexDict_.lookup("hotHExName"))
 {
	integVal_.open ("IntegralValues_"+ fluidName_, std::ios::trunc);

	integVal_ << std::setw(3) << " "<< std::left << std::setw(13) << "time [s]"<<"|"<< std::setw(3) << " "<< std::left << std::setw(13) <<"trel [s]" <<"|"<< std::setw(3) << " "<< std::left << std::setw(13) <<"Qcold [J]"<<"|"<< std::setw(3) << " "<< std::left << std::setw(13) <<"Qhot [J]" <<"|"<< std::setw(3) << " "<< std::left << std::setw(13) <<"Wmag [J]"<< "|"<< std::setw(3) << " "<< std::left << std::setw(13) <<"TcoldHEx [K]"<<"|"<< std::setw(3) << " "<< std::left << std::setw(13) <<"ThotHEx [K]"<<"\n";
	integVal_ <<"-------------------------------------------------------------------------------------------------------------------------\n";
	integVal_.close();
}
 
 // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
 
 Foam::therCycleControl::~therCycleControl()
 {}
 
 // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::therCycleControl::writeCycleVars(const double writeTime)
{
      	
	integVal_.open ("IntegralValues_"+ fluidName_, std::ios::app);
	integVal_ << std::setw(3) << " "<< std::left << std::setw(13) << writeTime<<"|"<< std::setw(3) << " "<< std::left << std::setw(13) << trel() <<"|"<< std::setw(3) << " "<< std::left << std::setw(13) <<coldSideQ<<"|"<< std::setw(3) << " "<< std::left << std::setw(13) <<hotSideQ <<"|"<< std::setw(3) << " "<< std::left << std::setw(13) << sumWmag << "|"<< std::setw(3) << " "<< std::left << std::setw(13) <<coldSideT<<"|"<< std::setw(3) << " "<< std::left << std::setw(13) <<hotSideT<<"\n";
	integVal_.close();

}

///////////////////


 void Foam::therCycleControl::initLinearTGrad(volScalarField& Tfluid, volScalarField& Tsolid)
{
        const double initTcold = readScalar(gradDict_.lookup("initialTcold"));
        const double initThot = readScalar(gradDict_.lookup("initialThot"));
        const vector gradDir = gradDict_.lookup("TgradVector");

	const vector minPoint = Tfluid.mesh().bounds().min();
	const vector maxPoint = Tfluid.mesh().bounds().max();

        const volVectorField cellCenters = Tfluid.mesh().C();

        forAll(Tfluid.internalField(),cellI)
	{
	  scalar& localT= Tfluid.ref()[cellI];
	  localT = initTcold + fabs((cellCenters[cellI] & gradDir)-(minPoint & gradDir))/fabs((maxPoint & gradDir)-(minPoint & gradDir))*(initThot-initTcold);
	}
	//Boundary
	forAll(Tfluid.boundaryField(),patchID)
	{
                const vectorField& faceCenters = Tfluid.boundaryFieldRef()[patchID].patch().Cf();

		forAll(Tfluid.boundaryField()[patchID],cellI)
		{
		  scalar& localT= Tfluid.boundaryFieldRef()[patchID][cellI];
	  	  localT = initTcold + fabs((faceCenters[cellI] & gradDir)-(minPoint & gradDir))/fabs((maxPoint & gradDir)-(minPoint & gradDir))*(initThot-initTcold);
		}
	}

       //Gradient on the assigned solid phase

        const volVectorField& cellCentersSolid = Tsolid.mesh().C();
	forAll(Tsolid.internalField(),cellI)
	{
		scalar& localT= Tsolid.ref()[cellI];
		localT = initTcold + fabs((cellCentersSolid[cellI] & gradDir)-(minPoint & gradDir))/fabs((maxPoint & gradDir)-(minPoint & gradDir))*(initThot-initTcold);
	}

	//Boundary
	forAll(Tsolid.boundaryField(),patchID)
	{
            const vectorField& faceCentersSolid = Tsolid.boundaryFieldRef()[patchID].patch().Cf();
	    forAll(Tsolid.boundaryField()[patchID],cellI)
	    {
		scalar& localT= Tsolid.boundaryFieldRef()[patchID][cellI];
		localT = initTcold + fabs((faceCentersSolid[cellI] & gradDir)-(minPoint & gradDir))/fabs((maxPoint & gradDir)-(minPoint & gradDir))*(initThot-initTcold);
	     }
	 }

}

////////////////////////////////////////////////////

void Foam::therCycleControl::updateHEx(const volScalarField& Tfluid, const surfaceScalarField& phiFlow)
{
        const double Isum = iphiSumRelStep();
        const double fixedQcold = readScalar(hexDict_.lookup("fixedQcold"));
        const double fixedQhot = readScalar(hexDict_.lookup("fixedQhot"));
        const double fixedTcold = readScalar(hexDict_.lookup("fixedTcold"));
        const double fixedThot = readScalar(hexDict_.lookup("fixedThot"));

        if (!(  (fixedTcold>0) && (fixedThot>0)   ))
             {
                 units ="[K]";
             }   
             else
             {
                 units ="[J]";
             }

        forAll(Tfluid.boundaryField(),patchID)
	{
		if (Tfluid.boundaryField()[patchID].type()=="heatExchangerBC")
		{
		        const fvsPatchField<double>& phiBound = phiFlow.boundaryField()[patchID];
			const fvPatchField<double>& TbcField = Tfluid.boundaryField()[patchID];

			scalar phiSum=gSum(phiBound);
		
		   	const word& sideName = TbcField.patch().name();

                        scalar Tavg=gSum(TbcField * TbcField.patch().magSf())/gSum(TbcField.patch().magSf());
                      
			if (sideName == coldName)
			{
				if (phiSum<0)
				{
					coldSideT = fixedQcold*dt()/(rhof_.value()*cf_.value()*phiSum)+coldSideT;
					coldSideQ = fixedQcold;
				}
				else if (phiSum>0)
				{

					if (!(Isum==0))
					{	
		                               coldSideQ=(coldSideT-(Tavg+coldSideT*(Isum-1))/Isum)*cf_.value()*fabs(phiSum)*rhof_.value()*dt();
						coldSideT = (Tavg+coldSideT*(Isum-1))/Isum;
					}
					else
					{
						coldSideQ = 0;
					}	
				}
				else
				{
					coldSideQ = 0;
				}

				if (fixedTcold>0)
				{
					coldSideQ = (fixedTcold-Tavg)*cf_.value()*fabs(phiSum)*rhof_.value()*dt();
					coldSideT = fixedTcold;
				}		
			}


			if (sideName == hotName)
			{

				if (phiSum<0)
				{
					hotSideT = fixedQhot*dt()/(rhof_.value()*cf_.value()*phiSum)+hotSideT;
					hotSideQ = fixedQhot;
				}		
				else if (phiSum>0)
				{

					if (!(Isum==0))
					{	
						hotSideQ = (hotSideT-(Tavg+hotSideT*(hotSideT-1))/Isum)*cf_.value()*fabs(phiSum)*rhof_.value()*dt();
						hotSideT = (Tavg+hotSideT*(Isum-1))/Isum;
					}
					else
					{
						hotSideQ = 0;
					}	

				}
				else
				{
					hotSideQ = 0;
				}

				if (fixedThot>0)
				{
					hotSideQ = (Tavg-fixedThot)*cf_.value()*fabs(phiSum)*rhof_.value()*dt();
					hotSideT = fixedThot;
				}

			}

		}
	}

}

////////////////////////////////////////////////////

bool Foam::therCycleControl::checkCycle()
{
   const double fixedTcold = readScalar(hexDict_.lookup("fixedTcold"));
   const double fixedThot = readScalar(hexDict_.lookup("fixedThot"));

   cycleAvgColdSideT = cycleAvgColdSideT + coldSideT*dt()/tau();
   cycleAvgHotSideT = cycleAvgHotSideT + hotSideT*dt()/tau();

   cycleColdSideQ = cycleColdSideQ + coldSideQ;
   cycleHotSideQ = cycleHotSideQ + hotSideQ;
   cycleSumWmag = cycleSumWmag + sumWmag;


        if( (trel()<(dt()/2)) || (trel()>(tau()-dt()/2)) )
	{
	     if (!(  (fixedTcold>0) && (fixedThot>0)   ))
             {
                 tolerance_ = max(fabs(cycleAvgColdSideT-OldCycleAvgColdSideT),fabs(cycleAvgHotSideT-OldCycleAvgHotSideT));
                // units ="[K]";
             }   
             else
             {
                 tolerance_ = max(fabs(cycleColdSideQ-OldCycleColdSideQ),fabs(cycleHotSideQ-OldCycleHotSideQ));
                 tolerance_ = max(tolerance_,fabs(cycleSumWmag-OldCycleSumWmag) );
              //   units ="[J]";
             }


             //Store and reset values for new cycle
             OldCycleAvgColdSideT=cycleAvgColdSideT;
             OldCycleAvgHotSideT=cycleAvgHotSideT;
             OldCycleColdSideQ=cycleColdSideQ;
             OldCycleHotSideQ=cycleHotSideQ;
             OldCycleSumWmag = cycleSumWmag;

             cycleAvgColdSideT = 0;
             cycleAvgHotSideT = 0;
             cycleColdSideQ = 0;
             cycleHotSideQ = 0;
             cycleSumWmag =0;

             // Check cycle tolerance
             if (tolerance_ <= cycleTol_ )
	     {
		return true;
             }
        }
        
        return false;
}

////////////////////////////////////////////////////

void Foam::therCycleControl::setCycleTol(double CycleSsTol)
{
       cycleTol_=CycleSsTol;
}

////////////////////////////////////////////////////

word Foam::therCycleControl::solidPhaseName()
{
       return gradDict_.lookup("solidPhaseName");
}


////////////////////////////////////////////////////

void Foam::therCycleControl::calcWmag(const magneticProp& magVars)
{
      	
      sumWmag=gSum(magVars.Wmag);
      reduce(sumWmag, sumOp<scalar>());//Sum over the processors
}

////////////////////////////////////////////////////
double Foam::therCycleControl::tolerance( )
{
       return tolerance_;
}
//////////////////////////////////////////////////////
 bool Foam::therCycleControl::writeData(Ostream& os) const
{
    return os.good();
}
 // ************************************************************************* //
