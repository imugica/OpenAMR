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
 
#include "cycleIDvars.H"   

 // * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
 
 namespace Foam
 {
     defineTypeNameAndDebug(cycleIDvars, 0);
 }
 
 
 // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
 Foam::cycleIDvars::cycleIDvars(const scalarField& inputTimeLine, const word firstTime, const double dt, const bool symH): 
timeLine_(inputTimeLine),
dt_(dt),
firstTime_(firstTime)
 {
     
        tau_ =gSum(timeLine_);

        Nsteps_=timeLine_.size();
   
        update(std::stod(firstTime), symH);

 }
 
 Foam::cycleIDvars::cycleIDvars()
 {}


 // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
 
 Foam::cycleIDvars::~cycleIDvars()
 {}
 
 
 // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

 void Foam::cycleIDvars::update(double timeValue, const bool symH)
 {
	trel_=(timeValue/tau_-floor(timeValue/tau_))*tau_;

	tcompareLow_=0;

	tcompareUp_=0;

        fileN_=0;

      if (!symH)
      {
	   for( int j = 0; j < Nsteps_ ; j++ ) 
	   {
	        tcompareUp_=tcompareUp_+timeLine_[j];

		if ((trel_>=tcompareLow_) && (trel_<tcompareUp_))
		{
			fileN_=j;
		}
		tcompareLow_=tcompareLow_+timeLine_[j];
	   }

       }
       else
       {
                for( int j = 0; j < 2*Nsteps_ ; j++ ) 
		{

			if (j<Nsteps_)
			{
				tcompareUp_=tcompareUp_+timeLine_[j];
				if ((trel_>tcompareLow_) && (trel_<=tcompareUp_))
				{
					fileN_=j;
				}
				tcompareLow_=tcompareLow_+timeLine_[j];
			}
			else
			{
				const int k=2*Nsteps_-j-1;
				tcompareUp_=tcompareUp_+timeLine_[k];
				if ((trel_>tcompareLow_) && (trel_<=tcompareUp_))
				{
					fileN_=k;
				}
				tcompareLow_=tcompareLow_+timeLine_[k];
			}
		}

       }

}

////////////////////////////////////////////////////

std::string Foam::cycleIDvars::fileN() const
{
       return std::to_string(fileN_);
}


////////////////////////////////////////////////////

void Foam::cycleIDvars::setDeltaT(scalar deltaT)
{
       dt_=deltaT;
}

////////////////////////////////////////////////////

int Foam::cycleIDvars::iphiSumRelStep()
{

        int itrel=int(round(trel_ /dt_));
        double tLow=0;
        double tUp=0;
        int itLow=0;
        int itUp=0;
        int itCompareLow=0;

	for( int j = 0; j < Nsteps_ ; j++ ) 
	{
		tUp=tUp+timeLine_[j];

		itUp=int(round((tUp/tau_-floor(tUp/tau_))*tau_/dt_));

		itLow=int(round((tLow/tau_-floor(tLow/tau_))*tau_/dt_)); 


		if ((itrel>itLow) && (itrel<=itUp))
		{
			itCompareLow=itLow;

		}
		tLow=tLow+timeLine_[j];
	}

        return itrel-itCompareLow;
}
 

///////////////////////////////////////////

double Foam::cycleIDvars::dt()
{
       return dt_;
}

///////////////////////////////////////////

double Foam::cycleIDvars::tau()
{
       return tau_;
}
///////////////////////////////////////////

double Foam::cycleIDvars::trel()
{
       return trel_;
}

///////////////////////////////////////////

word Foam::cycleIDvars::firstTime() 
{
       return firstTime_;
}
 // ************************************************************************* //
