if  (Pstream::master())
{
	forAll(fluidRegions, i)
	{
             caloricCycleControl[i].writeCycleVars(runTime.value());     
	}
}


