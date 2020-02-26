forAll(fluidRegions, i)
{
        cyclesTolReached = cyclesTolReached && caloricCycleControl[i].checkCycle() ;
	std::cout << "\n---------------- End of time step ----------------\n" ;
	std::cout << std::setw(30) << std::left << "Region: " << std::setw(15) << std::right << fluidNames[i] << "\n";
	std::cout << std::setw(30) << std::left << "Cold heat exchanger T: "<< std::setw(15) << std::right<<  caloricCycleControl[i].coldSideT <<"[K]\n";
	std::cout << std::setw(30) << std::left  <<"Hot heat exchanger T: "<<  std::setw(15) << std::right<< caloricCycleControl[i].hotSideT <<"[K]\n";
	std::cout << std::setw(30) << std::left  << "Heat exchanged at cold side: "<< std::setw(15) << std::right <<  caloricCycleControl[i].coldSideQ <<"[J]\n";
	std::cout << std::setw(30) << std::left <<  "Heat exchanged at hot side: "<< std::setw(15) << std::right << caloricCycleControl[i].hotSideQ <<"[J]\n";
	std::cout << std::setw(30) << std::left <<  "Magnetic work: "<< std::setw(15) << std::right<< caloricCycleControl[i].sumWmag <<"[J]\n";
	std::cout  << "--------------------------------------------------\n";
	std::cout << "\n";



        if( (caloricCycleControl[i].trel()<(caloricCycleControl[i].dt()/2)) || (caloricCycleControl[i].trel()>(caloricCycleControl[i].tau()-caloricCycleControl[i].dt()/2)) )
	{
 	std::cout << "\n########## End of cycle for region: "<< fluidNames[i]  <<" ###########\n" ;
	std::cout << std::setw(30) << std::left << "Tolerance reached towards SS: " << std::setw(15) << std::right << caloricCycleControl[i].tolerance() << caloricCycleControl[i].units << "\n";
	std::cout << std::setw(30) << std::left << "Average cold side T: "<< std::setw(15) << std::right<<  caloricCycleControl[i].cycleAvgColdSideT <<"[K]\n";
	std::cout << std::setw(30) << std::left  <<"Average hot side T: "<<  std::setw(15) << std::right<< caloricCycleControl[i].cycleAvgHotSideT <<"[K]\n";
	std::cout << std::setw(30) << std::left  << "Cooling power: "<< std::setw(15) << std::right <<  caloricCycleControl[i].cycleColdSideQ/caloricCycleControl[i].tau() << "[W]\n";
	std::cout << std::setw(30) << std::left <<  "Pumped heat power: "<< std::setw(15) << std::right << caloricCycleControl[i].cycleHotSideQ/caloricCycleControl[i].tau() << "[W]\n";
	std::cout << std::setw(30) << std::left <<  "Magnetic power: "<< std::setw(15) << std::right<< caloricCycleControl[i].cycleSumWmag/caloricCycleControl[i].tau() << "[W]\n";
	std::cout  << "####################################################\n";
	std::cout << "\n";
         }

}



