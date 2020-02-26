forAll(fluidRegions, i)
{
   caloricCycleControl[i].updateHEx(Tfluid[i], phiFluid[i]);
}

