
rhoM = (gT*rhoL)+((scalar(1)-gT)*rhoS);

//Update the rho and cp for material
rho == alphaM*rhoM + (scalar(1) - alphaM)*rhoG;

scalar rhoTMetal = 0.0;

forAll(rhoT,cellI)
{
    rhoTMetal = rhoM[cellI];
    if(T[cellI]>Tl.value())
    { 
        rhoTMetal  = rhoLCoeff.component(0).value() + rhoLCoeff.component(1).value()*T[cellI] + rhoLCoeff.component(2).value()*T[cellI];
    }
    else if(T[cellI]>=TVap.value()) // TODO this else if will never execute because the first if is satisfied 
    {
        rhoTMetal  = rhoLCoeff.component(0).value() + rhoLCoeff.component(1).value()*TVap.value() + rhoLCoeff.component(2).value()*TVap.value(); 
    }
    rhoT[cellI] = (alphaM[cellI]*rhoTMetal)+((1-alphaM[cellI])*rhoG.value());
}

rhoT.correctBoundaryConditions();

//Update cp based on Temperature dependent solid cp
forAll(cpS.internalField(), cellI)
{
    if (T[cellI]<Ts.value())
    {
        cpS[cellI] = cpSCoeff.component(0).value() + (cpSCoeff.component(1).value() * T[cellI]) + (cpSCoeff.component(2).value() * pow(T[cellI],2));
    }
    else
    {
       cpS[cellI] = cpSCoeff.component(0).value() + (cpSCoeff.component(1).value() * Ts.value()) + (cpSCoeff.component(2).value() * pow(Ts.value(),2));
    }
    //cM[cellI] = gT[cellI]*cpL.value() + (scalar(1.0)-gT[cellI])*cpS[cellI];
}

cpS.correctBoundaryConditions();

cM = gT*cpL + (scalar(1.0)-gT)*cpS;

cM.correctBoundaryConditions();

//rhoCp == alphaM * rhoM * cM+(scalar(1)-alphaM) * rhoG * cpG;
//rhoCp == ((alphaM*cM)+((scalar(1)-alphaM)*cpG))*rho;
//Update momentum sink term
//A = -1.0 * C * alphaM * (Foam::pow(scalar(1)-gT,2))/(Foam::pow(gT,3)+SMALL);

//Update Effective thermal conductivity
forAll(kM.internalField(), cellI)
{
    if(T[cellI]<Ts.value())
    {
          kM[cellI]=kSCoeff.component(0).value()+(kSCoeff.component(1).value()*T[cellI])+(kSCoeff.component(2).value()*pow(T[cellI],2));
    }
    else if (T[cellI]>=Ts.value() && T[cellI]<Tl.value())
    {
        kM[cellI]= gT[cellI]*(kLCoeff.component(0).value()+kLCoeff.component(1).value()*Tl.value() + kLCoeff.component(2).value()*pow(Tl.value(),2)) + (1-gT[cellI])*(kSCoeff.component(0).value()+kSCoeff.component(1).value()*Ts.value()+kSCoeff.component(2).value()*pow(Ts.value(),2));
    }
    else if (T[cellI]>=Tl.value() && T[cellI]<TVap.value())
    {
        kM[cellI]= kLCoeff.component(0).value()+(kLCoeff.component(1).value()*T[cellI]) + (kLCoeff.component(2).value()*pow(T[cellI],2));
    }
    else if (T[cellI]>=TVap.value())
    {
        kM[cellI]= (kLCoeff.component(0).value()+kLCoeff.component(1).value()*TVap.value() + kLCoeff.component(2).value()*pow(TVap.value(),2));
    }
}

kM.correctBoundaryConditions();

kEff = alphaM*(kM)+(scalar(1)-alphaM)*kG;

//kEff = (kM*kG)/((alphaM*kG)+((scalar(1.0)-alphaM)*kM)); //Harmonic interpolation
kEff.correctBoundaryConditions();
