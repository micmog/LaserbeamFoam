// NOTE: since we now define the properties using the new dimensionedScalar way, we dont define their dimension here but in /constant/thermalProperties

IOdictionary thermalProperties
(
    IOobject
    (
      	"thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
     )
);

dictionary metal = thermalProperties.subDict("metalProperties");
dictionary gas = thermalProperties.subDict("gasProperties");

Info<< "\nReading thermalProperties\n" << endl;
//const dimensionedScalar betaM (metal.lookup("betaM"));
const dimensionedScalar Tl ("Tl", metal);
const dimensionedScalar Ts ("Ts", metal);
const dimensionedScalar Hf ("Hf", metal);
const dimensionedScalar C ("C", metal);
//const dimensionedScalar betaM ("betaM", dimensionSet(0,0,0,-1,0,0,0),metal.lookup("betaM"));

const dimensionedVector kSCoeff ("kSCoeff",  metal); //Coefficients of the solid thermal conductivity polynomial
const dimensionedVector kLCoeff ("kLCoeff",  metal); //Coefficients of the solid thermal conductivity polynomial
const dimensionedVector cpSCoeff ("cpSCoeff",metal);
const dimensionedScalar cpL ("cpL", metal);

const dimensionedScalar kG ("kG",gas); //Constant gas thermal conductivity
const dimensionedScalar cpG ("cpG", gas);

const dimensionedVector rhoLCoeff ("rhoLCoeff",metal);

//Vapour Pressure parameters
Info<<"Reading vapour pressure parameters \n" <<endl;

const dimensionedScalar p0 ("p0", metal); // Atmospheric pressure
const dimensionedScalar TVap ("TVap",metal); // Vaporisation temperature
const dimensionedScalar Hv ("Hv",metal); // Enthalpy change due to vaporisation
const dimensionedScalar R ("R",metal); // Universal gas constant
const dimensionedScalar Mm ("Mm", metal); //Molar mass of metal

const dimensionedScalar TRef ("TRef",thermalProperties);
const dimensionedScalar hC ("hC",thermalProperties);

dimensionedScalar TSmooth("TSmooth",dimensionSet(0, 0, 0, 1, 0),50.0);
