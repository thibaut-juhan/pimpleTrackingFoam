Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField maxDifference
(
    IOobject
    (
        "maxDifference",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh

);

Info<< "Reading simulationProperties \n" << endl;
IOdictionary geometryProperties
(
    IOobject
    (
     "geometryProperties",
	 runTime.constant(),
	 mesh,
	 IOobject::MUST_READ_IF_MODIFIED,
	 IOobject::NO_WRITE   
    )
);

Info<<"Reading R\n" << endl;
dimensionedScalar R ("R",dimLength,geometryProperties);
Info<<"Reading Rinj \n" << endl;
dimensionedScalar rInj ("rInjection",dimLength,geometryProperties);
Info<<"Reading L \n" << endl;
dimensionedScalar L ("L",dimLength,geometryProperties);
Info<<"Reading velInitial\n" << endl;
dimensionedScalar velInitial ("velInitial",dimVelocity,geometryProperties);
Info << "Reading injector centers" << endl;
dimensionedScalar xC1 ("xC1",dimLength,geometryProperties);
dimensionedScalar xC2 ("xC2",dimLength,geometryProperties);
dimensionedScalar yC1 ("yC1",dimLength,geometryProperties);
dimensionedScalar yC2 ("yC2",dimLength,geometryProperties);


#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
singlePhaseTransportModel laminarTransport(U, phi);
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
#include "createMRF.H"
#include "createFvOptions.H"
