#include "readGravitationalAcceleration.H"
dimensionedScalar rhoInfValue ("rhoInf",dimDensity,laminarTransport); // fluid density == particles density
volScalarField rhoInf
(
	IOobject
	(
		"rho",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	rhoInfValue
);
volScalarField mu
(
	IOobject
	(
		"mu",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),

	laminarTransport.nu() * rhoInfValue
);
// add trackingParticles: 
const word kinematicCloudName // cloudName
(
	args.getOrDefault<word>("cloudName","kinematicCloud")
);
IOdictionary kinematicCloudProperties // update dictionary => change the SOI
(
    IOobject
    (
    	        "kinematicCloudProperties",
	 	runTime.constant(),
	 	mesh,
	 	IOobject::MUST_READ_IF_MODIFIED, // must read => modify 
	 	IOobject::AUTO_WRITE   
    )
);
dictionary& subDictKinematicCloud = kinematicCloudProperties.subDict("subModels").subDict("injectionModels").subDict("topInjection");
scalar tau = L.value()*R.value()*R.value()/( 2.*rInj.value()*rInj.value()*velInitial.value() );

if (runTime.value() > tau) // If we've already launched a simulation and the inital time > 2.tau => injected after 0.1tau
{
	subDictKinematicCloud.set("SOI", 0.01*tau + runTime.value()); // injection after 
}
else //If note : injected the particles after 5tau ! 
{
	subDictKinematicCloud.set("SOI", 0.1*tau + runTime.value() ); // injection after 	
}
kinematicCloudProperties.regIOobject::write();
Info << " Constructing kinematicCloud "  << kinematicCloudName << endl;
basicKinematicCollidingCloud kinematicCloud
(
	kinematicCloudName,
	rhoInf,
	U,
	mu,
	g
);



