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

if (runTime.value() > tau) // If we've already launched a simulation => injection if t > tau
{
	subDictKinematicCloud.set("SOI", runTime.value() + 0.01*tau ); 
	// injection after 10 dt => rentre dans la run loop
}
else  //If not (tau>t) => injection at t = tau 
{
	subDictKinematicCloud.set("SOI", tau ); // injection after 	
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

scalar SOI = subDictKinematicCloud.lookup("SOI"); 







