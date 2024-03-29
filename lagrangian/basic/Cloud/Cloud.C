/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017, 2020 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "Cloud.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "mapPolyMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "wallPolyPatch.H"
#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::checkPatches() const
{
    const polyBoundaryMesh& pbm = polyMesh_.boundaryMesh();
    bool ok = true;
    for (const polyPatch& pp : pbm)
    {
        if (isA<cyclicAMIPolyPatch>(pp))
        {
            const cyclicAMIPolyPatch& cami =
                refCast<const cyclicAMIPolyPatch>(pp);

            if (cami.owner())
            {
                ok = ok && (cami.AMI().singlePatchProc() != -1);
            }
        }
    }

    if (!ok)
    {
        FatalErrorInFunction
            << "Particle tracking across AMI patches is only currently "
            << "supported for cases where the AMI patches reside on a "
            << "single processor" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh, cloudName),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    labels_(),
    globalPositionsPtr_(),
    geometryType_(cloud::geometryType::COORDINATES)
{
    checkPatches();
    polyMesh_.oldCellCentres();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();

    if (particles.size())
    {
        IDLList<ParticleType>::operator=(particles);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::addParticle(ParticleType* pPtr)
{
    this->append(pPtr);
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::deleteParticle(ParticleType& p)
{
    delete(this->remove(&p));
}
/*
template <class ParticleType>
void Foam::Cloud<ParticleType>::reinjectionParticles()
{
    // create a new particle from this position : 
   // p.position().x() = 0.;
   // p.position().y() = 0.;
   // p.position().z() = 1e-5;
}
*/

template<class ParticleType>
  void Foam::Cloud<ParticleType>::deleteLostParticles()
  {
      for (ParticleType& p : *this)
      {
          if (p.cell() == -1)
          {
              WarningInFunction
                  << "deleting lost particle at position " << p.position()
                  << endl;
   
              deleteParticle(p);
          }
      }
  }

template<class ParticleType>
void Foam::Cloud<ParticleType>::cloudReset(const Cloud<ParticleType>& c)
{
    // Reset particle count and particles only
    // - not changing the cloud object registry or reference to the polyMesh
    ParticleType::particleCount_ = 0;
    IDLList<ParticleType>::operator=(c);
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::Cloud<ParticleType>::move
(
    TrackCloudType& cloud,
    typename ParticleType::trackingData& td,
    const scalar trackTime
)
{
    const polyBoundaryMesh& pbm = pMesh().boundaryMesh();
    const globalMeshData& pData = polyMesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Indexing of equivalent patch on neighbour processor into the
    // procPatches list on the neighbour
    const labelList& procPatchNeighbours = pData.processorPatchNeighbours();

    // Which processors this processor is connected to
    const labelList& neighbourProcs = pData[Pstream::myProcNo()];

    // Indexing from the processor number into the neighbourProcs list
    labelList neighbourProcIndices(Pstream::nProcs(), -1);

    forAll(neighbourProcs, i)
    {
        neighbourProcIndices[neighbourProcs[i]] = i;
    }

    // Initialise the stepFraction moved for the particles
    forAllIters(*this, pIter)
    {
        pIter().reset();
    }

    // List of lists of particles to be transferred for all of the
    // neighbour processors
    List<IDLList<ParticleType>> particleTransferLists
    (
        neighbourProcs.size()
    );

    // List of destination processorPatches indices for all of the
    // neighbour processors
    List<DynamicList<label>> patchIndexTransferLists
    (
        neighbourProcs.size()
    );

    // Allocate transfer buffers
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Clear the global positions as there are about to change
    globalPositionsPtr_.clear();

    // While there are particles to transfer
    while (true)
    {
        particleTransferLists = IDLList<ParticleType>();
        forAll(patchIndexTransferLists, i)
        {
            patchIndexTransferLists[i].clear();
        }

        // Loop over all particles
        for (ParticleType& p : *this)
        {
            // Move the particle
            bool keepParticle = p.move(cloud, td, trackTime);

            // If the particle is to be kept
            // (i.e. it hasn't passed through an inlet or outlet)
            if (keepParticle)
            {
                if (td.switchProcessor)
                {
                    #ifdef FULLDEBUG
                    if
                    (
                        !Pstream::parRun()
                     || !p.onBoundaryFace()
                     || procPatchNeighbours[p.patch()] < 0
                    )
                    {
                        FatalErrorInFunction
                            << "Switch processor flag is true when no parallel "
                            << "transfer is possible. This is a bug."
                            << exit(FatalError);
                    }
                    #endif

                    const label patchi = p.patch();

                    const label n = neighbourProcIndices
                    [
                        refCast<const processorPolyPatch>
                        (
                            pbm[patchi]
                        ).neighbProcNo()
                    ];

                    p.prepareForParallelTransfer();

                    particleTransferLists[n].append(this->remove(&p));

                    patchIndexTransferLists[n].append
                    (
                        procPatchNeighbours[patchi]
                    );
                }
            }
            else
            {
                p.position().x() = 0.;
                p.position().y() = 0.;
                p.position().z() = 0.0003;

            }
        }

        if (!Pstream::parRun())
        {
            break;
        }


        // Clear transfer buffers
        pBufs.clear();

        // Stream into send buffers
        forAll(particleTransferLists, i)
        {
            if (particleTransferLists[i].size())
            {
                UOPstream particleStream
                (
                    neighbourProcs[i],
                    pBufs
                );

                particleStream
                    << patchIndexTransferLists[i]
                    << particleTransferLists[i];
            }
        }

        // Start sending. Sets number of bytes transferred
        labelList allNTrans(Pstream::nProcs());
        pBufs.finishedSends(allNTrans);


        bool transferred = false;

        for (const label n : allNTrans)
        {
            if (n)
            {
                transferred = true;
                break;
            }
        }
        reduce(transferred, orOp<bool>());

        if (!transferred)
        {
            break;
        }

        // Retrieve from receive buffers
        for (const label neighbProci : neighbourProcs)
        {
            label nRec = allNTrans[neighbProci];

            if (nRec)
            {
                UIPstream particleStream(neighbProci, pBufs);

                labelList receivePatchIndex(particleStream);

                IDLList<ParticleType> newParticles
                (
                    particleStream,
                    typename ParticleType::iNew(polyMesh_)
                );

                label pI = 0;
                for (ParticleType& newp : newParticles)
                {
                    label patchi = procPatches[receivePatchIndex[pI++]];

                    newp.correctAfterParallelTransfer(patchi, td);

                    addParticle(newParticles.remove(&newp));
                }
            }
        }
    }
}

template<class ParticleType>
void Foam::Cloud<ParticleType>::autoMap(const mapPolyMesh& mapper)
{
    if (!globalPositionsPtr_)
    {
        FatalErrorInFunction
            << "Global positions are not available. "
            << "Cloud::storeGlobalPositions has not been called."
            << exit(FatalError);
    }

    // Reset stored data that relies on the mesh
    //    polyMesh_.clearCellTree();
    cellWallFacesPtr_.clear();

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();
    polyMesh_.oldCellCentres();

    const vectorField& positions = globalPositionsPtr_();

    label i = 0;
    forAllIters(*this, iter)
    {
        iter().autoMap(positions[i], mapper);
        ++i;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/this->name() + "_positions.obj"
    );

    forAllConstIters(*this, pIter)
    {
        const ParticleType& p = pIter();
        const point position(p.position());
        pObj<< "v " << position.x() << " " << position.y() << " "
            << position.z() << nl;
    }

    pObj.flush();
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::storeGlobalPositions() const
{
    // Store the global positions for later use by autoMap. It would be
    // preferable not to need this. If the mapPolyMesh object passed to autoMap
    // had a copy of the old mesh then the global positions could be recovered
    // within autoMap, and this pre-processing would not be necessary.

    globalPositionsPtr_.reset(new vectorField(this->size()));

    vectorField& positions = globalPositionsPtr_();

    label i = 0;
    forAllConstIters(*this, iter)
    {
        positions[i] = iter().position();
        ++i;
    }
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudIO.C"

// ************************************************************************* //
