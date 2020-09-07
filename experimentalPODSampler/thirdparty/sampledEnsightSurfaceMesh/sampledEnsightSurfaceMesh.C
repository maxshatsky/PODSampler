/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "sampledEnsightSurfaceMesh.H"
#include "meshSearch.H"
#include "Tuple2.H"
#include "globalIndex.H"
#include "treeDataCell.H"
#include "treeDataFace.H"
#include "meshTools.H"
#include "Xfer.H"
#include "surfZoneList.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(sampledEnsightSurfaceMesh, 0);
    addToRunTimeSelectionTable
    (
        sampledSurface,
        sampledEnsightSurfaceMesh,
        word
    );

    //- Private class for finding nearest
    //  Comprising:
    //  - global index
    //  - sqr(distance)
    typedef Tuple2<scalar, label> nearInfo;

    class nearestEqOp
    {
    public:

        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first() < x.first())
            {
                x = y;
            }
        }
    };
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledEnsightSurfaceMesh::sampledEnsightSurfaceMesh
(
    const word& name,
    const polyMesh& mesh,
    const word& surfaceName
)
:
    sampledSurface(name, mesh),
    MeshStorage(),
    ensName_(surfaceName),
    surface_
    (
        ensightSurfaceReader(ensName_).geometry()
    ),
    needsUpdate_(true),
    sampleElements_(0)
{
	deleteFaces(mesh);
}


Foam::sampledEnsightSurfaceMesh::sampledEnsightSurfaceMesh
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    MeshStorage(),
    ensName_(dict.get<fileName>("surfaceName")),
    surface_
    (
        ensightSurfaceReader(ensName_).geometry()
    ),
    needsUpdate_(true),
    sampleElements_(0)
{
	deleteFaces(mesh);
}


Foam::sampledEnsightSurfaceMesh::sampledEnsightSurfaceMesh
(
    const word& name,
    const polyMesh& mesh,
    const meshedSurface& surface
)
:
    sampledSurface(name, mesh),
    MeshStorage(),
    ensName_(""),
    surface_
    (
        surface
    ),
    needsUpdate_(true),
    sampleElements_(0)
{
	deleteFaces(mesh);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledEnsightSurfaceMesh::~sampledEnsightSurfaceMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::sampledEnsightSurfaceMesh::deleteFaces
(
	const polyMesh& curMesh
)
{
	//Pout<<nl<<nl<<"surface_ BEFORE = "<<nl<<surface_<<nl<<nl;
	faceList newFaces(0);
        label cellId = -1;

	forAll(surface_.surfFaces(), i)
	{
                cellId = curMesh.findCell(surface_.Cf()[i]);
                label maxCellId = cellId;
                if (Pstream::parRun)
            	    reduce(maxCellId,maxOp<scalar>());
		if ( cellId != -1 && maxCellId == cellId)
                {
                    sampleElements_.append(cellId);
		    newFaces.append(surface_.surfFaces()[i] );
                }
               
	}
	
	surface_.swapFaces(newFaces);	// since surface_ is not const anymore
	
	//forAll(surface_.Cf(), i)
	//{
	//    Pout << "face #" << i << ", center = " << surface_.Cf()[i] << "; cell #" << curMesh.cellCentres()[sampleElements_[i]] << endl;
	//}
	
	//Pout << "Cf() now: " << surface_.Cf() << endl;
	
	//Pout<<nl<<nl<<"surface_ AFTER = "<<nl<<surface_<<nl<<nl;

}


bool Foam::sampledEnsightSurfaceMesh::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledEnsightSurfaceMesh::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();
    MeshStorage::clear();

    sampleElements_.clear();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledEnsightSurfaceMesh::update()
{
    if (needsUpdate_)
    {
        needsUpdate_ = false;
        return true;
    }
    return needsUpdate_;
}







Foam::tmp<Foam::scalarField> Foam::sampledEnsightSurfaceMesh::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledEnsightSurfaceMesh::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledEnsightSurfaceMesh::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledEnsightSurfaceMesh::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledEnsightSurfaceMesh::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}





Foam::tmp<Foam::scalarField> Foam::sampledEnsightSurfaceMesh::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledEnsightSurfaceMesh::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledEnsightSurfaceMesh::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledEnsightSurfaceMesh::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledEnsightSurfaceMesh::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}



void Foam::sampledEnsightSurfaceMesh::print(Ostream& os) const
{
    os  << "sampledEnsightSurfaceMesh: " << name() << " :"
        << " faces:"   << faces().size()
        << " points:"  << points().size()
        << " zoneids:" << zoneIds().size();
}


// ************************************************************************* //
