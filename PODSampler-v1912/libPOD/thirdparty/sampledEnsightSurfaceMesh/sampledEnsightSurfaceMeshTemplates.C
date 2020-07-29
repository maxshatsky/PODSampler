/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledEnsightSurfaceMesh::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // One value per face
    tmp<Field<Type>> tvalues(new Field<Type>(sampleElements_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(sampleElements_, elemI)
    {
        values[elemI] = vField[sampleElements_[elemI]];
    }

    return tvalues;
}



template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledEnsightSurfaceMesh::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    const labelList& elements = sampleElements_;

    return sampledSurface::sampleOnFaces
    (
        sampler,
        elements,
        faces(),
        points()
    );
}



template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledEnsightSurfaceMesh::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    tmp<Field<Type>> tvalues(new Field<Type>(sampleElements_.size()));;
    
    FatalErrorInFunction
        << "Interpolation in ensightFoamReader is forbidden" << nl
        << exit(FatalError);

    return tvalues;
}


// ************************************************************************* //
