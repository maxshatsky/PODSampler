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

Class
    Foam::sampledEnsightSurfaceMesh

Description
    A sampledSurface from a mesh with ensight reader. 

SourceFiles
    sampledEnsightSurfaceMesh.C
    sampledEnsightSurfaceMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef sampledEnsightSurfaceMesh_H
#define sampledEnsightSurfaceMesh_H

#include "sampledSurface.H"
#include "ensightSurfaceReader.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class treeDataFace;
class meshSearch;

/*---------------------------------------------------------------------------*\
                       Class sampledEnsightSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class sampledEnsightSurfaceMesh
:
    public sampledSurface,
    public meshedSurface
{

private:

    //- Private typedefs for convenience
        typedef meshedSurface MeshStorage;


    // Private data

        //- For compatibility with the meshSurf interface
        labelList zoneIds_;

        //- Filename for ensight reader
        fileName ensName_;
        
        //- Ensight reader for mesh
        //ensightSurfaceReader ensReader_;

        //- Surface to sample on
        meshedSurface surface_; //ensight surface

        //- Track if the surface needs an update
        mutable bool needsUpdate_;

        //- From local surface triangle to mesh cell/face.
        labelList sampleElements_;

    // Private Member Functions

        //- Sample field on faces
        template<class Type>
        tmp<Field<Type>> sampleField
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        ) const;

        //- Sample volume field onto surface faces
        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;
        
        template<class Type>
        tmp<Field<Type>>
        interpolateField(const interpolation<Type>&) const;


	// MS edition
        void deleteFaces
        (
          const polyMesh& curMesh
        );

public:

    //- Runtime type information
    TypeName("sampledEnsightSurfaceMesh");


    // Constructors

        //- Construct from components
        sampledEnsightSurfaceMesh
        (
            const word& name,
            const polyMesh& mesh,
            const word& surfaceName
        );

        //- Construct from dictionary
        sampledEnsightSurfaceMesh
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from triSurface
        sampledEnsightSurfaceMesh
        (
            const word& name,
            const polyMesh& mesh,
            const meshedSurface& surface
        );


    //- Destructor
    virtual ~sampledEnsightSurfaceMesh();


    // Member Functions


        //- Does the surface need an update?
        virtual bool needsUpdate() const;

        //- Mark the surface as needing an update.
        //  May also free up unneeded data.
        //  Return false if surface was already marked as expired.
        virtual bool expire();

        //- Update the surface as required.
        //  Do nothing (and return false) if no update was needed
        virtual bool update();

        //- Points of surface
        virtual const pointField& points() const
        {
            return surface_.points();
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return surface_.surfFaces();
        }

        //- Face area vectors
        virtual const vectorField& Sf() const
        {
            return surface_.Sf();
        }

        //- Face area magnitudes
        virtual const scalarField& magSf() const
        {
            return surface_.magSf();
        }

        //- Face centres
        virtual const vectorField& Cf() const
        {
            return surface_.Cf();
        }

        //- Const access to per-face zone/region information
        virtual const labelList& zoneIds() const
        {
            return zoneIds_;
        }


        
        //- Interpolate field on surface
        virtual tmp<scalarField> sample
        (
            const interpolation<scalar>&
        ) const;


        //- Interpolate field on surface
        virtual tmp<vectorField> sample
        (
            const interpolation<vector>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<sphericalTensorField> sample
        (
            const interpolation<sphericalTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<symmTensorField> sample
        (
            const interpolation<symmTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<tensorField> sample
        (
            const interpolation<tensor>&
        ) const;





        //- Interpolate field on surface
        virtual tmp<scalarField> interpolate
        (
            const interpolation<scalar>&
        ) const;


        //- Interpolate field on surface
        virtual tmp<vectorField> interpolate
        (
            const interpolation<vector>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const interpolation<sphericalTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<symmTensorField> interpolate
        (
            const interpolation<symmTensor>&
        ) const;

        //- Interpolate field on surface
        virtual tmp<tensorField> interpolate
        (
            const interpolation<tensor>&
        ) const;

        //- Write information
        virtual void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledEnsightSurfaceMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
