/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::PODSampler

Description

    Return sound amplitudes of sound pressure in each point of STL-surface.
    Data are sorted by frequencies.

    Input data:
    -   path to triSurface;
    -   reference pressure.

    Output:
    -   file with pressure history.

SourceFiles
    PODSampler.C

\*---------------------------------------------------------------------------*/

#ifndef PODSampler_H
#define PODSampler_H

#include "IOmanip.H"
#include "OFstream.H"
#include "IFstream.H"
#include "PatchTools.H"
#include "surfaceWriter.H"
#include "foamFFTWDriver.H"
#include "sampledSurface.H"
#include <Eigen/Eigenvalues>
#include "addToRunTimeSelectionTable.H"
//#include "StringStream.H"
//#include "ISstream.H"


//trash?
/*
#include "Time.H"
#include <iomanip>
#include <fstream>
#include "Tuple2.H"
#include "Switch.H"
#include "IOmanip.H"
#include "argList.H"
#include "HashSet.H"

#include "polyMesh.H"
#include "OFstream.H"
#include "Eigen/Dense"
#include "volFields.H"
#include "PatchTools.H"
#include "dictionary.H"
#include "Foam2Eigen.H"
#include "wordReList.H"
#include "mapPolyMesh.H"
#include "ListListOps.H"
#include "volFieldsFwd.H"
#include "IOobjectList.H"
#include "pointFieldFwd.H"
#include "meshedSurfRef.H"
#include "surfaceFields.H"
#include "stringListOps.H"
#include "functionObject.H"
#include <Eigen/Eigenvalues>
#include "surfaceFieldsFwd.H"
#include "coordinateSystem.H"
#include "primitiveFieldsFwd.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                    Structure surfaceData declaration
\*---------------------------------------------------------------------------*/
/*
struct surfaceData
{
    DynamicList<Field<scalar>> values_;

    Field<Field<scalar>> cov_;

    Field<scalar> previousValues_;

    Field<scalar> currentValues_;

    Field<scalar> averageValue_;
}*/


namespace Foam
{

// Forward declaration of classes
class Time;
class fvMesh;
class dictionary;

struct surfaceData
{
    DynamicList<Field<scalar>> values_;

    Field<Field<scalar>> cov_;

    Field<scalar> previousValues_;

    Field<scalar> currentValues_;

    Field<scalar> averageValue_;
};

/*---------------------------------------------------------------------------*\
                    Class PODSampler Declaration
\*---------------------------------------------------------------------------*/

class PODSampler
:
    public functionObject
{


private:

//- Class used for surface merging information
    class mergeInfo
    {

    public:

	pointField points;
	faceList faces;
	labelList pointsMap;

	//- Clear all storage
	void clear()
	{
	    points.clear();
	    faces.clear();
    	    pointsMap.clear();
	}
    };
    
    //- Output verbosity
    static bool verbose_;
    
    //- Tolerance for merging points (fraction of mesh bounding box)
    static scalar mergeTol_;
    
protected:

    // Private data
    
	//- Const reference to fvMesh
	const fvMesh& mesh_;
        
        //- Name of this set of PODSampler,
        // Also used as the name of the probes directory.
        word name_;

        //- Interpolation scheme to use for sampling
        word interpolationScheme_;

        //- List of control surfaces
        PtrList<sampledSurface> controlSurfaces_;

        //- Writer for surface geometry
        autoPtr<surfaceWriter> formatter_;

	//- glogal string to keep directory name for any output
	std::string dirName_;

        //- Name of pressure field
        word fieldName_;

        //- Information for merging surfaces
        List<mergeInfo> mergeList_;

	//- number of snapshots to POD
	label nSnapshots_;

	//- number of covMatrix eigenvectors to be used in POD
	int nmodes_;

	//- Snapshoting starts after this time
	scalar PODStartTime_;

        //- time of the end of the current snapshot interval
        scalar endSnapshotTime_;

        //- time of the previous interpolated values
        scalar previousValuesTime_;

        //- time between snapshots
        scalar PODdeltaT_;

	//- average values vector
	//Field<scalar> averageValue_;

	//- current values to append later
	//Field<scalar> currentValues_;

        //- previous Values at the previsious time
        //Field<scalar> previousValues_;

	//- main data matrix
	//DynamicList<Field<scalar>> values_;

	//- covariance matrix
	//Field<Field<scalar>> cov_;

	//- list of time values where values are written
	Eigen::VectorXd timeList_;

	//- true if user wants to stop soving after POD
	bool stopAfterPOD_;

	//- true if user wants to write and check values before POD
	//- WARNING! It may be big data.
	bool writeInitValues_;

        //- flag to check reviving.
        //bool revived_;

	//- frequency of writing of initValues
	//- default 1
	label writeFreq_;

        //- number of snapshots being saved to file. It is important if you use continueWriting_.
        label nSavedSnapshots_;

	//- true if user wants to write and check values restored after POD
	//- WARNING! It may be big data.
	bool writeRestoredValues_;

        //- True if user wants to start writing from some time value
        //- also it is a flag to check if the case already was revived
        bool continueWriting_;

        List<surfaceData> dataBank_;

private:    

        //- Disallow default bitwise copy construct
        PODSampler(const PODSampler&);

        //- Disallow default bitwise assignment
        void operator=(const PODSampler&);


        //- Write geometry only
        void writeGeometry() const;


	void appendFieldToFile
	(
		std::string dirName,
		std::string fileName,
		const Field<scalar>& values
	) const;


	template<class Type>
	void appendValueToFile
	(
		std::string dirName,
		std::string fileName,
		Type& value
	) const;


	template<class Type>
	void valuesGather
	(
    		Field<Type>& values
	) const;


	void covCalc
	(
		DynamicList<Field<scalar>>& values_,
		Field<Field<scalar>>& cov_
	) ;


	void covGather
	(
		Field<Field<scalar>>& values
	) const;


        void writeVectorToFile
        (
            std::string dirName,
            std::string fileName,
            Eigen::VectorXd& vector
        ) const;


	void writeAcoeffs
	(
		std::string dirName,
		Eigen::MatrixXd& acoeffs
	) const;


	void writeBasis
	(
		std::string dirName,
		Eigen::VectorXd& eigenAver,
		Eigen::MatrixXd& basis
	) const;

        void writeEigenvalues
        (
            std::string dirName,
            std::vector<std::pair<double, Eigen::Matrix<double, -1, -1>>>& eigenpairs
        ) const;

	void writeFFTCoeffs
	(
		std::string dirName,
		List<autoPtr<Pair<List<scalar>>>>& fftcoeffs
	) const;


	void writeFrequencies
	(
		std::string dirName,
		List<autoPtr<Pair<List<scalar>>>>& fftcoeffs
	) const;


	void writeCoordinates
	(
		std::string fileName,
		const sampledSurface& surface,
		const label surfI
	) const;


	void writeRestoredValues
	(
		std::string dirName,
		Eigen::VectorXd& eigenAver,
		Eigen::MatrixXd& basis,
		Eigen::MatrixXd& acoeffs
	) const;


	void writeNewValues
	(
	    label n
	) ;

	void stopTheCase() ;

        void refreshPreviousValues
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vField
        ) ;

        void refreshPreviousValues
        (
            const List<Field<scalar>>& values,
            scalar time
        ) ;

        void initialize
        (
            const List<Field<scalar>>& values,
            scalar time
        ) ;

        void gainCurrentValues
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vField,
            const List<Field<scalar>>& values
        ) ;

        void getValues
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vField,
            const sampledSurface& surface,
            Field<scalar>& values
        ) ;

        void makeSnapshot
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vField,
            List<Field<scalar>>& values
        ) ;

	void prepareAndGatherData();


        void calcPOD
        (
            label surfI
        );


        void sampleAndCalc
        (
		const GeometricField<scalar, fvPatchField, volMesh>& vField
        ) ;

        bool readVectorFromFile
        (
            std::string fileName,
            Eigen::VectorXd& vector
        ) ;

        bool readValuesFromFile
        (
            label surfI
        ) ;

        bool writeValuesToFile
        (
            std::string fileName,
            const DynamicList<Field<scalar>>& values,
            label n
        ) ;

        bool writeAllValuesToFile
        (
            std::string fileName
        ) ;

        static bool comparePairs
        (
            std::pair< double, Eigen::Matrix<double, -1, -1> > a,
            std::pair< double, Eigen::Matrix<double, -1, -1> > b
        ) ;

        bool reviveTheCase
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vField
        ) ;

        bool readAndCheckData
        (
            const List<Field<scalar>>& values
        ) ;

        void checkSettingsPt1();

        void checkSettingsPt2();
        

public:

    //- Runtime type information
    TypeName("PODSampler");

    // Constructors
    
    //- Construct from Time and dictionary
    PODSampler
    (
        const word& name,
        const Time& time,
        const dictionary& dict
    );
        
    //- Construct for given objectRegistry and dictionary
    //  allow the possibility to load fields from files
    PODSampler
    (
        const word& name,
        const objectRegistry&,
        const dictionary&,
        const bool loadFromFiles = false
    );

    //- Destructor
    virtual ~PODSampler();

    // Member Functions

        //- Return name of the set of PODSampler
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the PODSampler data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        //virtual bool end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the PODSampler
        virtual bool write();

    // Copy-pasted functions from libSampling for correct working of surface sampling

        //- Mark the control surfaces as needing an update.
        // The function is copy-pasted from sampledSurfaces.C
        // May also free up unneeded data.
        // Return false if all surfaces were already marked as expired.
        virtual bool expire();

        virtual bool needsUpdate() const;

        virtual bool update();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {};

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {};

        //- Get merge tolerance
        static scalar mergeTol();

        //- Set tolerance (and return old tolerance)
        static scalar mergeTol(const scalar);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PODSamplerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
