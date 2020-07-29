/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/


#include "PODSampler.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PODSampler, 0);
        
    addToRunTimeSelectionTable
    (
        functionObject,
        PODSampler,
        dictionary
    );
}


Foam::scalar Foam::PODSampler::mergeTol_ = 1e-10;

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //


void Foam::PODSampler::writeGeometry() //const
{
    
    fileName outputDir = mesh_.time().path();

    if (Pstream::parRun())
        outputDir += "/..";
        
    outputDir += "/surfaceGeometryData";
    
    Info << "write geometry \n";

    Info << controlSurfaces_.size() << nl;


    forAll(controlSurfaces_, surfI)
    {
        Info << "in surfaces list\n";        
        
        const sampledSurface& s = controlSurfaces_.operator[](surfI);
        

        if (Pstream::parRun())
        {
            Info << "surfI size = " << mergeList_[surfI].points.size() << nl;

            if (Pstream::master() && mergeList_[surfI].faces.size())
            {
                /*formatter_->write
                (
                    outputDir,
                    s.name(),
                    meshedSurfRef
                    (
                        mergeList_[surfI].points,
                        mergeList_[surfI].faces
                    )
                );*/
                
                formatter_.open
                (
                    mergeList_[surfI].points,
                    mergeList_[surfI].faces,
                    outputDir,
                    true
                );
                formatter_.write();
            }
        }
        else if (s.points().size())
        {
            Info << s.points().size() << nl;

            /*formatter_->write
            (
                outputDir,
                s.name(),
                meshedSurfRef
                (
                    s.points(),
                    s.faces()
                )
            );*/
            
            formatter_.open
            (
                s.points(),
                s.faces(),
                outputDir,
                true
            );
            formatter_.write();
        }
        else
        {
            Info << "empty surface info" << nl;
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::PODSampler::PODSampler
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    mesh_
    (
        refCast<const fvMesh>
        (
            t.lookupObject<objectRegistry>
            (
                dict.lookupOrDefault("region", polyMesh::defaultRegion)
            )
        )                                                                
    ),
    name_(name),
    interpolationScheme_(word::null),
    controlSurfaces_(0),
    dirName_(word::null),
    fieldName_(word::null),
    mergeList_(),
    nSnapshots_(1),
    nmodes_(1),
    PODStartTime_(0),
    endSnapshotTime_(0),
    previousValuesTime_(0),
    PODdeltaT_(0),
    averageValue_(),
    currentValues_(),
    previousValues_(),
    values_(),
    cov_(),
    timeList_(0),
    stopAfterPOD_(false),
    writeInitValues_(false),
    writeFreq_(1),
    nSavedSnapshots_(0),
    writeRestoredValues_(false),
    continueWriting_(false)

{
    Info << "in constructor" << nl;
    read(dict);
}




Foam::PODSampler::PODSampler
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObject(name),
    mesh_(refCast<const fvMesh>(obr)),
    name_(name),
    interpolationScheme_(word::null),
    controlSurfaces_(0),
    dirName_(word::null),
    fieldName_(word::null),
    mergeList_(),
    nSnapshots_(1),
    nmodes_(1),
    PODStartTime_(0),
    endSnapshotTime_(0),
    previousValuesTime_(0),
    PODdeltaT_(0),
    averageValue_(),
    currentValues_(),
    previousValues_(),
    values_(),
    cov_(),
    timeList_(0),
    stopAfterPOD_(false),
    writeInitValues_(false),
    writeFreq_(1),
    nSavedSnapshots_(0),
    writeRestoredValues_(false),
    continueWriting_(false)
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::PODSampler::~PODSampler()
{}

// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

bool Foam::PODSampler::read(const dictionary& dict)
{
    dict.lookup("interpolationScheme") >> interpolationScheme_;
    dict.lookup("fieldName") >> fieldName_;

    const word writeType (dict.lookup("outputGeometryFormat"));

    // Define the surface formatter
    // Optionally defined extra controls for the output formats
    //formatter_ = surfaceWriter::New
    //formatter_ = Foam::surfaceWriters::foamWriter::New
    formatter_.New
    (
        writeType,
        dict.subOrEmptyDict("formatOptions").subOrEmptyDict(writeType)
    );

    //read surfaces for sampling
    PtrList<sampledSurface> newList
    (
        dict.lookup("surfaces"),
        sampledSurface::iNew(mesh_)
    );

    controlSurfaces_.transfer(newList);

    //Parallel fix as it was implemented in sampledSurfaces class
    if (Pstream::parRun())
    {
        mergeList_.setSize(controlSurfaces_.size());
    }
        
    // Ensure all surfaces and merge information are expired
    expire();

    if (controlSurfaces_.size())
    {
        Info<< "Function object "<< name_<<":" << nl;

        Info<< " Reading control surface description:" << nl;
            
        forAll(controlSurfaces_, surfI)
        {
            Info<< " " << controlSurfaces_.operator[](surfI).name() << nl;
            Info<< "needsUpdate: " << controlSurfaces_.operator[](surfI).needsUpdate() << nl;
            Info<< "size: " << controlSurfaces_.operator[](surfI).points().size() << nl;
        }

        Info<< endl;

        update();

        writeGeometry();

        Info << "writeGeometry OK" << nl;
    }

    if (Pstream::master() && debug)
    {
        Pout<< "PODSampler control surfaces additional info:" << nl << "(" << nl;
            
        forAll(controlSurfaces_, surfI)
        {
            Pout<< " " << controlSurfaces_.operator[](surfI) << endl;
        }
            
        Pout<< ")" << endl;
    }

    // Here begins reading for POD operations!

    dict.lookup("nSnapshots") >> nSnapshots_;
    dict.lookup("PODdeltaT") >> PODdeltaT_;
    dict.lookup("nmodes") >> nmodes_;
    stopAfterPOD_ = dict.lookupOrDefault("stopAfterPOD", true);
    writeInitValues_ = dict.lookupOrDefault("writeInitValues", false);
    writeRestoredValues_ = dict.lookupOrDefault("writeRestoredValues", false);
    continueWriting_ = dict.lookupOrDefault("continueWriting", false);
    writeFreq_ = dict.lookupOrDefault("writeFreq", 1);

    //Setting default adress for files to save.
    if (Pstream::parRun() )
    {
        dirName_ = mesh_.time().rootPath() + "/" + mesh_.time().caseName().path() + "/POD/";
    }
    else
    {
        dirName_ = mesh_.time().rootPath() + "/" + mesh_.time().caseName() + "/POD/";
    }

    checkSettingsPt1();

    if(!continueWriting_)
    {
        dict.lookup("PODStartTime") >> PODStartTime_;
        endSnapshotTime_ = PODStartTime_ + PODdeltaT_;
        nSavedSnapshots_ = 0;

        checkSettingsPt2();
    }

    return true;
}

void Foam::PODSampler::checkSettingsPt1()
{
    if(writeFreq_ < 1)
    {
        stopTheCase();

        FatalErrorInFunction
            << nl << "In field "<< fieldName_
            << nl << "writeFreq_ < 1"
            << nl << "It leads to crash and mean nothing."
            << nl << "Set writeFreq_ >= 1."
            << nl << exit(FatalError);
    }

    if ( nSnapshots_ < nmodes_ )
    {
        stopTheCase();

        FatalErrorInFunction
            << nl << "In field "<< fieldName_
            << nl << "nSnapshots is less than nmodes: "
            << nSnapshots_ << " < " << nmodes_
            << nl << "It leads to crash at the POD calculation."
            << nl << "Set nmodes <= nSnapshots."
            << nl << exit(FatalError);
    }


    if ( nSnapshots_ <= 1 )
    {
        stopTheCase();

        FatalErrorInFunction
            << nl << "In field "<< fieldName_
            << nl << "nSnapshots <= 1, that causes error at FFTW."
            << nl << "Set nSnapshots as any integer >= 2."
            << nl << exit(FatalError);
    }

    if ( PODdeltaT_ <= 0 )
    {
        stopTheCase();

        FatalErrorInFunction
            << nl << "In field "<< fieldName_
            << nl << "PODdeltaT <= 0, that means something wrong."
            << nl << "Set PODdeltaT as any scalar > 0. And > maxDeltaT."
            << nl << exit(FatalError);
    }

    if ( nmodes_ <= 0 )
    {
        stopTheCase();

        FatalErrorInFunction
            << nl << "In field "<< fieldName_
            << nl << "nmodes < 1, that causes errors at the POD calculation."
            << nl << "Set nmodes as any integer >= 1, but less than nSnapshots."
            << nl << exit(FatalError);
    }
}

bool Foam::PODSampler::execute()
{
    return true;
}

bool Foam::PODSampler::expire()
{
    bool justExpired = false;
    
    forAll(controlSurfaces_, surfI)
    {
        if (controlSurfaces_.operator[](surfI).expire())
        {
            justExpired = true;
        }
        //Clear merge information
        if (Pstream::parRun())
        {
            mergeList_[surfI].clear();
        }
    }
    
    // true if any surfaces just expired
    return justExpired;
}

bool Foam::PODSampler::needsUpdate() const
{
    forAll(controlSurfaces_, surfI)
    {
        if (controlSurfaces_.operator[](surfI).needsUpdate())
        {
            return true;
        }
    }

    return false;
}

bool Foam::PODSampler::update()
{
    //Actually the update() function copy-pasted from libSampling
    bool updated = false;
    if (!needsUpdate())
    {
        return updated;
    }
    
    // Serial: quick and easy, no merging required
    // Just like sampledSurfaces
    if (!Pstream::parRun())
    {
        forAll(controlSurfaces_, surfI)
        {
            if (controlSurfaces_.operator[](surfI).update())
            {
                updated = true;
            }
        }
	
        return updated;
    }

    //const fvMesh& mesh = refCast<const fvMesh>(obr_);

    // Dimension as fraction of mesh bounding box
    scalar mergeDim = mergeTol_ * mesh_.bounds().mag();

    if (Pstream::master() && debug)
    {
        Pout<< nl << "Merging all points within "
        << mergeDim << " metre" << endl;
    }
    
    forAll(controlSurfaces_, surfI)
    {
        sampledSurface& s = controlSurfaces_.operator[](surfI);
        
        if (s.update())
        {
            updated = true;
        }
        else
        {
            continue;
        }

        PatchTools::gatherAndMerge
        (
            mergeDim,
            primitivePatch
            (
                SubList<face>(s.faces(), s.faces().size()),
                s.points()
            ),
            mergeList_[surfI].points,
            mergeList_[surfI].faces,
            mergeList_[surfI].pointsMap
        );
    }

    return updated;
}

void Foam::PODSampler::timeSet()
{
// Do nothing - only valid on write
}

Foam::scalar Foam::PODSampler::mergeTol()
{
    return mergeTol_;
}

Foam::scalar Foam::PODSampler::mergeTol(const scalar tol)
{
    scalar oldTol = mergeTol_;
    mergeTol_ = tol;
    return oldTol;
}

/*==================================================================*\
|                                                                    |
|                                                                    |
|                            Main part                               |
|                                                                    |
|                                                                    |
\*==================================================================*/

bool Foam::PODSampler::write()
{
    update();

    //if(values_.size() == 15) stopTheCase();

    if(continueWriting_)
    {
        const volScalarField& tmp = mesh_.lookupObject<volScalarField>(fieldName_);

        forAll(controlSurfaces_, surfI)
        {
            reviveTheCase(tmp , surfI);
        }

        checkSettingsPt2();

        continueWriting_ = false;

        if ( values_.size() == nSnapshots_)
        {
            if(Pstream::master())
            {
                calcPOD(); //this function also writes results

                // Saving coordinates
                const sampledSurface& surface = controlSurfaces_.operator[](0);
                writeCoordinates(dirName_, surface, 0);

                if(stopAfterPOD_)
                {
                    Info
                        << nl << "In field "<< fieldName_
                        << nl << nl << "Everything is ok, the POD was calculated successfully,"
                        << nl << "case stops because you set the stopAfterPOD == true."<<nl<<nl;

                    stopTheCase();
                }
            }
            else
            {
                if(stopAfterPOD_)
                {
                    stopTheCase();
                }
            }
        }
    }

    if (values_.size() < nSnapshots_ )
    {
        const volScalarField& tmp = mesh_.lookupObject<volScalarField>(fieldName_);

        forAll(controlSurfaces_, surfI)
        {
            if(tmp.time().timeOutputValue() > PODStartTime_)
            {
                sampleAndCalc(tmp, surfI);
            }
            else
            {
                refreshPreviousValues(tmp, surfI);
            }
        }
    }
    
    return true;  
}

void Foam::PODSampler::checkSettingsPt2()
{
    if ( PODStartTime_ < mesh_.time().startTime().value() )
    {
        double tmp = mesh_.time().startTime().value();

        stopTheCase();

        FatalErrorInFunction
            << nl << "In field "<< fieldName_
            << nl<< "PODStartTime is less than startTime: "
            << PODStartTime_ << " < " << tmp
            << nl << "The program can't make snapshots out of case time range."
            << nl << "Set PODStartTime as any scalar >= startTime."
            << nl << exit(FatalError);
    }

    if
    (
         (nSnapshots_-values_.size()) * PODdeltaT_
         >
         mesh_.time().endTime().value() - PODStartTime_
    )
    {
        double tmp = mesh_.time().endTime().value();

        stopTheCase();

        FatalErrorInFunction
            << nl << "In field "<< fieldName_
            << nl << "(nSnapshots_-values_.size())*PODdeltaT > mesh_.time().endTime().value() - PODStartTime_: "
            << (nSnapshots_-values_.size()) * PODdeltaT_ << " > " << tmp - PODStartTime_
            << nl << "The program is set to calculate POD only after gaining nSnapshots."
            << nl << "So it won't output any result at this conditions because having not enough values case will stop at the endTime."
            << nl << "Set nSnapshots*PODdeltaT smaller than (endTime - PODStartTime)."
            << nl << exit(FatalError);
    }
}

void Foam::PODSampler::stopTheCase()
{
    Info<<nl<<nl<<"Field: "<< fieldName_ << ":  "
        << "Setting end time for now and finishing the case."<<nl<<nl;

    const Time* curtime = &mesh_.time();
    Time* endtime = const_cast<Time*> (curtime);
    endtime->writeAndEnd();
}

void Foam::PODSampler::refreshPreviousValues
(
    const GeometricField<scalar, fvPatchField, volMesh>& vField,
    const label surfI
)
{
    Field<scalar> values;
    const sampledSurface& surface = controlSurfaces_.operator[](surfI);
    getValues(vField, surface, values);
    previousValues_ = values;
    previousValuesTime_ = vField.time().timeOutputValue();
}


void Foam::PODSampler::sampleAndCalc
(
    const GeometricField<scalar, fvPatchField, volMesh>& vField,
    const label surfI
) 
{
    Field<scalar> values;
    const sampledSurface& surface = controlSurfaces_.operator[](surfI);

    getValues(vField, surface, values); //ok

    if // i.e. it is the very first call of the function
    (
        vField.time().timeOutputValue() - PODStartTime_ <= 1.1*vField.time().deltaTValue()
    )
    {
        initialize(values, vField.time().timeOutputValue());
    }

    if (!values.empty())
    {
        gainCurrentValues(vField, values);
    }

    refreshPreviousValues(values, vField.time().timeOutputValue());

    // Make a new snapshot.
    if
    ( vField.time().timeOutputValue() >= endSnapshotTime_ )
    {
        makeSnapshot(vField, values);

        // Main course. If this snapshot was the last.
        if (  ( values_.size() == nSnapshots_ ) && cov_.empty() )
        {
            prepareAndGatherData();

            if ( Pstream::master() )
            {
                calcPOD(); //this function also writes results

                // Saving coordinates
                writeCoordinates(dirName_, surface, surfI);

                Info << nl << "All files in field "<< fieldName_ <<" were saved."<< nl
                <<"POD for field "<< fieldName_<<" is finished"<< nl;
            }

            if (stopAfterPOD_)
            {
                stopTheCase();
            }
        }
    }
}

/*==================================================================*\
|                                                                    |
|                                                                    |
|                          Gaining Values                            |
|                                                                    |
|                                                                    |
\*==================================================================*/

void Foam::PODSampler::getValues
(
    const GeometricField<scalar, fvPatchField, volMesh>& vField,
    const sampledSurface& surface,
    Field<scalar>& values
) 
{
    autoPtr<interpolation<scalar> > interpolatorPtr;

    if (surface.interpolate())
    {
        if (interpolatorPtr.empty())
        {
            interpolatorPtr = interpolation<scalar>::New
            (
                interpolationScheme_,
                vField
            );
        }
        values = surface.interpolate(interpolatorPtr());
    }
    else
    {
        if (interpolatorPtr.empty())
        {
            interpolatorPtr = interpolation<scalar>::New
            (
                "cell",
                vField
            );
        }		
        values = surface.sample(interpolatorPtr());
    }
}

void Foam::PODSampler::initialize
(
    const Field<scalar>& values,
    scalar time
)
{
    currentValues_ = Field<scalar>(values.size(), 0);

    // Here we need to set previousValues_ exactly at PODStartTime_
    if
    (
        previousValues_.empty() != values.empty()
        ||
        previousValuesTime_ == time
    )
    {
        previousValues_ = values;		// Honestly, we may don't know the values at PODStartTime in this case,
        previousValuesTime_ = PODStartTime_;	// so this is a little lie.
    }
    else
    {
        if
        ( previousValuesTime_ != PODStartTime_ )
        {
            previousValues_ += (values - previousValues_)*(PODStartTime_ - previousValuesTime_)/(time - previousValuesTime_);
            previousValuesTime_ = PODStartTime_;
        }
    }
}

void Foam::PODSampler::gainCurrentValues
(
    const GeometricField<scalar, fvPatchField, volMesh>& vField,
    const Field<scalar>& values
) 
{
    if ( vField.time().timeOutputValue() <= endSnapshotTime_ ) // i.e. it is inside the current snapshot
    {
        currentValues_ += (previousValues_ + values)*(vField.time().timeOutputValue()-previousValuesTime_)/(2*PODdeltaT_);
    }
    else // i.e. current time is above the end of current snapshot. So we need to linearly interpolate values at the very endSnapshotTime_
    {
        Field<scalar> tmp = previousValues_ + (values - previousValues_)*(endSnapshotTime_ - previousValuesTime_)/vField.time().deltaTValue();
        currentValues_ += (previousValues_ + tmp)*(endSnapshotTime_ - previousValuesTime_)/(2*PODdeltaT_);
    }
}

void Foam::PODSampler::refreshPreviousValues
(
    const Field<scalar>& values,
    scalar time
)
{
    if(time <= endSnapshotTime_)
    {
        previousValues_ = values;
        previousValuesTime_ = time;
    }
    else
    {
        previousValues_ += (values - previousValues_)*(endSnapshotTime_ - previousValuesTime_)/(time - previousValuesTime_);
        previousValuesTime_ = endSnapshotTime_;
    }
}

void Foam::PODSampler::makeSnapshot
(
    const GeometricField<scalar, fvPatchField, volMesh>& vField,
    Field<scalar>& values
)
{
    values_.append(currentValues_);
    currentValues_ = Field<scalar>(values.size(), 0);
    endSnapshotTime_ += PODdeltaT_;

    // Add new time value to list of snapshot times.
    if (Pstream::master())
    {
        timeList_.conservativeResize(timeList_.rows()+1);
        timeList_(timeList_.rows()-1)=endSnapshotTime_ - 3*PODdeltaT_/2;
    }

    Info<<"Snapshot No "<<values_.size()<<" was made."<<nl;

    // write initial values before POD
    if
    (
        writeInitValues_
        &&
        (
            (values_.size() % writeFreq_ == 0)
            ||
            (values_.size() == nSnapshots_)
        )
    )
    {
        writeNewValues(values_.size() - nSavedSnapshots_);
        nSavedSnapshots_ = values_.size();
    }

    // If current time is inside the next snapshot, we must take it into account
    if ( vField.time().timeOutputValue() > endSnapshotTime_ - PODdeltaT_ )
    {
        if (!values.empty())
        {
            gainCurrentValues(vField, values);
        }

        refreshPreviousValues(values, vField.time().timeOutputValue()); 
    }
}

void Foam::PODSampler::writeNewValues
(
    label n
)
{
    if(Pstream::parRun()) // gathering values_ if parRun()
    {
       for(int i=0; i<n; i++)
       {
           valuesGather(values_[values_.size()-n+i]);
       }
    }

    if(Pstream::master()) // writing values_ and timeList to file
    {
        for(int i=0; i<n; i++)
        {
            appendFieldToFile(dirName_+fieldName_+"/", fieldName_+"_initValues.txt", values_[values_.size()-n+i] );
            appendValueToFile(dirName_+fieldName_+"/", "timeList.txt", timeList_(timeList_.rows()+i-n)); //appendValueToFile is a template function
        }
        Info<<"initValues from "<<values_.size()-n+1<<" to "<< values_.size() <<" were written ro file."<<nl;
    }
}

void Foam::PODSampler::appendFieldToFile
(
    std::string dirName,
    std::string fileName,
    const Field<scalar>& values
) const
{
    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fileName;

    Foam::OFstream theFile
    (
        dirName,
        Foam::IOstream::ASCII,
        Foam::IOstream::currentVersion,
        Foam::IOstream::UNCOMPRESSED,
        true
    );
    theFile << Foam::scientific << Foam::setprecision(10);

    for(int i=0; i<values.size()-1; i++)
    {
        theFile << values[i] << "\t";
    }
    theFile << values[values.size()-1] << nl;
}


/*==================================================================*\
|                                                                    |
|                                                                    |
|                          Preparing Data                            |
|                                                                    |
|                                                                    |
\*==================================================================*/


void Foam::PODSampler::prepareAndGatherData() 
{
    Info << nl <<"Field: "<< fieldName_ << ":  " << "POD begins.";

    if(Pstream::parRun() && !writeInitValues_)
    {
        forAll(values_, i)
        {
            valuesGather(values_[i]);
        }
    }
}

void Foam::PODSampler::covCalc
(
    DynamicList<Field<scalar>>& values,
    Field<Field<scalar>>& cov
)
{
    scalar s = 0.0; // "s" is sum
    cov.setSize(values.size());

    //It has no parallel implementation yet
    forAll ( cov , i )
    {
        cov[i].setSize(i+1); // <-- Memory econonomy used!
        
        forAll (cov[i], j)
        {
            forAll (values[i], k)
            {
                s += values[i][k]*values[j][k];
            }
            cov[i][j] = s;
            s = 0.0;
        }
    }
}

void Foam::PODSampler::covGather
(
	Field<Field<scalar> >& locCov
) const
{
    if (!Pstream::parRun())
    {
        return;
    }
    List<Field<Field<scalar>>> gatheredValues(Pstream::nProcs());
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    gatheredValues[Pstream::myProcNo()] = locCov;

    //Pstream::gatherList(gatheredValues);

    if (!Pstream::master())
    {
        UOPstream toMasterStream
        (
            Pstream::masterNo(),
            pBufs
        );
        forAll(locCov,i)
        {
            toMasterStream.write
            (
                reinterpret_cast<const char*>(locCov[i].cdata()),
                sizeof(scalar)*locCov[i].size()
            );
        }
    }
    
    pBufs.finishedSends(true);
    
    if (Pstream::master())
    {
        //read data from slaves
        for (label iproc=Pstream::firstSlave(); iproc<=Pstream::lastSlave(); iproc++)
        {
            UIPstream procStream(iproc, pBufs);
            gatheredValues[iproc].resize(locCov.size()); //allocate memory for gatherd values
            forAll(locCov,j)
            {
                gatheredValues[iproc][j].resize(locCov[j].size());
                procStream.read
                (
                    reinterpret_cast<char*>(gatheredValues[iproc][j].data()),
                    sizeof(scalar)*gatheredValues[iproc][j].size()
                );
            }
        }
    }
    
    if (Pstream::master())
    {
        for(label iproc=Pstream::firstSlave(); iproc<=Pstream::lastSlave(); iproc++)
        {
            locCov += gatheredValues[iproc];
        }
    }
}


/*==================================================================*\
|                                                                    |
|                                                                    |
|                         Calculating POD                            |
|                                                                    |
|                                                                    |
\*==================================================================*/


void Foam::PODSampler::calcPOD() 
{
    Info << nl <<"Field: "<< fieldName_ << ":  " << "Averaging of values.";

    averageValue_ = Field<scalar>(values_[0].size(), 0);

    forAll(values_, i)
    {
        averageValue_ += values_[i]/nSnapshots_;
    }

    forAll(values_, i)
    {
        values_[i] -= averageValue_;
    }

    covCalc(values_, cov_);

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Copying averagevalues to eigen format.";

    // EigenAver
    Eigen::VectorXd eigenAver= Eigen::VectorXd::Zero(values_[0].size());
    for (int i=0; i<eigenAver.rows(); i++)
    {
        eigenAver(i) = averageValue_[i];
    }

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Copying values to eigen format.";

    // Copying values matrix into eigen format
    Eigen::MatrixXd valuesMatrix;
    valuesMatrix.resize(values_[1].size(),values_.size());
    for (int i=0; i<values_[1].size();i++)
    {
        for (int j =0; j<values_.size(); j++)
        {
            //values were averaged in the beginning of final processing
            valuesMatrix(i, j)=values_[j][i];
        }
    }

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Copying cov matrix to eigen format.";

    // Making covariance matrix in eigen format
    Eigen::MatrixXd covMatrix= MatrixXd::Ones(3,3);
    covMatrix.resize(cov_.size(), cov_.size());

    // Copying covariance matrix into eigen format
    for (int i=0; i<cov_.size();i++)
    {
        for (int j =0; j<cov_[i].size(); j++)
        {
            covMatrix(i, j)=cov_[i][j];
        }
    }

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Solving eigen problem of cov matrix.";

    Eigen::SelfAdjointEigenSolver<MatrixXd> es(covMatrix);

    //eigenvalues and eigenvectors calculating and saving
    Eigen::VectorXd eivals = covMatrix.selfadjointView<Eigen::Lower>().eigenvalues();
    Eigen::MatrixXd eivecs = es.eigenvectors();

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Pairing and sorting eigenvalues and eigenvectors.";

    //pairing eivecs and eivals for SORTING
    std::vector<std::pair<double, Eigen::Matrix<double, -1, -1>>>  eigenpairs;
    for(int i=0; i< cov_.size();i++)
    {
        eigenpairs.push_back(std::make_pair(eivals[i] , eivecs.col(i)));
    }

    //sorting
    std::sort
    (
        eigenpairs.begin(),
        eigenpairs.end(),
        comparePairs
    );

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Collecting top valued eigenvectors into single matrix.";

    // Collecting top valued vectors into single matrix!
    Eigen::MatrixXd rawvecs;
    rawvecs.resize( cov_.size(), 1);
    rawvecs <<   eigenpairs[0].second;
    for(int i=1; i<nmodes_; i++)
    {
        rawvecs.conservativeResize(rawvecs.rows(), rawvecs.cols()+1);
        rawvecs.col(rawvecs.cols()-1) = eigenpairs[i].second;
    }

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Making basis matrix.";

    // Making basis matrix - matrix of basis vectors.
    Eigen::MatrixXd basis;
    basis.resize(values_[0].size(), nmodes_);
    basis=valuesMatrix*rawvecs;

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Making Phi matrix.";

    // Creating an intermediate matrix Phi.
    // Later it will be described in some paper.
    Eigen::MatrixXd phi=basis;
    for(int i=0; i<basis.cols();i++)
    {
        phi.col(i) = basis.col(i)/( basis.col(i).transpose() * basis.col(i) + VSMALL );
    }

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Making acoeffs matrix.";

    // Calculating acoeffs
    Eigen::MatrixXd acoeffs;
    acoeffs=phi.transpose()*valuesMatrix;

    // Desctruction of Phi matrix
    phi.resize(0, 0);

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Calculating fftw of acoeffs.";

    // Calculating fftw of acoeffs.
    List<autoPtr<Pair<List<scalar>>>> fftcoeffs (acoeffs.rows());
    autoPtr<Pair<List<scalar>>> curFFTCoeff;

    for (int i = 0; i < acoeffs.rows(); i++)
    {
        if ( (acoeffs.cols() > 0) )
        {
            List<scalar> acoeffsCur(acoeffs.cols());

            for (int j = 0; j < acoeffs.cols(); ++j)
            {
                acoeffsCur[j] = acoeffs(i,j);
            }

            FoamFftwDriver fftw (acoeffsCur, timeList_(timeList_.rows()-1) - timeList_(0));
            curFFTCoeff = fftw.simpleScalarForwardTransform();
            fftcoeffs[i] = curFFTCoeff;
        }
    }

    // Saving results to files
    Info << nl <<"Field: "<< fieldName_ << ":  " << "Saving results into files.";

    // Saving timeList
    if (!(isFile(dirName_ + "timeList.txt")))
    {
        writeVectorToFile(dirName_, "timeList.txt", timeList_);
    }

    //Saving eigenvalues of covariance matrix
    writeEigenvalues(dirName_, eigenpairs);

    // Saving acoeffs
    writeAcoeffs(dirName_, acoeffs);

    // Saving basis. It includes average values as a zero frequency basis.
    writeBasis (dirName_, eigenAver, basis);

    // Saving averageValues
    writeVectorToFile(dirName_+fieldName_+"/", fieldName_+"_averageValues.txt", eigenAver);

    // Saving fftw of acoeffs to file
    writeFFTCoeffs(dirName_, fftcoeffs);

    // Saving frequencies
    writeFrequencies(dirName_, fftcoeffs);

    if (writeRestoredValues_)
    {
        writeRestoredValues(dirName_, eigenAver, basis, acoeffs);
    }
}

bool Foam::PODSampler::comparePairs(std::pair<double, Eigen::Matrix<double, -1, -1>> a, std::pair<double, Eigen::Matrix<double, -1, -1>> b)
{
    return a.first > b.first;
}



/*==================================================================*\
|                                                                    |
|                                                                    |
|                          Saving Results                            |
|                                                                    |
|                                                                    |
\*==================================================================*/

void Foam::PODSampler::writeVectorToFile
(
    std::string dirName,
    std::string fileName,
    Eigen::VectorXd& vector
) const
{
    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fileName.c_str();

    Info << nl <<"Field: "<< fieldName_ << ":  Saving "<< fileName << " to file.";

    Foam::OFstream theFile(dirName);
    theFile << Foam::scientific << Foam::setprecision(10);

    for(int i=0; i<vector.rows(); i++)
    {
        theFile << vector(i) << nl;
    }
}

void Foam::PODSampler::writeEigenvalues
(
    std::string dirName,
    std::vector<std::pair<double, Eigen::Matrix<double, -1, -1>>>& eigenpairs
) const
{
    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "/";

    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "_eigenvalues.txt";

    Info << nl <<"Field: "<< fieldName_ << ":  Saving eigenvalues.";

    Foam::OFstream theFile(dirName);
    theFile << Foam::scientific << Foam::setprecision(10);

    for(unsigned int i=0; i<eigenpairs.size(); i++)
    {
        theFile << eigenpairs[i].first << nl;
    }
}

void Foam::PODSampler::writeAcoeffs
(
    std::string dirName,
    Eigen::MatrixXd& acoeffs
) const
{
    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "/";

    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "_acoeffs.txt";

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Saving acoeffs.";

    Foam::OFstream theFile(dirName);
    theFile << Foam::scientific << Foam::setprecision(10);

    for(int i=0; i<acoeffs.rows(); i++)
    {
        for (int j=0; j<acoeffs.cols()-1; j++)
        {
            theFile<< acoeffs(i, j) << "\t";
        }
        theFile << acoeffs(i, acoeffs.cols()-1) << nl;
    }
}

void Foam::PODSampler::writeBasis
(
    std::string dirName,
    Eigen::VectorXd& eigenAver,
    Eigen::MatrixXd& basis
) const
{
    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "/";

    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "_basis.txt";

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Saving basis.";

    Foam::OFstream theFile(dirName);
    theFile << Foam::scientific << Foam::setprecision(10);

    for(int i=0; i<basis.rows(); i++)
    {
        theFile << eigenAver(i) << "\t";
        for (int j=0; j<basis.cols()-1; j++)
        {
            theFile << basis(i, j) << "\t";
        }
        theFile << basis(i, basis.cols()-1) << nl;
    }
}

void Foam::PODSampler::writeFFTCoeffs
(
    std::string dirName,
    List<autoPtr<Pair<List<scalar>>>>& fftcoeffs
) const
{
    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "/";

    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "_fftw.txt";

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Saving fftw to acoeffs.";

    Foam::OFstream theFile(dirName);
    theFile << Foam::scientific << Foam::setprecision(10);

    for (int i = 0; i < fftcoeffs[0]().first().size(); ++i)
    {
        if (i==0)
            theFile << "1.0000000000e+00\t";
        else
            theFile << "0.0000000000e+00\t";

        for (int j=0; j<fftcoeffs.size()-1; j++)
            theFile <<  fftcoeffs[j]().second()[i] << "\t";

        theFile << fftcoeffs[fftcoeffs.size()-1]().second()[i] << nl;
    }
}

void Foam::PODSampler::writeFrequencies
(
    std::string dirName,
    List<autoPtr<Pair<List<scalar>>>>& fftcoeffs
) const
{
    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += "frequencies.txt";

    if (!(isFile(dirName)))
    {

        Info << nl <<"Field: "<< fieldName_ << ":  " << "Saving frequencies for all fields.";

        Foam::OFstream theFile(dirName);
        theFile << Foam::scientific << Foam::setprecision(10);

        for (int i = 0; i < fftcoeffs[0]().first().size(); ++i)
        {
            theFile << fftcoeffs[0]().first()[i] << nl;
        }
    }
}

void Foam::PODSampler::writeRestoredValues
(
    std::string dirName,
    Eigen::VectorXd& eigenAver,
    Eigen::MatrixXd& basis,
    Eigen::MatrixXd& acoeffs
) const
{
    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "/";

    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += fieldName_.c_str();
    dirName += "_restoredValues.txt";

    Info << nl <<"Field: "<< fieldName_ << ":  " << "Saving restoredValues.";

    Foam::OFstream theFile(dirName);
    theFile << Foam::scientific << Foam::setprecision(10);

    for(int i=0; i<basis.rows(); i++)
    {
        for (int j=0; j<acoeffs.cols()-1; j++)
        {
            scalar tmp= eigenAver(i);

            for (int k=0; k<nmodes_; k++)
                tmp+=acoeffs(k, j)*basis(i, k);

             theFile << tmp << "\t";
        }

        scalar tmp= eigenAver(i);

        for (int k=0; k<nmodes_; k++)
            tmp+=acoeffs(k, acoeffs.cols()-1)*basis(i, k);

        theFile << tmp << nl;
    }
}


void Foam::PODSampler::writeCoordinates
(
    std::string dirName,
    const sampledSurface& surface,
    const label surfI
) const
{

    if (!(isDir(dirName)))
    {
        mkDir(dirName);
    }

    dirName += "facecentres.txt";
    if (!(isFile(dirName)))
    {
        Info << nl <<"Field: "<< fieldName_ << ":  " << "Saving facecentres for all fields.";

        Foam::OFstream theFile(dirName);
        theFile << Foam::scientific << Foam::setprecision(10);

        if (Pstream::parRun())
        {
            for (int i=0; i<mergeList_[surfI].faces.size(); i++)
            {
                scalar x=0, y=0, z=0;
                for (int j=0 ;j<mergeList_[surfI].faces[i].size(); j++)
                {
                    x+=mergeList_[surfI].points[ mergeList_[surfI].faces[i][j] ].x();
                    y+=mergeList_[surfI].points[ mergeList_[surfI].faces[i][j] ].y();
                    z+=mergeList_[surfI].points[ mergeList_[surfI].faces[i][j] ].z();
                }
                x/=mergeList_[surfI].faces[i].size();
                y/=mergeList_[surfI].faces[i].size();
                z/=mergeList_[surfI].faces[i].size();
                theFile << x <<"\t"<<y<<"\t"<<z<<nl;
            }
        }
        else
        {
            for (int i=0; i<surface.faces().size(); i++)
            {
                scalar x=0, y=0, z=0;
                for (int j=0 ;j<surface.faces()[i].size(); j++)
                {
                    x+=surface.points()[ surface.faces()[i][j] ].x();
                    y+=surface.points()[ surface.faces()[i][j] ].y();
                    z+=surface.points()[ surface.faces()[i][j] ].z();
                }
                x/=surface.faces()[i].size();
                y/=surface.faces()[i].size();
                z/=surface.faces()[i].size();
                theFile << x <<"\t"<<y<<"\t"<<z<<nl;
            }
        }
    }
}



/*==================================================================*\
|                                                                    |
|                                                                    |
|                      Reviving Crashed Case                         |
|                                                                    |
|                                                                    |
\*==================================================================*/

bool Foam::PODSampler::reviveTheCase
(
    const GeometricField<scalar, fvPatchField, volMesh>& vField,
    const label surfI
)
{
    Field<scalar> values;
    bool stopNow = false;
    const sampledSurface& surface = controlSurfaces_.operator[](surfI);
    getValues(vField, surface, values);
    currentValues_ = Field<scalar>(values.size(), 0);

    if (Pstream::parRun())
    {
        valuesGather(values);
    }

    if (Pstream::master())
    {
        stopNow = !readAndCheckData(values);

        if(!stopNow)
        {
            PODStartTime_ = timeList_(timeList_.rows()-1) + PODdeltaT_/2;
            endSnapshotTime_ = PODStartTime_ + PODdeltaT_;
            nSavedSnapshots_ = timeList_.rows();
        }
    }
    else
    {
        PODStartTime_ = 0;
        endSnapshotTime_ = 0;
        nSavedSnapshots_ = 0;
    }

    reduce(stopNow, sumOp<bool>());
    reduce(PODStartTime_, sumOp<scalar>());
    reduce(endSnapshotTime_, sumOp<scalar>());
    reduce(nSavedSnapshots_, sumOp<label>());

    // This will stop the case at all the processors
    if (stopNow)
    {
        stopTheCase();

        FatalErrorInFunction
            << nl << "In field "<< fieldName_<<": reviving crashed and the case has been stopped, check log for more info."
            << nl << exit(FatalError);

        return false;
    }

    if ( PODStartTime_ < mesh_.time().timeOutputValue() - mesh_.time().deltaTValue() )
    {
        double tmp = mesh_.time().timeOutputValue() - mesh_.time().deltaTValue();

        stopTheCase();

        FatalErrorInFunction
            << nl << "In field "<< fieldName_
            << nl<< "Case revivng crashed because last saved snapshot ends earlier than you restarted the case: "
            << PODStartTime_ << " < " << tmp
            << nl << "The program can't make snapshots out of case time range."
            << nl << "Set case's startTime earlier to let the POD continue writing."
            << nl << exit(FatalError);

        return false;
    }

    if(!Pstream::master())
    {
        values_.setSize(nSavedSnapshots_);
    }

    Info<<nl<<"Field "<<fieldName_<<": The case was revived successfully and POD will contine writing from "<<PODStartTime_<<nl;

    return true;
}

bool Foam::PODSampler::readAndCheckData
(
    const Field<scalar>& values
)
{
    Info
        << nl <<"Field "<<fieldName_<<": Reading timeList_. If the case freezed here, it means that there is no timeList.txt file in POD/"<<fieldName_<<" folder."
        << nl <<"Maybe you have forgotten to set continueWriting = false."<<nl<<nl;
    readVectorFromFile(dirName_+fieldName_+"/timeList.txt", timeList_); //ok

    if( timeList_.rows() == 0 )
    {
        Info
            << nl << "In field "<< fieldName_
            << nl << "timeList_.rows() == 0 after reading from " << dirName_+fieldName_+"/timeList.txt"
            << nl << "Check your entry data, pls..."
            << nl << "Or refuse to use continueWriting and restart the whole case. Good luck!"
            << nl;

        return false;
    }

    if (timeList_(timeList_.rows()-1) + PODdeltaT_/2 < mesh_.time().timeOutputValue())
    {
        Info
            << nl << "In field "<< fieldName_
            << nl << "last timeList_(i) + PODdeltaT_ < mesh_.time().timeOutputValue():"
            << nl << timeList_(timeList_.rows()-1) + PODdeltaT_ << " < "<< mesh_.time().timeOutputValue()
            << nl << "It leads to unsurance of further results, because the function object looses values between this two times."
            << nl << "You should restart the case from earlier time."
            << nl << "Or refuse to use continueWriting and restart the whole case. Good luck!"
            << nl;

        return false;
    }

    Info
        << nl <<"Field "<<fieldName_<<": Reading values_. If the case freezed here, it means that there is no "<<fieldName_<<"_initValues.txt file in POD/"<<fieldName_<<" folder."
        << nl <<"Maybe you have forgotten to set continueWriting = false."<<nl<<nl;
    if(!( readValuesFromFile(dirName_+fieldName_+"/"+fieldName_+"_initValues.txt") ) )
    {
        return false;
    }

    if(values_.size() != timeList_.rows())
    {
        Info
            << nl << "In field "<< fieldName_
            << nl << "values_.size() != timeList_.rows():"
            << nl << values_.size() << " != " << timeList_.rows()
            << nl << "It means that initValues and timeList files has different number of strings!"
            << nl << "Check your entry data, pls..."
            << nl << "Or refuse to use continueWriting and restart the whole case. Good luck!"
            << nl;

        return false;
    }

    forAll(values_ ,i)
    {
        if (values_[i].size() != values.size())
        {
            Info
                << nl << "In field "<< fieldName_
                << nl << "values_["<<i<<"].size() != values.size() that was interpolated from surface:"
                << nl << values_[i].size() << " != " << values.size()
                << nl << "It means that your old case was set some other way"
                << nl << "and now you have different number of faces at the surface."
                << nl << "It is impossible to continue the case with such a data because this mistake can mean anything."
                << nl << "Maybe you should JUST delete last strings from initValues and timeList!"
                << nl << "Check your entry data, pls. Or refuse to use continueWriting and restart the whole case. Good luck!"
                << nl;

            return false;
        }
    }

    if (values_.size() > nSnapshots_)
    {
        Info
            << nl << "In field "<< fieldName_
            << nl << "You already have "<<values_.size()<<" snapshots, that is more than nSnapshots ("<<nSnapshots_<<")."
            << nl << "If you want to just calculate POD out of this data,"
            << nl << "just set nSnapshots = "<< values_.size()<<", and POD will do the calculations."
            << nl << nl << "Otherwise program doesn't understand, what do you mean when setting"
            << nl << "nSnapshots less than the number of gained snapshots."
            << nl << nl << "Set nSnapshots >= [number of gained snapshots]"
            << nl << nl;

        return false;
    }

    return true;
}

bool Foam::PODSampler::readVectorFromFile
(
    std::string fileName,
    Eigen::VectorXd& vector
)
{
    IFstream theFile(fileName);
    vector.resize(0);
    std::string s;
    while(!theFile.eof())
    {
        vector.conservativeResize(vector.rows()+1);
        theFile.getLine(s);
        std::stringstream ss;
        ss.str(s);

        if(!(ss >> vector(vector.rows()-1)))
        {
            vector.conservativeResize(vector.rows()-1);
        }
    }

    return true;
}



bool Foam::PODSampler::readValuesFromFile
(
    std::string fileName
)
{
    IFstream theFile(fileName);

    std::string s;
    std::stringstream ss1;
    double x; //buffer
    int i=0;

    if(theFile.eof())
    {

        Info<<nl<<nl<<"In field "<<fieldName_<<" initValues file seems to be empty. There's EOF at the very beginning."<<nl<<nl;
        return false;
    }

    theFile.getLine(s);

    values_.setSize(1);
    values_[0].setSize(0);
    ss1.str(s);
    while(ss1>>x) values_[0].append(x); //setting first snapshot

    while(!theFile.eof())
    {
        theFile.getLine (s);
        std::stringstream ss;
        ss.str(s);

        if (!( s.empty() ))
        {
            values_.append(Field<scalar>(values_[0].size())); //just resizing
            i=0;

            while (ss>>x)
            {

                if ( i > values_[0].size() )
                {
                    Info<<nl<<nl<<"In field "<<fieldName_<<"Size of old snapshots vary, check input data for initValues lines length sameness."
                        <<nl<<"Error was on line "<<values_.size()<<nl<<nl;
                    return false;
                }
                values_[values_.size()-1][i] = x;
                i++;
            }

            if ( i != values_[0].size() )
            {
                Info<<nl<<nl<<"In field "<<fieldName_<<" Size of old snapshots vary, check input data for initValues lines length sameness."
                    <<nl<<"Error was on line "<<values_.size()<<nl<<nl;
                return false;
            }

        }
        else
        {
            Info<<nl<<nl<<"Line "<<values_.size()+1<<" seems to be empty. If it was the last line, then it's OK! =)";
        }

    }
    return true;
}


bool Foam::PODSampler::writeValuesToFile
(
    std::string fileName,
    const DynamicList<Field<scalar>>& values,
    label n
)
{
    OFstream theFile(fileName);

    theFile<<Foam::scientific<<Foam::setprecision(10);
    for(label i=0; i<n; i++)
    {
        for(int j=0; j< values[i].size()-1; j++)
        {
            theFile<<values[i][j]<<"\t";
        }
        theFile<<values[i][values[i].size()-1]<<nl;
    }

    return true;
}

bool Foam::PODSampler::writeAllValuesToFile
(
    std::string fileName
)
{
    OFstream theFile(fileName);

    theFile<<Foam::scientific<<Foam::setprecision(10);
    forAll(values_, i)
    {
        for(int j=0; j< values_[i].size()-1; j++)
        {
            theFile<<values_[i][j]<<"\t";
        }
        theFile<<values_[i][values_[i].size()-1]<<nl;
    }

    return true;
}

// ************************************************************************* 
