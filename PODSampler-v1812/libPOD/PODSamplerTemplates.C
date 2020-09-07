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

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

using Eigen::MatrixXd;

template<class Type>
void Foam::PODSampler::valuesGather
(
    Field<Type>& values
) const
{

    List<Field<Type>> gatheredValues(Pstream::nProcs());
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    gatheredValues[Pstream::myProcNo()] = values;
    
    if (!Pstream::parRun())
    {
        return;
    }
    
    if (!Pstream::master())
    {
        UOPstream toMasterStream
        (
            Pstream::masterNo(),
            pBufs
        );
        toMasterStream << values;
    }
    
    pBufs.finishedSends();

    if (Pstream::master())
    {
        //read data from slaves
        for (label iproc=Pstream::firstSlave(); iproc<=Pstream::lastSlave(); iproc++)
        {
            UIPstream procStream(iproc, pBufs);
            Field<Type> recievedValues(procStream);
            gatheredValues[iproc] = recievedValues;
        }
    }
    
    if (Pstream::master())
    {
        // Combine values into single field
        Field<Type> allValues
        (
            ListListOps::combine<Field<Type> >
            (
                gatheredValues,
                accessOp<Field<Type> >()
            )
        );
        values = allValues;
    }
}

template<class Type>
void Foam::PODSampler::appendValueToFile
(
	std::string dirName,
	std::string fileName,
	Type& value
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
    theFile << Foam::scientific << Foam::setprecision(10) << value << nl;
}
