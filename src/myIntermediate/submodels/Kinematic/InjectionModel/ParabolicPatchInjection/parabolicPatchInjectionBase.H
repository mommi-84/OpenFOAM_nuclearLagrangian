/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::ParabolicPatchInjectionBase

Description
    Base class for patch-based injection models.

    Class handles injecting at a point adjacent to the patch faces corresponding 
    to a probability distribution function to provide a custom injection process.  
    Patch faces are triangulated, and area fractions accumulated.  The fractional 
    areas are then applied to determine across which face a parcel is to be injected.

SourceFiles
    ParabolicPatchInjectionBase.C

\*---------------------------------------------------------------------------*/

#ifndef parabolicPatchInjectionBase_H
#define parabolicPatchInjectionBase_H

#include "word.H"
#include "labelList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMesh;
class fvMesh;
class Random;

/*---------------------------------------------------------------------------*\
                     Class parabolicPatchInjectionBase Declaration
\*---------------------------------------------------------------------------*/

class parabolicPatchInjectionBase
{
protected:

    // Protected data

        //- Patch name
        const word patchName_;

        //- Patch ID
        const label patchId_;

        //- Patch area - total across all processors
        scalar patchArea_;

        //- Patch face normal directions
        vectorList patchNormal_;

        //- Patch face centroid
        vector patchCentroid_;

        //- Patch face center
        vector patchCenter_;

        //- List of cell labels corresponding to injector positions
        labelList cellOwners_;

        //- Decomposed patch faces as a list of triangles
        faceList triFace_;

        //- Addressing from per triangle to patch face
        labelList triToFace_;

        //- Cumulative triangle area per triangle face
        scalarList triCumulativeMagSf_;

        //- Cumulative area fractions per processor
        scalarList sumTriMagSf_;


public:

    // Constructors

        //- Construct from mesh and patch name
        parabolicPatchInjectionBase(const polyMesh& mesh, const word& patchName);

        //- Copy constructor
        parabolicPatchInjectionBase(const parabolicPatchInjectionBase& pib);


    //- Destructor
    virtual ~parabolicPatchInjectionBase() = default;


    // Member Functions

        //- Update patch geometry and derived info for injection locations
        virtual void updateMesh(const polyMesh& mesh);

        //- Return the processor that has the location specified by the fraction
        label whichProc(const scalar fraction01) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
