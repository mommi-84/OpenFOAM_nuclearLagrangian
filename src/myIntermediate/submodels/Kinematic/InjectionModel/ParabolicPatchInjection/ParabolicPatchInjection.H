/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::ParabolicPatchInjection

Group
    grpLagrangianIntermediateInjectionSubModels

Description
    Patch injection.

    User specifies:
      - Total mass to inject
      - Name of patch
      - Injection duration
      - Initial parcel velocity
      - Injection volume flow rate

    Properties:
      - Parcel diameters obtained by distribution model
      - Parcels injected according to probability
        density function specified in base file across the patch

SourceFiles
    ParabolicPatchInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ParabolicPatchInjection_H
#define ParabolicPatchInjection_H

#include "InjectionModel.H"
#include "parabolicPatchInjectionBase.H"

namespace Foam
{

class distributionModel;

/*---------------------------------------------------------------------------*\
                       Class ParabolicPatchInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParabolicPatchInjection
:
    public InjectionModel<CloudType>,
    public parabolicPatchInjectionBase
{
    // Private data

        //- Name of carrier (mass or volume) flux field
        const word phiName_;

        //- Name of carrier density field
        const word rhoName_;

        //- Injection duration [s]
        scalar duration_;

        //- Concentration profile of particle volume to carrier volume [-]
        autoPtr<Function1<scalar>> concentration_;

        //- Parcels to introduce per unit volume flow rate m3 [n/m3]
        const scalar parcelConcentration_;

        //- Center point of the circular patch
        const vector centerPoint_;

        //- Radius of the circular patch
        const scalar R_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;

public:

    //- Runtime type information
    TypeName("parabolicPatchInjection");


    // Constructors

        //- Construct from dictionary
        ParabolicPatchInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParabolicPatchInjection(const ParabolicPatchInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new ParabolicPatchInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParabolicPatchInjection();


    // Member Functions

        //- Inherit updateMesh from parabolicPatchInjectionBase
        using patchInjectionBase::updateMesh;

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Return the total volumetric flow rate across the patch [m3/s]
        virtual scalar flowRate() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParabolicPatchInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
