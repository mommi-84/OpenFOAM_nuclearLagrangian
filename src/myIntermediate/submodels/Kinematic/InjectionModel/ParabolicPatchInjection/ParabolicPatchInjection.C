/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "ParabolicPatchInjection.H"
#include "distributionModel.H"
#include "polyMesh.H"
#include "SubField.H"
#include "Random.H"
#include "triangle.H"
#include "volFields.H"
#include "polyMeshTetDecomposition.H"
#include "axisAngleRotation.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParabolicPatchInjection<CloudType>::ParabolicPatchInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    parabolicPatchInjectionBase(owner.mesh(), this->coeffDict().getWord("patch")),
    phiName_(this->coeffDict().template getOrDefault<word>("phi", "phi")),
    rhoName_(this->coeffDict().template getOrDefault<word>("rho", "rho")),
    duration_(this->coeffDict().getScalar("duration")),
    concentration_
    (
        Function1<scalar>::New
        (
            "concentration",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    parcelConcentration_
    (
        this->coeffDict().getScalar("parcelConcentration")
    ),
    centerPoint_
    (
       this->coeffDict().template get<vector>("centerPoint")
    ),
    R_
    (
       this->coeffDict().template getScalar("R")
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    ),
{
    // Convert from user time to reduce the number of time conversion calls
    const Time& time = owner.db().time();
    duration_ = time.userTimeToTime(duration_);
    concentration_->userTimeToTime(time);

    patchInjectionBase::updateMesh(owner.mesh());

    // Re-initialise total mass/volume to inject to zero
    // - will be reset during each injection
    this->volumeTotal_ = 0.0;
    this->massTotal_ = 0.0;
}


template<class CloudType>
Foam::ParabolicPatchInjection<CloudType>::ParabolicPatchInjection
(
    const ParabolicPatchInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    parabolicPatchInjectionBase(im),
    duration_(im.duration_),
    phiName_(im.phiName_),
    rhoName_(im.rhoName_),
    centerPoint_(im.centerPoint_),
    R_(im.R_),
    concentration_(im.concentration_.clone()),
    parcelConcentration_(im.parcelConcentration_),
    sizeDistribution_(im.sizeDistribution_.clone()),
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParabolicPatchInjection<CloudType>::~ParabolicPatchInjection()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParabolicPatchInjection<CloudType>::updateMesh()
{
    parabolicPatchInjectionBase::updateMesh(this->owner().mesh());
}


template<class CloudType>
Foam::scalar Foam::ParabolicPatchInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}

template<class CloudType>
Foam::scalar Foam::ParabolicPatchInjection<CloudType>::flowRate() const
{
   const polyMesh& mesh = this->owner().mesh();

    const auto& phi = mesh.lookupObject<surfaceScalarField>(phiName_);

    const scalarField& phip = phi.boundaryField()[patchId_];

    scalar flowRateIn = 0.0;
    if (phi.dimensions() == dimVolume/dimTime)
    {
        flowRateIn = max(0.0, -sum(phip));
    }
    else
    {
        const auto& rho = mesh.lookupObject<volScalarField>(rhoName_);
        const scalarField& rhop = rho.boundaryField()[patchId_];

        flowRateIn = max(0.0, -sum(phip/rhop));
    }

    reduce(flowRateIn, sumOp<scalar>());

    return flowRateIn;
}

template<class CloudType>
Foam::label Foam::ParabolicPatchInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar dt = time1 - time0;

        scalar c = concentration_->value(0.5*(time0 + time1));

        scalar nParcels = parcelConcentration_*c*flowRate()*dt;

        Random& rnd = this->owner().rndGen();

        label nParcelsToInject = floor(nParcels);

        // Inject an additional parcel with a probability based on the
        // remainder after the floor function
        if
        (
            nParcelsToInject > 0
         && (
               nParcels - scalar(nParcelsToInject)
             > rnd.globalPosition(scalar(0), scalar(1))
            )
        )
        {
            ++nParcelsToInject;
        }

        return nParcelsToInject;
    }

    return 0;
}


template<class CloudType>
Foam::scalar Foam::ParabolicPatchInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar volume = 0.0;

    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar c = concentration_->value(0.5*(time0 + time1));

        volume = c*(time1 - time0)*flowRate();
    }

    this->volumeTotal_ = volume;
    this->massTotal_ = volume*this->owner().constProps().rho0();

    return volume;
}


template<class CloudType>
void Foam::ParabolicPatchInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    const polyMesh& mesh = this->owner().mesh();
    Random& rnd = this->owner().rndGen();
    /*const label patchId = mesh.boundaryMesh().findPatchID(this->coeffDict().getWord("patch"))
    const polyPatch& patch = mesh.boundaryMesh()[patchId];*/

    // Assume the patch is circular with a known radius R
    const scalar R = R_;

    // Uniform angular distribution
    scalar theta = rnd.sample01<scalar>() * 2 * Foam::constant::mathematical::pi;

    // Radial distribution based on the parabolic profile 1/(pi*R^2)*(1 - r^2/R^2)
    scalar u = rnd.sample01<scalar>();

    // Invert the CDF to get the radial position
    scalar r = sqrt(u) * R; // Better initial guess
    scalar maxIter = 1000;  // Maximum number of iterations
    scalar tolerance = 1e-6; // Tolerance for convergence

    for (int iter = 0; iter < maxIter; ++iter)
    {
        // Function and derivative for Newton method
        scalar f = r * r / (R * R * Foam::constant::mathematical::pi) * (1 - (r * r / (2 * R * R))) - u;
        scalar df = 2 * r / (R * R * Foam::constant::mathematical::pi) * (1.0 - (r * r) / (R * R));

        // Handle zero derivative
        if (mag(df) < SMALL) 
        {
            r += SMALL; // Perturb slightly
            continue;
        }

        scalar delta = -f / df;
        r += delta;

        // Clamp r to valid range
        r = max(0.0, min(r, R));

        // Check for convergence
        if (mag(delta) < tolerance)
        {
            break;
        }

        // If the loop ends without convergence
        if (iter == maxIter - 1)
        {
            Warning << "Newton-Raphson did not converge for u = " << u
                    << " after " << maxIter << " iterations." << endl;
        }
    }

    // Set the arbitrary XY position
    vector position(0.0, r * Foam::cos(theta), r * Foam::sin(theta));

    this->findCellAtPosition
    (
        cellOwner,
        tetFacei,
        tetPti,
        position,
        false
    );

    /*Info << "Calculated position: " << position << endl;
    Info << "Assigned cellOwner: " << cellOwner << endl;
    Info << "Assigned tetFacei: " << tetFacei << endl;
    Info << "Assigned tetPti: " << tetPti << endl;
    Info << "Assigned cell center: " << mesh.cellCentres()[cellOwner] << endl;*/

}


template<class CloudType>
void Foam::ParabolicPatchInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // Set particle velocity to carrier velocity
    parcel.U() = this->owner().U()[parcel.cell()];

    // Set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::ParabolicPatchInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::ParabolicPatchInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
