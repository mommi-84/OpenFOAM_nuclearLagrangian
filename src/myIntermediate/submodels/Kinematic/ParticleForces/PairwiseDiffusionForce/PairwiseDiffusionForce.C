/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "PairwiseDiffusionForce.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"
#include "Cloud.H"

namespace Foam
{

template<class CloudType>
PairwiseDiffusionForce<CloudType>::PairwiseDiffusionForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    D_(this->coeffs().getScalar("D")),
    h_(this->coeffs().getScalar("h"))
{}

template<class CloudType>
PairwiseDiffusionForce<CloudType>::PairwiseDiffusionForce
(
    const PairwiseDiffusionForce& dff
)
:
    ParticleForce<CloudType>(dff),
    D_(dff.D_),
    h_(dff.h_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::PairwiseDiffusionForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{

    const label cellI = p.cell();
    const vector& pos = p.position();
    CloudType& cloud = const_cast<CloudType&>(this->owner());
    auto& cellOccupancy = cloud.cellOccupancy();

    const auto& mesh = cloud.mesh();

    const point& pCellCentre = mesh.C()[cellI];

    vector totalForce = vector::zero;
    scalar weightSum = 0.0;
    scalar sigma = h_/2.0;

    forAll(mesh.C(), otherCellI)
    {
        const point& otherCentre = mesh.C()[otherCellI];

        if (magSqr(pCellCentre - otherCentre) > sqr(h_))
        {
            continue; // Skip cells outside influence radius
        }

        const auto& parcelList = cellOccupancy[otherCellI];

        for (const auto* qPtr : parcelList)
        {
            // Skip self-comparison (compare pointers)
            if (qPtr == &p) continue;

            const vector dr = p.position() - qPtr->position();
            const scalar distSqr = magSqr(dr);

            if (distSqr > sqr(h_)) continue;
            if (mag(dr) < VSMALL) continue;

            const scalar w = 1/sqrt(2*pi*sqr(sigma))* exp(-0.5*sqr(dr/sigma));  // Gaussian Kernel

            totalForce += w * (dr / mag(dr));  // Normalized direction scaled by weight
            weightSum += w;
        }
    }

    vector averagedForce = (weightSum > VSMALL)
        ? D_ * totalForce / weightSum
        : vector::zero;

    return Foam::forceSuSp(averagedForce, 0.0);
}
    
    
// ************************************************************************* //
