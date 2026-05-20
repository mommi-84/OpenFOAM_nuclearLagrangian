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

\*---------------------------------------------------------------------------*/

#include "MorsiAlexanderDragForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::MorsiAlexanderDragForce<CloudType>::CdRe(const scalar Re) const
{
    scalar C1=24;
    scalar C2=0;
    scalar C3=0;
    scalar Cd=0;
    if (Re > 0.1)
    {
        C1 = 22.73;
	C2 = 0.09;
	C3 = 3.69;
    }
    if (Re > 1)
    {
        C1 = 29.17;
	C2 = -3.89;
	C3 = 1.22;
    }
    if (Re > 10)
    {
        C1 = 46.5;
	C2 = -116.67;
	C3 = 0.62;
    }
    if (Re > 100)
    {
        C1 = 98.33;
	C2 = -2778;
	C3 = 0.36;
    }
    if (Re > 1000)
    {
        C1 = 148.62;
	C2 = -47500;
	C3 = 0.36;
    }
    if (Re > 5000)
    {
        C1 = -490.546;
	C2 = 578800;
	C3 = 0.46;
    }
    if (Re > 10000)
    {
        C1 = -1662.5;
	C2 = 5416700;
	C3 = 0.52;
    }

    return (C1/(Re + SMALL) + C2/(pow(Re, 2.0) + SMALL) + C3)*Re;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MorsiAlexanderDragForce<CloudType>::MorsiAlexanderDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::MorsiAlexanderDragForce<CloudType>::MorsiAlexanderDragForce
(
    const MorsiAlexanderDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::MorsiAlexanderDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    // (AOB:Eq. 34)
    return forceSuSp (Zero, mass * 0.75 * muc * CdRe(Re) / (p.rho() * sqr(p.d())));

}
// ************************************************************************* //
