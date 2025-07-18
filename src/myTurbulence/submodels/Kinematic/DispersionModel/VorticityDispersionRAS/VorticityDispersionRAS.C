/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "VorticityDispersionRAS.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VorticityDispersionRAS<CloudType>::VorticityDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner)
{}


template<class CloudType>
Foam::VorticityDispersionRAS<CloudType>::VorticityDispersionRAS
(
    const VorticityDispersionRAS<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VorticityDispersionRAS<CloudType>::~VorticityDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::VorticityDispersionRAS<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const scalar d,
    const scalar rho,
    const vector& Uc,
    const scalar rhoc,
    const scalar muc,
    const vector& vortc,
    vector& UTurb,
    scalar& tTurb
)
{
    Random& rnd = this->owner().rndGen();

    const scalar cps = 0.16432;

    const scalar k = this->kPtr_->primitiveField()[celli];
    const scalar epsilon =
        this->epsilonPtr_->primitiveField()[celli] + ROOTVSMALL;

    const scalar UrelMag = mag(U - Uc - UTurb);

    //const scalar Cd = 18*muc/(rho*d*d)*(1 + 1/6*pow(rho*UrelMag*d/muc,2/3));

    //const scalar tRel = 4/3*rho*d/(rhoc*UrelMag*Cd + ROOTVSMALL);

    const scalar tTurbLoc =
        min(k/epsilon, cps*pow(k, 1.5)/epsilon/(UrelMag + SMALL));

    const vector H(
    vortc.x() * Uc.x(),
    vortc.y() * Uc.y(),
    vortc.z() * Uc.z()
    );
    Info<<"tTurbLoc value: "<<tTurbLoc<<" with dt = "<<dt<<" and tTurb: "<<tTurb<<endl;


    // Parcel is perturbed by the turbulence
    if (dt < tTurbLoc)
    {
        tTurb += dt;

        if (tTurb > tTurbLoc)
        {
            tTurb = 0;

            const scalar sigma = sqrt(2*k/3.0)*2*sqrt(3.0);

            scalar Gx = 0.5*abs(H[0])/(mag(H) + ROOTVSMALL);
            scalar Gy = 0.5*abs(H[1])/(mag(H) + ROOTVSMALL);
            scalar Gz = 0.5*abs(H[2])/(mag(H) + ROOTVSMALL);

            // Calculate a random direction dir distributed uniformly
            // in spherical coordinates

            const scalar theta = rnd.sample01<scalar>()*twoPi;
            const scalar u = 2*rnd.sample01<scalar>() - 1;

            const scalar a = sqrt(1 - sqr(u));
            const vector dir(a*cos(theta)*Gx, a*sin(theta)*Gy, u*Gz);

            UTurb = sigma*mag(rnd.GaussNormal<scalar>())*dir;
        }
    }
    else
    {
        tTurb = GREAT;
        UTurb = Zero;
    }

    return Uc + UTurb;
}


// ************************************************************************* //
