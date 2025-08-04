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

#include "StochasticDispersionRAS.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticDispersionRAS<CloudType>::StochasticDispersionRAS
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner)
{}


template<class CloudType>
Foam::StochasticDispersionRAS<CloudType>::StochasticDispersionRAS
(
    const StochasticDispersionRAS<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StochasticDispersionRAS<CloudType>::~StochasticDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::StochasticDispersionRAS<CloudType>::update
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

    const scalar eddyLength = cps*pow(k, 1.5)/epsilon;

    const scalar tTurbLoc =
      min(eddyLength/pow(2*k/3, 0.5), eddyLength/(UrelMag + SMALL)); //Interaction time taken as the minimum between the eddy lifetime and particle crossing time
     // Parcel is perturbed by the turbulence - REMOVED, the time discretisation is made so that it is always considered particle-eddy interaction
    //if (dt < tTurbLoc)
    //{
    tTurb += dt;  //Tracks the cumulative time a particle interacts with an eddy 

    if (tTurb > tTurbLoc)  //Interaction ended, reset interaction time for new edd
    {
        tTurb = dt;
    }

    if (tTurb == dt) //Computes the direction imposed by the specific eddy, it will remain constant throughout the whole interaction
    {
        const scalar sigma = sqrt(2*k/3.0);

        // Calculate a random direction dir distributed uniformly
        // in spherical coordinates

        const scalar theta = rnd.sample01<scalar>()*twoPi;
        const scalar u = 2*rnd.sample01<scalar>() - 1;

        const scalar a = sqrt(1 - sqr(u));
        const vector dir(a*cos(theta), a*sin(theta), u);

        UTurb = sigma*mag(rnd.GaussNormal<scalar>())*dir;
    }
	//}
    return Uc + UTurb;
}


// ************************************************************************* //
