/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Tommaso Pernatsch, Politecnico di Milano
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

#include "parabolicTurbulentInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parabolicTurbulentInletFvPatchVectorField::
parabolicTurbulentInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(0),
    n_(1, 0, 0),
    y_(0, 1, 0)
{
}


Foam::parabolicTurbulentInletFvPatchVectorField::
parabolicTurbulentInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(dict.get<scalar>("maxValue")),
    n_(dict.get<vector>("n")),
    y_(dict.get<vector>("y"))
{

    Info << "Using the parabolicVelocity boundary condition" << endl;
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("parabolicVelocityFvPatchVectorField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }
    n_ /= mag(n_);
    y_ /= mag(y_);
}


Foam::parabolicTurbulentInletFvPatchVectorField::
parabolicTurbulentInletFvPatchVectorField
(
    const parabolicTurbulentInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_)
{}


Foam::parabolicTurbulentInletFvPatchVectorField::
parabolicTurbulentInletFvPatchVectorField
(
    const parabolicTurbulentInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_)
{}


Foam::parabolicTurbulentInletFvPatchVectorField::
parabolicTurbulentInletFvPatchVectorField
(
    const parabolicTurbulentInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::parabolicTurbulentInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vector z_ = n_ ^ y_;

    // Use patch face centers to calculate the bounding box
    const vectorField& faceCentroids = patch().Cf();

    // Initialize min and max points for the bounding box
    point bbMin(GREAT, GREAT, GREAT);
    point bbMax(-GREAT, -GREAT, -GREAT);

    // Loop through the patch face centers to find the bounding box
    forAll(faceCentroids, i)
    {
        const point& pt = faceCentroids[i];
        bbMin = Foam::min(bbMin, pt);
        bbMax = Foam::max(bbMax, pt);
    }

    // Synchronize bounding box across processors for parallel cases
    reduce(bbMin, minOp<point>());
    reduce(bbMax, maxOp<point>());

    // Check for zero extents and handle gracefully
    vector bbExtent = bbMax - bbMin;
    if (bbExtent.x() == 0 || bbExtent.y() == 0 || bbExtent.z() == 0)
    {
        FatalErrorInFunction
            << "Bounding box has zero extent in one or more directions. Ensure patch geometry is valid."
            << nl << "Bounding box: min = " << bbMin << ", max = " << bbMax
            << exit(FatalError);
    }

    // Compute the center of the bounding box
    vector ctr = 0.5 * (bbMax + bbMin);

    // Calculate local 2D coordinates for the parabolic profile
    scalarField coordY = 2 * ((faceCentroids - ctr) & y_) / bbExtent.y();
    scalarField coordZ = 2 * ((faceCentroids - ctr) & z_) / bbExtent.z();

    // Calculate radial distance from the center
    scalarField radialCoord = sqrt(sqr(coordY) + sqr(coordZ));

    // Compute the parabolic velocity profile
    vectorField::operator=(n_ * maxValue_ * (1.0 - sqr(radialCoord)));

    // Apply the fixedValue boundary condition
    fixedValueFvPatchVectorField::updateCoeffs();
}



void Foam::parabolicTurbulentInletFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("maxValue") << maxValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y") << y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        parabolicTurbulentInletFvPatchVectorField
    );
}

// ************************************************************************* //
