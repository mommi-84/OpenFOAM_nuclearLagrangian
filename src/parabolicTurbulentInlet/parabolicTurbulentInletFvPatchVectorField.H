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

Class
    Foam::parabolicTurbulentInletFvPatchVectorField

Group
    grpGenericBoundaryConditions

 Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y.

Usage
    \table
        Property   | Description                  | Req'd | Default
        scalarData | single scalar value          | yes |
        data       | single vector value          | yes |
        fieldData  | vector field across patch    | yes |
        timeVsData | vector function of time      | yes |
        wordData   | word, eg name of data object | no  | wordDefault
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       parabolicTurbulentInlet;
        scalarData -1;
        data       (1 0 0);
        fieldData  uniform (3 0 0);
        timeVsData table (
                             (0 (0 0 0))
                             (1 (2 0 0))
                          );
        wordName   anotherName;
        value      uniform (4 0 0); // optional initial value
    }
    \endverbatim

SourceFiles
    parabolicTurbulentInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef parabolicTurbulentInletFvPatchVectorField_H
#define parabolicTurbulentInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class parabolicTurbulentInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class parabolicTurbulentInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Peak velocity magnitude
        scalar maxValue_;
        //- Flow direction
        vector n_;
        //- Direction of the y-coordinate
        vector y_;


    // Private Member Functions

        //- Return peak velocity magnitude
        scalar maxValue() const;

        //- Return flow direction   
        vector n() const;

        //- Return y direction
        vector y() const;


public:

    //- Runtime type information
    TypeName("parabolicTurbulentInlet");


    // Constructors

        //- Construct from patch and internal field
        parabolicTurbulentInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicTurbulentInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        parabolicTurbulentInletFvPatchVectorField
        (
            const parabolicTurbulentInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        parabolicTurbulentInletFvPatchVectorField
        (
            const parabolicTurbulentInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicTurbulentInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicTurbulentInletFvPatchVectorField
        (
            const parabolicTurbulentInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicTurbulentInletFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Return max value
            scalar& maxValue()
            {
                return maxValue_;
            }
            //- Return flow direction
            vector& n()
            {
                return n_;
            }
            //- Return y direction
            vector& y()
            {
                return y_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
