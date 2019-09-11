/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "robinFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
// #include "volFields.H"
// #include "EulerDdtScheme.H"
// #include "CrankNicolsonDdtScheme.H"
// #include "backwardDdtScheme.H"
// #include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class Type>
Foam::robinFvPatchField<Type>::robinFvPatchField(
    const fvPatch &p,
    const DimensionedField<Type, volMesh> &iF)
    : mixedFvPatchField<Type>(p, iF),
      beta_(0.0)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}

template <class Type>
Foam::robinFvPatchField<Type>::robinFvPatchField(
    const robinFvPatchField &ptf,
    const fvPatch &p,
    const DimensionedField<Type, volMesh> &iF,
    const fvPatchFieldMapper &mapper)
    : mixedFvPatchField<Type>(ptf, p, iF, mapper),
      beta_(ptf.beta_)
{
}

template <class Type>
Foam::robinFvPatchField<Type>::robinFvPatchField(
    const fvPatch &p,
    const DimensionedField<Type, volMesh> &iF,
    const dictionary &dict)
    : mixedFvPatchField<Type>(p, iF),
      beta_(readScalar(dict.lookup("beta"))) // exchange coeff
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=(
            Field<Type>("value", dict, p.size()));
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}

template <class Type>
Foam::robinFvPatchField<Type>::robinFvPatchField(
    const robinFvPatchField &ptpsf)
    : mixedFvPatchField<Type>(ptpsf),
      beta_(ptpsf.beta_)
{
}

template <class Type>
Foam::robinFvPatchField<Type>::robinFvPatchField(
    const robinFvPatchField &ptpsf,
    const DimensionedField<Type, volMesh> &iF)
    : mixedFvPatchField<Type>(ptpsf, iF),
      beta_(ptpsf.beta_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class Type>
void Foam::robinFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Ajit: (Sept 11, 19) just upadating the value fraction is enough I think.

     this->valueFraction() = (beta_ / this->patch().deltaCoeffs()) / (1 + beta_ / this->patch().deltaCoeffs());

    mixedFvPatchField<Type>::updateCoeffs();
}

template <class Type>
void Foam::robinFvPatchField<Type>::write(Ostream &os) const
{
    fvPatchField<Type>::write(os);

    // os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    // os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);

    // if (lInf_ > 0)
    // {
    //     os.writeEntry("fieldInf", fieldInf_);
    //     os.writeEntry("lInf", lInf_);
    // }

    os.writeEntry("beta", beta_);
    this->writeEntry("value", os);
}

// ************************************************************************* //
