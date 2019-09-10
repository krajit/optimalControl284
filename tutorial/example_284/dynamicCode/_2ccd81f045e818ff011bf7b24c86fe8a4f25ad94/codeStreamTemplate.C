/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR AUTHOR,AFFILIATION
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 0 "/home/ajit/Desktop/optimalControl283/tutorial/example_283/0/yd.#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_2ccd81f045e818ff011bf7b24c86fe8a4f25ad94
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 0 "/home/ajit/Desktop/optimalControl283/tutorial/example_283/0/yd.#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&>(dict);
	    const fvMesh& mesh = refCast<const fvMesh>(d.db());

	    scalarField yd(mesh.nCells(), 0.);
		scalar pi = constant::mathematical::pi;

	    forAll(yd, i)
	    {
		const scalar x = mesh.C()[i][0];
		const scalar y = mesh.C()[i][1];

		yd[i] = sin(pi*x)*sin(pi*y);
	    }

	    yd.writeEntry("", os);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

