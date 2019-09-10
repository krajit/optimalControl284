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
#line 0 "/home/ajit/Desktop/optimalControl283/tutorial/example_283/0/p.#codeStream"
#include "fvCFD.H"
//                #include "constants.H"
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
    void codeStream_9e399175bd53084aec4836217dcf0c3831cb94ed
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 0 "/home/ajit/Desktop/optimalControl283/tutorial/example_283/0/p.#codeStream"
const IOdictionary &d = static_cast<const IOdictionary &>(dict);
        const fvMesh &mesh = refCast<const fvMesh>(d.db());
        const pointField& points = mesh.boundary()["dirichletBoundary"].Cf();
        const volScalarField& y = d.db().lookupObject<volScalarField>("y");
        const volScalarField& yd = d.db().lookupObject<volScalarField>("yd");
        label patchi = mesh.boundaryMesh().findPatchID("dirichletBoundary");
        scalarField pb(points.size(), 0);
        scalar lambda2 = 0.01; // need to be updated everytime this is changed in physicalProp dict
        scalar beta = 0.1;
        pb = (lambda2/beta)*(y.boundaryField()[patchi] - yd.boundaryField()[patchi]);
        os << "nonuniform " <<  pb;
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

