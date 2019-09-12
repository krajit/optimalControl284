/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

Application
    laplaceAdjointFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"

    // Disable solvers performance output
    lduMatrix::debug = 0;
    solverPerformance::debug = 0;

    // Cost function value
    scalar J = 0;
    scalar Jold = 0;
    scalar Jk = 0;

    // Compute cost function value
#include "costFunctionValue.H"

    std::ofstream file("results.csv");
    file << 0 << "," << J << nl;
    file.close();

    while (runTime.loop() && (fabs(J - Jold) > tol) && (alpha > tol))
    {
        // save old cost value
        Jold = J;

        // Primal equation
        solve(fvm::laplacian(k, y));

        // Adjoint equation
        solve(fvm::laplacian(k, p) + (y - yd));

        // Save current control
        uk = u;

// calculate current cost
#include "costFunctionValue.H"
        Jk = J;

        bool alphaFound = false;

        // calculate derivative^2 integrate((lambda*u + beta*p)^2 dv). Why??
        scalar phip0 = 0;

        // integrate in the boundary patches
        forAll(patchesToIntegrate, i)
        {
            label patchi = mesh.boundaryMesh().findPatchID(patchesToIntegrate[i]);
            //phi0 += 0.5 * lambda.value() * gSum(mesh.magSf().boundaryField()[patchi] * Foam::pow(u.boundaryField()[patchi], 2));
            phip0 += gSum(
                mesh.magSf().boundaryField()[patchi] *
                Foam::pow(lambda * u.boundaryField()[patchi] + beta * p.boundaryField()[patchi], 2)
            );
        }

        while ((!alphaFound) && (alpha > tol))
        {
            u = uk - alpha * (lambda * uk + beta * p);

            // truncate u for constrained control set
            forAll(u, i)
            {
                u[i] = min(u[i], uMax[i]);
                u[i] = max(u[i], uMin[i]);
            }
            u.correctBoundaryConditions();

            // get new y
            solve(fvm::laplacian(k, y));
// get new cost
#include "costFunctionValue.H"

            // backtracking step to find alpha
            if (J <= Jk - c1 * alpha * phip0)
            {
                Info << "alpha found, alpha = " << alpha << ", J = " << J << ", phip0" << phip0 << endl;
                alphaFound = true;
            }
            else
            {
                Info << "alpha NOT found, alpha = " << alpha << ", J = " << J << ", phip0" << phip0 << endl;
                alpha = c2 * alpha;
            }
        }

        Info << "Iteration no. " << runTime.timeName() << " - "
             << "Cost value " << J
             << " - "
             << "Cost variation" << fabs(J - Jold) << endl;

        file.open("results.csv", std::ios::app);
        file << runTime.value() << "," << J << nl;
        file.close();
        runTime.write();
    }

    file.close();

    runTime++;
    y.write();
    yd.write();
    p.write();
    u.write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << endl;
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    Info << "\nEnd\n"
         << endl;
    return 0;
}

// ************************************************************************* //
