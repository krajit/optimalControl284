
// integrate in the interior
J = 0.5  * gSum(volField * (Foam::pow(y.internalField() - yd.internalField(), 2)));

// integrate in the boundary patches
forAll(patchesToIntegrate, i)
{
    label patchi = mesh.boundaryMesh().findPatchID(patchesToIntegrate[i]);
    J += 0.5 * lambda * gSum(mesh.magSf().boundaryField()[patchi] * Foam::pow(u.boundaryField()[patchi], 2));
}
