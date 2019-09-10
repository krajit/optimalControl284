
// integrate in the interior
J = 0.5 *lambda1.value()*gSum(volField * (Foam::pow(y.internalField() - yd.internalField(), 2) + lambda * Foam::pow(u.internalField(), 2)));


// integrate in the boundary patches
forAll(patchesToIntegrate, i)
{
    label patchi = mesh.boundaryMesh().findPatchID(patchesToIntegrate[i]);
    J += 0.5*lambda2.value() * gSum(
        mesh.magSf().boundaryField()[patchi] * (y.boundaryField()[patchi] - yd.boundaryField()[patchi]) * (y.boundaryField()[patchi] - yd.boundaryField()[patchi]));
}
