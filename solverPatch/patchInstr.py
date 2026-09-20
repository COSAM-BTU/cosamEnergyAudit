# patch a copy of pimpleFoam into cosamPimpleFoamInstr: writes rTrueV = -(UEqn == -grad p).residual(), phiUsed, phiHbyA
import re, sys, os
d = sys.argv[1]
c = open(d + "/pimpleFoam.C").read()
c = c.replace('            if (pimple.turbCorr())\n',
'''            // --- V4 instrumentation: residual of the momentum system actually assembled, with the final U and p
            {
                // as-solved residual of the assembled momentum system with the final U, p:  A U - H(U) + grad p  (x V)
                rTrueV.primitiveFieldRef() = (UEqn.A()().primitiveField()*U.primitiveField() - UEqn.H()().primitiveField() + fvc::grad(p)().primitiveField())*mesh.V();
                rTrueV.correctBoundaryConditions();
            }

            if (pimple.turbCorr())
''')
assert "V4 instrumentation" in c
open(d + "/pimpleFoam.C", "w").write(c)
u = open(d + "/UEqn.H").read()
u = "phiUsed = phi;   // V4: flux actually used in the momentum assembly\n" + u
u = u.replace("    solve(UEqn == -fvc::grad(p));", "    solve(UEqn == -fvc::grad(p));\n    rPredV.primitiveFieldRef() = (UEqn.A()().primitiveField()*U.primitiveField() - UEqn.H()().primitiveField() + fvc::grad(p)().primitiveField())*mesh.V();   // V4: right after the predictor solve")
assert "rPredV" in u
open(d + "/UEqn.H", "w").write(u)
pe = open(d + "/pEqn.H").read()
pe = pe.replace("MRF.makeRelative(phiHbyA);", "phiHbyAOut = phiHbyA;   // V4\nMRF.makeRelative(phiHbyA);", 1)
assert "V4" in pe
open(d + "/pEqn.H", "w").write(pe)
cf = open(d + "/createFields.H").read()
cf += '''
// --- V4 instrumentation fields
volVectorField rTrueV
(
    IOobject("rTrueV", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh, dimensionedVector(dimVelocity*dimVolume/dimTime, Zero)
);
volVectorField rPredV
(
    IOobject("rPredV", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh, dimensionedVector(dimVelocity*dimVolume/dimTime, Zero)
);
surfaceScalarField phiUsed(IOobject("phiUsed", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), phi);
surfaceScalarField phiHbyAOut(IOobject("phiHbyA", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), phi);
'''
open(d + "/createFields.H", "w").write(cf)
mf = open(d + "/Make/files").read()
mf = re.sub(r"EXE\s*=.*", "EXE = $(FOAM_USER_APPBIN)/cosamPimpleFoamInstr", mf)
open(d + "/Make/files", "w").write(mf)
print("patched", d)
