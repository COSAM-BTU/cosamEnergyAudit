/*---------------------------------------------------------------------------*\
  cosamEnergyAudit — runtime discrete kinetic-energy audit for OpenFOAM (v2406)
  Copyright (C) 2026 COSAM, Bursa Technical University. GPL-3.0-or-later.
\*---------------------------------------------------------------------------*/

#include "kineticEnergyAudit.H"
#include "fvCFD.H"
#include "linear.H"
#include "midPoint.H"
#include "convectionScheme.H"
#include "coupledFvPatchField.H"
#include "emptyFvPatch.H"
#include "turbulentTransportModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(kineticEnergyAudit, 0);
    addToRunTimeSelectionTable(functionObject, kineticEnergyAudit, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::volScalarField>
Foam::functionObjects::kineticEnergyAudit::newField(const word& name) const
{
    return autoPtr<volScalarField>::New
    (
        IOobject
        (
            name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );
}


void Foam::functionObjects::kineticEnergyAudit::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Discrete kinetic-energy audit (primary convention: U^{n+1}, midPoint ref, phi^n, half-half)");
    writeHeader(os, "All global values are per unit volume; e_* signed; dKdt = -epsNu - eTime - eConv - eMesh - eCont - ePres - eDev + eIter");
    writeCommented(os, "t");
    writeTabbed(os, "dt"); writeTabbed(os, "timeIndex"); writeTabbed(os, "CoMax");
    writeTabbed(os, "K"); writeTabbed(os, "dKdt");
    writeTabbed(os, "epsNu"); writeTabbed(os, "twoNuSS"); writeTabbed(os, "Omega"); writeTabbed(os, "twoNuOmega");
    writeTabbed(os, "eTime"); writeTabbed(os, "eTimeDiss"); writeTabbed(os, "eTimeStore");
    writeTabbed(os, "eConv"); writeTabbed(os, "eMesh"); writeTabbed(os, "eCont"); writeTabbed(os, "ePres");
    writeTabbed(os, "eDev"); writeTabbed(os, "ePhi"); writeTabbed(os, "eIter"); writeTabbed(os, "sumE");
    writeTabbed(os, "Rclosure"); writeTabbed(os, "cumClosure"); writeTabbed(os, "RcellMax");
    writeTabbed(os, "chkConv"); writeTabbed(os, "chkDiff"); writeTabbed(os, "chkTime"); writeTabbed(os, "chkPres");
    writeTabbed(os, "eConvOwner"); writeTabbed(os, "eMeshOwner"); writeTabbed(os, "epsNuOwner");
    writeTabbed(os, "eTimeHalf"); writeTabbed(os, "eConvHalf"); writeTabbed(os, "eDiffHalf"); writeTabbed(os, "ePresHalf"); writeTabbed(os, "eIterHalf");
    writeTabbed(os, "maskFrac"); writeTabbed(os, "sumFluxC"); writeTabbed(os, "sumFluxD"); writeTabbed(os, "sumFluxX"); writeTabbed(os, "eIterInstr"); writeTabbed(os, "chkInstr");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::kineticEnergyAudit::kineticEnergyAudit
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    UName_("U"), pName_("p"), phiName_("phi"),
    writeFields_(true), maskCoeff_(1e-3), offline_(false), instrName_(word::null),
    ddtType_("Euler"), psi_(1.0),
    K0_(0), cumInt_(0), started_(false)
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::kineticEnergyAudit::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    dict.readIfPresent("U", UName_);
    dict.readIfPresent("p", pName_);
    dict.readIfPresent("phi", phiName_);
    dict.readIfPresent("writeFields", writeFields_);
    dict.readIfPresent("maskCoeff", maskCoeff_);
    dict.readIfPresent("offline", offline_);

    // ddt scheme type (and CN off-centring coefficient)
    {
        ITstream& is = mesh_.ddtScheme("ddt(" + UName_ + ')');
        is.rewind();
        ddtType_ = word(is);
        if (ddtType_ == "CrankNicolson")
        {
            token t(is);
            psi_ = t.isNumber() ? t.number() : 1.0;
        }
        Info<< type() << ": ddt scheme = " << ddtType_;
        if (ddtType_ == "CrankNicolson") { Info<< " psi=" << psi_; }
        Info<< nl;
    }

    dict.readIfPresent("instrumentedResidual", instrName_);

    // Guarantee old-time storage of the fields we need (must happen before the first step)
    if (!offline_)
    {
        const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
        const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);
        U.oldTime();
        if (ddtType_ == "backward") { U.oldTime().oldTime(); }
        phi.oldTime();
    }

    // Patch restriction check
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& pp = mesh_.boundary()[patchi];
        if (!pp.coupled() && !isA<emptyFvPatch>(pp) && pp.size() > 0)
        {
            FatalErrorInFunction
                << "kineticEnergyAudit v0.1 supports only coupled (cyclic/processor) and empty patches; patch "
                << pp.name() << " is of type " << pp.type() << exit(FatalError);
        }
    }

    if (writeFields_)
    {
        eNumPtr_ = newField("eNum");   nuNumPtr_ = newField("nuNum");
        eTimePtr_ = newField("eTime"); eConvPtr_ = newField("eConv");
        eMeshPtr_ = newField("eMesh"); eContPtr_ = newField("eCont");
        ePresPtr_ = newField("ePres"); eDevPtr_ = newField("eDev");
        eIterPtr_ = newField("eIter"); epsNuPtr_ = newField("epsNu");
        SSPtr_ = newField("SS");
    }
    return true;
}


bool Foam::functionObjects::kineticEnergyAudit::execute()
{
    const Time& rt = mesh_.time();
    const scalar dt = rt.deltaTValue();
    const label nCells = mesh_.nCells();
    const scalarField& V = mesh_.V();
    const scalar Vtot = gSum(V);

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const volScalarField& p = mesh_.lookupObject<volScalarField>(pName_);
    const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);
    // old-time levels: runtime from the registry, offline from the written raw inputs
    autoPtr<volVectorField> UprevPtr, UprevprevPtr; autoPtr<surfaceScalarField> phiPrevPtr; autoPtr<volVectorField> ddt0Ptr;
    if (offline_)
    {
        UprevPtr.reset(new volVectorField(IOobject("U_prev", rt.timeName(), mesh_, IOobject::MUST_READ, IOobject::NO_WRITE), mesh_));
        phiPrevPtr.reset(new surfaceScalarField(IOobject("phi_prev", rt.timeName(), mesh_, IOobject::MUST_READ, IOobject::NO_WRITE), mesh_));
        if (ddtType_ == "backward")
        {
            UprevprevPtr.reset(new volVectorField(IOobject("U_prevprev", rt.timeName(), mesh_, IOobject::MUST_READ, IOobject::NO_WRITE), mesh_));
        }
        if (ddtType_ == "CrankNicolson")
        {
            ddt0Ptr.reset(new volVectorField(IOobject("ddt0(" + UName_ + ')', rt.timeName(), mesh_, IOobject::MUST_READ, IOobject::NO_WRITE), mesh_));
        }
    }
    const surfaceScalarField& phiOld = offline_ ? phiPrevPtr() : phi.oldTime();
    const volVectorField& Uold = offline_ ? UprevPtr() : U.oldTime();

    const auto& turb = mesh_.lookupObject<incompressible::turbulenceModel>
    (
        incompressible::turbulenceModel::propertiesName
    );
    const volScalarField nuEff(turb.nuEff());

    // ---- operators with the case's own schemes
    const volTensorField gradU(fvc::grad(U, "grad(" + UName_ + ')'));
    const volVectorField gradp(fvc::grad(p, "grad(" + pName_ + ')'));
    // explicit temporal operator T(a) mirroring the scheme (variable-dt backward coefficients as in OpenFOAM)
    const label stepInRun = rt.timeIndex() - rt.startTimeIndex();
    const bool firstStepEuler = !offline_ && ((ddtType_ == "backward" && rt.timeIndex() < 2)
                                           || (ddtType_ == "CrankNicolson" && stepInRun < 2));
    vectorField Tvec(nCells, Zero);
    {
        const vectorField& a_ = U.primitiveField(); const vectorField& b_ = Uold.primitiveField();
        if (ddtType_ == "Euler" || firstStepEuler) { Tvec = (a_ - b_)/dt; }
        else if (ddtType_ == "backward")
        {
            const vectorField& d_ = (offline_ ? UprevprevPtr() : Uold.oldTime()).primitiveField();
            // previous time step: Time::deltaT0 (also read from <time>/uniform/time in -postProcess mode); fall back to dt
            const scalar dt0 = (rt.deltaT0Value() > SMALL ? rt.deltaT0Value() : dt);
            const scalar c = 1.0 + dt/(dt + dt0), c00 = dt*dt/(dt0*(dt + dt0)), c0 = c + c00;
            Tvec = (c*a_ - c0*b_ + c00*d_)/dt;
        }
        else if (ddtType_ == "CrankNicolson")
        {
            const vectorField& d0 = (offline_ ? ddt0Ptr() : mesh_.lookupObject<volVectorField>("ddt0(" + UName_ + ')')).primitiveField();
            Tvec = ((1.0 + psi_)/dt)*(a_ - b_) - psi_*d0;
        }
        else { FatalErrorInFunction << "Unsupported ddt scheme " << ddtType_ << exit(FatalError); }
    }
    const volVectorField divConvOF(fvc::div(phiOld, U, "div(" + phiName_ + ',' + UName_ + ')'));
    const volVectorField divConvNew(fvc::div(phi, U, "div(" + phiName_ + ',' + UName_ + ')'));
    const volVectorField lapOF(fvc::laplacian(nuEff, U, "laplacian(nuEff," + UName_ + ')'));
    const volVectorField D2
    (
        fvc::div(nuEff*dev2(T(gradU)), "div((nuEff*dev2(T(grad(" + UName_ + ")))))")
    );

    // ---- face fields
    tmp<fv::convectionScheme<vector>> cs = fv::convectionScheme<vector>::New
    (
        mesh_, phiOld, mesh_.divScheme("div(" + phiName_ + ',' + UName_ + ')')
    );
    const surfaceVectorField Ufsch(cs().interpolate(phiOld, U));
    tmp<fv::convectionScheme<vector>> cs1 = fv::convectionScheme<vector>::New
    (
        mesh_, phi, mesh_.divScheme("div(" + phiName_ + ',' + UName_ + ')')
    );
    const surfaceVectorField UfschNew(cs1().interpolate(phi, U));
    const surfaceVectorField Ulin(linearInterpolate(U));
    const surfaceTensorField gradUf(linearInterpolate(gradU));
    const surfaceScalarField nuf(linearInterpolate(nuEff));
    const surfaceScalarField& w = mesh_.weights();
    const surfaceScalarField& delta = mesh_.nonOrthDeltaCoeffs();
    const surfaceVectorField& kv = mesh_.nonOrthCorrectionVectors();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    // ---- per-cell accumulators (x V_P)
    scalarField eCont(nCells, 0), eConv(nCells, 0), eMesh(nCells, 0), GTa(nCells, 0), epsNu(nCells, 0);
    scalarField eConvOwn(nCells, 0), eMeshOwn(nCells, 0), epsNuOwn(nCells, 0);
    scalarField fluxC(nCells, 0), fluxD(nCells, 0), fluxX(nCells, 0);
    vectorField convVec(nCells, Zero), convVecNew(nCells, Zero), diffVec(nCells, Zero);
    scalarField ePhi(nCells, 0);

    // generic face contribution: cell cP (this side), other cell values (o), s = +1 (owner side) or -1
    auto faceContrib = [&]
    (
        const label cP, const vector& aP, const vector& aN, const scalar pP, const scalar pN,
        const vector& S, const scalar mS, const scalar wf, const scalar phif, const scalar phi1f,
        const vector& Uf, const vector& Uf1, const vector& Ul, const scalar df, const vector& kf,
        const tensor& gUf, const scalar nf, const scalar s, const bool ownerSide
    )
    {
        // orient: phi is owner-outward; for the neighbour side the outward flux is -phi
        const scalar phiOut = s*phif;          // outward volumetric flux of THIS cell
        const scalar phiOut1 = s*phi1f;
        const vector Umid = 0.5*(aP + aN);
        const vector dsch = Uf - Umid;
        const vector dlin = Ul - Umid;
        // convective (as solved, phi^n)
        convVec[cP] += phiOut*Uf;
        convVecNew[cP] += phiOut1*Uf1;
        const scalar dfsch = phiOut*((aP - aN) & dsch);          // symmetric in P<->N
        const scalar mf = phiOut*((aP - aN) & dlin);
        eCont[cP] += 0.5*magSqr(aP)*phiOut;
        eConv[cP] += 0.5*(dfsch - mf);
        eMesh[cP] += 0.5*mf;
        if (ownerSide) { eConvOwn[cP] += (dfsch - mf); eMeshOwn[cP] += mf; }
        fluxC[cP] += phiOut*(0.5*(aP & aN) + (Umid & dsch));
        ePhi[cP] += (aP & (phiOut1*Uf1 - phiOut*Uf));
        // pressure: G^T a (owner-oriented S, a_O - a_N in owner orientation)
        const vector aO = ownerSide ? aP : aN;
        const vector aNN = ownerSide ? aN : aP;
        const scalar pO = ownerSide ? pP : pN;
        const scalar pNN = ownerSide ? pN : pP;
        const scalar omega = ownerSide ? wf : (1.0 - wf);
        GTa[cP] += omega*(S & (aO - aNN));
        fluxX[cP] += s*(wf*pO*(S & aNN) + (1.0 - wf)*pNN*(S & aO));
        // diffusion (owner-oriented J)
        const vector corr = (kf & gUf);
        const vector J = nf*mS*(df*(aNN - aO) + corr);
        diffVec[cP] += s*J;
        const scalar gf = nf*mS*df*magSqr(aNN - aO);
        const scalar epsf = 0.5*(gf - ((aO - aNN) & (nf*mS*corr)));
        epsNu[cP] += epsf;
        if (ownerSide) { epsNuOwn[cP] += 2.0*epsf; }
        fluxD[cP] += s*(0.5*(aO + aNN) & J);
    };

    // internal faces
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, facei)
    {
        const label O = own[facei], N = nei[facei];
        faceContrib(O, U[O], U[N], p[O], p[N], Sf[facei], magSf[facei], w[facei], phiOld[facei], phi[facei],
                    Ufsch[facei], UfschNew[facei], Ulin[facei], delta[facei], kv[facei], gradUf[facei], nuf[facei], 1.0, true);
        faceContrib(N, U[N], U[O], p[N], p[O], Sf[facei], magSf[facei], w[facei], phiOld[facei], phi[facei],
                    Ufsch[facei], UfschNew[facei], Ulin[facei], delta[facei], kv[facei], gradUf[facei], nuf[facei], -1.0, false);
    }
    // coupled boundary faces (this side is always "owner" of the patch face)
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& pp = mesh_.boundary()[patchi];
        if (!pp.coupled() || pp.size() == 0) continue;
        const labelUList& fc = pp.faceCells();
        const vectorField UN(refCast<const coupledFvPatchField<vector>>(U.boundaryField()[patchi]).patchNeighbourField());
        const scalarField pN(refCast<const coupledFvPatchField<scalar>>(p.boundaryField()[patchi]).patchNeighbourField());
        const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& mSp = magSf.boundaryField()[patchi];
        const fvsPatchScalarField& wp = w.boundaryField()[patchi];
        const fvsPatchScalarField& phip = phiOld.boundaryField()[patchi];
        const fvsPatchScalarField& phi1p = phi.boundaryField()[patchi];
        const fvsPatchVectorField& Ufp = Ufsch.boundaryField()[patchi];
        const fvsPatchVectorField& Uf1p = UfschNew.boundaryField()[patchi];
        const fvsPatchVectorField& Ulp = Ulin.boundaryField()[patchi];
        const fvsPatchScalarField& dp = delta.boundaryField()[patchi];
        const fvsPatchVectorField& kp = kv.boundaryField()[patchi];
        const fvsPatchTensorField& gp = gradUf.boundaryField()[patchi];
        const fvsPatchScalarField& np = nuf.boundaryField()[patchi];
        forAll(fc, i)
        {
            const label O = fc[i];
            faceContrib(O, U[O], UN[i], p[O], pN[i], Sfp[i], mSp[i], wp[i], phip[i], phi1p[i],
                        Ufp[i], Uf1p[i], Ulp[i], dp[i], kp[i], gp[i], np[i], 1.0, true);
        }
    }

    // ---- cell-based terms
    const vectorField& a = U.primitiveField();
    const vectorField& b = Uold.primitiveField();
    const scalarField K1(0.5*magSqr(a)*V);
    const scalarField K0(0.5*magSqr(b)*V);
    const scalarField dKdt((K1 - K0)/dt);

    scalarField eTime(nCells, 0), eTimeDiss(nCells, 0), eTimeStore(nCells, 0);
    if (ddtType_ == "Euler" || firstStepEuler)
    {
        eTimeDiss = 0.5*magSqr(a - b)*V/dt;
        eTime = eTimeDiss;
    }
    else if (ddtType_ == "backward")
    {
        const vectorField& d = (offline_ ? UprevprevPtr() : Uold.oldTime()).primitiveField();
        auto X = [](const vector& u, const vector& v){ return 0.25*((u - v) & (3.0*u - v)); };
        const scalar dt0 = (rt.deltaT0Value() > SMALL ? rt.deltaT0Value() : dt);
        const bool constDt = (mag(dt - dt0) <= 1e-12*dt);
        forAll(a, i)
        {
            // storage: change of the BDF2 G-norm (constant-step form, telescopes exactly for any step sequence)
            eTimeStore[i] = (X(a[i], b[i]) - X(b[i], d[i]))*V[i]/dt;
            if (constDt)
            {
                // constant step: closed form (v1.0 behaviour, bit-identical): dissipative part >= 0
                eTimeDiss[i] = magSqr(a[i] - 2.0*b[i] + d[i])*V[i]/(4.0*dt);
                eTime[i] = eTimeStore[i] + eTimeDiss[i];
            }
            else
            {
                // variable step (v1.1 fix): exact operator route with the variable-step BDF2 coefficients; the
                // remainder after the G-norm change is reported as the dissipative part (not sign-definite for dt != dt0)
                eTime[i] = (a[i] & Tvec[i])*V[i] - dKdt[i];
                eTimeDiss[i] = eTime[i] - eTimeStore[i];
            }
        }
    }
    else if (ddtType_ == "CrankNicolson")
    {
        const vectorField& d0 = (offline_ ? ddt0Ptr() : mesh_.lookupObject<volVectorField>("ddt0(" + UName_ + ')')).primitiveField();
        forAll(a, i)
        {
            eTimeDiss[i] = ((1.0 + psi_)/(2.0*dt))*magSqr(a[i] - b[i])*V[i];
            eTimeStore[i] = psi_*(dKdt[i] - (a[i] & d0[i])*V[i]);
            eTime[i] = eTimeStore[i] + eTimeDiss[i];
        }
    }
    else
    {
        FatalErrorInFunction << "Unsupported ddt scheme " << ddtType_ << exit(FatalError);
    }
    // check: closed form vs operator route
    scalarField eTimeOp((a & Tvec)*V - dKdt);
    if (!offline_) { const volVectorField ddtU(fvc::ddt(U)); eTimeOp = (a & ddtU.primitiveField())*V - dKdt; }
    const scalar chkTime = gMax(mag(eTimeOp - eTime))/(gMax(mag(eTime)) + SMALL);

    const scalarField eDev(-(a & D2.primitiveField())*V);
    const scalarField ePres(p.primitiveField()*GTa);
    // momentum residual as solved: T + C(phi^n) - D - D2 + Gp
    const vectorField rV
    (
        Tvec*V + convVec - diffVec - D2.primitiveField()*V + gradp.primitiveField()*V
    );
    const scalarField eIter(a & rV);
    scalar gIterInstr = 0, chkInstr = -1;
    if (instrName_ != word::null && mesh_.foundObject<volVectorField>(instrName_))
    {
        const volVectorField& rT = mesh_.lookupObject<volVectorField>(instrName_);
        const scalarField eIterInstr(a & rT.primitiveField());
        gIterInstr = gSum(eIterInstr)/Vtot;
        chkInstr = gMax(mag(rV - rT.primitiveField()))/(gMax(mag(rV)) + SMALL);
    }

    // consistency checks of the manual face loops against fvc operators
    const scalar chkConv = gMax(mag(convVec - divConvOF.primitiveField()*V))/(gMax(mag(convVec)) + SMALL);
    const scalar chkDiff = gMax(mag(diffVec - lapOF.primitiveField()*V))/(gMax(mag(diffVec)) + SMALL);
    const scalar aGp = gSum((a & gradp.primitiveField())*V);
    const scalar chkPres = mag(aGp - gSum(ePres))/(mag(aGp) + SMALL);

    // closure: cell form and global
    const scalarField rhsCell
    (
        -eTime - eConv - eMesh - eCont - epsNu - eDev - ePres + eIter - (fluxC - fluxD + fluxX)
    );
    const scalar RcellMax = gMax(mag(dKdt - rhsCell))/(gMax(mag(dKdt)) + SMALL);

    const scalar gK1 = gSum(K1)/Vtot, gdKdt = gSum(dKdt)/Vtot;
    const scalar gEps = gSum(epsNu)/Vtot, gTime = gSum(eTime)/Vtot, gTD = gSum(eTimeDiss)/Vtot, gTS = gSum(eTimeStore)/Vtot;
    const scalar gConv = gSum(eConv)/Vtot, gMesh = gSum(eMesh)/Vtot, gCont = gSum(eCont)/Vtot;
    const scalar gPres = gSum(ePres)/Vtot, gDev = gSum(eDev)/Vtot, gIter = gSum(eIter)/Vtot, gPhi = gSum(ePhi)/Vtot;
    const scalar sumE = gTime + gConv + gMesh + gCont + gPres + gDev;
    const scalar rhs = -gEps - sumE + gIter;
    const scalar Rclosure = mag(gdKdt - rhs)/max(max(mag(gdKdt), mag(rhs)), SMALL);

    if (!started_) { K0_ = gSum(K0)/Vtot; cumInt_ = 0; started_ = true; }
    cumInt_ += (gEps + sumE - gIter)*dt;
    const scalar cumClosure = mag(gK1 - K0_ + cumInt_)/max(K0_, SMALL);

    // diagnostics
    const volSymmTensorField S(symm(gradU));
    const scalarField SS(magSqr(S.primitiveField()));
    const scalarField twoNuSS(2.0*nuEff.primitiveField()*SS*V);
    const volVectorField vort(fvc::curl(U));
    const scalar Omega = gSum(0.5*magSqr(vort.primitiveField())*V)/Vtot;
    const scalar nuMean = gSum(nuEff.primitiveField()*V)/Vtot;
    const scalar CoMax = 0.5*gMax(fvc::surfaceSum(mag(phi))().primitiveField()/V)*dt;

    // alternative multiplier U^{n+1/2} (global only)
    const vectorField m(0.5*(a + b));
    const scalar eTimeHalf = gSum((m & Tvec)*V - dKdt)/Vtot;
    const scalar eConvHalf = gSum(m & convVec)/Vtot;
    const scalar eDiffHalf = gSum(-(m & diffVec))/Vtot;
    const scalar ePresHalf = gSum((m & gradp.primitiveField())*V)/Vtot;
    const scalar eIterHalf = gSum(m & rV)/Vtot;

    // nu_num indicator with mask
    const scalar SSmean = gSum(SS*V)/Vtot;
    label nMasked = 0;
    if (writeFields_)
    {
        volScalarField& eNum = eNumPtr_(); volScalarField& nuNum = nuNumPtr_();
        forAll(a, i)
        {
            const scalar en = (eTime[i] + eConv[i] + eMesh[i] + eCont[i] + ePres[i] + eDev[i] - eIter[i])/V[i];
            eNum[i] = en;
            if (SS[i] <= maskCoeff_*SSmean || SS[i] <= VSMALL) { nuNum[i] = 0; ++nMasked; }
            else { nuNum[i] = en/(2.0*SS[i]); }
        }
        eTimePtr_().primitiveFieldRef() = eTime/V; eConvPtr_().primitiveFieldRef() = eConv/V;
        eMeshPtr_().primitiveFieldRef() = eMesh/V; eContPtr_().primitiveFieldRef() = eCont/V;
        ePresPtr_().primitiveFieldRef() = ePres/V; eDevPtr_().primitiveFieldRef() = eDev/V;
        eIterPtr_().primitiveFieldRef() = eIter/V; epsNuPtr_().primitiveFieldRef() = epsNu/V;
        SSPtr_().primitiveFieldRef() = SS;
        eNum.correctBoundaryConditions(); nuNum.correctBoundaryConditions();
    }
    reduce(nMasked, sumOp<label>());
    const scalar maskFrac = scalar(nMasked)/scalar(returnReduce(nCells, sumOp<label>()));
    const scalar gTwoNuSS = gSum(twoNuSS)/Vtot;
    const scalar gFluxC = gSum(fluxC)/Vtot, gFluxD = gSum(fluxD)/Vtot, gFluxX = gSum(fluxX)/Vtot;
    const scalar gConvOwn = gSum(eConvOwn)/Vtot, gMeshOwn = gSum(eMeshOwn)/Vtot, gEpsOwn = gSum(epsNuOwn)/Vtot;

    if (Pstream::master())
    {
        Ostream& os = file();
        writeCurrentTime(os);
        os << tab << dt << tab << rt.timeIndex() << tab << CoMax
           << tab << gK1 << tab << gdKdt
           << tab << gEps << tab << gTwoNuSS << tab << Omega << tab << 2.0*nuMean*Omega
           << tab << gTime << tab << gTD << tab << gTS
           << tab << gConv << tab << gMesh << tab << gCont << tab << gPres
           << tab << gDev << tab << gPhi << tab << gIter << tab << sumE
           << tab << Rclosure << tab << cumClosure << tab << RcellMax
           << tab << chkConv << tab << chkDiff << tab << chkTime << tab << chkPres
           << tab << gConvOwn << tab << gMeshOwn << tab << gEpsOwn
           << tab << eTimeHalf << tab << eConvHalf << tab << eDiffHalf << tab << ePresHalf << tab << eIterHalf
           << tab << maskFrac << tab << gFluxC << tab << gFluxD << tab << gFluxX << tab << gIterInstr << tab << chkInstr << endl;
    }

    Log << type() << " " << name() << ": K=" << gK1 << " dKdt=" << gdKdt
        << " epsNu=" << gEps << " sumE=" << sumE << " eIter=" << gIter
        << " Rclosure=" << Rclosure << " RcellMax=" << RcellMax
        << " chk(conv,diff,time,pres)=(" << chkConv << "," << chkDiff << "," << chkTime << "," << chkPres << ")"
        << " timeIndex=" << rt.timeIndex() << " U.timeIndex=" << U.timeIndex() << nl;

    return true;
}


bool Foam::functionObjects::kineticEnergyAudit::write()
{
    if (!writeFields_) return true;

    eNumPtr_().write(); nuNumPtr_().write();
    eTimePtr_().write(); eConvPtr_().write(); eMeshPtr_().write(); eContPtr_().write();
    ePresPtr_().write(); eDevPtr_().write(); eIterPtr_().write(); epsNuPtr_().write(); SSPtr_().write();

    // raw inputs for offline recomputation
    if (offline_) return true;
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    const surfaceScalarField& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);
    volVectorField Uprev(IOobject("U_prev", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), U.oldTime());
    Uprev.write();
    surfaceScalarField phiPrev(IOobject("phi_prev", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), phi.oldTime());
    phiPrev.write();
    if (ddtType_ == "backward")
    {
        volVectorField Upp(IOobject("U_prevprev", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), U.oldTime().oldTime());
        Upp.write();
    }
    return true;
}

// ************************************************************************* //
