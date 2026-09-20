/*---------------------------------------------------------------------------*\
  cosamEnergyAudit — runtime discrete kinetic-energy audit for OpenFOAM (v2406)
  Copyright (C) 2026 COSAM, Bursa Technical University. GPL-3.0-or-later.

Class
    Foam::functionObjects::kineticEnergyAudit

Description
    Reconstructs the discrete momentum equation solved by pimpleFoam (laminar,
    PISO/PIMPLE) from the registered final fields U, p, phi (and their old-time
    levels) with the case's own fvSchemes operators, dots it with U^{n+1} and
    closes the discrete kinetic-energy budget exactly (see DERIVATION.md):

        dK/dt = -eps_nu - e_time - e_conv - e_mesh - e_cont - e_pres - e_dev + e_iter

    Primary convention: multiplier U^{n+1}, convective reference midPoint
    (e_conv = scheme - linear, e_mesh = linear - midPoint), convective flux
    phi^n ("as solved" in PISO), half-half face-to-cell distribution,
    e_pres = p_P (G^T U)_P V_P.  Alternatives (U^{n+1/2} multiplier, owner
    distribution, phi^{n+1} via e_phi) are written as global time series.

    Restrictions (v0.1): all non-empty boundary patches must be coupled
    (cyclic/processor); laminar incompressible; Gauss linear gradients.

Usage
    functions
    {
        kineticEnergyAudit
        {
            type            kineticEnergyAudit;
            libs            (cosamEnergyAuditFunctionObjects);
            U               U;      // optional
            p               p;      // optional
            phi             phi;    // optional
            writeFields     true;   // optional
            maskCoeff       1e-3;   // nu_num mask: S:S < maskCoeff*<S:S>
            writeControl    writeTime;
        }
    }

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_kineticEnergyAudit_H
#define functionObjects_kineticEnergyAudit_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
namespace functionObjects
{

class kineticEnergyAudit
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private data

        word UName_;
        word pName_;
        word phiName_;
        bool writeFields_;
        scalar maskCoeff_;
        bool offline_;
        word instrName_;

        //- ddt scheme type read from fvSchemes ("Euler", "backward", "CrankNicolson")
        word ddtType_;
        scalar psi_;

        //- cumulative closure bookkeeping
        scalar K0_;
        scalar cumInt_;
        bool started_;

        //- per-cell fields (primary convention), written at write times
        autoPtr<volScalarField> eNumPtr_;
        autoPtr<volScalarField> nuNumPtr_;
        autoPtr<volScalarField> eTimePtr_;
        autoPtr<volScalarField> eConvPtr_;
        autoPtr<volScalarField> eMeshPtr_;
        autoPtr<volScalarField> eContPtr_;
        autoPtr<volScalarField> ePresPtr_;
        autoPtr<volScalarField> eDevPtr_;
        autoPtr<volScalarField> eIterPtr_;
        autoPtr<volScalarField> epsNuPtr_;
        autoPtr<volScalarField> SSPtr_;


    // Private Member Functions

        void writeFileHeader(Ostream& os);

        autoPtr<volScalarField> newField(const word& name) const;

        //- No copy construct / assignment
        kineticEnergyAudit(const kineticEnergyAudit&) = delete;
        void operator=(const kineticEnergyAudit&) = delete;


public:

    //- Runtime type information
    TypeName("kineticEnergyAudit");


    // Constructors

        kineticEnergyAudit
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~kineticEnergyAudit() = default;


    // Member Functions

        virtual bool read(const dictionary& dict);
        virtual bool execute();
        virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif
