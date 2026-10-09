"""Write 0/, constant/ and system/ dictionaries for the simpleFoam case."""
import math
import os
import sys

import params as P

HDR = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
"""

N_PROCS = P.N_PROCS
END_ITER = P.N_ITER_MAX
WRITE_INTERVAL = P.WRITE_INTERVAL
# split streamwise in 2 for >= 8 processes (fewer, larger processor faces
# than N thin spanwise slabs); scotch is not built in the Ubuntu package
NX = 2 if N_PROCS >= 8 and N_PROCS % 2 == 0 else 1


def w(case, loc, obj, body, cls="dictionary"):
    os.makedirs(os.path.join(case, loc), exist_ok=True)
    with open(os.path.join(case, loc, obj), "w") as f:
        f.write(HDR.format(cls=cls, loc=loc, obj=obj) + body)


def freestream():
    a = math.radians(P.AOA_DEG)
    U = (P.U_INF * math.cos(a), 0.0, P.U_INF * math.sin(a))
    k = 1.5 * (P.TU * P.U_INF) ** 2
    omega = k / (P.NU * P.VISC_RATIO)
    tu = 100 * P.TU
    # Langtry-Menter correlation (Tu <= 1.3 %, zero pressure gradient)
    re_theta = 1173.51 - 589.428 * tu + 0.2196 / tu ** 2
    return U, k, omega, re_theta, a


def field(case, name, cls, dim, internal, bcs):
    body = f"dimensions      {dim};\n\ninternalField   {internal};\n\nboundaryField\n{{\n"
    for patch, bc in bcs.items():
        body += f"    {patch}\n    {{\n"
        for line in bc.strip().split("\n"):
            body += f"        {line.strip()}\n"
        body += "    }\n"
    body += "}\n"
    w(case, "0", name, body, cls)


def main(case="case"):
    U, k, omega, rt, a = freestream()
    Us = f"uniform ({U[0]:.6g} {U[1]:.6g} {U[2]:.6g})"
    sym = "type symmetryPlane;"

    def scalar_ff(v):
        return f"type inletOutlet;\ninletValue uniform {v:.6g};\nvalue uniform {v:.6g};"

    field(case, "U", "volVectorField", "[0 1 -1 0 0 0 0]", Us, {
        "wing": "type noSlip;", "root": sym,
        "side": f"type freestreamVelocity;\nfreestreamValue {Us};\nvalue {Us};",
        "farfield": f"type freestreamVelocity;\nfreestreamValue {Us};\nvalue {Us};",
        "outlet": f"type inletOutlet;\ninletValue {Us};\nvalue {Us};"})
    field(case, "p", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", {
        "wing": "type zeroGradient;", "root": sym,
        "side": "type freestreamPressure;\nfreestreamValue uniform 0;\nvalue uniform 0;",
        "farfield": "type freestreamPressure;\nfreestreamValue uniform 0;\nvalue uniform 0;",
        "outlet": "type fixedValue;\nvalue uniform 0;"})
    for name, v, wall in (
            ("k", k, "type fixedValue;\nvalue uniform 1e-14;"),
            ("omega", omega, "type omegaWallFunction;\nvalue uniform 1e3;"),
            ("ReThetat", rt, "type zeroGradient;"),
            ("gammaInt", 1.0, "type zeroGradient;")):
        dim = {"k": "[0 2 -2 0 0 0 0]", "omega": "[0 0 -1 0 0 0 0]"}.get(
            name, "[0 0 0 0 0 0 0]")
        ff = scalar_ff(v)
        field(case, name, "volScalarField", dim, f"uniform {v:.6g}", {
            "wing": wall, "root": sym, "side": ff, "farfield": ff, "outlet": ff})
    field(case, "nut", "volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", {
        "wing": "type nutLowReWallFunction;\nvalue uniform 0;", "root": sym,
        "side": "type calculated;\nvalue uniform 0;",
        "farfield": "type calculated;\nvalue uniform 0;",
        "outlet": "type calculated;\nvalue uniform 0;"})

    w(case, "constant", "transportProperties",
      f"transportModel  Newtonian;\n\nnu              {P.NU:.6g};\n")
    w(case, "constant", "turbulenceProperties",
      "simulationType  RAS;\n\nRAS\n{\n    RASModel        kOmegaSSTLM;\n"
      "    turbulence      on;\n    printCoeffs     on;\n}\n")

    lift = f"({-math.sin(a):.8f} 0 {math.cos(a):.8f})"
    drag = f"({math.cos(a):.8f} 0 {math.sin(a):.8f})"
    w(case, "system", "controlDict", f"""
application     simpleFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         {END_ITER};
deltaT          1;
writeControl    timeStep;
writeInterval   {WRITE_INTERVAL};
purgeWrite      2;
writeFormat     binary;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{{
    forceCoeffs
    {{
        type            forceCoeffs;
        libs            ("libforces.so");
        writeControl    timeStep;
        writeInterval   1;
        log             false;
        patches         (wing);
        rho             rhoInf;
        rhoInf          {P.RHO};
        liftDir         {lift};
        dragDir         {drag};
        CofR            ({0.25 * P.C} 0 0);
        pitchAxis       (0 1 0);
        magUInf         {P.U_INF:.6g};
        lRef            {P.C};
        Aref            {P.C * P.B};
    }}
    residuals
    {{
        type            solverInfo;
        libs            ("libutilityFunctionObjects.so");
        writeControl    timeStep;
        writeInterval   1;
        fields          (U p k omega ReThetat gammaInt);
    }}
    yPlus
    {{
        type            yPlus;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        log             false;
    }}
    wallShearStress
    {{
        type            wallShearStress;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        patches         (wing);
        log             false;
    }}
}}
""")
    w(case, "system", "fvSchemes", """
ddtSchemes      { default steadyState; }
gradSchemes
{
    default         cellLimited Gauss linear 1;
    grad(U)         cellLimited Gauss linear 1;
}
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss limitedLinear 1;
    div(phi,omega)  bounded Gauss limitedLinear 1;
    div(phi,ReThetat) bounded Gauss limitedLinear 1;
    div(phi,gammaInt) bounded Gauss limitedLinear 1;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited corrected 0.33; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited corrected 0.33; }
wallDist        { method meshWave; }
""")
    w(case, "system", "fvSolution", """
solvers
{
    p
    {
        solver          GAMG;
        smoother        DICGaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
        maxIter         50;
        nCellsInCoarsestLevel 50;
    }
    "(U|k|omega|ReThetat|gammaInt)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-9;
        relTol          0.1;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      no;
}
relaxationFactors
{
    fields          { p %s; }
    equations
    {
        U               %s;
        ".*"            %s;
    }
}
""" % (P.RELAX_P, P.RELAX_U, P.RELAX_TURB))
    w(case, "system", "decomposeParDict",
      f"numberOfSubdomains {N_PROCS};\n\nmethod          hierarchical;\n\nhierarchicalCoeffs\n{{\n    n           ({NX} {N_PROCS // NX} 1);\n    order       xyz;\n}}\n")
    print(f"U_inf={P.U_INF:.4f} m/s  k={k:.4g}  omega={omega:.4g}  ReThetat={rt:.1f}")
    print(f"first layer height (y+={P.YPLUS_TARGET}) = {P.H1:.4e} m")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "case")
