#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "surfaceFilmModel.H"
#include "basicNuclearCloud.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "radiationModel.H"
#include "localEulerDdtScheme.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " with nuclear particle clouds"
        " and surface film modelling. For the thermic part the" 
        " Boussinesq approximation is used."
        " The particles are modelled as nuclear particles,"
        " i.e. they undergo decay, and emit decay heat."
    );

    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRegionControls.H"
    #include "createUfIfPresent.H"


    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        FP.storeGlobalPositions();
        FP.evolve();
        
        surfaceFilm.evolve();

        if (solvePrimaryRegion)
        {
            while (pimple.loop())
            {
                #include "UEqn.H"
                
                vort = fvc::curl(U);
                #include "TEqn.H"

                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    laminarTransport.correct();
                    turbulence->correct();
                }
            }
        }

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;
    return 0;
}
