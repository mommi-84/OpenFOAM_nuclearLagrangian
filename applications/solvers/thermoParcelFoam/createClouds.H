/*const word kinematicCloudName
(
    args.getOrDefault<word>("cloud", "kinematicCloud")
);

Info<< "\nConstructing kinematicCloud " << kinematicCloudName << endl;

basicKinematicCloud parcels
(
    kinematicCloudName,
    rhoInf,
    U,
    muc,
    g
);

const word kinematicCloudName2
(
    args.getOrDefault<word>("cloud2", "kinematicCloud2")
);

Info<< "\nConstructing kinematicCloud " << kinematicCloudName2 << endl;

basicKinematicCloud parcels2
(
    kinematicCloudName2,
    rhoInf,
    U,
    muc,
    g
);

*/
Info<< "\nConstructing first fission product cloud" << endl;
basicKinematicCloud FP1
(
    "FPCloud1",
    rhoInf,
    U,
    muc,
    vort,
    g
);

/*Info<< "\nConstructing second fission product cloud" << endl;
basicKinematicCloud FP2
(
    "FPCloud2",
    rhoInf,
    U,
    muc,
    g
);*/
