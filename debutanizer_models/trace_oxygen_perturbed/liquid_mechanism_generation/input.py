
# Data sources
database(
    thermoLibraries = ['thermo_combined','multi_trays_combined','Conjugated_diene','thermo_DFT_CCSDTF12_BAC','primaryThermoLibrary', 'primaryNS', 'NitrogenCurran', 'NOx2018', 'FFCM1(-)', 'SulfurLibrary', 'SulfurGlarborgH2S'],
    reactionLibraries = ['kinetics_combined','Xu_cyclopentadiene','Conjugated_diene','Klippenstein_Glarborg2016'],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species

species(
    label='N-BUTANE',
    reactive=True,
    structure=SMILES("CCCC")
)

species(
    label='2-BUTENE',
    reactive=True,
    structure=SMILES("CC=CC"),
)

species(
    label='1,3-BUTADIENE',
    reactive=True,
    structure=SMILES("C=CC=C"),
)

species(
    label='CYCLOPENTADIENE',
    reactive=True,
    structure=SMILES("C1=CCC=C1"),
)

species(
    label='BENZENE',
    reactive=True,
    structure=SMILES("C1=CC=CC=C1"),
)

species(
    label='1,3-CYCLOHEXADIENE',
    reactive=True,
    structure=SMILES("C1=CCCC=C1"),
)

species(
    label='TOLUENE',
    reactive=True,
    structure=SMILES("CC1=CC=CC=C1"),
)

species(
    label='STYRENE',
    reactive=True,
    structure=SMILES("C=CC1=CC=CC=C1"),
)

species(
    label='OXYGEN',
    reactive=True,
    structure=SMILES("[O][O]"),
)

# tray 1
constantTVLiquidReactor(
    temperature=(312.454,'K'),
    initialConcentrations={
        'N-BUTANE': (296.5433727912449, 'mol/m^3'),
        '2-BUTENE': (281.9891302971655, 'mol/m^3'),
        '1,3-BUTADIENE': (8502.019971200936, 'mol/m^3'),
        'CYCLOPENTADIENE': (0.3963735099235698, 'mol/m^3'),
        'BENZENE': (0.002313093411301333, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (1.6668544921121383e-05, 'mol/m^3'),
        'TOLUENE': (4.412356030093125e-06, 'mol/m^3'),
        'STYRENE': (4.5568709469732974e-07, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.952524987968715,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.03555062251952099,
        '2-BUTENE': 0.028217828680617063,
        '1,3-BUTADIENE': 0.936209746697449,
        'CYCLOPENTADIENE': 2.1704539884802992e-05,
        'BENZENE': 9.672668402175344e-08,
        '1,3-CYCLOHEXADIENE': 6.622546099317859e-10,
        'TOLUENE': 1.580833858714708e-10,
        'STYRENE': 1.5506067072355674e-11,

    },
)

# tray 1
constantTVLiquidReactor(
    temperature=(312.454,'K'),
    initialConcentrations={
        'N-BUTANE': (296.5433727912449, 'mol/m^3'),
        '2-BUTENE': (281.9891302971655, 'mol/m^3'),
        '1,3-BUTADIENE': (8502.019971200936, 'mol/m^3'),
        'CYCLOPENTADIENE': (0.3963735099235698, 'mol/m^3'),
        'BENZENE': (0.002313093411301333, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (1.6668544921121383e-05, 'mol/m^3'),
        'TOLUENE': (4.412356030093125e-06, 'mol/m^3'),
        'STYRENE': (4.5568709469732974e-07, 'mol/m^3'),
        'OXYGEN': (0.014366199425993986, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.952524987968715,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.035547691666747615,
        '2-BUTENE': 0.02821550235563072,
        '1,3-BUTADIENE': 0.9361325640002667,
        'CYCLOPENTADIENE': 2.1702750526237435e-05,
        'BENZENE': 9.67187097121624e-08,
        '1,3-CYCLOHEXADIENE': 6.622000125542256e-10,
        'TOLUENE': 1.5807035321276977e-10,
        'STYRENE': 1.550478872625475e-11,
        'OXYGEN': 8.244167234384159e-05,

    },
)

# tray 2
constantTVLiquidReactor(
    temperature=(312.538,'K'),
    initialConcentrations={
        'N-BUTANE': (278.3455685282178, 'mol/m^3'),
        '2-BUTENE': (303.33492828220324, 'mol/m^3'),
        '1,3-BUTADIENE': (8504.175010230736, 'mol/m^3'),
        'CYCLOPENTADIENE': (0.6984654262421716, 'mol/m^3'),
        'BENZENE': (0.005161927706812445, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (3.893928202055236e-05, 'mol/m^3'),
        'TOLUENE': (1.1314682345781493e-05, 'mol/m^3'),
        'STYRENE': (1.2245440974854323e-06, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.95431617490664,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.033680198687136106,
        '2-BUTENE': 0.030049333954571662,
        '1,3-BUTADIENE': 0.9362343848852366,
        'CYCLOPENTADIENE': 3.588184635147815e-05,
        'BENZENE': 1.9879866272394598e-07,
        '1,3-CYCLOHEXADIENE': 1.420268177169796e-09,
        'TOLUENE': 3.69865784339409e-10,
        'STYRENE': 3.790772296423821e-11,

    },
)

# tray 2
constantTVLiquidReactor(
    temperature=(312.538,'K'),
    initialConcentrations={
        'N-BUTANE': (278.3455685282178, 'mol/m^3'),
        '2-BUTENE': (303.33492828220324, 'mol/m^3'),
        '1,3-BUTADIENE': (8504.175010230736, 'mol/m^3'),
        'CYCLOPENTADIENE': (0.6984654262421716, 'mol/m^3'),
        'BENZENE': (0.005161927706812445, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (3.893928202055236e-05, 'mol/m^3'),
        'TOLUENE': (1.1314682345781493e-05, 'mol/m^3'),
        'STYRENE': (1.2245440974854323e-06, 'mol/m^3'),
        'OXYGEN': (0.0037644440004335374, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.95431617490664,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.03367918143572464,
        '2-BUTENE': 0.030048426367070093,
        '1,3-BUTADIENE': 0.9362061075654282,
        'CYCLOPENTADIENE': 3.5880762603155106e-05,
        'BENZENE': 1.9879265835852843e-07,
        '1,3-CYCLOHEXADIENE': 1.4202252804570626e-09,
        'TOLUENE': 3.698546131912021e-10,
        'STYRENE': 3.790657802786034e-11,
        'OXYGEN': 3.020324852927674e-05,

    },
)

# tray 3
constantTVLiquidReactor(
    temperature=(312.573,'K'),
    initialConcentrations={
        'N-BUTANE': (265.5744617039428, 'mol/m^3'),
        '2-BUTENE': (320.8062025094457, 'mol/m^3'),
        '1,3-BUTADIENE': (8500.628743472838, 'mol/m^3'),
        'CYCLOPENTADIENE': (1.1555555695567925, 'mol/m^3'),
        'BENZENE': (0.010811294325993766, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (8.54339308337641e-05, 'mol/m^3'),
        'TOLUENE': (2.730207125963592e-05, 'mol/m^3'),
        'STYRENE': (3.100100934795691e-06, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.95858806976151,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.032372649269745,
        '2-BUTENE': 0.03155534936321892,
        '1,3-BUTADIENE': 0.9360142455655718,
        'CYCLOPENTADIENE': 5.7350532465257566e-05,
        'BENZENE': 4.0131238515996403e-07,
        '1,3-CYCLOHEXADIENE': 3.003474480615726e-09,
        'TOLUENE': 8.605700214449536e-10,
        'STYRENE': 9.256955712952872e-11,

    },
)

# tray 3
constantTVLiquidReactor(
    temperature=(312.573,'K'),
    initialConcentrations={
        'N-BUTANE': (265.5744617039428, 'mol/m^3'),
        '2-BUTENE': (320.8062025094457, 'mol/m^3'),
        '1,3-BUTADIENE': (8500.628743472838, 'mol/m^3'),
        'CYCLOPENTADIENE': (1.1555555695567925, 'mol/m^3'),
        'BENZENE': (0.010811294325993766, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (8.54339308337641e-05, 'mol/m^3'),
        'TOLUENE': (2.730207125963592e-05, 'mol/m^3'),
        'STYRENE': (3.100100934795691e-06, 'mol/m^3'),
        'OXYGEN': (0.003609581259013653, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.95858806976151,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.03237169628183002,
        '2-BUTENE': 0.031554420435025696,
        '1,3-BUTADIENE': 0.9359866911242644,
        'CYCLOPENTADIENE': 5.73488441769775e-05,
        'BENZENE': 4.0130057130284036e-07,
        '1,3-CYCLOHEXADIENE': 3.0033860641608625e-09,
        'TOLUENE': 8.60544687935064e-10,
        'STYRENE': 9.256683206157077e-11,
        'OXYGEN': 2.9438057634280254e-05,

    },
)

# tray 4
constantTVLiquidReactor(
    temperature=(312.594,'K'),
    initialConcentrations={
        'N-BUTANE': (256.66151125242675, 'mol/m^3'),
        '2-BUTENE': (335.14948759638867, 'mol/m^3'),
        '1,3-BUTADIENE': (8494.954716660202, 'mol/m^3'),
        'CYCLOPENTADIENE': (1.8471594242720164, 'mol/m^3'),
        'BENZENE': (0.022013678224441835, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.0001824954338572565, 'mol/m^3'),
        'TOLUENE': (6.433146130626823e-05, 'mol/m^3'),
        'STYRENE': (7.675239702068004e-06, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.96286081439831,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.031460983769913756,
        '2-BUTENE': 0.032793914179010064,
        '1,3-BUTADIENE': 0.9356544469724923,
        'CYCLOPENTADIENE': 8.984359610455241e-05,
        'BENZENE': 8.029503473617415e-07,
        '1,3-CYCLOHEXADIENE': 6.309003616219479e-09,
        'TOLUENE': 1.997212123794867e-09,
        'STYRENE': 2.259160793685877e-10,

    },
)

# tray 4
constantTVLiquidReactor(
    temperature=(312.594,'K'),
    initialConcentrations={
        'N-BUTANE': (256.66151125242675, 'mol/m^3'),
        '2-BUTENE': (335.14948759638867, 'mol/m^3'),
        '1,3-BUTADIENE': (8494.954716660202, 'mol/m^3'),
        'CYCLOPENTADIENE': (1.8471594242720164, 'mol/m^3'),
        'BENZENE': (0.022013678224441835, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.0001824954338572565, 'mol/m^3'),
        'TOLUENE': (6.433146130626823e-05, 'mol/m^3'),
        'STYRENE': (7.675239702068004e-06, 'mol/m^3'),
        'OXYGEN': (0.0036066987806489003, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.96286081439831,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.03146005801580793,
        '2-BUTENE': 0.032792949202805906,
        '1,3-BUTADIENE': 0.9356269149044473,
        'CYCLOPENTADIENE': 8.984095241487684e-05,
        'BENZENE': 8.029267201736574e-07,
        '1,3-CYCLOHEXADIENE': 6.3088179708485616e-09,
        'TOLUENE': 1.9971533548975777e-09,
        'STYRENE': 2.259094316826928e-10,
        'OXYGEN': 2.9425465922880613e-05,

    },
)

# tray 5
constantTVLiquidReactor(
    temperature=(312.609,'K'),
    initialConcentrations={
        'N-BUTANE': (250.45281652859975, 'mol/m^3'),
        '2-BUTENE': (346.9376420193079, 'mol/m^3'),
        '1,3-BUTADIENE': (8488.44413461237, 'mol/m^3'),
        'CYCLOPENTADIENE': (2.8935081636689355, 'mol/m^3'),
        'BENZENE': (0.04422594738732593, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.00038510911182181696, 'mol/m^3'),
        'TOLUENE': (0.00015009528587843368, 'mol/m^3'),
        'STYRENE': (1.8835495770030587e-05, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.967410154051485,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.0308263738476881,
        '2-BUTENE': 0.033812450018390554,
        '1,3-BUTADIENE': 0.9352205476790917,
        'CYCLOPENTADIENE': 0.00013901068959559475,
        'BENZENE': 1.599374607469522e-06,
        '1,3-CYCLOHEXADIENE': 1.3209586103098423e-08,
        'TOLUENE': 4.62985030746937e-09,
        'STYRENE': 5.511901592520867e-10,

    },
)

# tray 5
constantTVLiquidReactor(
    temperature=(312.609,'K'),
    initialConcentrations={
        'N-BUTANE': (250.45281652859975, 'mol/m^3'),
        '2-BUTENE': (346.9376420193079, 'mol/m^3'),
        '1,3-BUTADIENE': (8488.44413461237, 'mol/m^3'),
        'CYCLOPENTADIENE': (2.8935081636689355, 'mol/m^3'),
        'BENZENE': (0.04422594738732593, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.00038510911182181696, 'mol/m^3'),
        'TOLUENE': (0.00015009528587843368, 'mol/m^3'),
        'STYRENE': (1.8835495770030587e-05, 'mol/m^3'),
        'OXYGEN': (0.003606425990898293, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.967410154051485,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.030825466744132994,
        '2-BUTENE': 0.03381145504591114,
        '1,3-BUTADIENE': 0.9351930276766482,
        'CYCLOPENTADIENE': 0.00013900659903692727,
        'BENZENE': 1.599327543925832e-06,
        '1,3-CYCLOHEXADIENE': 1.3209197394956011e-08,
        'TOLUENE': 4.6297140684912055e-09,
        'STYRENE': 5.511739398110499e-10,
        'OXYGEN': 2.94262166413577e-05,

    },
)

# tray 6
constantTVLiquidReactor(
    temperature=(312.622,'K'),
    initialConcentrations={
        'N-BUTANE': (246.128189682302, 'mol/m^3'),
        '2-BUTENE': (356.6198595675372, 'mol/m^3'),
        '1,3-BUTADIENE': (8481.497088963566, 'mol/m^3'),
        'CYCLOPENTADIENE': (4.476397970544081, 'mol/m^3'),
        'BENZENE': (0.08826603402457128, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.0008080377946678594, 'mol/m^3'),
        'TOLUENE': (0.00034872805208246185, 'mol/m^3'),
        'STYRENE': (4.60588217216578e-05, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.972374168985446,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.030384879565865292,
        '2-BUTENE': 0.03465001621580587,
        '1,3-BUTADIENE': 0.934748489977608,
        'CYCLOPENTADIENE': 0.0002133960491212727,
        'BENZENE': 3.178505827588453e-06,
        '1,3-CYCLOHEXADIENE': 2.7614068120036804e-08,
        'TOLUENE': 1.0727091599277584e-08,
        'STYRENE': 1.3446122473355004e-09,

    },
)

# tray 6
constantTVLiquidReactor(
    temperature=(312.622,'K'),
    initialConcentrations={
        'N-BUTANE': (246.128189682302, 'mol/m^3'),
        '2-BUTENE': (356.6198595675372, 'mol/m^3'),
        '1,3-BUTADIENE': (8481.497088963566, 'mol/m^3'),
        'CYCLOPENTADIENE': (4.476397970544081, 'mol/m^3'),
        'BENZENE': (0.08826603402457128, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.0008080377946678594, 'mol/m^3'),
        'TOLUENE': (0.00034872805208246185, 'mol/m^3'),
        'STYRENE': (4.60588217216578e-05, 'mol/m^3'),
        'OXYGEN': (0.0036063350609814235, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.972374168985446,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.030383985404925713,
        '2-BUTENE': 0.034648996541168396,
        '1,3-BUTADIENE': 0.9347209823619778,
        'CYCLOPENTADIENE': 0.0002133897693395355,
        'BENZENE': 3.1784122910729932e-06,
        '1,3-CYCLOHEXADIENE': 2.7613255498053483e-08,
        'TOLUENE': 1.0726775924295802e-08,
        'STYRENE': 1.3445726783206606e-09,
        'OXYGEN': 2.942782569334396e-05,

    },
)

# tray 7
constantTVLiquidReactor(
    temperature=(312.637,'K'),
    initialConcentrations={
        'N-BUTANE': (243.1047699464021, 'mol/m^3'),
        '2-BUTENE': (364.5562227118784, 'mol/m^3'),
        '1,3-BUTADIENE': (8474.113579713789, 'mol/m^3'),
        'CYCLOPENTADIENE': (6.870582681709341, 'mol/m^3'),
        'BENZENE': (0.17557748578185003, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.001690787246232104, 'mol/m^3'),
        'TOLUENE': (0.0008087606875069692, 'mol/m^3'),
        'STYRENE': (0.00011246484908129616, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.978442856200395,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.03007751943546281,
        '2-BUTENE': 0.03533799052950231,
        '1,3-BUTADIENE': 0.9342521746646492,
        'CYCLOPENTADIENE': 0.0003259202521934712,
        'BENZENE': 6.309310256766132e-06,
        '1,3-CYCLOHEXADIENE': 5.767999438230342e-08,
        'TOLUENE': 2.4847996413088264e-08,
        'STYRENE': 3.2799447408402502e-09,

    },
)

# tray 7
constantTVLiquidReactor(
    temperature=(312.637,'K'),
    initialConcentrations={
        'N-BUTANE': (243.1047699464021, 'mol/m^3'),
        '2-BUTENE': (364.5562227118784, 'mol/m^3'),
        '1,3-BUTADIENE': (8474.113579713789, 'mol/m^3'),
        'CYCLOPENTADIENE': (6.870582681709341, 'mol/m^3'),
        'BENZENE': (0.17557748578185003, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.001690787246232104, 'mol/m^3'),
        'TOLUENE': (0.0008087606875069692, 'mol/m^3'),
        'STYRENE': (0.00011246484908129616, 'mol/m^3'),
        'OXYGEN': (0.0036062350380728675, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.978442856200395,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.030076634255606954,
        '2-BUTENE': 0.0353369505342517,
        '1,3-BUTADIENE': 0.9342246796710185,
        'CYCLOPENTADIENE': 0.00032591066037714885,
        'BENZENE': 6.3091245740885555e-06,
        '1,3-CYCLOHEXADIENE': 5.767829686302419e-08,
        'TOLUENE': 2.4847265137827003e-08,
        'STYRENE': 3.279848212234587e-09,
        'OXYGEN': 2.9429948761419477e-05,

    },
)

# tray 8
constantTVLiquidReactor(
    temperature=(312.654,'K'),
    initialConcentrations={
        'N-BUTANE': (240.96973549831398, 'mol/m^3'),
        '2-BUTENE': (371.0295234937946, 'mol/m^3'),
        '1,3-BUTADIENE': (8465.975352153999, 'mol/m^3'),
        'CYCLOPENTADIENE': (10.491130088697416, 'mol/m^3'),
        'BENZENE': (0.34866349184148476, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.003533181942860176, 'mol/m^3'),
        'TOLUENE': (0.0018741747025738471, 'mol/m^3'),
        'STYRENE': (0.00027445194948787997, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.98672110165269,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029862888497561864,
        '2-BUTENE': 0.03590186934658502,
        '1,3-BUTADIENE': 0.9337264345900084,
        'CYCLOPENTADIENE': 0.0004961056419747389,
        'BENZENE': 1.2515942049193627e-05,
        '1,3-CYCLOHEXADIENE': 1.204310755622165e-07,
        'TOLUENE': 5.7550197487320874e-08,
        'STYRENE': 8.00054791357617e-09,

    },
)

# tray 8
constantTVLiquidReactor(
    temperature=(312.654,'K'),
    initialConcentrations={
        'N-BUTANE': (240.96973549831398, 'mol/m^3'),
        '2-BUTENE': (371.0295234937946, 'mol/m^3'),
        '1,3-BUTADIENE': (8465.975352153999, 'mol/m^3'),
        'CYCLOPENTADIENE': (10.491130088697416, 'mol/m^3'),
        'BENZENE': (0.34866349184148476, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.003533181942860176, 'mol/m^3'),
        'TOLUENE': (0.0018741747025738471, 'mol/m^3'),
        'STYRENE': (0.00027445194948787997, 'mol/m^3'),
        'OXYGEN': (0.003605971341313947, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.98672110165269,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029862009556451612,
        '2-BUTENE': 0.035900812662838365,
        '1,3-BUTADIENE': 0.933698952635302,
        'CYCLOPENTADIENE': 0.000496091040318111,
        'BENZENE': 1.2515573673040014e-05,
        '1,3-CYCLOHEXADIENE': 1.204275309679529e-07,
        'TOLUENE': 5.7548503637964045e-08,
        'STYRENE': 8.000312437009115e-09,
        'OXYGEN': 2.9432555069915414e-05,

    },
)

# tray 9
constantTVLiquidReactor(
    temperature=(312.677,'K'),
    initialConcentrations={
        'N-BUTANE': (239.42756410821286, 'mol/m^3'),
        '2-BUTENE': (376.25617511543453, 'mol/m^3'),
        '1,3-BUTADIENE': (8456.47317584117, 'mol/m^3'),
        'CYCLOPENTADIENE': (15.964747364554002, 'mol/m^3'),
        'BENZENE': (0.6917466146847085, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.007378098919722929, 'mol/m^3'),
        'TOLUENE': (0.004341585275793849, 'mol/m^3'),
        'STYRENE': (0.0006695997241320533, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.998592149616336,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02971154640640123,
        '2-BUTENE': 0.03636147472423755,
        '1,3-BUTADIENE': 0.9331483167359153,
        'CYCLOPENTADIENE': 0.000753439202019762,
        'BENZENE': 2.4818747168752272e-05,
        '1,3-CYCLOHEXADIENE': 2.5138844326770904e-07,
        'TOLUENE': 1.332810050255664e-07,
        'STYRENE': 1.951480937238415e-08,

    },
)

# tray 9
constantTVLiquidReactor(
    temperature=(312.677,'K'),
    initialConcentrations={
        'N-BUTANE': (239.42756410821286, 'mol/m^3'),
        '2-BUTENE': (376.25617511543453, 'mol/m^3'),
        '1,3-BUTADIENE': (8456.47317584117, 'mol/m^3'),
        'CYCLOPENTADIENE': (15.964747364554002, 'mol/m^3'),
        'BENZENE': (0.6917466146847085, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.007378098919722929, 'mol/m^3'),
        'TOLUENE': (0.004341585275793849, 'mol/m^3'),
        'STYRENE': (0.0006695997241320533, 'mol/m^3'),
        'OXYGEN': (0.003605389389845984, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.998592149616336,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029710671812760076,
        '2-BUTENE': 0.03636040438228558,
        '1,3-BUTADIENE': 0.9331208484388126,
        'CYCLOPENTADIENE': 0.000753417023667748,
        'BENZENE': 2.4818016600300594e-05,
        '1,3-CYCLOHEXADIENE': 2.5138104335890146e-07,
        'TOLUENE': 1.3327708174544197e-07,
        'STYRENE': 1.9514234931459673e-08,
        'OXYGEN': 2.9436153513870548e-05,

    },
)

# tray 10
constantTVLiquidReactor(
    temperature=(312.709,'K'),
    initialConcentrations={
        'N-BUTANE': (238.25820537727532, 'mol/m^3'),
        '2-BUTENE': (380.39075843547573, 'mol/m^3'),
        '1,3-BUTADIENE': (8444.543170747935, 'mol/m^3'),
        'CYCLOPENTADIENE': (24.23600539272436, 'mol/m^3'),
        'BENZENE': (1.371686888963124, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.015401072809882016, 'mol/m^3'),
        'TOLUENE': (0.010055757646727828, 'mol/m^3'),
        'STYRENE': (0.0016335468635629612, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.01654920860817,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.0296024653783191,
        '2-BUTENE': 0.03673213601056928,
        '1,3-BUTADIENE': 0.9324729020299947,
        'CYCLOPENTADIENE': 0.0011424138889626952,
        'BENZENE': 4.920179834202551e-05,
        '1,3-CYCLOHEXADIENE': 5.246506539650216e-07,
        'TOLUENE': 3.0864424892357703e-07,
        'STYRENE': 4.759890944816251e-08,

    },
)

# tray 10
constantTVLiquidReactor(
    temperature=(312.709,'K'),
    initialConcentrations={
        'N-BUTANE': (238.25820537727532, 'mol/m^3'),
        '2-BUTENE': (380.39075843547573, 'mol/m^3'),
        '1,3-BUTADIENE': (8444.543170747935, 'mol/m^3'),
        'CYCLOPENTADIENE': (24.23600539272436, 'mol/m^3'),
        'BENZENE': (1.371686888963124, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.015401072809882016, 'mol/m^3'),
        'TOLUENE': (0.010055757646727828, 'mol/m^3'),
        'STYRENE': (0.0016335468635629612, 'mol/m^3'),
        'OXYGEN': (0.0036042527658851193, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.01654920860817,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02960159384178295,
        '2-BUTENE': 0.03673105456690656,
        '1,3-BUTADIENE': 0.9324454487691699,
        'CYCLOPENTADIENE': 0.0011423802547558509,
        'BENZENE': 4.9200349774673055e-05,
        '1,3-CYCLOHEXADIENE': 5.246352075416308e-07,
        'TOLUENE': 3.0863516202030154e-07,
        'STYRENE': 4.759750807202255e-08,
        'OXYGEN': 2.944134973254521e-05,

    },
)

# tray 11
constantTVLiquidReactor(
    temperature=(312.757,'K'),
    initialConcentrations={
        'N-BUTANE': (237.2861645659439, 'mol/m^3'),
        '2-BUTENE': (383.52147547328127, 'mol/m^3'),
        '1,3-BUTADIENE': (8428.521319395586, 'mol/m^3'),
        'CYCLOPENTADIENE': (36.726866213208574, 'mol/m^3'),
        'BENZENE': (2.7189136302884407, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.03213908818749144, 'mol/m^3'),
        'TOLUENE': (0.02328887937861508, 'mol/m^3'),
        'STYRENE': (0.0039851945596168226, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.04434312603619,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029520237054064532,
        '2-BUTENE': 0.03702468946696165,
        '1,3-BUTADIENE': 0.9316255660068633,
        'CYCLOPENTADIENE': 0.0017300663315436443,
        'BENZENE': 9.75156157404171e-05,
        '1,3-CYCLOHEXADIENE': 1.0947398389644808e-06,
        'TOLUENE': 7.146887209511412e-07,
        'STYRENE': 1.1609626653499118e-07,

    },
)

# tray 11
constantTVLiquidReactor(
    temperature=(312.757,'K'),
    initialConcentrations={
        'N-BUTANE': (237.2861645659439, 'mol/m^3'),
        '2-BUTENE': (383.52147547328127, 'mol/m^3'),
        '1,3-BUTADIENE': (8428.521319395586, 'mol/m^3'),
        'CYCLOPENTADIENE': (36.726866213208574, 'mol/m^3'),
        'BENZENE': (2.7189136302884407, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.03213908818749144, 'mol/m^3'),
        'TOLUENE': (0.02328887937861508, 'mol/m^3'),
        'STYRENE': (0.0039851945596168226, 'mol/m^3'),
        'OXYGEN': (0.0036022068427555634, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.04434312603619,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029519367710814304,
        '2-BUTENE': 0.03702359912463747,
        '1,3-BUTADIENE': 0.9315981305091042,
        'CYCLOPENTADIENE': 0.0017300153827153873,
        'BENZENE': 9.751274399713876e-05,
        '1,3-CYCLOHEXADIENE': 1.094707599904611e-06,
        'TOLUENE': 7.146676740396807e-07,
        'STYRENE': 1.1609284760900141e-07,
        'OXYGEN': 2.9449060609964105e-05,

    },
)

# tray 12
constantTVLiquidReactor(
    temperature=(312.83,'K'),
    initialConcentrations={
        'N-BUTANE': (236.3550422172036, 'mol/m^3'),
        '2-BUTENE': (385.66014711804416, 'mol/m^3'),
        '1,3-BUTADIENE': (8405.734282228172, 'mol/m^3'),
        'CYCLOPENTADIENE': (55.570181956089904, 'mol/m^3'),
        'BENZENE': (5.3874248076563, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.06704917466165398, 'mol/m^3'),
        'TOLUENE': (0.053934261821696666, 'mol/m^3'),
        'STYRENE': (0.009722863221069794, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.088389028665084,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02945253330212797,
        '2-BUTENE': 0.0372450333396661,
        '1,3-BUTADIENE': 0.9304878428325718,
        'CYCLOPENTADIENE': 0.0026171538639517926,
        'BENZENE': 0.00019321502040276457,
        '1,3-CYCLOHEXADIENE': 2.283736139606457e-06,
        'TOLUENE': 1.654749461726872e-06,
        'STYRENE': 2.831556784500323e-07,

    },
)

# tray 12
constantTVLiquidReactor(
    temperature=(312.83,'K'),
    initialConcentrations={
        'N-BUTANE': (236.3550422172036, 'mol/m^3'),
        '2-BUTENE': (385.66014711804416, 'mol/m^3'),
        '1,3-BUTADIENE': (8405.734282228172, 'mol/m^3'),
        'CYCLOPENTADIENE': (55.570181956089904, 'mol/m^3'),
        'BENZENE': (5.3874248076563, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.06704917466165398, 'mol/m^3'),
        'TOLUENE': (0.053934261821696666, 'mol/m^3'),
        'STYRENE': (0.009722863221069794, 'mol/m^3'),
        'OXYGEN': (0.0035986696689893527, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.088389028665084,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029451665601639577,
        '2-BUTENE': 0.037243936064491584,
        '1,3-BUTADIENE': 0.9304604297490388,
        'CYCLOPENTADIENE': 0.002617076760034698,
        'BENZENE': 0.00019320932809895327,
        '1,3-CYCLOHEXADIENE': 2.2836688585022013e-06,
        'TOLUENE': 1.6547007111864215e-06,
        'STYRENE': 2.8314733640631187e-07,
        'OXYGEN': 2.94609797905379e-05,

    },
)

# tray 13
constantTVLiquidReactor(
    temperature=(312.944,'K'),
    initialConcentrations={
        'N-BUTANE': (235.3002551815212, 'mol/m^3'),
        '2-BUTENE': (386.7185713504013, 'mol/m^3'),
        '1,3-BUTADIENE': (8371.926539136215, 'mol/m^3'),
        'CYCLOPENTADIENE': (83.95086388962793, 'mol/m^3'),
        'BENZENE': (10.67044388476342, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.13983020756307302, 'mol/m^3'),
        'TOLUENE': (0.12490315240982999, 'mol/m^3'),
        'STYRENE': (0.023724615540252146, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.15963452777455,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02938856188234794,
        '2-BUTENE': 0.03739369970298061,
        '1,3-BUTADIENE': 0.9288712148227402,
        'CYCLOPENTADIENE': 0.003954564324163605,
        'BENZENE': 0.0003826755744172996,
        '1,3-CYCLOHEXADIENE': 4.762389466104088e-06,
        'TOLUENE': 3.830731922078825e-06,
        'STYRENE': 6.905719621890482e-07,

    },
)

# tray 13
constantTVLiquidReactor(
    temperature=(312.944,'K'),
    initialConcentrations={
        'N-BUTANE': (235.3002551815212, 'mol/m^3'),
        '2-BUTENE': (386.7185713504013, 'mol/m^3'),
        '1,3-BUTADIENE': (8371.926539136215, 'mol/m^3'),
        'CYCLOPENTADIENE': (83.95086388962793, 'mol/m^3'),
        'BENZENE': (10.67044388476342, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.13983020756307302, 'mol/m^3'),
        'TOLUENE': (0.12490315240982999, 'mol/m^3'),
        'STYRENE': (0.023724615540252146, 'mol/m^3'),
        'OXYGEN': (0.0035926501084926135, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.15963452777455,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029387695509421243,
        '2-BUTENE': 0.03739259733910239,
        '1,3-BUTADIENE': 0.9288438317586716,
        'CYCLOPENTADIENE': 0.00395444774385998,
        'BENZENE': 0.00038266429316581444,
        '1,3-CYCLOHEXADIENE': 4.762249071166889e-06,
        'TOLUENE': 3.830618992346502e-06,
        'STYRENE': 6.905516041717226e-07,
        'OXYGEN': 2.9479936111284646e-05,

    },
)

# tray 14
constantTVLiquidReactor(
    temperature=(313.124,'K'),
    initialConcentrations={
        'N-BUTANE': (233.91357394926635, 'mol/m^3'),
        '2-BUTENE': (386.465786181505, 'mol/m^3'),
        '1,3-BUTADIENE': (8320.141951479221, 'mol/m^3'),
        'CYCLOPENTADIENE': (126.58444657274096, 'mol/m^3'),
        'BENZENE': (21.121746669872138, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.29146402763493306, 'mol/m^3'),
        'TOLUENE': (0.2892626143475304, 'mol/m^3'),
        'STYRENE': (0.05790408013237091, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.27666408131816,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029317251646030176,
        '2-BUTENE': 0.03746350888384523,
        '1,3-BUTADIENE': 0.9264745087024505,
        'CYCLOPENTADIENE': 0.005966800910052171,
        'BENZENE': 0.0007574538931153558,
        '1,3-CYCLOHEXADIENE': 9.925812219321084e-06,
        'TOLUENE': 8.866095685976073e-06,
        'STYRENE': 1.6840566010071493e-06,

    },
)

# tray 14
constantTVLiquidReactor(
    temperature=(313.124,'K'),
    initialConcentrations={
        'N-BUTANE': (233.91357394926635, 'mol/m^3'),
        '2-BUTENE': (386.465786181505, 'mol/m^3'),
        '1,3-BUTADIENE': (8320.141951479221, 'mol/m^3'),
        'CYCLOPENTADIENE': (126.58444657274096, 'mol/m^3'),
        'BENZENE': (21.121746669872138, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.29146402763493306, 'mol/m^3'),
        'TOLUENE': (0.2892626143475304, 'mol/m^3'),
        'STYRENE': (0.05790408013237091, 'mol/m^3'),
        'OXYGEN': (0.003582402306861458, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.27666408131816,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02931638647979331,
        '2-BUTENE': 0.03746240331762804,
        '1,3-BUTADIENE': 0.9264471679928259,
        'CYCLOPENTADIENE': 0.0059666248268793074,
        'BENZENE': 0.000757431540285622,
        '1,3-CYCLOHEXADIENE': 9.925519303814906e-06,
        'TOLUENE': 8.865834043215935e-06,
        'STYRENE': 1.684006903684571e-06,
        'OXYGEN': 2.9510482336876562e-05,

    },
)

# tray 15
constantTVLiquidReactor(
    temperature=(313.416,'K'),
    initialConcentrations={
        'N-BUTANE': (231.89856699144553, 'mol/m^3'),
        '2-BUTENE': (384.44350483033463, 'mol/m^3'),
        '1,3-BUTADIENE': (8238.996093665175, 'mol/m^3'),
        'CYCLOPENTADIENE': (190.357234469765, 'mol/m^3'),
        'BENZENE': (41.7747496182837, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.6070608552069439, 'mol/m^3'),
        'TOLUENE': (0.6699561694061805, 'mol/m^3'),
        'STYRENE': (0.1413878370390419, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.472659647847294,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02922494044599686,
        '2-BUTENE': 0.03743660837117322,
        '1,3-BUTADIENE': 0.9228108339200742,
        'CYCLOPENTADIENE': 0.008984449485837862,
        'BENZENE': 0.0014978789795579296,
        '1,3-CYCLOHEXADIENE': 2.0669437909017324e-05,
        'TOLUENE': 2.0513094818710082e-05,
        'STYRENE': 4.106264632215548e-06,

    },
)

# tray 15
constantTVLiquidReactor(
    temperature=(313.416,'K'),
    initialConcentrations={
        'N-BUTANE': (231.89856699144553, 'mol/m^3'),
        '2-BUTENE': (384.44350483033463, 'mol/m^3'),
        '1,3-BUTADIENE': (8238.996093665175, 'mol/m^3'),
        'CYCLOPENTADIENE': (190.357234469765, 'mol/m^3'),
        'BENZENE': (41.7747496182837, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.6070608552069439, 'mol/m^3'),
        'TOLUENE': (0.6699561694061805, 'mol/m^3'),
        'STYRENE': (0.1413878370390419, 'mol/m^3'),
        'OXYGEN': (0.0035648982978641426, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.472659647847294,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029224076536592498,
        '2-BUTENE': 0.037435501719199064,
        '1,3-BUTADIENE': 0.9227835549950976,
        'CYCLOPENTADIENE': 0.008984183899311923,
        'BENZENE': 0.0014978347012217794,
        '1,3-CYCLOHEXADIENE': 2.0668826906171126e-05,
        'TOLUENE': 2.0512488437473508e-05,
        'STYRENE': 4.106143248199758e-06,
        'OXYGEN': 2.9560689985210492e-05,

    },
)

# tray 16
constantTVLiquidReactor(
    temperature=(313.895,'K'),
    initialConcentrations={
        'N-BUTANE': (228.78785453535122, 'mol/m^3'),
        '2-BUTENE': (379.82881154922427, 'mol/m^3'),
        '1,3-BUTADIENE': (8109.566449993597, 'mol/m^3'),
        'CYCLOPENTADIENE': (285.07074517986166, 'mol/m^3'),
        'BENZENE': (82.51926327843971, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (1.262816499495666, 'mol/m^3'),
        'TOLUENE': (1.552019100098063, 'mol/m^3'),
        'STYRENE': (0.3455064051277872, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.8067332412389,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02909236253271466,
        '2-BUTENE': 0.03727849870207029,
        '1,3-BUTADIENE': 0.9170858221420719,
        'CYCLOPENTADIENE': 0.013485180538262918,
        'BENZENE': 0.0029577116723867444,
        '1,3-CYCLOHEXADIENE': 4.2980570109241796e-05,
        'TOLUENE': 4.743349807177332e-05,
        'STYRENE': 1.0010344312206744e-05,

    },
)

# tray 16
constantTVLiquidReactor(
    temperature=(313.895,'K'),
    initialConcentrations={
        'N-BUTANE': (228.78785453535122, 'mol/m^3'),
        '2-BUTENE': (379.82881154922427, 'mol/m^3'),
        '1,3-BUTADIENE': (8109.566449993597, 'mol/m^3'),
        'CYCLOPENTADIENE': (285.07074517986166, 'mol/m^3'),
        'BENZENE': (82.51926327843971, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (1.262816499495666, 'mol/m^3'),
        'TOLUENE': (1.552019100098063, 'mol/m^3'),
        'STYRENE': (0.3455064051277872, 'mol/m^3'),
        'OXYGEN': (0.003534718658455265, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.8067332412389,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02909150010871571,
        '2-BUTENE': 0.037277393605435756,
        '1,3-BUTADIENE': 0.9170586357346008,
        'CYCLOPENTADIENE': 0.013484780778933874,
        'BENZENE': 0.0029576239929648557,
        '1,3-CYCLOHEXADIENE': 4.297929597844169e-05,
        'TOLUENE': 4.7432091936845455e-05,
        'STYRENE': 1.0010047562117682e-05,
        'OXYGEN': 2.964434387143251e-05,

    },
)

# tray 17
constantTVLiquidReactor(
    temperature=(314.692,'K'),
    initialConcentrations={
        'N-BUTANE': (223.81853457845105, 'mol/m^3'),
        '2-BUTENE': (371.20229033584604, 'mol/m^3'),
        '1,3-BUTADIENE': (7900.227595377393, 'mol/m^3'),
        'CYCLOPENTADIENE': (423.9980186684238, 'mol/m^3'),
        'BENZENE': (162.69180726231988, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (2.6215367823132465, 'mol/m^3'),
        'TOLUENE': (3.597205697327759, 'mol/m^3'),
        'STYRENE': (0.8455272900938496, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(44.38714239796742,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.028890237129334063,
        '2-BUTENE': 0.036927867029165534,
        '1,3-BUTADIENE': 0.9079964800711661,
        'CYCLOPENTADIENE': 0.020135789779338806,
        'BENZENE': 0.0058264829322101416,
        '1,3-CYCLOHEXADIENE': 8.916399818475741e-05,
        'TOLUENE': 0.00010958383075977381,
        'STYRENE': 2.439522984089608e-05,

    },
)

# tray 17
constantTVLiquidReactor(
    temperature=(314.692,'K'),
    initialConcentrations={
        'N-BUTANE': (223.81853457845105, 'mol/m^3'),
        '2-BUTENE': (371.20229033584604, 'mol/m^3'),
        '1,3-BUTADIENE': (7900.227595377393, 'mol/m^3'),
        'CYCLOPENTADIENE': (423.9980186684238, 'mol/m^3'),
        'BENZENE': (162.69180726231988, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (2.6215367823132465, 'mol/m^3'),
        'TOLUENE': (3.597205697327759, 'mol/m^3'),
        'STYRENE': (0.8455272900938496, 'mol/m^3'),
        'OXYGEN': (0.0034823157473635622, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(44.38714239796742,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02888937659177766,
        '2-BUTENE': 0.03692676707916117,
        '1,3-BUTADIENE': 0.9079694340809,
        'CYCLOPENTADIENE': 0.020135190005676325,
        'BENZENE': 0.005826309381977152,
        '1,3-CYCLOHEXADIENE': 8.916134230592977e-05,
        'TOLUENE': 0.00010958056664665762,
        'STYRENE': 2.4394503193642157e-05,
        'OXYGEN': 2.9786448361543692e-05,

    },
)

# tray 18
constantTVLiquidReactor(
    temperature=(316.045,'K'),
    initialConcentrations={
        'N-BUTANE': (215.70576739538365, 'mol/m^3'),
        '2-BUTENE': (356.1570262906731, 'mol/m^3'),
        '1,3-BUTADIENE': (7558.322014957625, 'mol/m^3'),
        'CYCLOPENTADIENE': (623.2700221880423, 'mol/m^3'),
        'BENZENE': (319.82143150975253, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (5.423451240719866, 'mol/m^3'),
        'TOLUENE': (8.347348182606478, 'mol/m^3'),
        'STYRENE': (2.0750297959461577, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(45.441246275762765,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02857235707978077,
        '2-BUTENE': 0.0362803041925969,
        '1,3-BUTADIENE': 0.8934143592926131,
        'CYCLOPENTADIENE': 0.02980358141816677,
        'BENZENE': 0.011432969124794688,
        '1,3-CYCLOHEXADIENE': 0.00018422393224222617,
        'TOLUENE': 0.0002527870928283748,
        'STYRENE': 5.9417866976938e-05,

    },
)

# tray 18
constantTVLiquidReactor(
    temperature=(316.045,'K'),
    initialConcentrations={
        'N-BUTANE': (215.70576739538365, 'mol/m^3'),
        '2-BUTENE': (356.1570262906731, 'mol/m^3'),
        '1,3-BUTADIENE': (7558.322014957625, 'mol/m^3'),
        'CYCLOPENTADIENE': (623.2700221880423, 'mol/m^3'),
        'BENZENE': (319.82143150975253, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (5.423451240719866, 'mol/m^3'),
        'TOLUENE': (8.347348182606478, 'mol/m^3'),
        'STYRENE': (2.0750297959461577, 'mol/m^3'),
        'OXYGEN': (0.0033902673925169007, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(45.441246275762765,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.028571499013699587,
        '2-BUTENE': 0.036279214646559255,
        '1,3-BUTADIENE': 0.8933875288650076,
        'CYCLOPENTADIENE': 0.029802686376772698,
        'BENZENE': 0.011432625777447258,
        '1,3-CYCLOHEXADIENE': 0.00018421839975125393,
        'TOLUENE': 0.0002527795012940288,
        'STYRENE': 5.9416082578957576e-05,
        'OXYGEN': 3.0031336889251014e-05,

    },
)

# tray 19
constantTVLiquidReactor(
    temperature=(318.464,'K'),
    initialConcentrations={
        'N-BUTANE': (201.39067058267014, 'mol/m^3'),
        '2-BUTENE': (330.68573797728106, 'mol/m^3'),
        '1,3-BUTADIENE': (6978.789282783654, 'mol/m^3'),
        'CYCLOPENTADIENE': (906.6331035291009, 'mol/m^3'),
        'BENZENE': (635.5319214778499, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (11.334323207825234, 'mol/m^3'),
        'TOLUENE': (19.790350827073745, 'mol/m^3'),
        'STYRENE': (5.225333137845317, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(50.086537284921604,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.028066225101647864,
        '2-BUTENE': 0.03516318278238043,
        '1,3-BUTADIENE': 0.8699487142296402,
        'CYCLOPENTADIENE': 0.04343428642957027,
        'BENZENE': 0.022283533681416185,
        '1,3-CYCLOHEXADIENE': 0.0003778790931969842,
        'TOLUENE': 0.0005816013808958141,
        'STYRENE': 0.00014457730125232832,

    },
)

# tray 19
constantTVLiquidReactor(
    temperature=(318.464,'K'),
    initialConcentrations={
        'N-BUTANE': (201.39067058267014, 'mol/m^3'),
        '2-BUTENE': (330.68573797728106, 'mol/m^3'),
        '1,3-BUTADIENE': (6978.789282783654, 'mol/m^3'),
        'CYCLOPENTADIENE': (906.6331035291009, 'mol/m^3'),
        'BENZENE': (635.5319214778499, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (11.334323207825234, 'mol/m^3'),
        'TOLUENE': (19.790350827073745, 'mol/m^3'),
        'STYRENE': (5.225333137845317, 'mol/m^3'),
        'OXYGEN': (0.003149303112813599, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(50.086537284921604,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.028065369907542586,
        '2-BUTENE': 0.03516211133987154,
        '1,3-BUTADIENE': 0.8699222063893586,
        'CYCLOPENTADIENE': 0.04343296296175132,
        'BENZENE': 0.022282854689261736,
        '1,3-CYCLOHEXADIENE': 0.0003778675790025445,
        'TOLUENE': 0.0005815836591654865,
        'STYRENE': 0.00014457289589837207,
        'OXYGEN': 3.0470578147930116e-05,

    },
)

# tray 20
constantTVLiquidReactor(
    temperature=(329.86,'K'),
    initialConcentrations={
        'N-BUTANE': (146.4853681787276, 'mol/m^3'),
        '2-BUTENE': (223.2093041354276, 'mol/m^3'),
        '1,3-BUTADIENE': (4765.973383805642, 'mol/m^3'),
        'CYCLOPENTADIENE': (1229.836218647218, 'mol/m^3'),
        'BENZENE': (2295.4984723871485, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.398215875603285, 'mol/m^3'),
        'TOLUENE': (215.12836242326443, 'mol/m^3'),
        'STYRENE': (168.72773514409542, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.796955239488124,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.027183393892693395,
        '2-BUTENE': 0.0338349576327115,
        '1,3-BUTADIENE': 0.8411011570167976,
        'CYCLOPENTADIENE': 0.05719411783483731,
        'BENZENE': 0.03854236906197508,
        '1,3-CYCLOHEXADIENE': 0.0006754497695008717,
        'TOLUENE': 0.0011653601340782741,
        'STYRENE': 0.00030319465740618007,

    },
)

# tray 20
constantTVLiquidReactor(
    temperature=(329.86,'K'),
    initialConcentrations={
        'N-BUTANE': (146.4853681787276, 'mol/m^3'),
        '2-BUTENE': (223.2093041354276, 'mol/m^3'),
        '1,3-BUTADIENE': (4765.973383805642, 'mol/m^3'),
        'CYCLOPENTADIENE': (1229.836218647218, 'mol/m^3'),
        'BENZENE': (2295.4984723871485, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.398215875603285, 'mol/m^3'),
        'TOLUENE': (215.12836242326443, 'mol/m^3'),
        'STYRENE': (168.72773514409542, 'mol/m^3'),
        'OXYGEN': (0.0002837659008727891, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.796955239488124,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.027183335943528612,
        '2-BUTENE': 0.033834885503839625,
        '1,3-BUTADIENE': 0.841099363969552,
        'CYCLOPENTADIENE': 0.05719399590924656,
        'BENZENE': 0.0385422868979118,
        '1,3-CYCLOHEXADIENE': 0.000675448329586851,
        'TOLUENE': 0.0011653576497804448,
        'STYRENE': 0.0003031940110602075,
        'OXYGEN': 2.1317854941472357e-06,

    },
)

# tray 21
constantTVLiquidReactor(
    temperature=(329.894,'K'),
    initialConcentrations={
        'N-BUTANE': (142.7363277062117, 'mol/m^3'),
        '2-BUTENE': (227.9367505134558, 'mol/m^3'),
        '1,3-BUTADIENE': (4766.264359539623, 'mol/m^3'),
        'CYCLOPENTADIENE': (1230.054450447704, 'mol/m^3'),
        'BENZENE': (2295.5075653788354, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.39803401576954, 'mol/m^3'),
        'TOLUENE': (215.1183601324088, 'mol/m^3'),
        'STYRENE': (168.7177328532398, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.795469548899575,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02647539853489381,
        '2-BUTENE': 0.03460624406274787,
        '1,3-BUTADIENE': 0.8409997854225972,
        'CYCLOPENTADIENE': 0.05721717228972626,
        'BENZENE': 0.03855645808726871,
        '1,3-CYCLOHEXADIENE': 0.0006756875006235679,
        'TOLUENE': 0.0011658844285076435,
        'STYRENE': 0.0003033696736348421,

    },
)

# tray 21
constantTVLiquidReactor(
    temperature=(329.894,'K'),
    initialConcentrations={
        'N-BUTANE': (142.7363277062117, 'mol/m^3'),
        '2-BUTENE': (227.9367505134558, 'mol/m^3'),
        '1,3-BUTADIENE': (4766.264359539623, 'mol/m^3'),
        'CYCLOPENTADIENE': (1230.054450447704, 'mol/m^3'),
        'BENZENE': (2295.5075653788354, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.39803401576954, 'mol/m^3'),
        'TOLUENE': (215.1183601324088, 'mol/m^3'),
        'STYRENE': (168.7177328532398, 'mol/m^3'),
        'OXYGEN': (7.047977856529406e-06, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.795469548899575,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02647539713372358,
        '2-BUTENE': 0.034606242231264996,
        '1,3-BUTADIENE': 0.8409997409139617,
        'CYCLOPENTADIENE': 0.05721716926159437,
        'BENZENE': 0.03855645604672676,
        '1,3-CYCLOHEXADIENE': 0.0006756874648638337,
        'TOLUENE': 0.0011658843668049818,
        'STYRENE': 0.0003033696575794631,
        'OXYGEN': 5.292348034283288e-08,

    },
)

# tray 22
constantTVLiquidReactor(
    temperature=(329.915,'K'),
    initialConcentrations={
        'N-BUTANE': (138.96637335281588, 'mol/m^3'),
        '2-BUTENE': (233.236146068591, 'mol/m^3'),
        '1,3-BUTADIENE': (4765.045898653577, 'mol/m^3'),
        'CYCLOPENTADIENE': (1230.3363331899984, 'mol/m^3'),
        'BENZENE': (2295.5984952957047, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.399397964522585, 'mol/m^3'),
        'TOLUENE': (215.11745083324013, 'mol/m^3'),
        'STYRENE': (168.7150049557337, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.79498513047084,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02576784900255335,
        '2-BUTENE': 0.035476928208897036,
        '1,3-BUTADIENE': 0.8407843741318273,
        'CYCLOPENTADIENE': 0.05725163982244832,
        'BENZENE': 0.03857325507162806,
        '1,3-CYCLOHEXADIENE': 0.0006759796862313427,
        'TOLUENE': 0.0011664338203393246,
        'STYRENE': 0.0003035402560751136,

    },
)

# tray 22
constantTVLiquidReactor(
    temperature=(329.915,'K'),
    initialConcentrations={
        'N-BUTANE': (138.96637335281588, 'mol/m^3'),
        '2-BUTENE': (233.236146068591, 'mol/m^3'),
        '1,3-BUTADIENE': (4765.045898653577, 'mol/m^3'),
        'CYCLOPENTADIENE': (1230.3363331899984, 'mol/m^3'),
        'BENZENE': (2295.5984952957047, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.399397964522585, 'mol/m^3'),
        'TOLUENE': (215.11745083324013, 'mol/m^3'),
        'STYRENE': (168.7150049557337, 'mol/m^3'),
        'OXYGEN': (1.750219039897794e-07, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.79498513047084,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02576784896869315,
        '2-BUTENE': 0.03547692816227863,
        '1,3-BUTADIENE': 0.8407843730269959,
        'CYCLOPENTADIENE': 0.05725163974721689,
        'BENZENE': 0.03857325502094094,
        '1,3-CYCLOHEXADIENE': 0.0006759796853430727,
        'TOLUENE': 0.0011664338188065741,
        'STYRENE': 0.00030354025567624703,
        'OXYGEN': 1.314048471152285e-09,

    },
)

# tray 23
constantTVLiquidReactor(
    temperature=(329.932,'K'),
    initialConcentrations={
        'N-BUTANE': (135.19096320440792, 'mol/m^3'),
        '2-BUTENE': (239.19842071770253, 'mol/m^3'),
        '1,3-BUTADIENE': (4762.6271628648565, 'mol/m^3'),
        'CYCLOPENTADIENE': (1230.7273318325358, 'mol/m^3'),
        'BENZENE': (2295.7439831626953, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.402125862028655, 'mol/m^3'),
        'TOLUENE': (215.12199732908357, 'mol/m^3'),
        'STYRENE': (168.71591425490243, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.794952836711364,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02506140377460763,
        '2-BUTENE': 0.03645973349480782,
        '1,3-BUTADIENE': 0.840437671458623,
        'CYCLOPENTADIENE': 0.057300281898212874,
        'BENZENE': 0.038593763726793434,
        '1,3-CYCLOHEXADIENE': 0.0006763419678295413,
        'TOLUENE': 0.0011670759111517507,
        'STYRENE': 0.00030372776797397565,

    },
)

# tray 23
constantTVLiquidReactor(
    temperature=(329.932,'K'),
    initialConcentrations={
        'N-BUTANE': (135.19096320440792, 'mol/m^3'),
        '2-BUTENE': (239.19842071770253, 'mol/m^3'),
        '1,3-BUTADIENE': (4762.6271628648565, 'mol/m^3'),
        'CYCLOPENTADIENE': (1230.7273318325358, 'mol/m^3'),
        'BENZENE': (2295.7439831626953, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.402125862028655, 'mol/m^3'),
        'TOLUENE': (215.12199732908357, 'mol/m^3'),
        'STYRENE': (168.71591425490243, 'mol/m^3'),
        'OXYGEN': (4.345749865986457e-09, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.794952836711364,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02506140377378995,
        '2-BUTENE': 0.03645973349361825,
        '1,3-BUTADIENE': 0.840437671431202,
        'CYCLOPENTADIENE': 0.057300281896343334,
        'BENZENE': 0.03859376372553423,
        '1,3-CYCLOHEXADIENE': 0.0006763419678074743,
        'TOLUENE': 0.0011670759111136724,
        'STYRENE': 0.00030372776796406593,
        'OXYGEN': 3.262701999761505e-11,

    },
)

# tray 24
constantTVLiquidReactor(
    temperature=(329.948,'K'),
    initialConcentrations={
        'N-BUTANE': (131.41646235516868, 'mol/m^3'),
        '2-BUTENE': (245.91814157433424, 'mol/m^3'),
        '1,3-BUTADIENE': (4759.01724516515, 'mol/m^3'),
        'CYCLOPENTADIENE': (1231.3001903088116, 'mol/m^3'),
        'BENZENE': (2295.9622149631814, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.406126778370904, 'mol/m^3'),
        'TOLUENE': (215.1301810216018, 'mol/m^3'),
        'STYRENE': (168.71955145157716, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.795243484157606,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.024356102704351547,
        '2-BUTENE': 0.037569155164882045,
        '1,3-BUTADIENE': 0.8399364760631881,
        'CYCLOPENTADIENE': 0.0573700329322262,
        'BENZENE': 0.03861960703619545,
        '1,3-CYCLOHEXADIENE': 0.0006768044206022894,
        'TOLUENE': 0.0011678676693055105,
        'STYRENE': 0.00030395400924897674,

    },
)

# tray 24
constantTVLiquidReactor(
    temperature=(329.948,'K'),
    initialConcentrations={
        'N-BUTANE': (131.41646235516868, 'mol/m^3'),
        '2-BUTENE': (245.91814157433424, 'mol/m^3'),
        '1,3-BUTADIENE': (4759.01724516515, 'mol/m^3'),
        'CYCLOPENTADIENE': (1231.3001903088116, 'mol/m^3'),
        'BENZENE': (2295.9622149631814, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.406126778370904, 'mol/m^3'),
        'TOLUENE': (215.1301810216018, 'mol/m^3'),
        'STYRENE': (168.71955145157716, 'mol/m^3'),
        'OXYGEN': (1.0789198356194558e-10, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.795243484157606,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.024356102704331813,
        '2-BUTENE': 0.03756915516485161,
        '1,3-BUTADIENE': 0.8399364760625077,
        'CYCLOPENTADIENE': 0.057370032932179724,
        'BENZENE': 0.03861960703616416,
        '1,3-CYCLOHEXADIENE': 0.0006768044206017411,
        'TOLUENE': 0.0011678676693045644,
        'STYRENE': 0.0003039540092487305,
        'OXYGEN': 8.100720013547708e-13,

    },
)

# tray 25
constantTVLiquidReactor(
    temperature=(329.967,'K'),
    initialConcentrations={
        'N-BUTANE': (127.6437801042668, 'mol/m^3'),
        '2-BUTENE': (253.48987575203006, 'mol/m^3'),
        '1,3-BUTADIENE': (4754.061564695781, 'mol/m^3'),
        'CYCLOPENTADIENE': (1232.1731175107554, 'mol/m^3'),
        'BENZENE': (2296.2804696722237, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.411946293050526, 'mol/m^3'),
        'TOLUENE': (215.14382050913218, 'mol/m^3'),
        'STYRENE': (168.7268258449267, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.795824803424757,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02365203884074887,
        '2-BUTENE': 0.0388208844173888,
        '1,3-BUTADIENE': 0.8392491197734466,
        'CYCLOPENTADIENE': 0.05747340935921202,
        'BENZENE': 0.03865398061498637,
        '1,3-CYCLOHEXADIENE': 0.0006774264331289996,
        'TOLUENE': 0.0011688966298030203,
        'STYRENE': 0.00030424393128524255,

    },
)

# tray 25
constantTVLiquidReactor(
    temperature=(329.967,'K'),
    initialConcentrations={
        'N-BUTANE': (127.6437801042668, 'mol/m^3'),
        '2-BUTENE': (253.48987575203006, 'mol/m^3'),
        '1,3-BUTADIENE': (4754.061564695781, 'mol/m^3'),
        'CYCLOPENTADIENE': (1232.1731175107554, 'mol/m^3'),
        'BENZENE': (2296.2804696722237, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.411946293050526, 'mol/m^3'),
        'TOLUENE': (215.14382050913218, 'mol/m^3'),
        'STYRENE': (168.7268258449267, 'mol/m^3'),
        'OXYGEN': (2.678340701381417e-12, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.795824803424757,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.023652038840748395,
        '2-BUTENE': 0.03882088441738802,
        '1,3-BUTADIENE': 0.8392491197734298,
        'CYCLOPENTADIENE': 0.057473409359210866,
        'BENZENE': 0.038653980614985595,
        '1,3-CYCLOHEXADIENE': 0.000677426433128986,
        'TOLUENE': 0.0011688966298029969,
        'STYRENE': 0.0003042439312852365,
        'OXYGEN': 2.0111258174683333e-14,

    },
)

# tray 26
constantTVLiquidReactor(
    temperature=(329.989,'K'),
    initialConcentrations={
        'N-BUTANE': (123.8701885541962, 'mol/m^3'),
        '2-BUTENE': (262.0145554585148, 'mol/m^3'),
        '1,3-BUTADIENE': (4747.450959739392, 'mol/m^3'),
        'CYCLOPENTADIENE': (1233.5552522471671, 'mol/m^3'),
        'BENZENE': (2296.753305239943, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.42085742490371, 'mol/m^3'),
        'TOLUENE': (215.16382509084337, 'mol/m^3'),
        'STYRENE': (168.7368281357823, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.796761441720484,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.022948942409678673,
        '2-BUTENE': 0.040232324974048696,
        '1,3-BUTADIENE': 0.8383308552888243,
        'CYCLOPENTADIENE': 0.05763220876902122,
        'BENZENE': 0.03870240881354429,
        '1,3-CYCLOHEXADIENE': 0.0006783101649902727,
        'TOLUENE': 0.0011703123620094796,
        'STYRENE': 0.00030463721788293247,

    },
)

# tray 26
constantTVLiquidReactor(
    temperature=(329.989,'K'),
    initialConcentrations={
        'N-BUTANE': (123.8701885541962, 'mol/m^3'),
        '2-BUTENE': (262.0145554585148, 'mol/m^3'),
        '1,3-BUTADIENE': (4747.450959739392, 'mol/m^3'),
        'CYCLOPENTADIENE': (1233.5552522471671, 'mol/m^3'),
        'BENZENE': (2296.753305239943, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.42085742490371, 'mol/m^3'),
        'TOLUENE': (215.16382509084337, 'mol/m^3'),
        'STYRENE': (168.7368281357823, 'mol/m^3'),
        'OXYGEN': (6.648077175130481e-14, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.796761441720484,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02294894240967866,
        '2-BUTENE': 0.040232324974048675,
        '1,3-BUTADIENE': 0.838330855288824,
        'CYCLOPENTADIENE': 0.05763220876902119,
        'BENZENE': 0.038702408813544274,
        '1,3-CYCLOHEXADIENE': 0.0006783101649902724,
        'TOLUENE': 0.0011703123620094792,
        'STYRENE': 0.0003046372178829323,
        'OXYGEN': 4.992592736757076e-16,

    },
)

# tray 27
constantTVLiquidReactor(
    temperature=(330.018,'K'),
    initialConcentrations={
        'N-BUTANE': (120.0893226107761, 'mol/m^3'),
        '2-BUTENE': (271.5931129015132, 'mol/m^3'),
        '1,3-BUTADIENE': (4738.621664811395, 'mol/m^3'),
        'CYCLOPENTADIENE': (1235.7830352104618, 'mol/m^3'),
        'BENZENE': (2297.49893055827, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.43504249193529, 'mol/m^3'),
        'TOLUENE': (215.19565056174758, 'mol/m^3'),
        'STYRENE': (168.75410481998748, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.798279620907447,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02224621819415727,
        '2-BUTENE': 0.04182212220093525,
        '1,3-BUTADIENE': 0.8371164475996861,
        'CYCLOPENTADIENE': 0.057883407534544226,
        'BENZENE': 0.038774606009481724,
        '1,3-CYCLOHEXADIENE': 0.0006796376109146268,
        'TOLUENE': 0.0011723562287286433,
        'STYRENE': 0.00030520462155215426,

    },
)

# tray 27
constantTVLiquidReactor(
    temperature=(330.018,'K'),
    initialConcentrations={
        'N-BUTANE': (120.0893226107761, 'mol/m^3'),
        '2-BUTENE': (271.5931129015132, 'mol/m^3'),
        '1,3-BUTADIENE': (4738.621664811395, 'mol/m^3'),
        'CYCLOPENTADIENE': (1235.7830352104618, 'mol/m^3'),
        'BENZENE': (2297.49893055827, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.43504249193529, 'mol/m^3'),
        'TOLUENE': (215.19565056174758, 'mol/m^3'),
        'STYRENE': (168.75410481998748, 'mol/m^3'),
        'OXYGEN': (1.6505325720341222e-15, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.798279620907447,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02224621819415727,
        '2-BUTENE': 0.04182212220093525,
        '1,3-BUTADIENE': 0.8371164475996861,
        'CYCLOPENTADIENE': 0.057883407534544226,
        'BENZENE': 0.038774606009481724,
        '1,3-CYCLOHEXADIENE': 0.0006796376109146268,
        'TOLUENE': 0.0011723562287286433,
        'STYRENE': 0.00030520462155215426,
        'OXYGEN': 1.2397465867736156e-17,

    },
)

# tray 28
constantTVLiquidReactor(
    temperature=(330.059,'K'),
    initialConcentrations={
        'N-BUTANE': (116.28845208564476, 'mol/m^3'),
        '2-BUTENE': (282.31920589540056, 'mol/m^3'),
        '1,3-BUTADIENE': (4726.7189386932205, 'mol/m^3'),
        'CYCLOPENTADIENE': (1239.4475108602894, 'mol/m^3'),
        'BENZENE': (2298.7264844360043, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.45822962073694, 'mol/m^3'),
        'TOLUENE': (215.247480614363, 'mol/m^3'),
        'STYRENE': (168.7813837950482, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.800767329731357,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.021542451680341582,
        '2-BUTENE': 0.04360920461848188,
        '1,3-BUTADIENE': 0.8355080043837025,
        'CYCLOPENTADIENE': 0.05828943983636617,
        'BENZENE': 0.03888759329133631,
        '1,3-CYCLOHEXADIENE': 0.0006817356354837896,
        'TOLUENE': 0.001175502820031265,
        'STYRENE': 0.00030606773425649445,

    },
)

# tray 28
constantTVLiquidReactor(
    temperature=(330.059,'K'),
    initialConcentrations={
        'N-BUTANE': (116.28845208564476, 'mol/m^3'),
        '2-BUTENE': (282.31920589540056, 'mol/m^3'),
        '1,3-BUTADIENE': (4726.7189386932205, 'mol/m^3'),
        'CYCLOPENTADIENE': (1239.4475108602894, 'mol/m^3'),
        'BENZENE': (2298.7264844360043, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.45822962073694, 'mol/m^3'),
        'TOLUENE': (215.247480614363, 'mol/m^3'),
        'STYRENE': (168.7813837950482, 'mol/m^3'),
        'OXYGEN': (4.0975202859253083e-17, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.800767329731357,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.021542451680341582,
        '2-BUTENE': 0.04360920461848188,
        '1,3-BUTADIENE': 0.8355080043837025,
        'CYCLOPENTADIENE': 0.05828943983636617,
        'BENZENE': 0.03888759329133631,
        '1,3-CYCLOHEXADIENE': 0.0006817356354837896,
        'TOLUENE': 0.001175502820031265,
        'STYRENE': 0.00030606773425649445,
        'OXYGEN': 3.0793673874023627e-19,

    },
)

# tray 29
constantTVLiquidReactor(
    temperature=(330.119,'K'),
    initialConcentrations={
        'N-BUTANE': (112.44302590124755, 'mol/m^3'),
        '2-BUTENE': (294.2655783736722, 'mol/m^3'),
        '1,3-BUTADIENE': (4710.3242746817095, 'mol/m^3'),
        'CYCLOPENTADIENE': (1245.512536315463, 'mol/m^3'),
        'BENZENE': (2300.7905935489343, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.49778413457503, 'mol/m^3'),
        'TOLUENE': (215.33295473622005, 'mol/m^3'),
        'STYRENE': (168.8259394543141, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.804904058968,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.020835312605364135,
        '2-BUTENE': 0.0456117275950952,
        '1,3-BUTADIENE': 0.8333515041776602,
        'CYCLOPENTADIENE': 0.05895593566834109,
        'BENZENE': 0.03907232363875581,
        '1,3-CYCLOHEXADIENE': 0.0006851984145450411,
        'TOLUENE': 0.0011805507142331826,
        'STYRENE': 0.00030744718600554763,

    },
)

# tray 29
constantTVLiquidReactor(
    temperature=(330.119,'K'),
    initialConcentrations={
        'N-BUTANE': (112.44302590124755, 'mol/m^3'),
        '2-BUTENE': (294.2655783736722, 'mol/m^3'),
        '1,3-BUTADIENE': (4710.3242746817095, 'mol/m^3'),
        'CYCLOPENTADIENE': (1245.512536315463, 'mol/m^3'),
        'BENZENE': (2300.7905935489343, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.49778413457503, 'mol/m^3'),
        'TOLUENE': (215.33295473622005, 'mol/m^3'),
        'STYRENE': (168.8259394543141, 'mol/m^3'),
        'OXYGEN': (1.0228524488779251e-18, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.804904058968,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.020835312605364135,
        '2-BUTENE': 0.0456117275950952,
        '1,3-BUTADIENE': 0.8333515041776602,
        'CYCLOPENTADIENE': 0.05895593566834109,
        'BENZENE': 0.03907232363875581,
        '1,3-CYCLOHEXADIENE': 0.0006851984145450411,
        'TOLUENE': 0.0011805507142331826,
        'STYRENE': 0.00030744718600554763,
        'OXYGEN': 7.828574736287717e-21,

    },
)

# tray 30
constantTVLiquidReactor(
    temperature=(330.21,'K'),
    initialConcentrations={
        'N-BUTANE': (108.5148534924994, 'mol/m^3'),
        '2-BUTENE': (307.44586982385846, 'mol/m^3'),
        '1,3-BUTADIENE': (4687.109866905011, 'mol/m^3'),
        'CYCLOPENTADIENE': (1255.6239430713147, 'mol/m^3'),
        'BENZENE': (2304.355046290206, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.56661808164498, 'mol/m^3'),
        'TOLUENE': (215.4802612015481, 'mol/m^3'),
        'STYRENE': (168.9032298836529, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.81185617554187,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.020121106639965188,
        '2-BUTENE': 0.047844115788558204,
        '1,3-BUTADIENE': 0.8303992740317603,
        'CYCLOPENTADIENE': 0.060060719820037536,
        'BENZENE': 0.039384912997021286,
        '1,3-CYCLOHEXADIENE': 0.0006911132280673653,
        'TOLUENE': 0.0011890103923734292,
        'STYRENE': 0.00030974710221654363,

    },
)

# tray 30
constantTVLiquidReactor(
    temperature=(330.21,'K'),
    initialConcentrations={
        'N-BUTANE': (108.5148534924994, 'mol/m^3'),
        '2-BUTENE': (307.44586982385846, 'mol/m^3'),
        '1,3-BUTADIENE': (4687.109866905011, 'mol/m^3'),
        'CYCLOPENTADIENE': (1255.6239430713147, 'mol/m^3'),
        'BENZENE': (2304.355046290206, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.56661808164498, 'mol/m^3'),
        'TOLUENE': (215.4802612015481, 'mol/m^3'),
        'STYRENE': (168.9032298836529, 'mol/m^3'),
        'OXYGEN': (0.0, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.81185617554187,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.020121106639965188,
        '2-BUTENE': 0.047844115788558204,
        '1,3-BUTADIENE': 0.8303992740317603,
        'CYCLOPENTADIENE': 0.060060719820037536,
        'BENZENE': 0.039384912997021286,
        '1,3-CYCLOHEXADIENE': 0.0006911132280673653,
        'TOLUENE': 0.0011890103923734292,
        'STYRENE': 0.00030974710221654363,
        'OXYGEN': 5.400711782234887e-22,

    },
)

# tray 31
constantTVLiquidReactor(
    temperature=(330.356,'K'),
    initialConcentrations={
        'N-BUTANE': (104.43664672091708, 'mol/m^3'),
        '2-BUTENE': (321.75914803823457, 'mol/m^3'),
        '1,3-BUTADIENE': (4653.19300791281, 'mol/m^3'),
        'CYCLOPENTADIENE': (1272.5005356422328, 'mol/m^3'),
        'BENZENE': (2310.6292105541784, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.68891881983402, 'mol/m^3'),
        'TOLUENE': (215.73850216545654, 'mol/m^3'),
        'STYRENE': (169.03689686145057, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.823831167392985,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.019392798177076974,
        '2-BUTENE': 0.05031099527076644,
        '1,3-BUTADIENE': 0.8262469223327892,
        'CYCLOPENTADIENE': 0.06190229418118434,
        'BENZENE': 0.03992809624675895,
        '1,3-CYCLOHEXADIENE': 0.0007014939340595701,
        'TOLUENE': 0.0012036698868550305,
        'STYRENE': 0.00031372997050938277,

    },
)

# tray 31
constantTVLiquidReactor(
    temperature=(330.356,'K'),
    initialConcentrations={
        'N-BUTANE': (104.43664672091708, 'mol/m^3'),
        '2-BUTENE': (321.75914803823457, 'mol/m^3'),
        '1,3-BUTADIENE': (4653.19300791281, 'mol/m^3'),
        'CYCLOPENTADIENE': (1272.5005356422328, 'mol/m^3'),
        'BENZENE': (2310.6292105541784, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.68891881983402, 'mol/m^3'),
        'TOLUENE': (215.73850216545654, 'mol/m^3'),
        'STYRENE': (169.03689686145057, 'mol/m^3'),
        'OXYGEN': (1.0726274853721091e-19, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.823831167392985,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.019392798177076974,
        '2-BUTENE': 0.05031099527076644,
        '1,3-BUTADIENE': 0.8262469223327892,
        'CYCLOPENTADIENE': 0.06190229418118434,
        'BENZENE': 0.03992809624675895,
        '1,3-CYCLOHEXADIENE': 0.0007014939340595701,
        'TOLUENE': 0.0012036698868550305,
        'STYRENE': 0.00031372997050938277,
        'OXYGEN': 1.350179873083092e-21,

    },
)

# tray 32
constantTVLiquidReactor(
    temperature=(330.595,'K'),
    initialConcentrations={
        'N-BUTANE': (100.09838038708897, 'mol/m^3'),
        '2-BUTENE': (336.87715601690286, 'mol/m^3'),
        '1,3-BUTADIENE': (4602.117673607397, 'mol/m^3'),
        'CYCLOPENTADIENE': (1300.6524379049279, 'mol/m^3'),
        'BENZENE': (2321.859055287521, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.909787587909236, 'mol/m^3'),
        'TOLUENE': (216.19860754481454, 'mol/m^3'),
        'STYRENE': (169.2724053461417, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.844479889616185,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.01863919528428359,
        '2-BUTENE': 0.052997686591585276,
        '1,3-BUTADIENE': 0.8202227924836333,
        'CYCLOPENTADIENE': 0.06497788356059545,
        'BENZENE': 0.04089158965442781,
        '1,3-CYCLOHEXADIENE': 0.0007201168178104449,
        'TOLUENE': 0.0012298696888429683,
        'STYRENE': 0.00032086591882092243,

    },
)

# tray 32
constantTVLiquidReactor(
    temperature=(330.595,'K'),
    initialConcentrations={
        'N-BUTANE': (100.09838038708897, 'mol/m^3'),
        '2-BUTENE': (336.87715601690286, 'mol/m^3'),
        '1,3-BUTADIENE': (4602.117673607397, 'mol/m^3'),
        'CYCLOPENTADIENE': (1300.6524379049279, 'mol/m^3'),
        'BENZENE': (2321.859055287521, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.909787587909236, 'mol/m^3'),
        'TOLUENE': (216.19860754481454, 'mol/m^3'),
        'STYRENE': (169.2724053461417, 'mol/m^3'),
        'OXYGEN': (1.5929648416642504e-19, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.844479889616185,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.01863919528428359,
        '2-BUTENE': 0.052997686591585276,
        '1,3-BUTADIENE': 0.8202227924836333,
        'CYCLOPENTADIENE': 0.06497788356059545,
        'BENZENE': 0.04089158965442781,
        '1,3-CYCLOHEXADIENE': 0.0007201168178104449,
        'TOLUENE': 0.0012298696888429683,
        'STYRENE': 0.00032086591882092243,
        'OXYGEN': 1.3749096521478575e-21,

    },
)

# tray 33
constantTVLiquidReactor(
    temperature=(330.992,'K'),
    initialConcentrations={
        'N-BUTANE': (95.32183185395144, 'mol/m^3'),
        '2-BUTENE': (352.0388103556683, 'mol/m^3'),
        '1,3-BUTADIENE': (4523.090482856401, 'mol/m^3'),
        'CYCLOPENTADIENE': (1347.3358572255602, 'mol/m^3'),
        'BENZENE': (2342.163705724407, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (46.31360734872522, 'mol/m^3'),
        'TOLUENE': (217.0342534808422, 'mol/m^3'),
        'STYRENE': (169.6915922629086, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.880266162144167,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.017841408314096273,
        '2-BUTENE': 0.05584982602601892,
        '1,3-BUTADIENE': 0.8112053780217061,
        'CYCLOPENTADIENE': 0.07010783267025002,
        'BENZENE': 0.04262941986530965,
        '1,3-CYCLOHEXADIENE': 0.0007541153514177536,
        'TOLUENE': 0.0012779705955342975,
        'STYRENE': 0.0003340491556669065,

    },
)

# tray 33
constantTVLiquidReactor(
    temperature=(330.992,'K'),
    initialConcentrations={
        'N-BUTANE': (95.32183185395144, 'mol/m^3'),
        '2-BUTENE': (352.0388103556683, 'mol/m^3'),
        '1,3-BUTADIENE': (4523.090482856401, 'mol/m^3'),
        'CYCLOPENTADIENE': (1347.3358572255602, 'mol/m^3'),
        'BENZENE': (2342.163705724407, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (46.31360734872522, 'mol/m^3'),
        'TOLUENE': (217.0342534808422, 'mol/m^3'),
        'STYRENE': (169.6915922629086, 'mol/m^3'),
        'OXYGEN': (0.0, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.880266162144167,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.017841408314096273,
        '2-BUTENE': 0.05584982602601892,
        '1,3-BUTADIENE': 0.8112053780217061,
        'CYCLOPENTADIENE': 0.07010783267025002,
        'BENZENE': 0.04262941986530965,
        '1,3-CYCLOHEXADIENE': 0.0007541153514177536,
        'TOLUENE': 0.0012779705955342975,
        'STYRENE': 0.0003340491556669065,
        'OXYGEN': 4.581242134858834e-22,

    },
)

# tray 34
constantTVLiquidReactor(
    temperature=(331.664,'K'),
    initialConcentrations={
        'N-BUTANE': (89.81557073793881, 'mol/m^3'),
        '2-BUTENE': (365.6983024677551, 'mol/m^3'),
        '1,3-BUTADIENE': (4398.425566828768, 'mol/m^3'),
        'CYCLOPENTADIENE': (1423.8533822709671, 'mol/m^3'),
        'BENZENE': (2379.126716931725, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (47.057232208881295, 'mol/m^3'),
        'TOLUENE': (218.5691504775938, 'mol/m^3'),
        'STYRENE': (170.44631057292273, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.94218151257688,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.016968292975126708,
        '2-BUTENE': 0.058738775682146874,
        '1,3-BUTADIENE': 0.7973236699080007,
        'CYCLOPENTADIENE': 0.07862036745116788,
        'BENZENE': 0.04580428103702765,
        '1,3-CYCLOHEXADIENE': 0.0008170566617385421,
        'TOLUENE': 0.0013684394334660746,
        'STYRENE': 0.0003591168513256236,

    },
)

# tray 34
constantTVLiquidReactor(
    temperature=(331.664,'K'),
    initialConcentrations={
        'N-BUTANE': (89.81557073793881, 'mol/m^3'),
        '2-BUTENE': (365.6983024677551, 'mol/m^3'),
        '1,3-BUTADIENE': (4398.425566828768, 'mol/m^3'),
        'CYCLOPENTADIENE': (1423.8533822709671, 'mol/m^3'),
        'BENZENE': (2379.126716931725, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (47.057232208881295, 'mol/m^3'),
        'TOLUENE': (218.5691504775938, 'mol/m^3'),
        'STYRENE': (170.44631057292273, 'mol/m^3'),
        'OXYGEN': (9.855348109948175e-20, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.94218151257688,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.016968292975126708,
        '2-BUTENE': 0.058738775682146874,
        '1,3-BUTADIENE': 0.7973236699080007,
        'CYCLOPENTADIENE': 0.07862036745116788,
        'BENZENE': 0.04580428103702765,
        '1,3-CYCLOHEXADIENE': 0.0008170566617385421,
        'TOLUENE': 0.0013684394334660746,
        'STYRENE': 0.0003591168513256236,
        'OXYGEN': 1.1453095258418563e-21,

    },
)

# tray 35
constantTVLiquidReactor(
    temperature=(332.812,'K'),
    initialConcentrations={
        'N-BUTANE': (83.13731392339946, 'mol/m^3'),
        '2-BUTENE': (374.95133080836194, 'mol/m^3'),
        '1,3-BUTADIENE': (4200.107418137107, 'mol/m^3'),
        'CYCLOPENTADIENE': (1546.6996999612163, 'mol/m^3'),
        'BENZENE': (2446.442134389971, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (48.42918279460334, 'mol/m^3'),
        'TOLUENE': (221.41525687559886, 'mol/m^3'),
        'STYRENE': (171.81116862512897, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.048178705956946,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.01596969386763755,
        '2-BUTENE': 0.061398676422908235,
        '1,3-BUTADIENE': 0.7754867022131062,
        'CYCLOPENTADIENE': 0.09260216444076882,
        'BENZENE': 0.0516575801634892,
        '1,3-CYCLOHEXADIENE': 0.0009347606410519136,
        'TOLUENE': 0.001542349407737827,
        'STYRENE': 0.00040807284330002805,

    },
)

# tray 35
constantTVLiquidReactor(
    temperature=(332.812,'K'),
    initialConcentrations={
        'N-BUTANE': (83.13731392339946, 'mol/m^3'),
        '2-BUTENE': (374.95133080836194, 'mol/m^3'),
        '1,3-BUTADIENE': (4200.107418137107, 'mol/m^3'),
        'CYCLOPENTADIENE': (1546.6996999612163, 'mol/m^3'),
        'BENZENE': (2446.442134389971, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (48.42918279460334, 'mol/m^3'),
        'TOLUENE': (221.41525687559886, 'mol/m^3'),
        'STYRENE': (171.81116862512897, 'mol/m^3'),
        'OXYGEN': (6.086284869737684e-20, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.048178705956946,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.01596969386763755,
        '2-BUTENE': 0.061398676422908235,
        '1,3-BUTADIENE': 0.7754867022131062,
        'CYCLOPENTADIENE': 0.09260216444076882,
        'BENZENE': 0.0516575801634892,
        '1,3-CYCLOHEXADIENE': 0.0009347606410519136,
        'TOLUENE': 0.001542349407737827,
        'STYRENE': 0.00040807284330002805,
        'OXYGEN': 9.950066179174586e-22,

    },
)

# tray 36
constantTVLiquidReactor(
    temperature=(334.782,'K'),
    initialConcentrations={
        'N-BUTANE': (74.66046335318791, 'mol/m^3'),
        '2-BUTENE': (374.7549221879246, 'mol/m^3'),
        '1,3-BUTADIENE': (3887.026621364866, 'mol/m^3'),
        'CYCLOPENTADIENE': (1736.6977612593482, 'mol/m^3'),
        'BENZENE': (2567.9517823022443, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (50.939212219859925, 'mol/m^3'),
        'TOLUENE': (226.6955571481916, 'mol/m^3'),
        'STYRENE': (174.26809497893396, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.224598579560872,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.01476699295614436,
        '2-BUTENE': 0.0633214697956589,
        '1,3-BUTADIENE': 0.7407286466724355,
        'CYCLOPENTADIENE': 0.1151429450768152,
        'BENZENE': 0.06249527018975611,
        '1,3-CYCLOHEXADIENE': 0.001155979448597803,
        'TOLUENE': 0.0018827291019382184,
        'STYRENE': 0.000505966758653856,

    },
)

# tray 36
constantTVLiquidReactor(
    temperature=(334.782,'K'),
    initialConcentrations={
        'N-BUTANE': (74.66046335318791, 'mol/m^3'),
        '2-BUTENE': (374.7549221879246, 'mol/m^3'),
        '1,3-BUTADIENE': (3887.026621364866, 'mol/m^3'),
        'CYCLOPENTADIENE': (1736.6977612593482, 'mol/m^3'),
        'BENZENE': (2567.9517823022443, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (50.939212219859925, 'mol/m^3'),
        'TOLUENE': (226.6955571481916, 'mol/m^3'),
        'STYRENE': (174.26809497893396, 'mol/m^3'),
        'OXYGEN': (1.5191479351498577e-19, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.224598579560872,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.01476699295614436,
        '2-BUTENE': 0.0633214697956589,
        '1,3-BUTADIENE': 0.7407286466724355,
        'CYCLOPENTADIENE': 0.1151429450768152,
        'BENZENE': 0.06249527018975611,
        '1,3-CYCLOHEXADIENE': 0.001155979448597803,
        'TOLUENE': 0.0018827291019382184,
        'STYRENE': 0.000505966758653856,
        'OXYGEN': 1.3160593722396795e-21,

    },
)

# tray 37
constantTVLiquidReactor(
    temperature=(338.152,'K'),
    initialConcentrations={
        'N-BUTANE': (63.66121688902557, 'mol/m^3'),
        '2-BUTENE': (357.38276157006953, 'mol/m^3'),
        '1,3-BUTADIENE': (3408.6534217078465, 'mol/m^3'),
        'CYCLOPENTADIENE': (2010.9878554952008, 'mol/m^3'),
        'BENZENE': (2781.9007837037157, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (55.420329453089515, 'mol/m^3'),
        'TOLUENE': (236.38141189309562, 'mol/m^3'),
        'STYRENE': (178.60817991109948, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.500814991279945,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.013245006119192828,
        '2-BUTENE': 0.06360262938441477,
        '1,3-BUTADIENE': 0.6855413167200882,
        'CYCLOPENTADIENE': 0.15028106942985406,
        'BENZENE': 0.0824964381133544,
        '1,3-CYCLOHEXADIENE': 0.0015705507255944352,
        'TOLUENE': 0.002557521181574786,
        'STYRENE': 0.0007054683259263666,

    },
)

# tray 37
constantTVLiquidReactor(
    temperature=(338.152,'K'),
    initialConcentrations={
        'N-BUTANE': (63.66121688902557, 'mol/m^3'),
        '2-BUTENE': (357.38276157006953, 'mol/m^3'),
        '1,3-BUTADIENE': (3408.6534217078465, 'mol/m^3'),
        'CYCLOPENTADIENE': (2010.9878554952008, 'mol/m^3'),
        'BENZENE': (2781.9007837037157, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (55.420329453089515, 'mol/m^3'),
        'TOLUENE': (236.38141189309562, 'mol/m^3'),
        'STYRENE': (178.60817991109948, 'mol/m^3'),
        'OXYGEN': (0.0, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.500814991279945,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.013245006119192828,
        '2-BUTENE': 0.06360262938441477,
        '1,3-BUTADIENE': 0.6855413167200882,
        'CYCLOPENTADIENE': 0.15028106942985406,
        'BENZENE': 0.0824964381133544,
        '1,3-CYCLOHEXADIENE': 0.0015705507255944352,
        'TOLUENE': 0.002557521181574786,
        'STYRENE': 0.0007054683259263666,
        'OXYGEN': 2.8844313326072757e-22,

    },
)

# tray 38
constantTVLiquidReactor(
    temperature=(343.811,'K'),
    initialConcentrations={
        'N-BUTANE': (49.72493410998902, 'mol/m^3'),
        '2-BUTENE': (313.7682269437719, 'mol/m^3'),
        '1,3-BUTADIENE': (2728.815898235494, 'mol/m^3'),
        'CYCLOPENTADIENE': (2357.4581177418017, 'mol/m^3'),
        'BENZENE': (3140.4647248939154, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (63.0278900180318, 'mol/m^3'),
        'TOLUENE': (253.70174245833525, 'mol/m^3'),
        'STYRENE': (186.03169832429853, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.887291609197888,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.011258999211870053,
        '2-BUTENE': 0.06079199574456029,
        '1,3-BUTADIENE': 0.6000489579965729,
        'CYCLOPENTADIENE': 0.20172398587932097,
        'BENZENE': 0.11882299168239058,
        '1,3-CYCLOHEXADIENE': 0.0023349098365563112,
        'TOLUENE': 0.003901199726916019,
        'STYRENE': 0.0011169599218128054,

    },
)

# tray 38
constantTVLiquidReactor(
    temperature=(343.811,'K'),
    initialConcentrations={
        'N-BUTANE': (49.72493410998902, 'mol/m^3'),
        '2-BUTENE': (313.7682269437719, 'mol/m^3'),
        '1,3-BUTADIENE': (2728.815898235494, 'mol/m^3'),
        'CYCLOPENTADIENE': (2357.4581177418017, 'mol/m^3'),
        'BENZENE': (3140.4647248939154, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (63.0278900180318, 'mol/m^3'),
        'TOLUENE': (253.70174245833525, 'mol/m^3'),
        'STYRENE': (186.03169832429853, 'mol/m^3'),
        'OXYGEN': (0.0, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.887291609197888,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.011258999211870053,
        '2-BUTENE': 0.06079199574456029,
        '1,3-BUTADIENE': 0.6000489579965729,
        'CYCLOPENTADIENE': 0.20172398587932097,
        'BENZENE': 0.11882299168239058,
        '1,3-CYCLOHEXADIENE': 0.0023349098365563112,
        'TOLUENE': 0.003901199726916019,
        'STYRENE': 0.0011169599218128054,
        'OXYGEN': 7.211069495225135e-22,

    },
)

# tray 39
constantTVLiquidReactor(
    temperature=(352.989,'K'),
    initialConcentrations={
        'N-BUTANE': (33.501763921692834, 'mol/m^3'),
        '2-BUTENE': (238.7437711333567, 'mol/m^3'),
        '1,3-BUTADIENE': (1874.6839101083144, 'mol/m^3'),
        'CYCLOPENTADIENE': (2675.5127809667856, 'mol/m^3'),
        'BENZENE': (3707.0763158807777, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (75.10902063310353, 'mol/m^3'),
        'TOLUENE': (286.79659500203854, 'mol/m^3'),
        'STYRENE': (201.56707462139633, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(22.792439765938695,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.00869749147857355,
        '2-BUTENE': 0.05303390901576453,
        '1,3-BUTADIENE': 0.47564908086034374,
        'CYCLOPENTADIENE': 0.26813604558312776,
        'BENZENE': 0.18227503098675524,
        '1,3-CYCLOHEXADIENE': 0.003687780626922706,
        'TOLUENE': 0.006553171114039088,
        'STYRENE': 0.0019674903344733567,

    },
)

# tray 39
constantTVLiquidReactor(
    temperature=(352.989,'K'),
    initialConcentrations={
        'N-BUTANE': (33.501763921692834, 'mol/m^3'),
        '2-BUTENE': (238.7437711333567, 'mol/m^3'),
        '1,3-BUTADIENE': (1874.6839101083144, 'mol/m^3'),
        'CYCLOPENTADIENE': (2675.5127809667856, 'mol/m^3'),
        'BENZENE': (3707.0763158807777, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (75.10902063310353, 'mol/m^3'),
        'TOLUENE': (286.79659500203854, 'mol/m^3'),
        'STYRENE': (201.56707462139633, 'mol/m^3'),
        'OXYGEN': (1.28038415943479e-19, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(22.792439765938695,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.00869749147857355,
        '2-BUTENE': 0.05303390901576453,
        '1,3-BUTADIENE': 0.47564908086034374,
        'CYCLOPENTADIENE': 0.26813604558312776,
        'BENZENE': 0.18227503098675524,
        '1,3-CYCLOHEXADIENE': 0.003687780626922706,
        'TOLUENE': 0.006553171114039088,
        'STYRENE': 0.0019674903344733567,
        'OXYGEN': 1.8027703064709523e-21,

    },
)

# tray 40
constantTVLiquidReactor(
    temperature=(376.855,'K'),
    initialConcentrations={
        'N-BUTANE': (11.04634816110057, 'mol/m^3'),
        '2-BUTENE': (91.8392160378622, 'mol/m^3'),
        '1,3-BUTADIENE': (626.0142870791539, 'mol/m^3'),
        'CYCLOPENTADIENE': (2251.142858938367, 'mol/m^3'),
        'BENZENE': (5094.121267803104, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (101.88242535606209, 'mol/m^3'),
        'TOLUENE': (509.4121267803105, 'mol/m^3'),
        'STYRENE': (407.52970142424834, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(22.792439765938695,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.005699799772008009,
        '2-BUTENE': 0.03944089842236406,
        '1,3-BUTADIENE': 0.31824098727036054,
        'CYCLOPENTADIENE': 0.3323269867069205,
        'BENZENE': 0.28319298867228043,
        '1,3-CYCLOHEXADIENE': 0.00585708976571641,
        'TOLUENE': 0.011559899537604019,
        'STYRENE': 0.0036813498527460057,

    },
)

# tray 40
constantTVLiquidReactor(
    temperature=(376.855,'K'),
    initialConcentrations={
        'N-BUTANE': (11.04634816110057, 'mol/m^3'),
        '2-BUTENE': (91.8392160378622, 'mol/m^3'),
        '1,3-BUTADIENE': (626.0142870791539, 'mol/m^3'),
        'CYCLOPENTADIENE': (2251.142858938367, 'mol/m^3'),
        'BENZENE': (5094.121267803104, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (101.88242535606209, 'mol/m^3'),
        'TOLUENE': (509.4121267803105, 'mol/m^3'),
        'STYRENE': (407.52970142424834, 'mol/m^3'),
        'OXYGEN': (6.939653046563163e-20, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(22.792439765938695,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.005699799772008009,
        '2-BUTENE': 0.03944089842236406,
        '1,3-BUTADIENE': 0.31824098727036054,
        'CYCLOPENTADIENE': 0.3323269867069205,
        'BENZENE': 0.28319298867228043,
        '1,3-CYCLOHEXADIENE': 0.00585708976571641,
        'TOLUENE': 0.011559899537604019,
        'STYRENE': 0.0036813498527460057,
        'OXYGEN': 1.3909499443620023e-21,

    },
)


liquidVolumetricMassTransferCoefficientPowerLaw(
    prefactor=(14.5, "1/s"),
    diffusionCoefficientPower=1/2,
    solventViscosityPower=-1/6,
    solventDensityPower=-1/6,
)

solvation(
    solvent='benzene'
)

simulator(
    atol=1e-20,
    rtol=1e-6,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=0.01,
    toleranceInterruptSimulation=0.01,
    maximumEdgeSpecies=100000,
    filterReactions=True,
    toleranceBranchReactionToCore=0.01,
    branchingIndex=0.5,
    branchingRatioMax=1.0,
    maxNumObjsPerIter=10,
    toleranceReactionToCoreDeadendRadical=0.01,
)

options(
    units='si',
)

generatedSpeciesConstraints(
    #allows exceptions to the following restrictions
    allowed=['input species','seed mechanisms','reaction libraries'],
    # Constraints on generated species
    maximumRadicalElectrons = 1,
    maximumCarbonAtoms=16,
)
