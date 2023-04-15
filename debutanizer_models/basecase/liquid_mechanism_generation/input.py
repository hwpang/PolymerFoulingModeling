
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

# tray 1
constantTVLiquidReactor(
    temperature=(312.458,'K'),
    initialConcentrations={
        'N-BUTANE': (296.53703130726603, 'mol/m^3'),
        '2-BUTENE': (282.000153491539, 'mol/m^3'),
        '1,3-BUTADIENE': (8501.977600303433, 'mol/m^3'),
        'CYCLOPENTADIENE': (0.39646385857688005, 'mol/m^3'),
        'BENZENE': (0.002313434034374856, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (1.667035354094377e-05, 'mol/m^3'),
        'TOLUENE': (4.412611202563451e-06, 'mol/m^3'),
        'STYRENE': (4.5569979435629687e-07, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.95295388694515,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.035550012136590196,
        '2-BUTENE': 0.028219024422026593,
        '1,3-BUTADIENE': 0.9362091566870978,
        'CYCLOPENTADIENE': 2.1709180384481632e-05,
        'BENZENE': 9.673801038204599e-08,
        '1,3-CYCLOHEXADIENE': 6.623015959693516e-10,
        'TOLUENE': 1.5808345387120728e-10,
        'STYRENE': 1.5505496600949123e-11,

    },
)

# tray 2
constantTVLiquidReactor(
    temperature=(312.539,'K'),
    initialConcentrations={
        'N-BUTANE': (278.3393371979181, 'mol/m^3'),
        '2-BUTENE': (303.3467315613462, 'mol/m^3'),
        '1,3-BUTADIENE': (8504.123533351854, 'mol/m^3'),
        'CYCLOPENTADIENE': (0.6986221433186132, 'mol/m^3'),
        'BENZENE': (0.005162696639416604, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (3.894359278432994e-05, 'mol/m^3'),
        'TOLUENE': (1.131534129146735e-05, 'mol/m^3'),
        'STYRENE': (1.2245821498437603e-06, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.95557189406627,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.03367969703691793,
        '2-BUTENE': 0.03005064394437299,
        '1,3-BUTADIENE': 0.9362335685678159,
        'CYCLOPENTADIENE': 3.588979688383436e-05,
        'BENZENE': 1.9882581311556697e-07,
        '1,3-CYCLOHEXADIENE': 1.4204062697853981e-09,
        'TOLUENE': 3.6988182937883374e-10,
        'STYRENE': 3.790819775060338e-11,

    },
)

# tray 3
constantTVLiquidReactor(
    temperature=(312.574,'K'),
    initialConcentrations={
        'N-BUTANE': (265.56921697248055, 'mol/m^3'),
        '2-BUTENE': (320.81699073902575, 'mol/m^3'),
        '1,3-BUTADIENE': (8500.568195123325, 'mol/m^3'),
        'CYCLOPENTADIENE': (1.1558031770545196, 'mol/m^3'),
        'BENZENE': (0.01081286561241745, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (8.544287034516592e-05, 'mol/m^3'),
        'TOLUENE': (2.7303542725239142e-05, 'mol/m^3'),
        'STYRENE': (3.1001912847207396e-06, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.95984406456748,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.03237221900122265,
        '2-BUTENE': 0.031556645469975206,
        '1,3-BUTADIENE': 0.9360133674997134,
        'CYCLOPENTADIENE': 5.736270621194937e-05,
        'BENZENE': 4.0136595214386863e-07,
        '1,3-CYCLOHEXADIENE': 3.003748701076387e-09,
        'TOLUENE': 8.606053706858609e-10,
        'STYRENE': 9.257060085749577e-11,

    },
)

# tray 4
constantTVLiquidReactor(
    temperature=(312.595,'K'),
    initialConcentrations={
        'N-BUTANE': (256.65632047221436, 'mol/m^3'),
        '2-BUTENE': (335.16018900419607, 'mol/m^3'),
        '1,3-BUTADIENE': (8494.894202656313, 'mol/m^3'),
        'CYCLOPENTADIENE': (1.8475392394512273, 'mol/m^3'),
        'BENZENE': (0.02201663657124017, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.00018251342435558426, 'mol/m^3'),
        'TOLUENE': (6.433489093294596e-05, 'mol/m^3'),
        'STYRENE': (7.675447845068583e-06, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.964254938454715,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.031460676556310706,
        '2-BUTENE': 0.03279517100396579,
        '1,3-BUTADIENE': 0.9356534791223895,
        'CYCLOPENTADIENE': 8.9861735206941e-05,
        'BENZENE': 8.030493610709459e-07,
        '1,3-CYCLOHEXADIENE': 6.30955264723086e-09,
        'TOLUENE': 1.997294692086369e-09,
        'STYRENE': 2.259185721936904e-10,

    },
)

# tray 5
constantTVLiquidReactor(
    temperature=(312.61,'K'),
    initialConcentrations={
        'N-BUTANE': (250.44857262408973, 'mol/m^3'),
        '2-BUTENE': (346.9482720718813, 'mol/m^3'),
        '1,3-BUTADIENE': (8488.38366001788, 'mol/m^3'),
        'CYCLOPENTADIENE': (2.8940725968325176, 'mol/m^3'),
        'BENZENE': (0.044231499161359954, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.0003851440617464563, 'mol/m^3'),
        'TOLUENE': (0.0001501025609751018, 'mol/m^3'),
        'STYRENE': (1.883601826188685e-05, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.96866671820215,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.030826053797853824,
        '2-BUTENE': 0.03381364321777727,
        '1,3-BUTADIENE': 0.9352196484426352,
        'CYCLOPENTADIENE': 0.00013903658264902225,
        'BENZENE': 1.5995671935816107e-06,
        '1,3-CYCLOHEXADIENE': 1.3210693759489875e-08,
        'TOLUENE': 4.630002431602389e-09,
        'STYRENE': 5.511947241287402e-10,

    },
)

# tray 6
constantTVLiquidReactor(
    temperature=(312.624,'K'),
    initialConcentrations={
        'N-BUTANE': (246.1239719553215, 'mol/m^3'),
        '2-BUTENE': (356.63043101238685, 'mol/m^3'),
        '1,3-BUTADIENE': (8481.436656420448, 'mol/m^3'),
        'CYCLOPENTADIENE': (4.477216517432244, 'mol/m^3'),
        'BENZENE': (0.08827632033338052, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.0008081056470030199, 'mol/m^3'),
        'TOLUENE': (0.0003487441270594376, 'mol/m^3'),
        'STYRENE': (4.605999779147477e-05, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.97363105355151,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.030384602677682667,
        '2-BUTENE': 0.0346512185569107,
        '1,3-BUTADIENE': 0.9347475252696202,
        'CYCLOPENTADIENE': 0.0002134349530307756,
        'BENZENE': 3.178854674930114e-06,
        '1,3-CYCLOHEXADIENE': 2.7616074644314896e-08,
        'TOLUENE': 1.0727382615997624e-08,
        'STYRENE': 1.3446234785510156e-09,

    },
)

# tray 7
constantTVLiquidReactor(
    temperature=(312.638,'K'),
    initialConcentrations={
        'N-BUTANE': (243.10147981423958, 'mol/m^3'),
        '2-BUTENE': (364.5658368232237, 'mol/m^3'),
        '1,3-BUTADIENE': (8474.053191864015, 'mol/m^3'),
        'CYCLOPENTADIENE': (6.871759546688119, 'mol/m^3'),
        'BENZENE': (0.17559642744908302, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.0016909224986630356, 'mol/m^3'),
        'TOLUENE': (0.000808793982307111, 'mol/m^3'),
        'STYRENE': (0.00011246780549157609, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.97970013254082,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.030077250844218795,
        '2-BUTENE': 0.03533917230676937,
        '1,3-BUTADIENE': 0.934251206179535,
        'CYCLOPENTADIENE': 0.0003259749279179341,
        'BENZENE': 6.309929083104356e-06,
        '1,3-CYCLOHEXADIENE': 5.768383531565433e-08,
        'TOLUENE': 2.48486802228915e-08,
        'STYRENE': 3.2799600892659948e-09,

    },
)

# tray 8
constantTVLiquidReactor(
    temperature=(312.655,'K'),
    initialConcentrations={
        'N-BUTANE': (240.96736758345762, 'mol/m^3'),
        '2-BUTENE': (371.0390984214073, 'mol/m^3'),
        '1,3-BUTADIENE': (8465.924106502622, 'mol/m^3'),
        'CYCLOPENTADIENE': (10.492794242483383, 'mol/m^3'),
        'BENZENE': (0.3486977530825437, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.003533433344149106, 'mol/m^3'),
        'TOLUENE': (0.0018742451943801009, 'mol/m^3'),
        'STYRENE': (0.00027445756254380297, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.98797891250683,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02986259836498997,
        '2-BUTENE': 0.03590290131192615,
        '1,3-BUTADIENE': 0.9337256152836877,
        'CYCLOPENTADIENE': 0.0004961819780452667,
        'BENZENE': 1.2517070887700167e-05,
        '1,3-CYCLOHEXADIENE': 1.2043849892785912e-07,
        'TOLUENE': 5.7551384731119864e-08,
        'STYRENE': 8.00057960542448e-09,

    },
)

# tray 9
constantTVLiquidReactor(
    temperature=(312.678,'K'),
    initialConcentrations={
        'N-BUTANE': (239.42520552832076, 'mol/m^3'),
        '2-BUTENE': (376.26480911177964, 'mol/m^3'),
        '1,3-BUTADIENE': (8456.412894771058, 'mol/m^3'),
        'CYCLOPENTADIENE': (15.96710582075341, 'mol/m^3'),
        'BENZENE': (0.6918088058898746, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.0073785543706448, 'mol/m^3'),
        'TOLUENE': (0.004341731761409017, 'mol/m^3'),
        'STYRENE': (0.0006696129475373447, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(42.999850727162496,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02971125962278032,
        '2-BUTENE': 0.0363624915475712,
        '1,3-BUTADIENE': 0.9331474780319091,
        'CYCLOPENTADIENE': 0.0007535458881433515,
        'BENZENE': 2.4820709772311813e-05,
        '1,3-CYCLOHEXADIENE': 2.5140196451327086e-07,
        'TOLUENE': 1.3328303211082063e-07,
        'STYRENE': 1.9514827239123467e-08,

    },
)

# tray 10
constantTVLiquidReactor(
    temperature=(312.71,'K'),
    initialConcentrations={
        'N-BUTANE': (238.25585387566392, 'mol/m^3'),
        '2-BUTENE': (380.39845811098496, 'mol/m^3'),
        '1,3-BUTADIENE': (8444.482961891697, 'mol/m^3'),
        'CYCLOPENTADIENE': (24.239314004952753, 'mol/m^3'),
        'BENZENE': (1.371796794140139, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.015401888878847116, 'mol/m^3'),
        'TOLUENE': (0.01005606049538301, 'mol/m^3'),
        'STYRENE': (0.0016335733472376335, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.0178089463516,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029602307111076677,
        '2-BUTENE': 0.03673317291131709,
        '1,3-BUTADIENE': 0.9324718727780389,
        'CYCLOPENTADIENE': 0.0011425610003132402,
        'BENZENE': 4.9205277531678595e-05,
        '1,3-CYCLOHEXADIENE': 5.246734817745992e-07,
        'TOLUENE': 3.0864931059723185e-07,
        'STYRENE': 4.7598929827042985e-08,

    },
)

# tray 11
constantTVLiquidReactor(
    temperature=(312.758,'K'),
    initialConcentrations={
        'N-BUTANE': (237.28381894822212, 'mol/m^3'),
        '2-BUTENE': (383.52824690448796, 'mol/m^3'),
        '1,3-BUTADIENE': (8428.461207521701, 'mol/m^3'),
        'CYCLOPENTADIENE': (36.73128129839363, 'mol/m^3'),
        'BENZENE': (2.719115402829638, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.03214062130336771, 'mol/m^3'),
        'TOLUENE': (0.023289556772293288, 'mol/m^3'),
        'STYRENE': (0.003985261366081598, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.04560466057648,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029520081938117424,
        '2-BUTENE': 0.03702561868107829,
        '1,3-BUTADIENE': 0.9316245854720262,
        'CYCLOPENTADIENE': 0.0017302665354400769,
        'BENZENE': 9.752180120090163e-05,
        '1,3-CYCLOHEXADIENE': 1.0947781947262448e-06,
        'TOLUENE': 7.146970368784534e-07,
        'STYRENE': 1.1609690546483547e-07,

    },
)

# tray 12
constantTVLiquidReactor(
    temperature=(312.831,'K'),
    initialConcentrations={
        'N-BUTANE': (236.35270223568662, 'mol/m^3'),
        '2-BUTENE': (385.6659963099283, 'mol/m^3'),
        '1,3-BUTADIENE': (8405.67430828719, 'mol/m^3'),
        'CYCLOPENTADIENE': (55.576210638193416, 'mol/m^3'),
        'BENZENE': (5.387765007274786, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.06705195133209363, 'mol/m^3'),
        'TOLUENE': (0.0539354811498155, 'mol/m^3'),
        'STYRENE': (0.009723077155411831, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.089792073528486,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029452377285762905,
        '2-BUTENE': 0.03724595254121206,
        '1,3-BUTADIENE': 0.9304868048119173,
        'CYCLOPENTADIENE': 0.0026174182222198873,
        'BENZENE': 0.0001932254221582767,
        '1,3-CYCLOHEXADIENE': 2.2837998800398153e-06,
        'TOLUENE': 1.6547613991882648e-06,
        'STYRENE': 2.83155450675509e-07,

    },
)

# tray 13
constantTVLiquidReactor(
    temperature=(312.945,'K'),
    initialConcentrations={
        'N-BUTANE': (235.29792158476755, 'mol/m^3'),
        '2-BUTENE': (386.72532342917026, 'mol/m^3'),
        '1,3-BUTADIENE': (8371.85767690126, 'mol/m^3'),
        'CYCLOPENTADIENE': (83.95890308456254, 'mol/m^3'),
        'BENZENE': (10.671015800742124, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.13983481691460328, 'mol/m^3'),
        'TOLUENE': (0.12490603352943996, 'mol/m^3'),
        'STYRENE': (0.023724926578896785, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.160903529125335,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02938847426166059,
        '2-BUTENE': 0.03739447164961311,
        '1,3-BUTADIENE': 0.9288701628066134,
        'CYCLOPENTADIENE': 0.003954915437628141,
        'BENZENE': 0.0003826920433227189,
        '1,3-CYCLOHEXADIENE': 4.762495811639939e-06,
        'TOLUENE': 3.830736246738672e-06,
        'STYRENE': 6.905691036739546e-07,

    },
)

# tray 14
constantTVLiquidReactor(
    temperature=(313.126,'K'),
    initialConcentrations={
        'N-BUTANE': (233.9112487462748, 'mol/m^3'),
        '2-BUTENE': (386.471630496751, 'mol/m^3'),
        '1,3-BUTADIENE': (8320.073402703121, 'mol/m^3'),
        'CYCLOPENTADIENE': (126.59459186457508, 'mol/m^3'),
        'BENZENE': (21.12270996958625, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.29147226559552525, 'mol/m^3'),
        'TOLUENE': (0.28926631916524104, 'mol/m^3'),
        'STYRENE': (0.0579048207838302, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.27808055911667,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029317193155650866,
        '2-BUTENE': 0.037464213035878303,
        '1,3-BUTADIENE': 0.9264733800584011,
        'CYCLOPENTADIENE': 0.00596725689887344,
        'BENZENE': 0.0007574807560215492,
        '1,3-CYCLOHEXADIENE': 9.925951807016406e-06,
        'TOLUENE': 8.866096089703705e-06,
        'STYRENE': 1.684047278072041e-06,

    },
)

# tray 15
constantTVLiquidReactor(
    temperature=(313.417,'K'),
    initialConcentrations={
        'N-BUTANE': (231.8962539855535, 'mol/m^3'),
        '2-BUTENE': (384.44845209304845, 'mol/m^3'),
        '1,3-BUTADIENE': (8238.928036075518, 'mol/m^3'),
        'CYCLOPENTADIENE': (190.37063091126683, 'mol/m^3'),
        'BENZENE': (41.77622440824642, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (0.6070717292871599, 'mol/m^3'),
        'TOLUENE': (0.6699639348918746, 'mol/m^3'),
        'STYRENE': (0.14138879978737975, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.47409017692117,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029224780658760028,
        '2-BUTENE': 0.03743730837967463,
        '1,3-BUTADIENE': 0.9228096911105067,
        'CYCLOPENTADIENE': 0.008985016105843702,
        'BENZENE': 0.0014979149799062806,
        '1,3-CYCLOHEXADIENE': 2.0669530134662868e-05,
        'TOLUENE': 2.051299134408007e-05,
        'STYRENE': 4.106243829901883e-06,

    },
)

# tray 16
constantTVLiquidReactor(
    temperature=(313.896,'K'),
    initialConcentrations={
        'N-BUTANE': (228.78556035900684, 'mol/m^3'),
        '2-BUTENE': (379.83378674527745, 'mol/m^3'),
        '1,3-BUTADIENE': (8109.499175858431, 'mol/m^3'),
        'CYCLOPENTADIENE': (285.0862961891419, 'mol/m^3'),
        'BENZENE': (82.52149165958589, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (1.2628270413762228, 'mol/m^3'),
        'TOLUENE': (1.5520369843342596, 'mol/m^3'),
        'STYRENE': (0.3455097694945072, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(43.808044557038464,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.029092197003129056,
        '2-BUTENE': 0.03727919100618878,
        '1,3-BUTADIENE': 0.9170845438487676,
        'CYCLOPENTADIENE': 0.013485880524729366,
        'BENZENE': 0.002957763782437599,
        '1,3-CYCLOHEXADIENE': 4.2980491643759014e-05,
        'TOLUENE': 4.743308665737489e-05,
        'STYRENE': 1.0010256446351699e-05,

    },
)

# tray 17
constantTVLiquidReactor(
    temperature=(314.694,'K'),
    initialConcentrations={
        'N-BUTANE': (223.8162704820475, 'mol/m^3'),
        '2-BUTENE': (371.206408455691, 'mol/m^3'),
        '1,3-BUTADIENE': (7900.15249546093, 'mol/m^3'),
        'CYCLOPENTADIENE': (424.0145473201535, 'mol/m^3'),
        'BENZENE': (162.69445964228058, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (2.6215390996829773, 'mol/m^3'),
        'TOLUENE': (3.5972202947197727, 'mol/m^3'),
        'STYRENE': (0.8455339928246847, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(44.3884922526428,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.028890166720635244,
        '2-BUTENE': 0.03692855565208122,
        '1,3-BUTADIENE': 0.907995125150626,
        'CYCLOPENTADIENE': 0.02013648426246537,
        'BENZENE': 0.005826527481013035,
        '1,3-CYCLOHEXADIENE': 8.916321096933101e-05,
        'TOLUENE': 0.00010958252762040855,
        'STYRENE': 2.4394994589586668e-05,

    },
)

# tray 18
constantTVLiquidReactor(
    temperature=(316.046,'K'),
    initialConcentrations={
        'N-BUTANE': (215.70264311295188, 'mol/m^3'),
        '2-BUTENE': (356.16032618779616, 'mol/m^3'),
        '1,3-BUTADIENE': (7558.2489846401695, 'mol/m^3'),
        'CYCLOPENTADIENE': (623.286253914047, 'mol/m^3'),
        'BENZENE': (319.8231327636339, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (5.423427504799826, 'mol/m^3'),
        'TOLUENE': (8.347343119704037, 'mol/m^3'),
        'STYRENE': (2.075035421393852, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(45.44266753791055,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.028572196602983692,
        '2-BUTENE': 0.03628091283079546,
        '1,3-BUTADIENE': 0.8934134597983757,
        'CYCLOPENTADIENE': 0.02980408489050314,
        'BENZENE': 0.011432924843373621,
        '1,3-CYCLOHEXADIENE': 0.00018422038034939268,
        'TOLUENE': 0.0002527834765989916,
        'STYRENE': 5.941717701968112e-05,

    },
)

# tray 19
constantTVLiquidReactor(
    temperature=(318.466,'K'),
    initialConcentrations={
        'N-BUTANE': (201.38854224504828, 'mol/m^3'),
        '2-BUTENE': (330.688282761765, 'mol/m^3'),
        '1,3-BUTADIENE': (6978.728853389947, 'mol/m^3'),
        'CYCLOPENTADIENE': (906.6421642594829, 'mol/m^3'),
        'BENZENE': (635.5298931076896, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (11.334163670324266, 'mol/m^3'),
        'TOLUENE': (19.79023103350278, 'mol/m^3'),
        'STYRENE': (5.225356065843185, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(50.088107621787735,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02806617155324948,
        '2-BUTENE': 0.0351638144790856,
        '1,3-BUTADIENE': 0.8699482592577275,
        'CYCLOPENTADIENE': 0.043434465913076335,
        'BENZENE': 0.02228325563275365,
        '1,3-CYCLOHEXADIENE': 0.00037786958761622946,
        'TOLUENE': 0.000581588675390468,
        'STYRENE': 0.00014457490110069508,

    },
)

# tray 20
constantTVLiquidReactor(
    temperature=(329.861,'K'),
    initialConcentrations={
        'N-BUTANE': (146.48357219005186, 'mol/m^3'),
        '2-BUTENE': (223.2079530204399, 'mol/m^3'),
        '1,3-BUTADIENE': (4765.871791254884, 'mol/m^3'),
        'CYCLOPENTADIENE': (1229.856053098356, 'mol/m^3'),
        'BENZENE': (2295.530042119169, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.39875943858825, 'mol/m^3'),
        'TOLUENE': (215.13069739792766, 'mol/m^3'),
        'STYRENE': (168.72944169413398, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.797184655947174,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.027183345643206052,
        '2-BUTENE': 0.03383539813021173,
        '1,3-BUTADIENE': 0.8410976618812868,
        'CYCLOPENTADIENE': 0.05719601853049768,
        'BENZENE': 0.0385435019785629,
        '1,3-CYCLOHEXADIENE': 0.0006754715960532313,
        'TOLUENE': 0.0011653981098692463,
        'STYRENE': 0.00030320413031241776,

    },
)

# tray 21
constantTVLiquidReactor(
    temperature=(329.894,'K'),
    initialConcentrations={
        'N-BUTANE': (142.7336451173016, 'mol/m^3'),
        '2-BUTENE': (227.93628007628405, 'mol/m^3'),
        '1,3-BUTADIENE': (4766.162765227551, 'mol/m^3'),
        'CYCLOPENTADIENE': (1230.0742835778562, 'mol/m^3'),
        'BENZENE': (2295.5391350558148, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.398577579855335, 'mol/m^3'),
        'TOLUENE': (215.1206951676172, 'mol/m^3'),
        'STYRENE': (168.71943946382353, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.79569892358887,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.026475283650160757,
        '2-BUTENE': 0.034606824979107255,
        '1,3-BUTADIENE': 0.8409963126940846,
        'CYCLOPENTADIENE': 0.05721904119993079,
        'BENZENE': 0.0385575370451568,
        '1,3-CYCLOHEXADIENE': 0.0006757081327463426,
        'TOLUENE': 0.001165913790566225,
        'STYRENE': 0.0003033785082472222,

    },
)

# tray 22
constantTVLiquidReactor(
    temperature=(329.915,'K'),
    initialConcentrations={
        'N-BUTANE': (138.96371358393046, 'mol/m^3'),
        '2-BUTENE': (233.23655284715244, 'mol/m^3'),
        '1,3-BUTADIENE': (4764.944311717007, 'mol/m^3'),
        'CYCLOPENTADIENE': (1230.3561646138778, 'mol/m^3'),
        'BENZENE': (2295.6209714856277, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.40003244971867, 'mol/m^3'),
        'TOLUENE': (215.11978587395262, 'mol/m^3'),
        'STYRENE': (168.7167115828298, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.795182196874293,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.025767754907048055,
        '2-BUTENE': 0.035477736558536525,
        '1,3-BUTADIENE': 0.8407805666175437,
        'CYCLOPENTADIENE': 0.05725355343998733,
        'BENZENE': 0.038574364204100674,
        '1,3-CYCLOHEXADIENE': 0.0006760008466789526,
        'TOLUENE': 0.0011664740980835678,
        'STYRENE': 0.0003035493280210984,

    },
)

# tray 23
constantTVLiquidReactor(
    temperature=(329.932,'K'),
    initialConcentrations={
        'N-BUTANE': (135.18741699490724, 'mol/m^3'),
        '2-BUTENE': (239.20060999316814, 'mol/m^3'),
        '1,3-BUTADIENE': (4762.516497632564, 'mol/m^3'),
        'CYCLOPENTADIENE': (1230.7471608896494, 'mol/m^3'),
        'BENZENE': (2295.775551408607, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.40266940134597, 'mol/m^3'),
        'TOLUENE': (215.12433234227555, 'mol/m^3'),
        'STYRENE': (168.71762087649438, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.795182196874293,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02506120447988145,
        '2-BUTENE': 0.03646083467097352,
        '1,3-BUTADIENE': 0.8404336954022393,
        'CYCLOPENTADIENE': 0.05730218376516102,
        'BENZENE': 0.03859486496485268,
        '1,3-CYCLOHEXADIENE': 0.0006763629897120348,
        'TOLUENE': 0.0011671159493176707,
        'STYRENE': 0.000303737777862298,

    },
)

# tray 24
constantTVLiquidReactor(
    temperature=(329.949,'K'),
    initialConcentrations={
        'N-BUTANE': (131.41293899321317, 'mol/m^3'),
        '2-BUTENE': (245.92210876178353, 'mol/m^3'),
        '1,3-BUTADIENE': (4758.90660178416, 'mol/m^3'),
        'CYCLOPENTADIENE': (1231.3200158983384, 'mol/m^3'),
        'BENZENE': (2295.9937818881076, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.40667029347014, 'mol/m^3'),
        'TOLUENE': (215.13251598525684, 'mol/m^3'),
        'STYRENE': (168.72125805115272, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.795472852491475,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02435582288151105,
        '2-BUTENE': 0.037570586381571625,
        '1,3-BUTADIENE': 0.8399321721544316,
        'CYCLOPENTADIENE': 0.05737198058977342,
        'BENZENE': 0.0386207391086899,
        '1,3-CYCLOHEXADIENE': 0.0006768269828059796,
        'TOLUENE': 0.0011679086395248817,
        'STYRENE': 0.0003039632616914062,

    },
)

# tray 25
constantTVLiquidReactor(
    temperature=(329.967,'K'),
    initialConcentrations={
        'N-BUTANE': (127.64027957884831, 'mol/m^3'),
        '2-BUTENE': (253.49561569411546, 'mol/m^3'),
        '1,3-BUTADIENE': (4753.960044248815, 'mol/m^3'),
        'CYCLOPENTADIENE': (1232.1929378163402, 'mol/m^3'),
        'BENZENE': (2296.3029417340667, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.41248977292349, 'mol/m^3'),
        'TOLUENE': (215.1461553902256, 'mol/m^3'),
        'STYRENE': (168.7285324004694, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.7960541881016,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.023651746615516032,
        '2-BUTENE': 0.03882259722904456,
        '1,3-BUTADIENE': 0.8392443957772308,
        'CYCLOPENTADIENE': 0.057475428314623773,
        'BENZENE': 0.03865519336017374,
        '1,3-CYCLOHEXADIENE': 0.0006774486565159007,
        'TOLUENE': 0.0011689370153033605,
        'STYRENE': 0.0003042530315918131,

    },
)

# tray 26
constantTVLiquidReactor(
    temperature=(329.989,'K'),
    initialConcentrations={
        'N-BUTANE': (123.86671087081884, 'mol/m^3'),
        '2-BUTENE': (262.0220623869329, 'mol/m^3'),
        '1,3-BUTADIENE': (4747.331293433988, 'mol/m^3'),
        'CYCLOPENTADIENE': (1233.5750641865102, 'mol/m^3'),
        'BENZENE': (2296.784867376297, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.421491780202885, 'mol/m^3'),
        'TOLUENE': (215.1661598508465, 'mol/m^3'),
        'STYRENE': (168.73853463077984, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.796990852730758,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.022948638688791004,
        '2-BUTENE': 0.040234418478186595,
        '1,3-BUTADIENE': 0.8383257194258078,
        'CYCLOPENTADIENE': 0.057634299454348635,
        'BENZENE': 0.03870360255617019,
        '1,3-CYCLOHEXADIENE': 0.0006783320553317076,
        'TOLUENE': 0.0011703421727289386,
        'STYRENE': 0.00030464716863520103,

    },
)

# tray 27
constantTVLiquidReactor(
    temperature=(330.019,'K'),
    initialConcentrations={
        'N-BUTANE': (120.08586781347269, 'mol/m^3'),
        '2-BUTENE': (271.6023804370048, 'mol/m^3'),
        '1,3-BUTADIENE': (4738.5111448875095, 'mol/m^3'),
        'CYCLOPENTADIENE': (1235.8119266013903, 'mol/m^3'),
        'BENZENE': (2297.530488181257, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.43558583200396, 'mol/m^3'),
        'TOLUENE': (215.19798512910697, 'mol/m^3'),
        'STYRENE': (168.75581121040696, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.798509074604144,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.022245809916546802,
        '2-BUTENE': 0.041824606655381485,
        '1,3-BUTADIENE': 0.8371111361661239,
        'CYCLOPENTADIENE': 0.05788548601218639,
        'BENZENE': 0.03877569159069215,
        '1,3-CYCLOHEXADIENE': 0.000679659358196889,
        'TOLUENE': 0.0011723957928076122,
        'STYRENE': 0.0003052145080648844,

    },
)

# tray 28
constantTVLiquidReactor(
    temperature=(330.059,'K'),
    initialConcentrations={
        'N-BUTANE': (116.28411100184107, 'mol/m^3'),
        '2-BUTENE': (282.3320456791127, 'mol/m^3'),
        '1,3-BUTADIENE': (4726.59939788144, 'mol/m^3'),
        'CYCLOPENTADIENE': (1239.4672871330233, 'mol/m^3'),
        'BENZENE': (2298.758034628447, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.45886374981735, 'mol/m^3'),
        'TOLUENE': (215.24981486798836, 'mol/m^3'),
        'STYRENE': (168.7830900203445, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.80099685338279,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02154195463039694,
        '2-BUTENE': 0.04361211060031248,
        '1,3-BUTADIENE': 0.8355021188333733,
        'CYCLOPENTADIENE': 0.0582916478276189,
        'BENZENE': 0.038888798621993306,
        '1,3-CYCLOHEXADIENE': 0.0006817587289401365,
        'TOLUENE': 0.0011755329811516402,
        'STYRENE': 0.00030607777621324047,

    },
)

# tray 29
constantTVLiquidReactor(
    temperature=(330.119,'K'),
    initialConcentrations={
        'N-BUTANE': (112.43870809430935, 'mol/m^3'),
        '2-BUTENE': (294.280164431765, 'mol/m^3'),
        '1,3-BUTADIENE': (4710.195740172318, 'mol/m^3'),
        'CYCLOPENTADIENE': (1245.5322758758082, 'mol/m^3'),
        'BENZENE': (2300.822131247056, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.49832709486035, 'mol/m^3'),
        'TOLUENE': (215.33528847245938, 'mol/m^3'),
        'STYRENE': (168.8285547035738, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.805133698967378,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.02083482356418545,
        '2-BUTENE': 0.04561505159062335,
        '1,3-BUTADIENE': 0.83334494251313,
        'CYCLOPENTADIENE': 0.05895836668191271,
        'BENZENE': 0.039073544192178476,
        '1,3-CYCLOHEXADIENE': 0.0006852227749869584,
        'TOLUENE': 0.0011805913352488002,
        'STYRENE': 0.0003074573477342603,

    },
)

# tray 30
constantTVLiquidReactor(
    temperature=(330.211,'K'),
    initialConcentrations={
        'N-BUTANE': (108.50965016963579, 'mol/m^3'),
        '2-BUTENE': (307.46310398092436, 'mol/m^3'),
        '1,3-BUTADIENE': (4686.972379978808, 'mol/m^3'),
        'CYCLOPENTADIENE': (1255.643621425998, 'mol/m^3'),
        'BENZENE': (2304.386562412231, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.56716062526947, 'mol/m^3'),
        'TOLUENE': (215.48259404612222, 'mol/m^3'),
        'STYRENE': (168.90493537139895, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.812086011136945,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.020120517162801138,
        '2-BUTENE': 0.04784794081429351,
        '1,3-BUTADIENE': 0.8303917083241272,
        'CYCLOPENTADIENE': 0.0600635512342092,
        'BENZENE': 0.03938633359654256,
        '1,3-CYCLOHEXADIENE': 0.0006911395895420699,
        'TOLUENE': 0.0011890510142605153,
        'STYRENE': 0.00030975826422379935,

    },
)

# tray 31
constantTVLiquidReactor(
    temperature=(330.357,'K'),
    initialConcentrations={
        'N-BUTANE': (104.43146808397019, 'mol/m^3'),
        '2-BUTENE': (321.77902343615705, 'mol/m^3'),
        '1,3-BUTADIENE': (4653.055726289772, 'mol/m^3'),
        'CYCLOPENTADIENE': (1272.538297713995, 'mol/m^3'),
        'BENZENE': (2310.6606886978707, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.689551552522666, 'mol/m^3'),
        'TOLUENE': (215.74083344686446, 'mol/m^3'),
        'STYRENE': (169.03860154009303, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.82406134008704,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.019392206283074837,
        '2-BUTENE': 0.05031531630216249,
        '1,3-BUTADIENE': 0.8262382677011989,
        'CYCLOPENTADIENE': 0.0619056200574209,
        'BENZENE': 0.0399296129371946,
        '1,3-CYCLOHEXADIENE': 0.0007015232272935257,
        'TOLUENE': 0.0012037103900021666,
        'STYRENE': 0.000313743101652765,

    },
)

# tray 32
constantTVLiquidReactor(
    temperature=(330.595,'K'),
    initialConcentrations={
        'N-BUTANE': (100.09413730389777, 'mol/m^3'),
        '2-BUTENE': (336.8987584908834, 'mol/m^3'),
        '1,3-BUTADIENE': (4601.962515276717, 'mol/m^3'),
        'CYCLOPENTADIENE': (1300.6900295695596, 'mol/m^3'),
        'BENZENE': (2321.8995583921464, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (45.910509913016924, 'mol/m^3'),
        'TOLUENE': (216.2018453348093, 'mol/m^3'),
        'STYRENE': (169.27410859922065, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.844710644130274,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.018638593439215114,
        '2-BUTENE': 0.05300238134316177,
        '1,3-BUTADIENE': 0.8202127112851256,
        'CYCLOPENTADIENE': 0.06498207712630887,
        'BENZENE': 0.040893285605563476,
        '1,3-CYCLOHEXADIENE': 0.0007201507465069373,
        'TOLUENE': 0.0012299195670683126,
        'STYRENE': 0.0003208808870499278,

    },
)

# tray 33
constantTVLiquidReactor(
    temperature=(330.993,'K'),
    initialConcentrations={
        'N-BUTANE': (95.31579909650142, 'mol/m^3'),
        '2-BUTENE': (352.06213964150976, 'mol/m^3'),
        '1,3-BUTADIENE': (4522.917617014307, 'mol/m^3'),
        'CYCLOPENTADIENE': (1347.3822592460206, 'mol/m^3'),
        'BENZENE': (2342.2040859223384, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (46.31441815882577, 'mol/m^3'),
        'TOLUENE': (217.03657691889867, 'mol/m^3'),
        'STYRENE': (169.6942022722591, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.880497926668884,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.01784070080283154,
        '2-BUTENE': 0.0558551025134796,
        '1,3-BUTADIENE': 0.8111930365036867,
        'CYCLOPENTADIENE': 0.07011330315509864,
        'BENZENE': 0.04263160191842208,
        '1,3-CYCLOHEXADIENE': 0.0007541580339371115,
        'TOLUENE': 0.0012780300575113525,
        'STYRENE': 0.0003340670150330157,

    },
)

# tray 34
constantTVLiquidReactor(
    temperature=(331.665,'K'),
    initialConcentrations={
        'N-BUTANE': (89.80984409870382, 'mol/m^3'),
        '2-BUTENE': (365.723367658241, 'mol/m^3'),
        '1,3-BUTADIENE': (4398.226176789709, 'mol/m^3'),
        'CYCLOPENTADIENE': (1423.9175069941748, 'mol/m^3'),
        'BENZENE': (2379.1759663243815, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (47.05822037645663, 'mol/m^3'),
        'TOLUENE': (218.57237391838345, 'mol/m^3'),
        'STYRENE': (170.44891601386502, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(20.94244778312163,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.016967595130300196,
        '2-BUTENE': 0.058744383140362046,
        '1,3-BUTADIENE': 0.7973077711726697,
        'CYCLOPENTADIENE': 0.07862817743371307,
        'BENZENE': 0.045807286853308676,
        '1,3-CYCLOHEXADIENE': 0.0008171147654880623,
        'TOLUENE': 0.0013685296072320028,
        'STYRENE': 0.0003591418969262756,

    },
)

# tray 35
constantTVLiquidReactor(
    temperature=(332.813,'K'),
    initialConcentrations={
        'N-BUTANE': (83.1309912029575, 'mol/m^3'),
        '2-BUTENE': (374.97543069539756, 'mol/m^3'),
        '1,3-BUTADIENE': (4199.863763860401, 'mol/m^3'),
        'CYCLOPENTADIENE': (1546.7994528262477, 'mol/m^3'),
        'BENZENE': (2446.5000692502854, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (48.43043544568248, 'mol/m^3'),
        'TOLUENE': (221.41937238220038, 'mol/m^3'),
        'STYRENE': (171.81376580440772, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.048481409734308,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.015968896231340495,
        '2-BUTENE': 0.061404485508541426,
        '1,3-BUTADIENE': 0.7754658169900673,
        'CYCLOPENTADIENE': 0.09261337814324277,
        'BENZENE': 0.051661987807770884,
        '1,3-CYCLOHEXADIENE': 0.0009348467793761602,
        'TOLUENE': 0.001542479635974806,
        'STYRENE': 0.00040810890368629877,

    },
)

# tray 36
constantTVLiquidReactor(
    temperature=(334.784,'K'),
    initialConcentrations={
        'N-BUTANE': (74.65346450929478, 'mol/m^3'),
        '2-BUTENE': (374.7744767955242, 'mol/m^3'),
        '1,3-BUTADIENE': (3886.7303045871226, 'mol/m^3'),
        'CYCLOPENTADIENE': (1736.8327357879534, 'mol/m^3'),
        'BENZENE': (2568.0362604587717, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (50.940995253602814, 'mol/m^3'),
        'TOLUENE': (226.70054998611246, 'mol/m^3'),
        'STYRENE': (174.2706772861175, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.224907455876977,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.014766097430699048,
        '2-BUTENE': 0.06332708898108652,
        '1,3-BUTADIENE': 0.7407008711180484,
        'CYCLOPENTADIENE': 0.1151589799623375,
        'BENZENE': 0.06250188912467129,
        '1,3-CYCLOHEXADIENE': 0.0011561097988368949,
        'TOLUENE': 0.001882939672368497,
        'STYRENE': 0.0005060239119518393,

    },
)

# tray 37
constantTVLiquidReactor(
    temperature=(338.154,'K'),
    initialConcentrations={
        'N-BUTANE': (63.65346626070565, 'mol/m^3'),
        '2-BUTENE': (357.39241910330946, 'mol/m^3'),
        '1,3-BUTADIENE': (3408.287257092215, 'mol/m^3'),
        'CYCLOPENTADIENE': (2011.1666343933725, 'mol/m^3'),
        'BENZENE': (2782.011245609011, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (55.42272186761139, 'mol/m^3'),
        'TOLUENE': (236.38816468860546, 'mol/m^3'),
        'STYRENE': (178.61164524084828, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.501168177735067,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.013243992252264535,
        '2-BUTENE': 0.06360746278963428,
        '1,3-BUTADIENE': 0.6855025989809796,
        'CYCLOPENTADIENE': 0.15030491207162647,
        'BENZENE': 0.08250685173349175,
        '1,3-CYCLOHEXADIENE': 0.0015707690811000877,
        'TOLUENE': 0.002557848503658626,
        'STYRENE': 0.0007055645872447166,

    },
)

# tray 38
constantTVLiquidReactor(
    temperature=(343.815,'K'),
    initialConcentrations={
        'N-BUTANE': (49.71663133423855, 'mol/m^3'),
        '2-BUTENE': (313.76450907650127, 'mol/m^3'),
        '1,3-BUTADIENE': (2728.3992911450555, 'mol/m^3'),
        'CYCLOPENTADIENE': (2357.6529852836447, 'mol/m^3'),
        'BENZENE': (3140.618481048207, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (63.031145676663385, 'mol/m^3'),
        'TOLUENE': (253.70929970529303, 'mol/m^3'),
        'STYRENE': (186.03602801218815, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(21.887624211554915,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.011257798423908222,
        '2-BUTENE': 0.060795191488673204,
        '1,3-BUTADIENE': 0.5999969160004318,
        'CYCLOPENTADIENE': 0.20175597175416396,
        'BENZENE': 0.11883998336240235,
        '1,3-CYCLOHEXADIENE': 0.002335249673065046,
        'TOLUENE': 0.0039017594537536773,
        'STYRENE': 0.0011171298436018222,

    },
)

# tray 39
constantTVLiquidReactor(
    temperature=(352.994,'K'),
    initialConcentrations={
        'N-BUTANE': (33.49428678183601, 'mol/m^3'),
        '2-BUTENE': (238.72777728758373, 'mol/m^3'),
        '1,3-BUTADIENE': (1874.2906590633304, 'mol/m^3'),
        'CYCLOPENTADIENE': (2675.687537409039, 'mol/m^3'),
        'BENZENE': (3707.2539210712534, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (75.11265780991042, 'mol/m^3'),
        'TOLUENE': (286.8057705088729, 'mol/m^3'),
        'STYRENE': (201.57131027162998, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(22.792650964941515,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.00869621191316662,
        '2-BUTENE': 0.05303441166757057,
        '1,3-BUTADIENE': 0.475586104628943,
        'CYCLOPENTADIENE': 0.268174058998293,
        'BENZENE': 0.1822990401057888,
        '1,3-CYCLOHEXADIENE': 0.0036882908114239782,
        'TOLUENE': 0.006554091441900118,
        'STYRENE': 0.0019677904329138953,

    },
)

# tray 40
constantTVLiquidReactor(
    temperature=(376.859,'K'),
    initialConcentrations={
        'N-BUTANE': (11.043280627023393, 'mol/m^3'),
        '2-BUTENE': (91.82774859915178, 'mol/m^3'),
        '1,3-BUTADIENE': (625.8459155798549, 'mol/m^3'),
        'CYCLOPENTADIENE': (2251.1746971606926, 'mol/m^3'),
        'BENZENE': (5094.199547663385, 'mol/m^3'),
        '1,3-CYCLOHEXADIENE': (101.88362723580187, 'mol/m^3'),
        'TOLUENE': (509.41995476633855, 'mol/m^3'),
        'STYRENE': (407.53632753053665, 'mol/m^3'),

    },
    terminationTime=(8000,'hr'),
    residenceTime=(22.792650964941515,'s'),
    vaporPressure = (430000,'Pa'),
    vaporMoleFractions={
        'N-BUTANE': 0.005698587720564913,
        '2-BUTENE': 0.03943888422444632,
        '1,3-BUTADIENE': 0.318177872728851,
        'CYCLOPENTADIENE': 0.3323598670560533,
        'BENZENE': 0.28322388671044535,
        '1,3-CYCLOHEXADIENE': 0.005857747656900939,
        'TOLUENE': 0.011561295375481855,
        'STYRENE': 0.0036818585272565905,

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
