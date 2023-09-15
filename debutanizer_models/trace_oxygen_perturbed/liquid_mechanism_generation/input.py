# Data sources
database(
    thermoLibraries=[
        "hwpang_fouling",
        "Conjugated_diene",
        "thermo_DFT_CCSDTF12_BAC",
        "primaryThermoLibrary",
        "primaryNS",
        "NitrogenCurran",
        "NOx2018",
        "FFCM1(-)",
        "SulfurLibrary",
        "SulfurGlarborgH2S",
    ],
    reactionLibraries=[
        "hwpang_fouling",
        "Xu_cyclopentadiene",
        "Conjugated_diene",
        "Klippenstein_Glarborg2016",
    ],
    seedMechanisms=[],
    kineticsDepositories=["training"],
    kineticsFamilies="default",
    kineticsEstimator="rate rules",
)

# List of species

species(label="N-BUTANE", reactive=True, structure=SMILES("CCCC"))

species(
    label="2-BUTENE",
    reactive=True,
    structure=SMILES("CC=CC"),
)

species(
    label="1,3-BUTADIENE",
    reactive=True,
    structure=SMILES("C=CC=C"),
)

species(
    label="CYCLOPENTADIENE",
    reactive=True,
    structure=SMILES("C1=CCC=C1"),
)

species(
    label="BENZENE",
    reactive=True,
    structure=SMILES("C1=CC=CC=C1"),
)

species(
    label="1,3-CYCLOHEXADIENE",
    reactive=True,
    structure=SMILES("C1=CCCC=C1"),
)

species(
    label="TOLUENE",
    reactive=True,
    structure=SMILES("CC1=CC=CC=C1"),
)

species(
    label="STYRENE",
    reactive=True,
    structure=SMILES("C=CC1=CC=CC=C1"),
)

species(
    label="OXYGEN",
    reactive=True,
    structure=SMILES("[O][O]"),
)

# tray 1
constantTVLiquidReactor(
    temperature=(312.429, "K"),
    initialConcentrations={
        "N-BUTANE": (296.5765844815041, "mol/m^3"),
        "2-BUTENE": (281.92912082891564, "mol/m^3"),
        "1,3-BUTADIENE": (8502.257949481502, "mol/m^3"),
        "CYCLOPENTADIENE": (0.39589344443529784, "mol/m^3"),
        "BENZENE": (0.0023112402759985127, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (1.665862552001021e-05, "mol/m^3"),
        "TOLUENE": (4.4109889269461795e-06, "mol/m^3"),
        "STYRENE": (4.556226876697324e-07, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.95032898033134, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.035553883003067326,
        "2-BUTENE": 0.028211142674341633,
        "1,3-BUTADIENE": 0.9362131968868228,
        "CYCLOPENTADIENE": 2.1679934034644303e-05,
        "BENZENE": 9.666614481959447e-08,
        "1,3-CYCLOHEXADIENE": 6.620012274355941e-10,
        "TOLUENE": 1.580785016773818e-10,
        "STYRENE": 1.5508959924234848e-11,
    },
)

# tray 1
constantTVLiquidReactor(
    temperature=(312.429, "K"),
    initialConcentrations={
        "N-BUTANE": (296.5765844815041, "mol/m^3"),
        "2-BUTENE": (281.92912082891564, "mol/m^3"),
        "1,3-BUTADIENE": (8502.257949481502, "mol/m^3"),
        "CYCLOPENTADIENE": (0.39589344443529784, "mol/m^3"),
        "BENZENE": (0.0023112402759985127, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (1.665862552001021e-05, "mol/m^3"),
        "TOLUENE": (4.4109889269461795e-06, "mol/m^3"),
        "STYRENE": (4.556226876697324e-07, "mol/m^3"),
        "OXYGEN": (0.09096013093918687, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.95032898033134, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.03553532411959674,
        "2-BUTENE": 0.028196416651042857,
        "1,3-BUTADIENE": 0.9357245000088154,
        "CYCLOPENTADIENE": 2.1668617257532667e-05,
        "BENZENE": 9.661568575392497e-08,
        "1,3-CYCLOHEXADIENE": 6.616556673279597e-10,
        "TOLUENE": 1.5799598578195676e-10,
        "STYRENE": 1.5500864353352787e-11,
        "OXYGEN": 0.0005219931524491927,
    },
)

# tray 2
constantTVLiquidReactor(
    temperature=(312.531, "K"),
    initialConcentrations={
        "N-BUTANE": (278.37819299409085, "mol/m^3"),
        "2-BUTENE": (303.27378896383055, "mol/m^3"),
        "1,3-BUTADIENE": (8504.476711047948, "mol/m^3"),
        "CYCLOPENTADIENE": (0.6976331962013436, "mol/m^3"),
        "BENZENE": (0.0051578568060383695, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (3.89164413455092e-05, "mol/m^3"),
        "TOLUENE": (1.1311046413470471e-05, "mol/m^3"),
        "STYRENE": (1.2243471868484392e-06, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.94729848887025, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.03368313622931902,
        "2-BUTENE": 0.030042463741195305,
        "1,3-BUTADIENE": 0.9362383599145827,
        "CYCLOPENTADIENE": 3.58396363189146e-05,
        "BENZENE": 1.9865138397371436e-07,
        "1,3-CYCLOHEXADIENE": 1.4195154823820984e-09,
        "TOLUENE": 3.6977922780967813e-10,
        "STYRENE": 3.7905430302306235e-11,
    },
)

# tray 2
constantTVLiquidReactor(
    temperature=(312.531, "K"),
    initialConcentrations={
        "N-BUTANE": (278.37819299409085, "mol/m^3"),
        "2-BUTENE": (303.27378896383055, "mol/m^3"),
        "1,3-BUTADIENE": (8504.476711047948, "mol/m^3"),
        "CYCLOPENTADIENE": (0.6976331962013436, "mol/m^3"),
        "BENZENE": (0.0051578568060383695, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (3.89164413455092e-05, "mol/m^3"),
        "TOLUENE": (1.1311046413470471e-05, "mol/m^3"),
        "STYRENE": (1.2243471868484392e-06, "mol/m^3"),
        "OXYGEN": (0.023839047172996265, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.94729848887025, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.03367669409106582,
        "2-BUTENE": 0.030036717907328477,
        "1,3-BUTADIENE": 0.9360592977004393,
        "CYCLOPENTADIENE": 3.583278173475299e-05,
        "BENZENE": 1.986133904902384e-07,
        "1,3-CYCLOHEXADIENE": 1.4192439899970723e-09,
        "TOLUENE": 3.6970850491462195e-10,
        "STYRENE": 3.7898180620421153e-11,
        "OXYGEN": 0.00019125707919059743,
    },
)

# tray 3
constantTVLiquidReactor(
    temperature=(312.567, "K"),
    initialConcentrations={
        "N-BUTANE": (265.6030367453191, "mol/m^3"),
        "2-BUTENE": (320.7447176425728, "mol/m^3"),
        "1,3-BUTADIENE": (8500.930329855677, "mol/m^3"),
        "CYCLOPENTADIENE": (1.1542288653776733, "mol/m^3"),
        "BENZENE": (0.010803186440167105, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (8.538649276485263e-05, "mol/m^3"),
        "TOLUENE": (2.7293858984225422e-05, "mol/m^3"),
        "STYRENE": (3.0996280948955933e-06, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.9514310832657, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.03237511936331659,
        "2-BUTENE": 0.0315484741012161,
        "1,3-BUTADIENE": 0.9360187157645903,
        "CYCLOPENTADIENE": 5.728578529640691e-05,
        "BENZENE": 4.01030653446669e-07,
        "1,3-CYCLOHEXADIENE": 3.0019816555569598e-09,
        "TOLUENE": 8.603813254997271e-10,
        "STYRENE": 9.256413171161396e-11,
    },
)

# tray 3
constantTVLiquidReactor(
    temperature=(312.567, "K"),
    initialConcentrations={
        "N-BUTANE": (265.6030367453191, "mol/m^3"),
        "2-BUTENE": (320.7447176425728, "mol/m^3"),
        "1,3-BUTADIENE": (8500.930329855677, "mol/m^3"),
        "CYCLOPENTADIENE": (1.1542288653776733, "mol/m^3"),
        "BENZENE": (0.010803186440167105, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (8.538649276485263e-05, "mol/m^3"),
        "TOLUENE": (2.7293858984225422e-05, "mol/m^3"),
        "STYRENE": (3.0996280948955933e-06, "mol/m^3"),
        "OXYGEN": (0.02285606305277972, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.9514310832657, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.032369084649036164,
        "2-BUTENE": 0.031542593473409056,
        "1,3-BUTADIENE": 0.935844242106366,
        "CYCLOPENTADIENE": 5.727510724012813e-05,
        "BENZENE": 4.0095590143157723e-07,
        "1,3-CYCLOHEXADIENE": 3.00142208691528e-09,
        "TOLUENE": 8.602209506324305e-10,
        "STYRENE": 9.254687778024688e-11,
        "OXYGEN": 0.00018639975385740642,
    },
)

# tray 4
constantTVLiquidReactor(
    temperature=(312.587, "K"),
    initialConcentrations={
        "N-BUTANE": (256.68616136803564, "mol/m^3"),
        "2-BUTENE": (335.08937490105205, "mol/m^3"),
        "1,3-BUTADIENE": (8495.265213233153, "mol/m^3"),
        "CYCLOPENTADIENE": (1.8451366679128545, "mol/m^3"),
        "BENZENE": (0.021998020669952403, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.0001824003873320965, "mol/m^3"),
        "TOLUENE": (6.431344105612528e-05, "mol/m^3"),
        "STYRENE": (7.674141567945883e-06, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.95584006065909, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.031463124099745346,
        "2-BUTENE": 0.03278704785173842,
        "1,3-BUTADIENE": 0.9356592698738803,
        "CYCLOPENTADIENE": 8.97472243976501e-05,
        "BENZENE": 8.024214106804119e-07,
        "1,3-CYCLOHEXADIENE": 6.306108486491007e-09,
        "TOLUENE": 1.9968148108999133e-09,
        "STYRENE": 2.2590450776668837e-10,
    },
)

# tray 4
constantTVLiquidReactor(
    temperature=(312.587, "K"),
    initialConcentrations={
        "N-BUTANE": (256.68616136803564, "mol/m^3"),
        "2-BUTENE": (335.08937490105205, "mol/m^3"),
        "1,3-BUTADIENE": (8495.265213233153, "mol/m^3"),
        "CYCLOPENTADIENE": (1.8451366679128545, "mol/m^3"),
        "BENZENE": (0.021998020669952403, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.0001824003873320965, "mol/m^3"),
        "TOLUENE": (6.431344105612528e-05, "mol/m^3"),
        "STYRENE": (7.674141567945883e-06, "mol/m^3"),
        "OXYGEN": (0.022837785549711866, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.95584006065909, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.03145726187885791,
        "2-BUTENE": 0.03278093895688919,
        "1,3-BUTADIENE": 0.9354849374935995,
        "CYCLOPENTADIENE": 8.973050266169708e-05,
        "BENZENE": 8.02271903227201e-07,
        "1,3-CYCLOHEXADIENE": 6.304933530031796e-09,
        "TOLUENE": 1.9964427636278215e-09,
        "STYRENE": 2.258624171554764e-10,
        "OXYGEN": 0.00018632036884988105,
    },
)

# tray 5
constantTVLiquidReactor(
    temperature=(312.602, "K"),
    initialConcentrations={
        "N-BUTANE": (250.47453831049705, "mol/m^3"),
        "2-BUTENE": (346.8788190440747, "mol/m^3"),
        "1,3-BUTADIENE": (8488.75442109555, "mol/m^3"),
        "CYCLOPENTADIENE": (2.8905098172580366, "mol/m^3"),
        "BENZENE": (0.04419654827654472, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.0003849223959520299, "mol/m^3"),
        "TOLUENE": (0.00015005648153007772, "mol/m^3"),
        "STYRENE": (1.8832920922276428e-05, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.960249943311155, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.030828187835231734,
        "2-BUTENE": 0.03380581547381704,
        "1,3-BUTADIENE": 0.9352255097721821,
        "CYCLOPENTADIENE": 0.0001388701420234544,
        "BENZENE": 1.598392519842638e-06,
        "1,3-CYCLOHEXADIENE": 1.3204038079893814e-08,
        "TOLUENE": 4.6290214295322395e-09,
        "STYRENE": 5.511661384085738e-10,
    },
)

# tray 5
constantTVLiquidReactor(
    temperature=(312.602, "K"),
    initialConcentrations={
        "N-BUTANE": (250.47453831049705, "mol/m^3"),
        "2-BUTENE": (346.8788190440747, "mol/m^3"),
        "1,3-BUTADIENE": (8488.75442109555, "mol/m^3"),
        "CYCLOPENTADIENE": (2.8905098172580366, "mol/m^3"),
        "BENZENE": (0.04419654827654472, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.0003849223959520299, "mol/m^3"),
        "TOLUENE": (0.00015005648153007772, "mol/m^3"),
        "STYRENE": (1.8832920922276428e-05, "mol/m^3"),
        "OXYGEN": (0.02283596689269019, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.960249943311155, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.03082244376701433,
        "2-BUTENE": 0.03379951659853885,
        "1,3-BUTADIENE": 0.9350512536934432,
        "CYCLOPENTADIENE": 0.0001388442670166779,
        "BENZENE": 1.5980946990391176e-06,
        "1,3-CYCLOHEXADIENE": 1.3201577834877744e-08,
        "TOLUENE": 4.6281589261955795e-09,
        "STYRENE": 5.510634422684353e-10,
        "OXYGEN": 0.00018632519848746487,
    },
)

# tray 6
constantTVLiquidReactor(
    temperature=(312.616, "K"),
    initialConcentrations={
        "N-BUTANE": (246.14704392741606, "mol/m^3"),
        "2-BUTENE": (356.5640770130168, "mol/m^3"),
        "1,3-BUTADIENE": (8481.807151272742, "mol/m^3"),
        "CYCLOPENTADIENE": (4.4720230565936685, "mol/m^3"),
        "BENZENE": (0.08821123085092172, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.0008076764952689896, "mol/m^3"),
        "TOLUENE": (0.00034864382568361245, "mol/m^3"),
        "STYRENE": (4.605294243142175e-05, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.965212143708875, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.030386332374278104,
        "2-BUTENE": 0.03464362239220111,
        "1,3-BUTADIENE": 0.9347536357870375,
        "CYCLOPENTADIENE": 0.00021319306351400425,
        "BENZENE": 3.1767096741502824e-06,
        "1,3-CYCLOHEXADIENE": 2.760336693244429e-08,
        "TOLUENE": 1.0725365699270977e-08,
        "STYRENE": 1.344562376342903e-09,
    },
)

# tray 6
constantTVLiquidReactor(
    temperature=(312.616, "K"),
    initialConcentrations={
        "N-BUTANE": (246.14704392741606, "mol/m^3"),
        "2-BUTENE": (356.5640770130168, "mol/m^3"),
        "1,3-BUTADIENE": (8481.807151272742, "mol/m^3"),
        "CYCLOPENTADIENE": (4.4720230565936685, "mol/m^3"),
        "BENZENE": (0.08821123085092172, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.0008076764952689896, "mol/m^3"),
        "TOLUENE": (0.00034864382568361245, "mol/m^3"),
        "STYRENE": (4.605294243142175e-05, "mol/m^3"),
        "OXYGEN": (0.022835421295583684, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.965212143708875, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.030380670324709952,
        "2-BUTENE": 0.034637167058770645,
        "1,3-BUTADIENE": 0.9345794580891619,
        "CYCLOPENTADIENE": 0.0002131533380980398,
        "BENZENE": 3.1761177406646248e-06,
        "1,3-CYCLOHEXADIENE": 2.759822344787056e-08,
        "TOLUENE": 1.072336718390302e-08,
        "STYRENE": 1.3443118367671309e-09,
        "OXYGEN": 0.00018633540561617082,
    },
)

# tray 7
constantTVLiquidReactor(
    temperature=(312.631, "K"),
    initialConcentrations={
        "N-BUTANE": (243.12079864334515, "mol/m^3"),
        "2-BUTENE": (364.50342424114876, "mol/m^3"),
        "1,3-BUTADIENE": (8474.423403764731, "mol/m^3"),
        "CYCLOPENTADIENE": (6.8643026894781185, "mol/m^3"),
        "BENZENE": (0.17547676005057125, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.0016900961568148702, "mol/m^3"),
        "TOLUENE": (0.0008085821864657848, "mol/m^3"),
        "STYRENE": (0.00011245120096434672, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.971278612278006, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.03007882411493111,
        "2-BUTENE": 0.03533190501907566,
        "1,3-BUTADIENE": 0.9342572483475643,
        "CYCLOPENTADIENE": 0.000325630712960067,
        "BENZENE": 6.306021058142861e-06,
        "1,3-CYCLOHEXADIENE": 5.7659978100308243e-08,
        "TOLUENE": 2.4844575407159224e-08,
        "STYRENE": 3.2798573921922375e-09,
    },
)

# tray 7
constantTVLiquidReactor(
    temperature=(312.631, "K"),
    initialConcentrations={
        "N-BUTANE": (243.12079864334515, "mol/m^3"),
        "2-BUTENE": (364.50342424114876, "mol/m^3"),
        "1,3-BUTADIENE": (8474.423403764731, "mol/m^3"),
        "CYCLOPENTADIENE": (6.8643026894781185, "mol/m^3"),
        "BENZENE": (0.17547676005057125, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.0016900961568148702, "mol/m^3"),
        "TOLUENE": (0.0008085821864657848, "mol/m^3"),
        "STYRENE": (0.00011245120096434672, "mol/m^3"),
        "OXYGEN": (0.022834784765626094, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.971278612278006, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.030073218958617424,
        "2-BUTENE": 0.035325320956821914,
        "1,3-BUTADIENE": 0.934083150520662,
        "CYCLOPENTADIENE": 0.0003255700320292,
        "BENZENE": 6.304845937944903e-06,
        "1,3-CYCLOHEXADIENE": 5.7649233225805796e-08,
        "TOLUENE": 2.483994564742619e-08,
        "STYRENE": 3.279246194317627e-09,
        "OXYGEN": 0.00018634891750651575,
    },
)

# tray 8
constantTVLiquidReactor(
    temperature=(312.648, "K"),
    initialConcentrations={
        "N-BUTANE": (240.98478597138524, "mol/m^3"),
        "2-BUTENE": (370.9805712238522, "mol/m^3"),
        "1,3-BUTADIENE": (8466.294006877835, "mol/m^3"),
        "CYCLOPENTADIENE": (10.482284408692156, "mol/m^3"),
        "BENZENE": (0.34847923722315066, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.00353187740252273, "mol/m^3"),
        "TOLUENE": (0.0018737987025744878, "mol/m^3"),
        "STYRENE": (0.00027442170464343204, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.979553830744514, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029863936714688268,
        "2-BUTENE": 0.035895945265176095,
        "1,3-BUTADIENE": 0.9337317239679019,
        "CYCLOPENTADIENE": 0.0004956981356124432,
        "BENZENE": 1.2509979239022567e-05,
        "1,3-CYCLOHEXADIENE": 1.2039374988226417e-07,
        "TOLUENE": 5.754323994584882e-08,
        "STYRENE": 8.0003924689626e-09,
    },
)

# tray 8
constantTVLiquidReactor(
    temperature=(312.648, "K"),
    initialConcentrations={
        "N-BUTANE": (240.98478597138524, "mol/m^3"),
        "2-BUTENE": (370.9805712238522, "mol/m^3"),
        "1,3-BUTADIENE": (8466.294006877835, "mol/m^3"),
        "CYCLOPENTADIENE": (10.482284408692156, "mol/m^3"),
        "BENZENE": (0.34847923722315066, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.00353187740252273, "mol/m^3"),
        "TOLUENE": (0.0018737987025744878, "mol/m^3"),
        "STYRENE": (0.00027442170464343204, "mol/m^3"),
        "OXYGEN": (0.02283314797430659, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.979553830744514, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02985837112768263,
        "2-BUTENE": 0.03588925552402,
        "1,3-BUTADIENE": 0.9335577092290107,
        "CYCLOPENTADIENE": 0.0004956057549216909,
        "BENZENE": 1.2507647819070847e-05,
        "1,3-CYCLOHEXADIENE": 1.203713127234828e-07,
        "TOLUENE": 5.753251591064995e-08,
        "STYRENE": 7.99890147730967e-09,
        "OXYGEN": 0.00018636481381584184,
    },
)

# tray 9
constantTVLiquidReactor(
    temperature=(312.67, "K"),
    initialConcentrations={
        "N-BUTANE": (239.4416554884921, "mol/m^3"),
        "2-BUTENE": (376.21011948968516, "mol/m^3"),
        "1,3-BUTADIENE": (8456.791523939572, "mol/m^3"),
        "CYCLOPENTADIENE": (15.952259132790774, "mol/m^3"),
        "BENZENE": (0.6914142985157182, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.007375636297692889, "mol/m^3"),
        "TOLUENE": (0.0043407978591945475, "mol/m^3"),
        "STYRENE": (0.0006695294892453653, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.99142053688759, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029712474065865038,
        "2-BUTENE": 0.03635582653739624,
        "1,3-BUTADIENE": 0.9331536165601828,
        "CYCLOPENTADIENE": 0.0007528705404045995,
        "BENZENE": 2.480819340160545e-05,
        "1,3-CYCLOHEXADIENE": 2.5132046817474883e-07,
        "TOLUENE": 1.332677341394015e-07,
        "STYRENE": 1.9514547371591313e-08,
    },
)

# tray 9
constantTVLiquidReactor(
    temperature=(312.67, "K"),
    initialConcentrations={
        "N-BUTANE": (239.4416554884921, "mol/m^3"),
        "2-BUTENE": (376.21011948968516, "mol/m^3"),
        "1,3-BUTADIENE": (8456.791523939572, "mol/m^3"),
        "CYCLOPENTADIENE": (15.952259132790774, "mol/m^3"),
        "BENZENE": (0.6914142985157182, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.007375636297692889, "mol/m^3"),
        "TOLUENE": (0.0043407978591945475, "mol/m^3"),
        "STYRENE": (0.0006695294892453653, "mol/m^3"),
        "OXYGEN": (0.02282941972741215, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(42.99142053688759, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02970693602257601,
        "2-BUTENE": 0.03634905025411783,
        "1,3-BUTADIENE": 0.9329796880912047,
        "CYCLOPENTADIENE": 0.0007527302145052969,
        "BENZENE": 2.480356945650158e-05,
        "1,3-CYCLOHEXADIENE": 2.5127362510038597e-07,
        "TOLUENE": 1.3324289465686432e-07,
        "STYRENE": 1.9510910097634586e-08,
        "OXYGEN": 0.0001863878207097195,
    },
)

# tray 10
constantTVLiquidReactor(
    temperature=(312.703, "K"),
    initialConcentrations={
        "N-BUTANE": (238.2713496950428, "mol/m^3"),
        "2-BUTENE": (380.34756421400084, "mol/m^3"),
        "1,3-BUTADIENE": (8444.870227162477, "mol/m^3"),
        "CYCLOPENTADIENE": (24.21869182627139, "mol/m^3"),
        "BENZENE": (1.371094621927575, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.01539647754696621, "mol/m^3"),
        "TOLUENE": (0.010054172545789337, "mol/m^3"),
        "STYRENE": (0.0016333904308789742, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.009371025629505, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029603308088409504,
        "2-BUTENE": 0.03672688614185015,
        "1,3-BUTADIENE": 0.9324781101750769,
        "CYCLOPENTADIENE": 0.0011416316668933309,
        "BENZENE": 4.918318069630277e-05,
        "1,3-CYCLOHEXADIENE": 5.245287036700553e-07,
        "TOLUENE": 3.086195980595432e-07,
        "STYRENE": 4.759877208377542e-08,
    },
)

# tray 10
constantTVLiquidReactor(
    temperature=(312.703, "K"),
    initialConcentrations={
        "N-BUTANE": (238.2713496950428, "mol/m^3"),
        "2-BUTENE": (380.34756421400084, "mol/m^3"),
        "1,3-BUTADIENE": (8444.870227162477, "mol/m^3"),
        "CYCLOPENTADIENE": (24.21869182627139, "mol/m^3"),
        "BENZENE": (1.371094621927575, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.01539647754696621, "mol/m^3"),
        "TOLUENE": (0.010054172545789337, "mol/m^3"),
        "STYRENE": (0.0016333904308789742, "mol/m^3"),
        "OXYGEN": (0.022822236032176525, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.009371025629505, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029597789439914812,
        "2-BUTENE": 0.036720039515982655,
        "1,3-BUTADIENE": 0.932304277612051,
        "CYCLOPENTADIENE": 0.0011414188439256689,
        "BENZENE": 4.917401196808824e-05,
        "1,3-CYCLOHEXADIENE": 5.244309210326459e-07,
        "TOLUENE": 3.085620651961491e-07,
        "STYRENE": 4.7589898721003774e-08,
        "OXYGEN": 0.00018641999327262767,
    },
)

# tray 11
constantTVLiquidReactor(
    temperature=(312.75, "K"),
    initialConcentrations={
        "N-BUTANE": (237.29927751695632, "mol/m^3"),
        "2-BUTENE": (383.4811102623507, "mol/m^3"),
        "1,3-BUTADIENE": (8428.8478588015, "mol/m^3"),
        "CYCLOPENTADIENE": (36.703044817276925, "mol/m^3"),
        "BENZENE": (2.717864706190217, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.03213085013192928, "mol/m^3"),
        "TOLUENE": (0.023285629841299883, "mol/m^3"),
        "STYRENE": (0.003984841213626771, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.03715476774712, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029521072807420518,
        "2-BUTENE": 0.03701969335881033,
        "1,3-BUTADIENE": 0.9316308290861999,
        "CYCLOPENTADIENE": 0.0017289961014882666,
        "BENZENE": 9.748336939451182e-05,
        "1,3-CYCLOHEXADIENE": 1.0945313435680342e-06,
        "TOLUENE": 7.146478486595723e-07,
        "STYRENE": 1.1609749436072664e-07,
    },
)

# tray 11
constantTVLiquidReactor(
    temperature=(312.75, "K"),
    initialConcentrations={
        "N-BUTANE": (237.29927751695632, "mol/m^3"),
        "2-BUTENE": (383.4811102623507, "mol/m^3"),
        "1,3-BUTADIENE": (8428.8478588015, "mol/m^3"),
        "CYCLOPENTADIENE": (36.703044817276925, "mol/m^3"),
        "BENZENE": (2.717864706190217, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.03213085013192928, "mol/m^3"),
        "TOLUENE": (0.023285629841299883, "mol/m^3"),
        "STYRENE": (0.003984841213626771, "mol/m^3"),
        "OXYGEN": (0.022809323567322615, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.03715476774712, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02951556803976644,
        "2-BUTENE": 0.03701279033018753,
        "1,3-BUTADIENE": 0.9314571087310167,
        "CYCLOPENTADIENE": 0.0017286736971544002,
        "BENZENE": 9.746519175909343e-05,
        "1,3-CYCLOHEXADIENE": 1.0943272473017586e-06,
        "TOLUENE": 7.145145889238229e-07,
        "STYRENE": 1.1607584576631933e-07,
        "OXYGEN": 0.00018646909243394825,
    },
)

# tray 12
constantTVLiquidReactor(
    temperature=(312.823, "K"),
    initialConcentrations={
        "N-BUTANE": (236.36721579334673, "mol/m^3"),
        "2-BUTENE": (385.62257890537563, "mol/m^3"),
        "1,3-BUTADIENE": (8406.069179604992, "mol/m^3"),
        "CYCLOPENTADIENE": (55.538148127979376, "mol/m^3"),
        "BENZENE": (5.3855981315779005, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.06703460662480984, "mol/m^3"),
        "TOLUENE": (0.053927818268135266, "mol/m^3"),
        "STYRENE": (0.009722085773631073, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.0811845305962, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02945324165688815,
        "2-BUTENE": 0.03724025126029976,
        "1,3-BUTADIENE": 0.9304934048957203,
        "CYCLOPENTADIENE": 0.002615720805653496,
        "BENZENE": 0.00019316015048507942,
        "1,3-CYCLOHEXADIENE": 2.2833827714968046e-06,
        "TOLUENE": 1.6546882029781846e-06,
        "STYRENE": 2.831599787300747e-07,
    },
)

# tray 12
constantTVLiquidReactor(
    temperature=(312.823, "K"),
    initialConcentrations={
        "N-BUTANE": (236.36721579334673, "mol/m^3"),
        "2-BUTENE": (385.62257890537563, "mol/m^3"),
        "1,3-BUTADIENE": (8406.069179604992, "mol/m^3"),
        "CYCLOPENTADIENE": (55.538148127979376, "mol/m^3"),
        "BENZENE": (5.3855981315779005, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.06703460662480984, "mol/m^3"),
        "TOLUENE": (0.053927818268135266, "mol/m^3"),
        "STYRENE": (0.009722085773631073, "mol/m^3"),
        "OXYGEN": (0.022786954085955984, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.0811845305962, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029447747299039168,
        "2-BUTENE": 0.037233304273981804,
        "1,3-BUTADIENE": 0.930319825912389,
        "CYCLOPENTADIENE": 0.002615232855759671,
        "BENZENE": 0.0001931241174058927,
        "1,3-CYCLOHEXADIENE": 2.282956817634104e-06,
        "TOLUENE": 1.654379529005331e-06,
        "STYRENE": 2.8310715662411516e-07,
        "OXYGEN": 0.00018654509792116882,
    },
)

# tray 13
constantTVLiquidReactor(
    temperature=(312.937, "K"),
    initialConcentrations={
        "N-BUTANE": (235.31239472077394, "mol/m^3"),
        "2-BUTENE": (386.6846746060351, "mol/m^3"),
        "1,3-BUTADIENE": (8372.260345572013, "mol/m^3"),
        "CYCLOPENTADIENE": (83.90828833763338, "mol/m^3"),
        "BENZENE": (10.667332760647817, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.13980471189888316, "mol/m^3"),
        "TOLUENE": (0.124890814992619, "mol/m^3"),
        "STYRENE": (0.023722925922162175, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.15240388514281, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029389343780664573,
        "2-BUTENE": 0.037389120132589666,
        "1,3-BUTADIENE": 0.9288770004883815,
        "CYCLOPENTADIENE": 0.003952667447793632,
        "BENZENE": 0.00038258506428257857,
        "1,3-CYCLOHEXADIENE": 4.761840803915088e-06,
        "TOLUENE": 3.830658862919353e-06,
        "STYRENE": 6.905866212137674e-07,
    },
)

# tray 13
constantTVLiquidReactor(
    temperature=(312.937, "K"),
    initialConcentrations={
        "N-BUTANE": (235.31239472077394, "mol/m^3"),
        "2-BUTENE": (386.6846746060351, "mol/m^3"),
        "1,3-BUTADIENE": (8372.260345572013, "mol/m^3"),
        "CYCLOPENTADIENE": (83.90828833763338, "mol/m^3"),
        "BENZENE": (10.667332760647817, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.13980471189888316, "mol/m^3"),
        "TOLUENE": (0.124890814992619, "mol/m^3"),
        "STYRENE": (0.023722925922162175, "mol/m^3"),
        "OXYGEN": (0.022748853221351846, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.15240388514281, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029383857809090325,
        "2-BUTENE": 0.037382140880117605,
        "1,3-BUTADIENE": 0.9287036113559571,
        "CYCLOPENTADIENE": 0.003951929621817562,
        "BENZENE": 0.00038251364891505457,
        "1,3-CYCLOHEXADIENE": 4.760951933326949e-06,
        "TOLUENE": 3.82994381171612e-06,
        "STYRENE": 6.904577126337797e-07,
        "OXYGEN": 0.00018666533064463684,
    },
)

# tray 14
constantTVLiquidReactor(
    temperature=(313.118, "K"),
    initialConcentrations={
        "N-BUTANE": (233.92657807025589, "mol/m^3"),
        "2-BUTENE": (386.4355185940654, "mol/m^3"),
        "1,3-BUTADIENE": (8320.49227344997, "mol/m^3"),
        "CYCLOPENTADIENE": (126.52942496915018, "mol/m^3"),
        "BENZENE": (21.11660854439653, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.2914243291390925, "mol/m^3"),
        "TOLUENE": (0.28923921272754727, "mol/m^3"),
        "STYRENE": (0.057900310800584845, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.26939039104348, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029317950486895554,
        "2-BUTENE": 0.03745907647508385,
        "1,3-BUTADIENE": 0.9264808128225474,
        "CYCLOPENTADIENE": 0.005964372492997727,
        "BENZENE": 0.0007573124380121265,
        "1,3-CYCLOHEXADIENE": 9.925057532858639e-06,
        "TOLUENE": 8.866126153707932e-06,
        "STYRENE": 1.6841007768803442e-06,
    },
)

# tray 14
constantTVLiquidReactor(
    temperature=(313.118, "K"),
    initialConcentrations={
        "N-BUTANE": (233.92657807025589, "mol/m^3"),
        "2-BUTENE": (386.4355185940654, "mol/m^3"),
        "1,3-BUTADIENE": (8320.49227344997, "mol/m^3"),
        "CYCLOPENTADIENE": (126.52942496915018, "mol/m^3"),
        "BENZENE": (21.11660854439653, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.2914243291390925, "mol/m^3"),
        "TOLUENE": (0.28923921272754727, "mol/m^3"),
        "STYRENE": (0.057900310800584845, "mol/m^3"),
        "OXYGEN": (0.022684018098529053, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.26939039104348, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02931247217065436,
        "2-BUTENE": 0.03745207691803346,
        "1,3-BUTADIENE": 0.9263076917550871,
        "CYCLOPENTADIENE": 0.0059632579976747506,
        "BENZENE": 0.0007571709275395349,
        "1,3-CYCLOHEXADIENE": 9.92320294879072e-06,
        "TOLUENE": 8.864469440257719e-06,
        "STYRENE": 1.683786087876353e-06,
        "OXYGEN": 0.00018685877253389962,
    },
)

# tray 15
constantTVLiquidReactor(
    temperature=(313.41, "K"),
    initialConcentrations={
        "N-BUTANE": (231.91059676172665, "mol/m^3"),
        "2-BUTENE": (384.4158999714928, "mol/m^3"),
        "1,3-BUTADIENE": (8239.371076998057, "mol/m^3"),
        "CYCLOPENTADIENE": (190.28790283511066, "mol/m^3"),
        "BENZENE": (41.76673156273285, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.6070004235260587, "mol/m^3"),
        "TOLUENE": (0.6699168631909825, "mol/m^3"),
        "STYRENE": (0.14137966887965564, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.46531357950995, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029225622836447158,
        "2-BUTENE": 0.037432469748745256,
        "1,3-BUTADIENE": 0.9228175026598565,
        "CYCLOPENTADIENE": 0.008981442717379719,
        "BENZENE": 0.0014976734311890935,
        "1,3-CYCLOHEXADIENE": 2.0668676691182536e-05,
        "TOLUENE": 2.05135111083885e-05,
        "STYRENE": 4.106418582634226e-06,
    },
)

# tray 15
constantTVLiquidReactor(
    temperature=(313.41, "K"),
    initialConcentrations={
        "N-BUTANE": (231.91059676172665, "mol/m^3"),
        "2-BUTENE": (384.4158999714928, "mol/m^3"),
        "1,3-BUTADIENE": (8239.371076998057, "mol/m^3"),
        "CYCLOPENTADIENE": (190.28790283511066, "mol/m^3"),
        "BENZENE": (41.76673156273285, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (0.6070004235260587, "mol/m^3"),
        "TOLUENE": (0.6699168631909825, "mol/m^3"),
        "STYRENE": (0.14137966887965564, "mol/m^3"),
        "OXYGEN": (0.02257308002020674, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.46531357950995, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029220152519098947,
        "2-BUTENE": 0.037425463311626744,
        "1,3-BUTADIENE": 0.9226447739340299,
        "CYCLOPENTADIENE": 0.008979761612337621,
        "BENZENE": 0.001497393103580733,
        "1,3-CYCLOHEXADIENE": 2.0664808023565047e-05,
        "TOLUENE": 2.0509671483949464e-05,
        "STYRENE": 4.1056499621351076e-06,
        "OXYGEN": 0.00018717538985640062,
    },
)

# tray 16
constantTVLiquidReactor(
    temperature=(313.888, "K"),
    initialConcentrations={
        "N-BUTANE": (228.8006932546586, "mol/m^3"),
        "2-BUTENE": (379.8046950930302, "mol/m^3"),
        "1,3-BUTADIENE": (8109.946350050398, "mol/m^3"),
        "CYCLOPENTADIENE": (284.98810193937675, "mol/m^3"),
        "BENZENE": (82.50755869953758, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (1.2627390365760456, "mol/m^3"),
        "TOLUENE": (1.5519691560184843, "mol/m^3"),
        "STYRENE": (0.34548936507951333, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.79926298536231, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.029093061280977515,
        "2-BUTENE": 0.03727458322113507,
        "1,3-BUTADIENE": 0.9170927637802245,
        "CYCLOPENTADIENE": 0.013481705685191851,
        "BENZENE": 0.0029574586545334266,
        "1,3-CYCLOHEXADIENE": 4.2980862955705936e-05,
        "TOLUENE": 4.743569789568374e-05,
        "STYRENE": 1.0010817086314654e-05,
    },
)

# tray 16
constantTVLiquidReactor(
    temperature=(313.888, "K"),
    initialConcentrations={
        "N-BUTANE": (228.8006932546586, "mol/m^3"),
        "2-BUTENE": (379.8046950930302, "mol/m^3"),
        "1,3-BUTADIENE": (8109.946350050398, "mol/m^3"),
        "CYCLOPENTADIENE": (284.98810193937675, "mol/m^3"),
        "BENZENE": (82.50755869953758, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (1.2627390365760456, "mol/m^3"),
        "TOLUENE": (1.5519691560184843, "mol/m^3"),
        "STYRENE": (0.34548936507951333, "mol/m^3"),
        "OXYGEN": (0.022382030100079553, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(43.79926298536231, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02908760035276931,
        "2-BUTENE": 0.037267586576093296,
        "1,3-BUTADIENE": 0.9169206204057303,
        "CYCLOPENTADIENE": 0.01347917509461013,
        "BENZENE": 0.002956903523217572,
        "1,3-CYCLOHEXADIENE": 4.297279521045695e-05,
        "TOLUENE": 4.742679395332398e-05,
        "STYRENE": 1.0008938000683688e-05,
        "OXYGEN": 0.0001877055204149567,
    },
)

# tray 17
constantTVLiquidReactor(
    temperature=(314.686, "K"),
    initialConcentrations={
        "N-BUTANE": (223.8312129429256, "mol/m^3"),
        "2-BUTENE": (371.18062349623676, "mol/m^3"),
        "1,3-BUTADIENE": (7900.618926855348, "mol/m^3"),
        "CYCLOPENTADIENE": (423.9071278687031, "mol/m^3"),
        "BENZENE": (162.67705193200692, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (2.6214940706115404, "mol/m^3"),
        "TOLUENE": (3.5971490030307263, "mol/m^3"),
        "STYRENE": (0.8454918307207243, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(44.37945392635248, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.028890987253988486,
        "2-BUTENE": 0.036924074901811064,
        "1,3-BUTADIENE": 0.908003458273232,
        "CYCLOPENTADIENE": 0.020132053795482663,
        "BENZENE": 0.005826269960395669,
        "1,3-CYCLOHEXADIENE": 8.91681443866286e-05,
        "TOLUENE": 0.00010959109746919714,
        "STYRENE": 2.4396573234125716e-05,
    },
)

# tray 17
constantTVLiquidReactor(
    temperature=(314.686, "K"),
    initialConcentrations={
        "N-BUTANE": (223.8312129429256, "mol/m^3"),
        "2-BUTENE": (371.18062349623676, "mol/m^3"),
        "1,3-BUTADIENE": (7900.618926855348, "mol/m^3"),
        "CYCLOPENTADIENE": (423.9071278687031, "mol/m^3"),
        "BENZENE": (162.67705193200692, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (2.6214940706115404, "mol/m^3"),
        "TOLUENE": (3.5971490030307263, "mol/m^3"),
        "STYRENE": (0.8454918307207243, "mol/m^3"),
        "OXYGEN": (0.02205021612647454, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(44.37945392635248, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02888553829259526,
        "2-BUTENE": 0.03691711086638881,
        "1,3-BUTADIENE": 0.907832205011945,
        "CYCLOPENTADIENE": 0.02012825680568326,
        "BENZENE": 0.005825171101440049,
        "1,3-CYCLOHEXADIENE": 8.915132690053866e-05,
        "TOLUENE": 0.00010957042812848194,
        "STYRENE": 2.4391971938070612e-05,
        "OXYGEN": 0.00018860419498042288,
    },
)

# tray 18
constantTVLiquidReactor(
    temperature=(316.038, "K"),
    initialConcentrations={
        "N-BUTANE": (215.71818396922342, "mol/m^3"),
        "2-BUTENE": (356.13851126994433, "mol/m^3"),
        "1,3-BUTADIENE": (7558.6932202098105, "mol/m^3"),
        "CYCLOPENTADIENE": (623.186471018986, "mol/m^3"),
        "BENZENE": (319.8090186049199, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (5.423608063330955, "mol/m^3"),
        "TOLUENE": (8.347408397370794, "mol/m^3"),
        "STYRENE": (2.0749603557422063, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(45.43330766824121, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02857303083209566,
        "2-BUTENE": 0.0362767041263674,
        "1,3-BUTADIENE": 0.8934203476241,
        "CYCLOPENTADIENE": 0.029800345473745212,
        "BENZENE": 0.011433099029548678,
        "1,3-CYCLOHEXADIENE": 0.0001842404430591832,
        "TOLUENE": 0.00025281058821931653,
        "STYRENE": 5.942188286429222e-05,
    },
)

# tray 18
constantTVLiquidReactor(
    temperature=(316.038, "K"),
    initialConcentrations={
        "N-BUTANE": (215.71818396922342, "mol/m^3"),
        "2-BUTENE": (356.13851126994433, "mol/m^3"),
        "1,3-BUTADIENE": (7558.6932202098105, "mol/m^3"),
        "CYCLOPENTADIENE": (623.186471018986, "mol/m^3"),
        "BENZENE": (319.8090186049199, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (5.423608063330955, "mol/m^3"),
        "TOLUENE": (8.347408397370794, "mol/m^3"),
        "STYRENE": (2.0749603557422063, "mol/m^3"),
        "OXYGEN": (0.021467336551026985, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(45.43330766824121, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02856759751224658,
        "2-BUTENE": 0.03626980591043269,
        "1,3-BUTADIENE": 0.893250458803525,
        "CYCLOPENTADIENE": 0.02979467877322879,
        "BENZENE": 0.011430924962532169,
        "1,3-CYCLOHEXADIENE": 0.00018420540871990856,
        "TOLUENE": 0.0002527625148876808,
        "STYRENE": 5.9410583464605375e-05,
        "OXYGEN": 0.00019015553096239575,
    },
)

# tray 19
constantTVLiquidReactor(
    temperature=(318.458, "K"),
    initialConcentrations={
        "N-BUTANE": (201.40171589458024, "mol/m^3"),
        "2-BUTENE": (330.6682196813547, "mol/m^3"),
        "1,3-BUTADIENE": (6979.105413971493, "mol/m^3"),
        "CYCLOPENTADIENE": (906.5850667214169, "mol/m^3"),
        "BENZENE": (635.5442454812828, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (11.335234551858779, "mol/m^3"),
        "TOLUENE": (19.79108037414708, "mol/m^3"),
        "STYRENE": (5.225228955406403, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(50.078038566763695, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02806685462222585,
        "2-BUTENE": 0.03515965016771664,
        "1,3-BUTADIENE": 0.869951082820798,
        "CYCLOPENTADIENE": 0.043433092394927014,
        "BENZENE": 0.022285131855484912,
        "1,3-CYCLOHEXADIENE": 0.00037793115348076897,
        "TOLUENE": 0.0005816683371858885,
        "STYRENE": 0.0001445886481809265,
    },
)

# tray 19
constantTVLiquidReactor(
    temperature=(318.458, "K"),
    initialConcentrations={
        "N-BUTANE": (201.40171589458024, "mol/m^3"),
        "2-BUTENE": (330.6682196813547, "mol/m^3"),
        "1,3-BUTADIENE": (6979.105413971493, "mol/m^3"),
        "CYCLOPENTADIENE": (906.5850667214169, "mol/m^3"),
        "BENZENE": (635.5442454812828, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (11.335234551858779, "mol/m^3"),
        "TOLUENE": (19.79108037414708, "mol/m^3"),
        "STYRENE": (5.225228955406403, "mol/m^3"),
        "OXYGEN": (0.01994148330983979, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(50.078038566763695, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02806143950733814,
        "2-BUTENE": 0.0351528665951492,
        "1,3-BUTADIENE": 0.8697832376837648,
        "CYCLOPENTADIENE": 0.043424712575085685,
        "BENZENE": 0.022280832244755123,
        "1,3-CYCLOHEXADIENE": 0.0003778582368449974,
        "TOLUENE": 0.0005815561122531414,
        "STYRENE": 0.00014456075178313286,
        "OXYGEN": 0.00019293629302575158,
    },
)

# tray 20
constantTVLiquidReactor(
    temperature=(329.857, "K"),
    initialConcentrations={
        "N-BUTANE": (146.49646041014455, "mol/m^3"),
        "2-BUTENE": (223.21923485513116, "mol/m^3"),
        "1,3-BUTADIENE": (4766.509094828753, "mol/m^3"),
        "CYCLOPENTADIENE": (1229.7395049177119, "mol/m^3"),
        "BENZENE": (2295.3634001992705, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.39549791808543, "mol/m^3"),
        "TOLUENE": (215.11620849504826, "mol/m^3"),
        "STYRENE": (168.71863055801867, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.79588222962608, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.027183270046343688,
        "2-BUTENE": 0.0338326990078106,
        "1,3-BUTADIENE": 0.8411197267156377,
        "CYCLOPENTADIENE": 0.05718418806198105,
        "BENZENE": 0.03853645419833336,
        "1,3-CYCLOHEXADIENE": 0.0006753375093515854,
        "TOLUENE": 0.0011651817795350955,
        "STYRENE": 0.0003031426810069294,
    },
)

# tray 20
constantTVLiquidReactor(
    temperature=(329.857, "K"),
    initialConcentrations={
        "N-BUTANE": (146.49646041014455, "mol/m^3"),
        "2-BUTENE": (223.21923485513116, "mol/m^3"),
        "1,3-BUTADIENE": (4766.509094828753, "mol/m^3"),
        "CYCLOPENTADIENE": (1229.7395049177119, "mol/m^3"),
        "BENZENE": (2295.3634001992705, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.39549791808543, "mol/m^3"),
        "TOLUENE": (215.11620849504826, "mol/m^3"),
        "STYRENE": (168.71863055801867, "mol/m^3"),
        "OXYGEN": (0.001790240505713527, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.79588222962608, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.027182904481812153,
        "2-BUTENE": 0.03383224402080049,
        "1,3-BUTADIENE": 0.8411084152163832,
        "CYCLOPENTADIENE": 0.05718341904078219,
        "BENZENE": 0.03853593595454584,
        "1,3-CYCLOHEXADIENE": 0.0006753284273154723,
        "TOLUENE": 0.001165166110002671,
        "STYRENE": 0.0003031386043004853,
        "OXYGEN": 1.3448144057376089e-05,
    },
)

# tray 21
constantTVLiquidReactor(
    temperature=(329.891, "K"),
    initialConcentrations={
        "N-BUTANE": (142.75002694548945, "mol/m^3"),
        "2-BUTENE": (227.94046848340528, "mol/m^3"),
        "1,3-BUTADIENE": (4766.818266522439, "mol/m^3"),
        "CYCLOPENTADIENE": (1229.9486504752049, "mol/m^3"),
        "BENZENE": (2295.3634001992705, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.3951341866811, "mol/m^3"),
        "TOLUENE": (215.1052965529182, "mol/m^3"),
        "STYRENE": (168.7077186158886, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.79429986789779, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02647577889846454,
        "2-BUTENE": 0.034602716346700244,
        "1,3-BUTADIENE": 0.8410191653801523,
        "CYCLOPENTADIENE": 0.05720712431861196,
        "BENZENE": 0.038550626239470386,
        "1,3-CYCLOHEXADIENE": 0.0006755749372272671,
        "TOLUENE": 0.0011656934580365666,
        "STYRENE": 0.0003033204213364636,
    },
)

# tray 21
constantTVLiquidReactor(
    temperature=(329.891, "K"),
    initialConcentrations={
        "N-BUTANE": (142.75002694548945, "mol/m^3"),
        "2-BUTENE": (227.94046848340528, "mol/m^3"),
        "1,3-BUTADIENE": (4766.818266522439, "mol/m^3"),
        "CYCLOPENTADIENE": (1229.9486504752049, "mol/m^3"),
        "BENZENE": (2295.3634001992705, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.3951341866811, "mol/m^3"),
        "TOLUENE": (215.1052965529182, "mol/m^3"),
        "STYRENE": (168.7077186158886, "mol/m^3"),
        "OXYGEN": (4.446743723992355e-05, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.79429986789779, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.026475770058655514,
        "2-BUTENE": 0.03460270479344602,
        "1,3-BUTADIENE": 0.8410188845782602,
        "CYCLOPENTADIENE": 0.057207105218133145,
        "BENZENE": 0.03855061336807668,
        "1,3-CYCLOHEXADIENE": 0.0006755747116643686,
        "TOLUENE": 0.00116569306883149,
        "STYRENE": 0.00030332032006297096,
        "OXYGEN": 3.3388286934938034e-07,
    },
)

# tray 22
constantTVLiquidReactor(
    temperature=(329.913, "K"),
    initialConcentrations={
        "N-BUTANE": (138.98176959657422, "mol/m^3"),
        "2-BUTENE": (233.23185108797517, "mol/m^3"),
        "1,3-BUTADIENE": (4765.599766317915, "mol/m^3"),
        "CYCLOPENTADIENE": (1230.230542313565, "mol/m^3"),
        "BENZENE": (2295.4452397652462, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.396498179447356, "mol/m^3"),
        "TOLUENE": (215.10438722440736, "mol/m^3"),
        "STYRENE": (168.70499063035606, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.793783230505095, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.025768655427876305,
        "2-BUTENE": 0.035472135668104914,
        "1,3-BUTADIENE": 0.8408045823504157,
        "CYCLOPENTADIENE": 0.05724155130662463,
        "BENZENE": 0.03856746304527593,
        "1,3-CYCLOHEXADIENE": 0.0006758678235315874,
        "TOLUENE": 0.0011662540603741554,
        "STYRENE": 0.0003034903177967596,
    },
)

# tray 22
constantTVLiquidReactor(
    temperature=(329.913, "K"),
    initialConcentrations={
        "N-BUTANE": (138.98176959657422, "mol/m^3"),
        "2-BUTENE": (233.23185108797517, "mol/m^3"),
        "1,3-BUTADIENE": (4765.599766317915, "mol/m^3"),
        "CYCLOPENTADIENE": (1230.230542313565, "mol/m^3"),
        "BENZENE": (2295.4452397652462, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.396498179447356, "mol/m^3"),
        "TOLUENE": (215.10438722440736, "mol/m^3"),
        "STYRENE": (168.70499063035606, "mol/m^3"),
        "OXYGEN": (1.10432491670285e-06, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.793783230505095, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02576865521424028,
        "2-BUTENE": 0.035472135374021815,
        "1,3-BUTADIENE": 0.8408045753796931,
        "CYCLOPENTADIENE": 0.057241550832061366,
        "BENZENE": 0.0385674627255309,
        "1,3-CYCLOHEXADIENE": 0.0006758678179282793,
        "TOLUENE": 0.0011662540507052814,
        "STYRENE": 0.00030349031528066156,
        "OXYGEN": 8.290538141422056e-09,
    },
)

# tray 23
constantTVLiquidReactor(
    temperature=(329.93, "K"),
    initialConcentrations={
        "N-BUTANE": (135.20805627659396, "mol/m^3"),
        "2-BUTENE": (239.18613417694644, "mol/m^3"),
        "1,3-BUTADIENE": (4763.190045764193, "mol/m^3"),
        "CYCLOPENTADIENE": (1230.6215535732256, "mol/m^3"),
        "BENZENE": (2295.599825612089, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.39913523212879, "mol/m^3"),
        "TOLUENE": (215.10893386696156, "mol/m^3"),
        "STYRENE": (168.70589995886695, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.793783230505095, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.025062512465599127,
        "2-BUTENE": 0.03645354529430011,
        "1,3-BUTADIENE": 0.8404589611924519,
        "CYCLOPENTADIENE": 0.057290200403809305,
        "BENZENE": 0.03858797629949894,
        "1,3-CYCLOHEXADIENE": 0.0006762301870406831,
        "TOLUENE": 0.0011668962908517662,
        "STYRENE": 0.00030367786644812996,
    },
)

# tray 23
constantTVLiquidReactor(
    temperature=(329.93, "K"),
    initialConcentrations={
        "N-BUTANE": (135.20805627659396, "mol/m^3"),
        "2-BUTENE": (239.18613417694644, "mol/m^3"),
        "1,3-BUTADIENE": (4763.190045764193, "mol/m^3"),
        "CYCLOPENTADIENE": (1230.6215535732256, "mol/m^3"),
        "BENZENE": (2295.599825612089, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.39913523212879, "mol/m^3"),
        "TOLUENE": (215.10893386696156, "mol/m^3"),
        "STYRENE": (168.70589995886695, "mol/m^3"),
        "OXYGEN": (2.742161963999825e-08, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.793783230505095, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02506251246043973,
        "2-BUTENE": 0.03645354528679574,
        "1,3-BUTADIENE": 0.8404589610194341,
        "CYCLOPENTADIENE": 0.057290200392015475,
        "BENZENE": 0.03858797629155518,
        "1,3-CYCLOHEXADIENE": 0.0006762301869014736,
        "TOLUENE": 0.0011668962906115476,
        "STYRENE": 0.0003036778663856145,
        "OXYGEN": 2.058612799122281e-10,
    },
)

# tray 24
constantTVLiquidReactor(
    temperature=(329.946, "K"),
    initialConcentrations={
        "N-BUTANE": (131.4343429566137, "mol/m^3"),
        "2-BUTENE": (245.8960692584246, "mol/m^3"),
        "1,3-BUTADIENE": (4759.589104861273, "mol/m^3"),
        "CYCLOPENTADIENE": (1231.1944305350537, "mol/m^3"),
        "BENZENE": (2295.8180644546896, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.40313627757648, "mol/m^3"),
        "TOLUENE": (215.11711782355908, "mol/m^3"),
        "STYRENE": (168.71044660142113, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.7940738358796, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.024357623312915647,
        "2-BUTENE": 0.037561386445511456,
        "1,3-BUTADIENE": 0.8399591852701618,
        "CYCLOPENTADIENE": 0.05735978078027013,
        "BENZENE": 0.03861373929782876,
        "1,3-CYCLOHEXADIENE": 0.0006766919855015879,
        "TOLUENE": 0.0011676886448268092,
        "STYRENE": 0.00030390426298379985,
    },
)

# tray 24
constantTVLiquidReactor(
    temperature=(329.946, "K"),
    initialConcentrations={
        "N-BUTANE": (131.4343429566137, "mol/m^3"),
        "2-BUTENE": (245.8960692584246, "mol/m^3"),
        "1,3-BUTADIENE": (4759.589104861273, "mol/m^3"),
        "CYCLOPENTADIENE": (1231.1944305350537, "mol/m^3"),
        "BENZENE": (2295.8180644546896, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.40313627757648, "mol/m^3"),
        "TOLUENE": (215.11711782355908, "mol/m^3"),
        "STYRENE": (168.71044660142113, "mol/m^3"),
        "OXYGEN": (6.808333519636002e-10, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.7940738358796, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.024357623312791146,
        "2-BUTENE": 0.037561386445319464,
        "1,3-BUTADIENE": 0.8399591852658684,
        "CYCLOPENTADIENE": 0.05735978077997694,
        "BENZENE": 0.038613739297631386,
        "1,3-CYCLOHEXADIENE": 0.0006766919854981291,
        "TOLUENE": 0.0011676886448208407,
        "STYRENE": 0.00030390426298224646,
        "OXYGEN": 5.111466713166611e-12,
    },
)

# tray 25
constantTVLiquidReactor(
    temperature=(329.964, "K"),
    initialConcentrations={
        "N-BUTANE": (127.66335762216595, "mol/m^3"),
        "2-BUTENE": (253.45804515455856, "mol/m^3"),
        "1,3-BUTADIENE": (4754.651451047418, "mol/m^3"),
        "CYCLOPENTADIENE": (1232.0673859054589, "mol/m^3"),
        "BENZENE": (2296.1272361483752, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.40904691289693, "mol/m^3"),
        "TOLUENE": (215.1307577512217, "mol/m^3"),
        "STYRENE": (168.71681190099702, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.794655070997585, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02365384491362528,
        "2-BUTENE": 0.038811494100389594,
        "1,3-BUTADIENE": 0.8392733343398032,
        "CYCLOPENTADIENE": 0.057463023775678725,
        "BENZENE": 0.038648090336057074,
        "1,3-CYCLOHEXADIENE": 0.0006773126032504411,
        "TOLUENE": 0.001168706923501397,
        "STYRENE": 0.00030419300769431835,
    },
)

# tray 25
constantTVLiquidReactor(
    temperature=(329.964, "K"),
    initialConcentrations={
        "N-BUTANE": (127.66335762216595, "mol/m^3"),
        "2-BUTENE": (253.45804515455856, "mol/m^3"),
        "1,3-BUTADIENE": (4754.651451047418, "mol/m^3"),
        "CYCLOPENTADIENE": (1232.0673859054589, "mol/m^3"),
        "BENZENE": (2296.1272361483752, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.40904691289693, "mol/m^3"),
        "TOLUENE": (215.1307577512217, "mol/m^3"),
        "STYRENE": (168.71681190099702, "mol/m^3"),
        "OXYGEN": (1.6902234628063875e-11, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.794655070997585, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.023653844913622275,
        "2-BUTENE": 0.03881149410038467,
        "1,3-BUTADIENE": 0.8392733343396966,
        "CYCLOPENTADIENE": 0.057463023775671425,
        "BENZENE": 0.038648090336052175,
        "1,3-CYCLOHEXADIENE": 0.0006773126032503552,
        "TOLUENE": 0.0011687069235012486,
        "STYRENE": 0.0003041930076942797,
        "OXYGEN": 1.2690692355478186e-13,
    },
)

# tray 26
constantTVLiquidReactor(
    temperature=(329.987, "K"),
    initialConcentrations={
        "N-BUTANE": (123.89146295920736, "mol/m^3"),
        "2-BUTENE": (261.9720880015405, "mol/m^3"),
        "1,3-BUTADIENE": (4748.040632773621, "mol/m^3"),
        "CYCLOPENTADIENE": (1233.4495652419334, "mol/m^3"),
        "BENZENE": (2296.609180259119, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.417958332303144, "mol/m^3"),
        "TOLUENE": (215.15076297846014, "mol/m^3"),
        "STYRENE": (168.72772384312705, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.7955915737037, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.02295104603521349,
        "2-BUTENE": 0.04022103116945413,
        "1,3-BUTADIENE": 0.8383569871869024,
        "CYCLOPENTADIENE": 0.057621617703741766,
        "BENZENE": 0.03869641482556692,
        "1,3-CYCLOHEXADIENE": 0.0006781942702425768,
        "TOLUENE": 0.0011701225437136683,
        "STYRENE": 0.00030458626516476246,
    },
)

# tray 26
constantTVLiquidReactor(
    temperature=(329.987, "K"),
    initialConcentrations={
        "N-BUTANE": (123.89146295920736, "mol/m^3"),
        "2-BUTENE": (261.9720880015405, "mol/m^3"),
        "1,3-BUTADIENE": (4748.040632773621, "mol/m^3"),
        "CYCLOPENTADIENE": (1233.4495652419334, "mol/m^3"),
        "BENZENE": (2296.609180259119, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.417958332303144, "mol/m^3"),
        "TOLUENE": (215.15076297846014, "mol/m^3"),
        "STYRENE": (168.72772384312705, "mol/m^3"),
        "OXYGEN": (4.1956053758689454e-13, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.7955915737037, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.022951046035213418,
        "2-BUTENE": 0.040221031169454004,
        "1,3-BUTADIENE": 0.8383569871868998,
        "CYCLOPENTADIENE": 0.05762161770374159,
        "BENZENE": 0.038696414825566794,
        "1,3-CYCLOHEXADIENE": 0.0006781942702425746,
        "TOLUENE": 0.0011701225437136646,
        "STYRENE": 0.0003045862651647615,
        "OXYGEN": 3.1505937743652505e-15,
    },
)

# tray 27
constantTVLiquidReactor(
    temperature=(330.016, "K"),
    initialConcentrations={
        "N-BUTANE": (120.11229366816208, "mol/m^3"),
        "2-BUTENE": (271.5391332640736, "mol/m^3"),
        "1,3-BUTADIENE": (4739.229239503595, "mol/m^3"),
        "CYCLOPENTADIENE": (1235.6774200934883, "mol/m^3"),
        "BENZENE": (2297.3548296380072, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.432052924221146, "mol/m^3"),
        "TOLUENE": (215.1825894763395, "mol/m^3"),
        "STYRENE": (168.74409175632215, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.797109533100187, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.022248508416981515,
        "2-BUTENE": 0.04180880373430059,
        "1,3-BUTADIENE": 0.8371450794079719,
        "CYCLOPENTADIENE": 0.05787238201211754,
        "BENZENE": 0.0387683889183592,
        "1,3-CYCLOHEXADIENE": 0.000679519311297604,
        "TOLUENE": 0.0011721657119635142,
        "STYRENE": 0.0003051524870080691,
    },
)

# tray 27
constantTVLiquidReactor(
    temperature=(330.016, "K"),
    initialConcentrations={
        "N-BUTANE": (120.11229366816208, "mol/m^3"),
        "2-BUTENE": (271.5391332640736, "mol/m^3"),
        "1,3-BUTADIENE": (4739.229239503595, "mol/m^3"),
        "CYCLOPENTADIENE": (1235.6774200934883, "mol/m^3"),
        "BENZENE": (2297.3548296380072, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.432052924221146, "mol/m^3"),
        "TOLUENE": (215.1825894763395, "mol/m^3"),
        "STYRENE": (168.74409175632215, "mol/m^3"),
        "OXYGEN": (1.041290264331517e-14, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.797109533100187, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.022248508416981512,
        "2-BUTENE": 0.041808803734300584,
        "1,3-BUTADIENE": 0.8371450794079718,
        "CYCLOPENTADIENE": 0.05787238201211753,
        "BENZENE": 0.038768388918359195,
        "1,3-CYCLOHEXADIENE": 0.0006795193112976039,
        "TOLUENE": 0.0011721657119635142,
        "STYRENE": 0.00030515248700806907,
        "OXYGEN": 7.82071811035935e-17,
    },
)

# tray 28
constantTVLiquidReactor(
    temperature=(330.056, "K"),
    initialConcentrations={
        "N-BUTANE": (116.31130049285665, "mol/m^3"),
        "2-BUTENE": (282.25375110728504, "mol/m^3"),
        "1,3-BUTADIENE": (4727.344315866935, "mol/m^3"),
        "CYCLOPENTADIENE": (1239.3329207070594, "mol/m^3"),
        "BENZENE": (2298.5733298425307, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.45524080124753, "mol/m^3"),
        "TOLUENE": (215.23442120145728, "mol/m^3"),
        "STYRENE": (168.7713716116473, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.799564574687793, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.0215449382207204,
        "2-BUTENE": 0.04359367733518359,
        "1,3-BUTADIENE": 0.8355394822470414,
        "CYCLOPENTADIENE": 0.058277703384645795,
        "BENZENE": 0.038881268975371155,
        "1,3-CYCLOHEXADIENE": 0.0006816142091836069,
        "TOLUENE": 0.0011753020849858986,
        "STYRENE": 0.000306013542868025,
    },
)

# tray 28
constantTVLiquidReactor(
    temperature=(330.056, "K"),
    initialConcentrations={
        "N-BUTANE": (116.31130049285665, "mol/m^3"),
        "2-BUTENE": (282.25375110728504, "mol/m^3"),
        "1,3-BUTADIENE": (4727.344315866935, "mol/m^3"),
        "CYCLOPENTADIENE": (1239.3329207070594, "mol/m^3"),
        "BENZENE": (2298.5733298425307, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.45524080124753, "mol/m^3"),
        "TOLUENE": (215.23442120145728, "mol/m^3"),
        "STYRENE": (168.7713716116473, "mol/m^3"),
        "OXYGEN": (2.5830840343137173e-16, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.799564574687793, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.0215449382207204,
        "2-BUTENE": 0.04359367733518359,
        "1,3-BUTADIENE": 0.8355394822470414,
        "CYCLOPENTADIENE": 0.058277703384645795,
        "BENZENE": 0.038881268975371155,
        "1,3-CYCLOHEXADIENE": 0.0006816142091836069,
        "TOLUENE": 0.0011753020849858986,
        "STYRENE": 0.000306013542868025,
        "OXYGEN": 1.9403334421515258e-18,
    },
)

# tray 29
constantTVLiquidReactor(
    temperature=(330.116, "K"),
    initialConcentrations={
        "N-BUTANE": (112.46756887754182, "mol/m^3"),
        "2-BUTENE": (294.18595982650936, "mol/m^3"),
        "1,3-BUTADIENE": (4710.976402671839, "mol/m^3"),
        "CYCLOPENTADIENE": (1245.3890485892446, "mol/m^3"),
        "BENZENE": (2300.6375055621347, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.49470565861793, "mol/m^3"),
        "TOLUENE": (215.3198980814761, "mol/m^3"),
        "STYRENE": (168.8159287086784, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.8037006920028, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.020838129110866365,
        "2-BUTENE": 0.04559396369476728,
        "1,3-BUTADIENE": 0.8333861642404715,
        "CYCLOPENTADIENE": 0.05894318234362573,
        "BENZENE": 0.03906575457485914,
        "1,3-CYCLOHEXADIENE": 0.0006850729570469209,
        "TOLUENE": 0.0011803416489372836,
        "STYRENE": 0.0003073914294258269,
    },
)

# tray 29
constantTVLiquidReactor(
    temperature=(330.116, "K"),
    initialConcentrations={
        "N-BUTANE": (112.46756887754182, "mol/m^3"),
        "2-BUTENE": (294.18595982650936, "mol/m^3"),
        "1,3-BUTADIENE": (4710.976402671839, "mol/m^3"),
        "CYCLOPENTADIENE": (1245.3890485892446, "mol/m^3"),
        "BENZENE": (2300.6375055621347, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.49470565861793, "mol/m^3"),
        "TOLUENE": (215.3198980814761, "mol/m^3"),
        "STYRENE": (168.8159287086784, "mol/m^3"),
        "OXYGEN": (6.33475523119125e-18, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.8037006920028, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.020838129110866365,
        "2-BUTENE": 0.04559396369476728,
        "1,3-BUTADIENE": 0.8333861642404715,
        "CYCLOPENTADIENE": 0.05894318234362573,
        "BENZENE": 0.03906575457485914,
        "1,3-CYCLOHEXADIENE": 0.0006850729570469209,
        "TOLUENE": 0.0011803416489372836,
        "STYRENE": 0.0003073914294258269,
        "OXYGEN": 4.74486662857868e-20,
    },
)

# tray 30
constantTVLiquidReactor(
    temperature=(330.208, "K"),
    initialConcentrations={
        "N-BUTANE": (108.539269710719, "mol/m^3"),
        "2-BUTENE": (307.3530366634525, "mol/m^3"),
        "1,3-BUTADIENE": (4687.797618930563, "mol/m^3"),
        "CYCLOPENTADIENE": (1255.4825950595532, "mol/m^3"),
        "BENZENE": (2304.1929800395133, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.56345089403733, "mol/m^3"),
        "TOLUENE": (215.46720930023196, "mol/m^3"),
        "STYRENE": (168.89322163209968, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.810651779857782, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.020124014992391172,
        "2-BUTENE": 0.04782373562868305,
        "1,3-BUTADIENE": 0.830438618676771,
        "CYCLOPENTADIENE": 0.06004624473445234,
        "BENZENE": 0.039377929336557364,
        "1,3-CYCLOHEXADIENE": 0.0006909785147789936,
        "TOLUENE": 0.00118879088564921,
        "STYRENE": 0.0003096872307169869,
    },
)

# tray 30
constantTVLiquidReactor(
    temperature=(330.208, "K"),
    initialConcentrations={
        "N-BUTANE": (108.539269710719, "mol/m^3"),
        "2-BUTENE": (307.3530366634525, "mol/m^3"),
        "1,3-BUTADIENE": (4687.797618930563, "mol/m^3"),
        "CYCLOPENTADIENE": (1255.4825950595532, "mol/m^3"),
        "BENZENE": (2304.1929800395133, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.56345089403733, "mol/m^3"),
        "TOLUENE": (215.46720930023196, "mol/m^3"),
        "STYRENE": (168.89322163209968, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.810651779857782, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.020124014992391172,
        "2-BUTENE": 0.04782373562868305,
        "1,3-BUTADIENE": 0.830438618676771,
        "CYCLOPENTADIENE": 0.06004624473445234,
        "BENZENE": 0.039377929336557364,
        "1,3-CYCLOHEXADIENE": 0.0006909785147789936,
        "TOLUENE": 0.00118879088564921,
        "STYRENE": 0.0003096872307169869,
        "OXYGEN": 0.0,
    },
)

# tray 31
constantTVLiquidReactor(
    temperature=(330.353, "K"),
    initialConcentrations={
        "N-BUTANE": (104.46274999662948, "mol/m^3"),
        "2-BUTENE": (321.65313682490057, "mol/m^3"),
        "1,3-BUTADIENE": (4653.934225186933, "mol/m^3"),
        "CYCLOPENTADIENE": (1272.341545650501, "mol/m^3"),
        "BENZENE": (2310.4673467643, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.685573713042956, "mol/m^3"),
        "TOLUENE": (215.7245492687993, "mol/m^3"),
        "STYRENE": (169.02689292319295, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.822592620144594, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.019395907215277484,
        "2-BUTENE": 0.05028811870718016,
        "1,3-BUTADIENE": 0.8262923073807384,
        "CYCLOPENTADIENE": 0.06188502302122857,
        "BENZENE": 0.03992021485031993,
        "1,3-CYCLOHEXADIENE": 0.0007013432608996932,
        "TOLUENE": 0.0012034204476724065,
        "STYRENE": 0.0003136651166834235,
    },
)

# tray 31
constantTVLiquidReactor(
    temperature=(330.353, "K"),
    initialConcentrations={
        "N-BUTANE": (104.46274999662948, "mol/m^3"),
        "2-BUTENE": (321.65313682490057, "mol/m^3"),
        "1,3-BUTADIENE": (4653.934225186933, "mol/m^3"),
        "CYCLOPENTADIENE": (1272.341545650501, "mol/m^3"),
        "BENZENE": (2310.4673467643, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.685573713042956, "mol/m^3"),
        "TOLUENE": (215.7245492687993, "mol/m^3"),
        "STYRENE": (169.02689292319295, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.822592620144594, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.019395907215277484,
        "2-BUTENE": 0.05028811870718016,
        "1,3-BUTADIENE": 0.8262923073807384,
        "CYCLOPENTADIENE": 0.06188502302122857,
        "BENZENE": 0.03992021485031993,
        "1,3-CYCLOHEXADIENE": 0.0007013432608996932,
        "TOLUENE": 0.0012034204476724065,
        "STYRENE": 0.0003136651166834235,
        "OXYGEN": 0.0,
    },
)

# tray 32
constantTVLiquidReactor(
    temperature=(330.591, "K"),
    initialConcentrations={
        "N-BUTANE": (100.12616232844012, "mol/m^3"),
        "2-BUTENE": (336.7579927184409, "mol/m^3"),
        "1,3-BUTADIENE": (4602.9208957288865, "mol/m^3"),
        "CYCLOPENTADIENE": (1300.4488899205228, "mol/m^3"),
        "BENZENE": (2321.67936730294, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.906176809772404, "mol/m^3"),
        "TOLUENE": (216.18466949528363, "mol/m^3"),
        "STYRENE": (169.2614996789893, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.843205775329505, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.018642397781554663,
        "2-BUTENE": 0.05297209369632084,
        "1,3-BUTADIENE": 0.8202769023870485,
        "CYCLOPENTADIENE": 0.06495599227023692,
        "BENZENE": 0.04088229513500688,
        "1,3-CYCLOHEXADIENE": 0.000719938914327269,
        "TOLUENE": 0.0012295898536788072,
        "STYRENE": 0.0003207899618259945,
    },
)

# tray 32
constantTVLiquidReactor(
    temperature=(330.591, "K"),
    initialConcentrations={
        "N-BUTANE": (100.12616232844012, "mol/m^3"),
        "2-BUTENE": (336.7579927184409, "mol/m^3"),
        "1,3-BUTADIENE": (4602.9208957288865, "mol/m^3"),
        "CYCLOPENTADIENE": (1300.4488899205228, "mol/m^3"),
        "BENZENE": (2321.67936730294, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (45.906176809772404, "mol/m^3"),
        "TOLUENE": (216.18466949528363, "mol/m^3"),
        "STYRENE": (169.2614996789893, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.843205775329505, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.018642397781554663,
        "2-BUTENE": 0.05297209369632084,
        "1,3-BUTADIENE": 0.8202769023870485,
        "CYCLOPENTADIENE": 0.06495599227023692,
        "BENZENE": 0.04088229513500688,
        "1,3-CYCLOHEXADIENE": 0.000719938914327269,
        "TOLUENE": 0.0012295898536788072,
        "STYRENE": 0.0003207899618259945,
        "OXYGEN": 0.0,
    },
)

# tray 33
constantTVLiquidReactor(
    temperature=(330.988, "K"),
    initialConcentrations={
        "N-BUTANE": (95.3503689895157, "mol/m^3"),
        "2-BUTENE": (351.9092243660339, "mol/m^3"),
        "1,3-BUTADIENE": (4524.000274273203, "mol/m^3"),
        "CYCLOPENTADIENE": (1347.0610693861104, "mol/m^3"),
        "BENZENE": (2341.9573930946412, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (46.30937307147825, "mol/m^3"),
        "TOLUENE": (217.01852373972264, "mol/m^3"),
        "STYRENE": (169.68070012248592, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.87895395865322, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.017844795253284463,
        "2-BUTENE": 0.05582168515143175,
        "1,3-BUTADIENE": 0.8112717842017053,
        "CYCLOPENTADIENE": 0.07007848135912395,
        "BENZENE": 0.042617788663668214,
        "1,3-CYCLOHEXADIENE": 0.0007538907994650472,
        "TOLUENE": 0.0012776196601531702,
        "STYRENE": 0.0003339549111679936,
    },
)

# tray 33
constantTVLiquidReactor(
    temperature=(330.988, "K"),
    initialConcentrations={
        "N-BUTANE": (95.3503689895157, "mol/m^3"),
        "2-BUTENE": (351.9092243660339, "mol/m^3"),
        "1,3-BUTADIENE": (4524.000274273203, "mol/m^3"),
        "CYCLOPENTADIENE": (1347.0610693861104, "mol/m^3"),
        "BENZENE": (2341.9573930946412, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (46.30937307147825, "mol/m^3"),
        "TOLUENE": (217.01852373972264, "mol/m^3"),
        "STYRENE": (169.68070012248592, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.87895395865322, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.017844795253284463,
        "2-BUTENE": 0.05582168515143175,
        "1,3-BUTADIENE": 0.8112717842017053,
        "CYCLOPENTADIENE": 0.07007848135912395,
        "BENZENE": 0.042617788663668214,
        "1,3-CYCLOHEXADIENE": 0.0007538907994650472,
        "TOLUENE": 0.0012776196601531702,
        "STYRENE": 0.0003339549111679936,
        "OXYGEN": 0.0,
    },
)

# tray 34
constantTVLiquidReactor(
    temperature=(331.659, "K"),
    initialConcentrations={
        "N-BUTANE": (89.84620351340955, "mol/m^3"),
        "2-BUTENE": (365.5646106132974, "mol/m^3"),
        "1,3-BUTADIENE": (4399.476827998963, "mol/m^3"),
        "CYCLOPENTADIENE": (1423.4810374369877, "mol/m^3"),
        "BENZENE": (2378.876130634689, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (47.052203531982315, "mol/m^3"),
        "TOLUENE": (218.55165160899656, "mol/m^3"),
        "STYRENE": (170.43453345797113, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.94082676483657, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.01697199460290572,
        "2-BUTENE": 0.05870868133063935,
        "1,3-BUTADIENE": 0.7974077464243368,
        "CYCLOPENTADIENE": 0.07857937501175875,
        "BENZENE": 0.04578848543926164,
        "1,3-CYCLOHEXADIENE": 0.0008167487402739007,
        "TOLUENE": 0.0013679795649824986,
        "STYRENE": 0.0003589888858415344,
    },
)

# tray 34
constantTVLiquidReactor(
    temperature=(331.659, "K"),
    initialConcentrations={
        "N-BUTANE": (89.84620351340955, "mol/m^3"),
        "2-BUTENE": (365.5646106132974, "mol/m^3"),
        "1,3-BUTADIENE": (4399.476827998963, "mol/m^3"),
        "CYCLOPENTADIENE": (1423.4810374369877, "mol/m^3"),
        "BENZENE": (2378.876130634689, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (47.052203531982315, "mol/m^3"),
        "TOLUENE": (218.55165160899656, "mol/m^3"),
        "STYRENE": (170.43453345797113, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(20.94082676483657, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.01697199460290572,
        "2-BUTENE": 0.05870868133063935,
        "1,3-BUTADIENE": 0.7974077464243368,
        "CYCLOPENTADIENE": 0.07857937501175875,
        "BENZENE": 0.04578848543926164,
        "1,3-CYCLOHEXADIENE": 0.0008167487402739007,
        "TOLUENE": 0.0013679795649824986,
        "STYRENE": 0.0003589888858415344,
        "OXYGEN": 0.0,
    },
)

# tray 35
constantTVLiquidReactor(
    temperature=(332.806, "K"),
    initialConcentrations={
        "N-BUTANE": (83.17073198249219, "mol/m^3"),
        "2-BUTENE": (374.82521216767793, "mol/m^3"),
        "1,3-BUTADIENE": (4201.379611912771, "mol/m^3"),
        "CYCLOPENTADIENE": (1546.1858266895506, "mol/m^3"),
        "BENZENE": (2446.1118807260964, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (48.42256159781611, "mol/m^3"),
        "TOLUENE": (221.39512186238892, "mol/m^3"),
        "STYRENE": (171.79761689571825, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(21.046707534506478, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.015973805415120036,
        "2-BUTENE": 0.06136762080362346,
        "1,3-BUTADIENE": 0.7755972629274721,
        "CYCLOPENTADIENE": 0.09254273137198596,
        "BENZENE": 0.05163471750416923,
        "1,3-CYCLOHEXADIENE": 0.0009343073167301805,
        "TOLUENE": 0.0015416705226263074,
        "STYRENE": 0.0004078841382727229,
    },
)

# tray 35
constantTVLiquidReactor(
    temperature=(332.806, "K"),
    initialConcentrations={
        "N-BUTANE": (83.17073198249219, "mol/m^3"),
        "2-BUTENE": (374.82521216767793, "mol/m^3"),
        "1,3-BUTADIENE": (4201.379611912771, "mol/m^3"),
        "CYCLOPENTADIENE": (1546.1858266895506, "mol/m^3"),
        "BENZENE": (2446.1118807260964, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (48.42256159781611, "mol/m^3"),
        "TOLUENE": (221.39512186238892, "mol/m^3"),
        "STYRENE": (171.79761689571825, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(21.046707534506478, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.015973805415120036,
        "2-BUTENE": 0.06136762080362346,
        "1,3-BUTADIENE": 0.7755972629274721,
        "CYCLOPENTADIENE": 0.09254273137198596,
        "BENZENE": 0.05163471750416923,
        "1,3-CYCLOHEXADIENE": 0.0009343073167301805,
        "TOLUENE": 0.0015416705226263074,
        "STYRENE": 0.0004078841382727229,
        "OXYGEN": 0.0,
    },
)

# tray 36
constantTVLiquidReactor(
    temperature=(334.773, "K"),
    initialConcentrations={
        "N-BUTANE": (74.69779078420012, "mol/m^3"),
        "2-BUTENE": (374.6533490791294, "mol/m^3"),
        "1,3-BUTADIENE": (3888.5978840396115, "mol/m^3"),
        "CYCLOPENTADIENE": (1735.980873471787, "mol/m^3"),
        "BENZENE": (2567.4981436379435, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (50.93012589930468, "mol/m^3"),
        "TOLUENE": (226.67013655376374, "mol/m^3"),
        "STYRENE": (174.25189454647168, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(21.222962360618077, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.014771592555117355,
        "2-BUTENE": 0.06329126810120088,
        "1,3-BUTADIENE": 0.7408786265971722,
        "CYCLOPENTADIENE": 0.11505594201180525,
        "BENZENE": 0.06245996852017588,
        "1,3-CYCLOHEXADIENE": 0.0011552694177442136,
        "TOLUENE": 0.0018816690516387981,
        "STYRENE": 0.0005056637451454725,
    },
)

# tray 36
constantTVLiquidReactor(
    temperature=(334.773, "K"),
    initialConcentrations={
        "N-BUTANE": (74.69779078420012, "mol/m^3"),
        "2-BUTENE": (374.6533490791294, "mol/m^3"),
        "1,3-BUTADIENE": (3888.5978840396115, "mol/m^3"),
        "CYCLOPENTADIENE": (1735.980873471787, "mol/m^3"),
        "BENZENE": (2567.4981436379435, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (50.93012589930468, "mol/m^3"),
        "TOLUENE": (226.67013655376374, "mol/m^3"),
        "STYRENE": (174.25189454647168, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(21.222962360618077, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.014771592555117355,
        "2-BUTENE": 0.06329126810120088,
        "1,3-BUTADIENE": 0.7408786265971722,
        "CYCLOPENTADIENE": 0.11505594201180525,
        "BENZENE": 0.06245996852017588,
        "1,3-CYCLOHEXADIENE": 0.0011552694177442136,
        "TOLUENE": 0.0018816690516387981,
        "STYRENE": 0.0005056637451454725,
        "OXYGEN": 0.0,
    },
)

# tray 37
constantTVLiquidReactor(
    temperature=(338.139, "K"),
    initialConcentrations={
        "N-BUTANE": (63.70273602824606, "mol/m^3"),
        "2-BUTENE": (357.3297316191429, "mol/m^3"),
        "1,3-BUTADIENE": (3410.5638858917505, "mol/m^3"),
        "CYCLOPENTADIENE": (2010.0615799236543, "mol/m^3"),
        "BENZENE": (2781.2812765361023, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (55.40756855382294, "mol/m^3"),
        "TOLUENE": (236.34902922312997, "mol/m^3"),
        "STYRENE": (178.588482214661, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(21.499023332267956, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.01325030349808012,
        "2-BUTENE": 0.06357661678422681,
        "1,3-BUTADIENE": 0.6857461810369917,
        "CYCLOPENTADIENE": 0.15015603964119445,
        "BENZENE": 0.08244072176435054,
        "1,3-CYCLOHEXADIENE": 0.0015694204143269894,
        "TOLUENE": 0.0025557606747208176,
        "STYRENE": 0.000704956186108433,
    },
)

# tray 37
constantTVLiquidReactor(
    temperature=(338.139, "K"),
    initialConcentrations={
        "N-BUTANE": (63.70273602824606, "mol/m^3"),
        "2-BUTENE": (357.3297316191429, "mol/m^3"),
        "1,3-BUTADIENE": (3410.5638858917505, "mol/m^3"),
        "CYCLOPENTADIENE": (2010.0615799236543, "mol/m^3"),
        "BENZENE": (2781.2812765361023, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (55.40756855382294, "mol/m^3"),
        "TOLUENE": (236.34902922312997, "mol/m^3"),
        "STYRENE": (178.588482214661, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(21.499023332267956, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.01325030349808012,
        "2-BUTENE": 0.06357661678422681,
        "1,3-BUTADIENE": 0.6857461810369917,
        "CYCLOPENTADIENE": 0.15015603964119445,
        "BENZENE": 0.08244072176435054,
        "1,3-CYCLOHEXADIENE": 0.0015694204143269894,
        "TOLUENE": 0.0025557606747208176,
        "STYRENE": 0.000704956186108433,
        "OXYGEN": 0.0,
    },
)

# tray 38
constantTVLiquidReactor(
    temperature=(343.793, "K"),
    initialConcentrations={
        "N-BUTANE": (49.76873152526159, "mol/m^3"),
        "2-BUTENE": (313.7865358776358, "mol/m^3"),
        "1,3-BUTADIENE": (2731.004503171835, "mol/m^3"),
        "CYCLOPENTADIENE": (2356.38843656165, "mol/m^3"),
        "BENZENE": (3139.65673594271, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (63.010918897200035, "mol/m^3"),
        "TOLUENE": (253.65900675545387, "mol/m^3"),
        "STYRENE": (186.01042152012582, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(21.885529576555818, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.0112652045060818,
        "2-BUTENE": 0.060775024310009716,
        "1,3-BUTADIENE": 0.600322240128896,
        "CYCLOPENTADIENE": 0.20155408062163221,
        "BENZENE": 0.11873604749441898,
        "1,3-CYCLOHEXADIENE": 0.002333100933240373,
        "TOLUENE": 0.003898241559296623,
        "STYRENE": 0.0011160604464241783,
    },
)

# tray 38
constantTVLiquidReactor(
    temperature=(343.793, "K"),
    initialConcentrations={
        "N-BUTANE": (49.76873152526159, "mol/m^3"),
        "2-BUTENE": (313.7865358776358, "mol/m^3"),
        "1,3-BUTADIENE": (2731.004503171835, "mol/m^3"),
        "CYCLOPENTADIENE": (2356.38843656165, "mol/m^3"),
        "BENZENE": (3139.65673594271, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (63.010918897200035, "mol/m^3"),
        "TOLUENE": (253.65900675545387, "mol/m^3"),
        "STYRENE": (186.01042152012582, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(21.885529576555818, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.0112652045060818,
        "2-BUTENE": 0.060775024310009716,
        "1,3-BUTADIENE": 0.600322240128896,
        "CYCLOPENTADIENE": 0.20155408062163221,
        "BENZENE": 0.11873604749441898,
        "1,3-CYCLOHEXADIENE": 0.002333100933240373,
        "TOLUENE": 0.003898241559296623,
        "STYRENE": 0.0011160604464241783,
        "OXYGEN": 0.0,
    },
)

# tray 39
constantTVLiquidReactor(
    temperature=(352.966, "K"),
    initialConcentrations={
        "N-BUTANE": (33.54167331649539, "mol/m^3"),
        "2-BUTENE": (238.82785874367605, "mol/m^3"),
        "1,3-BUTADIENE": (1876.799486660255, "mol/m^3"),
        "CYCLOPENTADIENE": (2674.5897623594074, "mol/m^3"),
        "BENZENE": (3706.150211624951, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (75.08934764096722, "mol/m^3"),
        "TOLUENE": (286.7476526078495, "mol/m^3"),
        "STYRENE": (201.5444804707819, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(22.79135194321763, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.008704316866445927,
        "2-BUTENE": 0.05303078090891888,
        "1,3-BUTADIENE": 0.4759828286461817,
        "CYCLOPENTADIENE": 0.2679359035430747,
        "BENZENE": 0.1821469344271036,
        "1,3-CYCLOHEXADIENE": 0.0036850586733788776,
        "TOLUENE": 0.006548257642627249,
        "STYRENE": 0.0019659192922690544,
    },
)

# tray 39
constantTVLiquidReactor(
    temperature=(352.966, "K"),
    initialConcentrations={
        "N-BUTANE": (33.54167331649539, "mol/m^3"),
        "2-BUTENE": (238.82785874367605, "mol/m^3"),
        "1,3-BUTADIENE": (1876.799486660255, "mol/m^3"),
        "CYCLOPENTADIENE": (2674.5897623594074, "mol/m^3"),
        "BENZENE": (3706.150211624951, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (75.08934764096722, "mol/m^3"),
        "TOLUENE": (286.7476526078495, "mol/m^3"),
        "STYRENE": (201.5444804707819, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(22.79135194321763, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.008704316866445927,
        "2-BUTENE": 0.05303078090891888,
        "1,3-BUTADIENE": 0.4759828286461817,
        "CYCLOPENTADIENE": 0.2679359035430747,
        "BENZENE": 0.1821469344271036,
        "1,3-CYCLOHEXADIENE": 0.0036850586733788776,
        "TOLUENE": 0.006548257642627249,
        "STYRENE": 0.0019659192922690544,
        "OXYGEN": 0.0,
    },
)

# tray 40
constantTVLiquidReactor(
    temperature=(376.836, "K"),
    initialConcentrations={
        "N-BUTANE": (11.062890662862614, "mol/m^3"),
        "2-BUTENE": (91.90037661939398, "mol/m^3"),
        "1,3-BUTADIENE": (626.9147179134244, "mol/m^3"),
        "CYCLOPENTADIENE": (2250.9608890150207, "mol/m^3"),
        "BENZENE": (5093.703679598705, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (101.87389172627194, "mol/m^3"),
        "TOLUENE": (509.37127728838135, "mol/m^3"),
        "STYRENE": (407.4964762335986, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(22.79135194321763, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.00570622279604917,
        "2-BUTENE": 0.039451919331440474,
        "1,3-BUTADIENE": 0.3185781561032965,
        "CYCLOPENTADIENE": 0.33215116275406975,
        "BENZENE": 0.28302813868378796,
        "1,3-CYCLOHEXADIENE": 0.005853572868250705,
        "TOLUENE": 0.011552205660580775,
        "STYRENE": 0.003678621802524684,
    },
)

# tray 40
constantTVLiquidReactor(
    temperature=(376.836, "K"),
    initialConcentrations={
        "N-BUTANE": (11.062890662862614, "mol/m^3"),
        "2-BUTENE": (91.90037661939398, "mol/m^3"),
        "1,3-BUTADIENE": (626.9147179134244, "mol/m^3"),
        "CYCLOPENTADIENE": (2250.9608890150207, "mol/m^3"),
        "BENZENE": (5093.703679598705, "mol/m^3"),
        "1,3-CYCLOHEXADIENE": (101.87389172627194, "mol/m^3"),
        "TOLUENE": (509.37127728838135, "mol/m^3"),
        "STYRENE": (407.4964762335986, "mol/m^3"),
        "OXYGEN": (0.0, "mol/m^3"),
    },
    terminationTime=(8000, "hr"),
    residenceTime=(22.79135194321763, "s"),
    vaporPressure=(430000, "Pa"),
    vaporMoleFractions={
        "N-BUTANE": 0.00570622279604917,
        "2-BUTENE": 0.039451919331440474,
        "1,3-BUTADIENE": 0.3185781561032965,
        "CYCLOPENTADIENE": 0.33215116275406975,
        "BENZENE": 0.28302813868378796,
        "1,3-CYCLOHEXADIENE": 0.005853572868250705,
        "TOLUENE": 0.011552205660580775,
        "STYRENE": 0.003678621802524684,
        "OXYGEN": 0.0,
    },
)


liquidVolumetricMassTransferCoefficientPowerLaw(
    prefactor=(14.5, "1/s"),
    diffusionCoefficientPower=1 / 2,
    solventViscosityPower=-1 / 6,
    solventDensityPower=-1 / 6,
)

solvation(solvent="benzene")

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
    units="si",
)

generatedSpeciesConstraints(
    # allows exceptions to the following restrictions
    allowed=["input species", "seed mechanisms", "reaction libraries"],
    # Constraints on generated species
    maximumRadicalElectrons=1,
    maximumCarbonAtoms=16,
)
