names.radiomics_pre <- c(
    "MORPHOLOGICAL_Volume(IBSI:RNU0)" = "AAA",
    "MORPHOLOGICAL_ApproximateVolume(IBSI:YEKZ)" = "AAB",
    "MORPHOLOGICAL_Compactness1(IBSI:SKGS)" = "AAC",
    "MORPHOLOGICAL_CentreOfMassShift(IBSI:KLMA)" = "AAD",
    "INTENSITY-BASED_Mean(HU)IBSI:Q4LE" = "AAE",
    "INTENSITY-BASED_Variance(HU)IBSI:ECT3" = "AAF",
    "INTENSITY-BASED_Skewness(HU)IBSI:KE2A" = "AAG",
    "INTENSITY-BASED_Kurtosis(HU)IBSI:IPH6" = "AAH",
    "INTENSITY-BASED_MinimumGreyLevel(HU)IBSI:1GSF" = "AAI",
    "INTENSITY-BASED_90thPercentile(HU)IBSI:8DWT" = "AAJ",
    "INTENSITY-BASED_MaximumGreyLevel(HU)IBSI:84IY" = "AAK",
    "INTENSITY-BASED_InterquartileRange(HU)IBSI:SALO" = "AAL",
    "INTENSITY-BASED_Range(HU)IBSI:2OJQ" = "AAM",
    "INTENSITY-BASED_CoefficientOfVariation(HU)IBSI:7TET" = "AAN",
    "GLCM_JointMaximum(IBSI:GYBY)" = "AAO",
    "GLCM_DifferenceAverage(IBSI:TF7R)" = "AAP",
    "GLCM_DifferenceEntropy(IBSI:NTRS)" = "AAQ",
    "GLCM_NormalisedInverseDifference(IBSI:NDRX)" = "AAR",
    "GLCM_Correlation(IBSI:NI2N)" = "AAS",
    "GLCM_ClusterShade(IBSI:7NFM)" = "AAT",
    "GLRLM_ShortRunsEmphasis(IBSI:22OV)" = "AAU",
    "GLRLM_ShortRunLowGreyLevelEmphasis(IBSI:HTZT)" = "AAV",
    "GLRLM_ShortRunHighGreyLevelEmphasis(IBSI:GD3A)" = "AAW",
    "NGTDM_Contrast(IBSI:65HE)" = "AAX",
    "NGTDM_Strength(IBSI:1X9X)" = "AAY",
    "GLSZM_SmallZoneEmphasis(IBSI:5QRC)" = "AAZ",
    "GLSZM_ZoneSizeEntropy(IBSI:GU8N)" = "ABA"
)

names.radiomics_post <- c(
    "MORPHOLOGICAL_Volume(IBSI:RNU0)" = "BAA",
    "MORPHOLOGICAL_ApproximateVolume(IBSI:YEKZ)" = "BAB",
    "MORPHOLOGICAL_Compactness1(IBSI:SKGS)" = "BAC",
    "MORPHOLOGICAL_CentreOfMassShift(IBSI:KLMA)" = "BAD",
    "INTENSITY-BASED_Mean(HU)IBSI:Q4LE" = "BAE",
    "INTENSITY-BASED_Variance(HU)IBSI:ECT3" = "BAF",
    "INTENSITY-BASED_Skewness(HU)IBSI:KE2A" = "BAG",
    "INTENSITY-BASED_Kurtosis(HU)IBSI:IPH6" = "BAH",
    "INTENSITY-BASED_MinimumGreyLevel(HU)IBSI:1GSF" = "BAI",
    "INTENSITY-BASED_90thPercentile(HU)IBSI:8DWT" = "BAJ",
    "INTENSITY-BASED_MaximumGreyLevel(HU)IBSI:84IY" = "BAK",
    "INTENSITY-BASED_InterquartileRange(HU)IBSI:SALO" = "BAL",
    "INTENSITY-BASED_Range(HU)IBSI:2OJQ" = "BAM",
    "INTENSITY-BASED_CoefficientOfVariation(HU)IBSI:7TET" = "BAN",
    "GLCM_JointMaximum(IBSI:GYBY)" = "BAO",
    "GLCM_DifferenceAverage(IBSI:TF7R)" = "BAP",
    "GLCM_DifferenceEntropy(IBSI:NTRS)" = "BAQ",
    "GLCM_NormalisedInverseDifference(IBSI:NDRX)" = "BAR",
    "GLCM_Correlation(IBSI:NI2N)" = "BAS",
    "GLCM_ClusterShade(IBSI:7NFM)" = "BAT",
    "GLRLM_ShortRunsEmphasis(IBSI:22OV)" = "BAU",
    "GLRLM_ShortRunLowGreyLevelEmphasis(IBSI:HTZT)" = "BAV",
    "GLRLM_ShortRunHighGreyLevelEmphasis(IBSI:GD3A)" = "BAW",
    "NGTDM_Contrast(IBSI:65HE)" = "BAX",
    "NGTDM_Strength(IBSI:1X9X)" = "BAY",
    "GLSZM_SmallZoneEmphasis(IBSI:5QRC)" = "BAZ",
    "GLSZM_ZoneSizeEntropy(IBSI:GU8N)" = "BBA"
)

names.pre_operative <- c(
    "AGE" = "CA",
    "N mtx" = "CB",
    "SEX_1" = "CC",
    "SEX_2" = "CD",
    "Sinc 1/Meta 2_1" = "CE",
    "Sinc 1/Meta 2_2" = "CF",
    "extrahep_1" = "CG",
    "extrahep_2" = "CH"
)

names.surgery <- c(
    "Morb severa_1" = "DA",
    "Morb severa_2" = "DB",
    "Infective morbidity_1" = "DC",
    "Infective morbidity_2" = "DD",
    "r0 = Margine almeno 1 mm_1" = "DE",
    "r0 = Margine almeno 1 mm_2" = "DF"
)

names.relapse_indicator <- "Recidiva (0/1)"

names.dead_indicator <- "Morto =1"

names.post_time <- "tPost"

names.surgery_time <- "tSurgery"

names.relapse_time <- "tRelapse"

names.dead_time <- "tDead"