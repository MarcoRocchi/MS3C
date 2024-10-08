names.radiomics <- c(
    "MORPHOLOGICAL_Volume(IBSI:RNU0)" = "AA",
    "MORPHOLOGICAL_ApproximateVolume(IBSI:YEKZ)" = "AB",
    "MORPHOLOGICAL_Compactness1(IBSI:SKGS)" = "AC",
    "MORPHOLOGICAL_CentreOfMassShift(IBSI:KLMA)" = "AD",
    "INTENSITY-BASED_Mean(HU)IBSI:Q4LE" = "AE",
    "INTENSITY-BASED_Variance(HU)IBSI:ECT3" = "AF",
    "INTENSITY-BASED_Skewness(HU)IBSI:KE2A" = "AG",
    "INTENSITY-BASED_Kurtosis(HU)IBSI:IPH6" = "AH",
    "INTENSITY-BASED_MinimumGreyLevel(HU)IBSI:1GSF" = "AI",
    "INTENSITY-BASED_90thPercentile(HU)IBSI:8DWT" = "AJ",
    "INTENSITY-BASED_MaximumGreyLevel(HU)IBSI:84IY" = "AK",
    "INTENSITY-BASED_InterquartileRange(HU)IBSI:SALO" = "AL",
    "INTENSITY-BASED_Range(HU)IBSI:2OJQ" = "AM",
    "INTENSITY-BASED_CoefficientOfVariation(HU)IBSI:7TET" = "AN",
    "GLCM_JointMaximum(IBSI:GYBY)" = "AO",
    "GLCM_DifferenceAverage(IBSI:TF7R)" = "AP",
    "GLCM_DifferenceEntropy(IBSI:NTRS)" = "AQ",
    "GLCM_NormalisedInverseDifference(IBSI:NDRX)" = "AR",
    "GLCM_Correlation(IBSI:NI2N)" = "AS",
    "GLCM_ClusterShade(IBSI:7NFM)" = "AT",
    "GLRLM_ShortRunsEmphasis(IBSI:22OV)" = "AU",
    "GLRLM_ShortRunLowGreyLevelEmphasis(IBSI:HTZT)" = "AV",
    "GLRLM_ShortRunHighGreyLevelEmphasis(IBSI:GD3A)" = "AW",
    "NGTDM_Contrast(IBSI:65HE)" = "AX",
    "NGTDM_Strength(IBSI:1X9X)" = "AY",
    "GLSZM_SmallZoneEmphasis(IBSI:5QRC)" = "AZ",
    "GLSZM_ZoneSizeEntropy(IBSI:GU8N)" = "BA"
)

names.pre_operative <- c(
    "SEX" = "CA",
    "AGE" = "CB",
    "Sinc 1/Meta 2" = "CC",
    "N mtx" = "CD",
    "Malattia extrahep sinc fegato (0/1)" = "CE"
)

names.surgery <- c(
    "Morb severa" = "DA",
    "Infective morbidity (0/1)" = "DB",
    "r0 = Margine almeno 1 mm" = "DC"
)

names.relapse_indicator <- "Recidiva (0/1)"

names.dead_indicator <- "Morto =1"

names.post_time <- "tPost"

names.surgery_time <- "tSurgery"

names.relapse_time <- "tRelapse"

names.dead_time <- "tDead"