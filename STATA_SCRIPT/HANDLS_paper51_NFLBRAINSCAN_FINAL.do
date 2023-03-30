
*****DIRECTORY***

cd "E:\HANDLS_PAPER51_NFLSMRI\DATA"

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\DATA_MANAGEMENT.smcl", replace

//////////////////////Preliminary data management//////////


**Create NfL_w1w3 longitudinal dataset**

use 2021-10-29_nfl13,clear
capture rename HNDid HNDID
destring HNDID, replace
destring HNDwave, replace
sort HNDID
capture rename DilCorAvgConc NFL
save NFL_w1w3, replace

**Create NFLw1 datase**
use NFL_w1w3,clear
keep if HNDwave==1
keep HNDwave HNDID NFL
rename NFL NFLw1
save NFLw1, replace


**Create NFLw3 dataset**
use NFL_w1w3,clear
keep if HNDwave==3
keep HNDwave HNDID NFL
rename NFL NFLw3
save NFLw3, replace


**Create wave 1 covariates file and merge with HANDLS-SCAN file******
use 2021-10-29_covar13,clear
capture rename HNDid HNDID
destring HNDID, replace
destring HNDwave, replace
save COVARIATESW1W3, replace
keep if HNDwave==1
sort HNDID
save COVARIATESw1, replace


use COVARIATESW1W3
keep if HNDwave==3
keep HNDID DOV
capture rename DOV DOVw3
save DOVw3, replace

use 2021-10-29_smri,clear
capture rename HNDid HNDID
destring HNDID, replace
destring HNDwave, replace
sort HNDID
save SMRI_SCAN, replace

merge HNDID using COVARIATESw1
save SMRI_SCAN_COVARIATESw1, replace
capture drop rownames
destring HNDID-rxNSAID, replace

save, replace

***Merge final merged dataset with NfLw1 and NfLw3**

use SMRI_SCAN_COVARIATESw1,clear
capture drop _merge
sort HNDID
save HANDLS_paper51_NFLBRAINSCANFINALIZED, replace


use NFLw1,clear
capture drop _merge
sort HNDID
save, replace

use NFLw3,clear
capture drop _merge
sort HNDID
save, replace

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear
merge HNDID using NFLw1
capture drop _merge
sort HNDID
save, replace
merge HNDID using NFLw3
capture drop _merge
sort HNDID
save, replace


////////////Drop wave 1 variables with >20% missing/////////

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

su

foreach z of var BasoAbs BCratio BiliIndCalc  CRPinf EosinAbs  Hcy LDL LucAbs LucPct MonoAbs MPV  NeutAbs T3 T7 UBili2 Urobili PSA UKetones2 Basophils Ca2 CBratio T4FeeeDial HDLpct LDLHDLrat T4FreeCalc UIBC HIVwestern FTA RPRtiter LymphAbs {
	
capture drop `z'
}

save HANDLS_paper51_NFLBRAINSCANFINALIZED,replace

su 

///////////////////////SELECTION PROCESS: /////////////////

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear



**Complete data on socio-demographics: SAMPLE 1**
capture drop sample*

capture drop sample1
gen sample1=1 if PovStat~=. & Age~=. 
replace sample1=0 if sample1~=1
tab sample1

save, replace

**Create days between visit 1 and date of scan**

capture drop TIME_V1SCAN
gen TIME_V1SCAN=.
replace TIME_V1SCAN=DOScan-DOV


su TIME_V1SCAN if sample1==1


sort HNDID
save, replace

use DOVw3,clear
sort HNDID
save, replace

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

merge HNDID using DOVw3
tab _merge
capture drop _merge
sort HNDID
save, replace

capture drop TIME_V2SCAN 
gen TIME_V2SCAN=DOScan-DOVw3

su TIME_V2SCAN if sample1==1

save, replace

**Complete data on total brain volume and white matter lesion: SAMPLE 2**
capture drop sample2
gen sample2=.
replace sample2=1 if sample1==1 & TOTALBrain_volM2~=. & TOTALBrain_wmlM2~=.
replace sample2=0 if sample2~=1

tab sample2


**Complete data on NfL: SAMPLES 1A-1C**


**NfL at visit 1**

capture drop sample1A
gen sample1A=.
replace sample1A=1 if NFLw1~=.
replace sample1A=0 if sample1A~=1

tab sample1A


**NFL at visit 2**


capture drop sample1B
gen sample1B=.
replace sample1B=1 if NFLw3~=.
replace sample1B=0 if sample1B~=1

tab sample1B



**NFL at visit 1 OR 2**

capture drop sample1C
gen sample1C=.
replace sample1C=1 if NFLw1~=. |  NFLw3~=.
replace sample1C=0 if sample1C~=1

tab sample1C

save, replace


**NFL at visit 1 AND 2**

capture drop sample1D
gen sample1D=.
replace sample1D=1 if NFLw1~=. &  NFLw3~=.
replace sample1D=0 if sample1D~=1

tab sample1C

save, replace


****NFL and brain volume/WML data: SAMPLES 2A-2C***


**NFL at visit 1 and brain volume/WML data**

capture drop sample2A
gen sample2A=.
replace sample2A=1 if sample1A==1 & sample2==1
replace sample2A=0 if sample2A~=1

tab sample2A


**NFL at visit 2 and brain volume/WML data**

capture drop sample2B
gen sample2B=.
replace sample2B=1 if sample1B==1 & sample2==1
replace sample2B=0 if sample2B~=1

tab sample2B

**NFL at visit 1 OR 2 and brain volume/WML data**

capture drop sample2C
gen sample2C=.
replace sample2C=1 if sample1C==1 & sample2==1
replace sample2C=0 if sample2C~=1

tab sample2C

**NFL at visit 1 and 2 and brain volume/WML data**

capture drop sample2D
gen sample2D=.
replace sample2D=1 if sample1D==1 & sample2==1
replace sample2D=0 if sample2D~=1

tab sample2D



save, replace

/////////////FINAL SAMPLE: NFL at wave 1 AND wave 3, plus data on brain volume and WML**

capture drop sample_final
gen sample_final=.
replace sample_final=sample2D

tab sample_final

********************************Z-SCORING IN SAMPLE 1**********************************************

*****LIST OF KEY VOLUMETRIC VARIABLES: 
su ICV_volM2 TOTALBrain_volM2 GM_volM2 WM_volM2 FRONTAL_GM_L_volM2 FRONTAL_WM_L_volM2 OCCIPITAL_GM_L_volM2 OCCIPITAL_WM_L_volM2 PARIETAL_GM_L_volM2 PARIETAL_WM_L_volM2 TEMPORAL_GM_L_volM2 TEMPORAL_WM_L_volM2 FRONTAL_GM_R_volM2 FRONTAL_WM_R_volM2 OCCIPITAL_GM_R_volM2 OCCIPITAL_WM_R_volM2 PARIETAL_GM_R_volM2 PARIETAL_WM_R_volM2 TEMPORAL_GM_R_volM2 TEMPORAL_WM_R_volM2 Left_Hippocampus_volM2 Right_Hippocampus_volM2

*****LIST OF KEY WML VARIABLES**
su TOTALBrain_wmlM2


****LIST OF VOLUMETRIC SMALL ROIs: 
su _3rd_Ventricle_volM2 _4th_Ventricle_volM2 Right_Accumbens_Area_volM2 Left_Accumbens_Area_volM2 Right_Amygdala_volM2 Left_Amygdala_volM2 Brain_Stem_volM2 Right_Caudate_volM2 Left_Caudate_volM2 Right_Cerebellum_Exterior_volM2 Left_Cerebellum_Exterior_volM2 Right_Cerebellum_White_Matter_vo Left_Cerebellum_White_Matter_vol Right_Hippocampus_volM2 Left_Hippocampus_volM2 Right_Inf_Lat_Vent_volM2 Left_Inf_Lat_Vent_volM2 Right_Lateral_Ventricle_volM2 Left_Lateral_Ventricle_volM2 Right_Pallidum_volM2 Left_Pallidum_volM2 Right_Putamen_volM2 Left_Putamen_volM2 Right_Thalamus_Proper_volM2 Left_Thalamus_Proper_volM2 Right_Ventral_DC_volM2 Left_Ventral_DC_volM2 Cerebellar_Vermal_Lobules_I_V_vo Cerebellar_Vermal_Lobules_VI_VII Cerebellar_Vermal_Lobules_VIII_X Left_Basal_Forebrain_volM2 Right_Basal_Forebrain_volM2 frontal_lobe_WM_right_volM2 frontal_lobe_WM_left_volM2 occipital_lobe_WM_right_volM2 occipital_lobe_WM_left_volM2 parietal_lobe_WM_right_volM2 parietal_lobe_WM_left_volM2 temporal_lobe_WM_right_volM2 temporal_lobe_WM_left_volM2 fornix_right_volM2 fornix_left_volM2 anterior_limb_of_internal_capsul anterior_limb_of_internal_capsu0 posterior_limb_of_internal_capsu posterior_limb_of_internal_caps0 corpus_callosum_volM2 Right_ACgG_anterior_cingulate_gy Left_ACgG_anterior_cingulate_gyr Right_AIns_anterior_insula_volM2 Left_AIns_anterior_insula_volM2 Right_AOrG_anterior_orbital_gyru Left_AOrG_anterior_orbital_gyrus Right_AnG_angular_gyrus_volM2 Left_AnG_angular_gyrus_volM2 Right_Calc_calcarine_cortex_volM Left_Calc_calcarine_cortex_volM2 Right_CO_central_operculum_volM2 Left_CO_central_operculum_volM2 Right_Cun_cuneus_volM2 Left_Cun_cuneus_volM2 Right_Ent_entorhinal_area_volM2 Left_Ent_entorhinal_area_volM2 Right_FO_frontal_operculum_volM2 Left_FO_frontal_operculum_volM2 Right_FRP_frontal_pole_volM2 Left_FRP_frontal_pole_volM2 Right_FuG_fusiform_gyrus_volM2 Left_FuG_fusiform_gyrus_volM2 Right_GRe_gyrus_rectus_volM2 Left_GRe_gyrus_rectus_volM2 Right_IOG_inferior_occipital_gyr Left_IOG_inferior_occipital_gyru Right_ITG_inferior_temporal_gyru Left_ITG_inferior_temporal_gyrus Right_LiG_lingual_gyrus_volM2 Left_LiG_lingual_gyrus_volM2 Right_LOrG_lateral_orbital_gyrus Left_LOrG_lateral_orbital_gyrus_ Right_MCgG_middle_cingulate_gyru Left_MCgG_middle_cingulate_gyrus Right_MFC_medial_frontal_cortex_ Left_MFC_medial_frontal_cortex_v Right_MFG_middle_frontal_gyrus_v Left_MFG_middle_frontal_gyrus_vo Right_MOG_middle_occipital_gyrus Left_MOG_middle_occipital_gyrus_ Right_MOrG_medial_orbital_gyrus_ Left_MOrG_medial_orbital_gyrus_v Right_MPoG_postcentral_gyrus_med Left_MPoG_postcentral_gyrus_medi Right_MPrG_precentral_gyrus_medi Left_MPrG_precentral_gyrus_media Right_MSFG_superior_frontal_gyru Left_MSFG_superior_frontal_gyrus Right_MTG_middle_temporal_gyrus_ Left_MTG_middle_temporal_gyrus_v Right_OCP_occipital_pole_volM2 Left_OCP_occipital_pole_volM2 Right_OFuG_occipital_fusiform_gy Left_OFuG_occipital_fusiform_gyr Right_OpIFG_opercular_part_of_th Left_OpIFG_opercular_part_of_the Right_OrIFG_orbital_part_of_the_ Left_OrIFG_orbital_part_of_the_i Right_PCgG_posterior_cingulate_g Left_PCgG_posterior_cingulate_gy Right_PCu_precuneus_volM2 Left_PCu_precuneus_volM2 Right_PHG_parahippocampal_gyrus_ Left_PHG_parahippocampal_gyrus_v Right_PIns_posterior_insula_volM Left_PIns_posterior_insula_volM2 Right_PO_parietal_operculum_volM Left_PO_parietal_operculum_volM2 Right_PoG_postcentral_gyrus_volM Left_PoG_postcentral_gyrus_volM2 Right_POrG_posterior_orbital_gyr Left_POrG_posterior_orbital_gyru Right_PP_planum_polare_volM2 Left_PP_planum_polare_volM2 Right_PrG_precentral_gyrus_volM2 Left_PrG_precentral_gyrus_volM2 Right_PT_planum_temporale_volM2 Left_PT_planum_temporale_volM2 Right_SCA_subcallosal_area_volM2 Left_SCA_subcallosal_area_volM2 Right_SFG_superior_frontal_gyrus Left_SFG_superior_frontal_gyrus_ Right_SMC_supplementary_motor_co Left_SMC_supplementary_motor_cor Right_SMG_supramarginal_gyrus_vo Left_SMG_supramarginal_gyrus_vol Right_SOG_superior_occipital_gyr Left_SOG_superior_occipital_gyru Right_SPL_superior_parietal_lobu Left_SPL_superior_parietal_lobul Right_STG_superior_temporal_gyru Left_STG_superior_temporal_gyrus Right_TMP_temporal_pole_volM2 Left_TMP_temporal_pole_volM2 Right_TrIFG_triangular_part_of_t Left_TrIFG_triangular_part_of_th Right_TTG_transverse_temporal_gy Left_TTG_transverse_temporal_gyr

****LIST OF WML SMALL ROIs:

su _3rd_Ventricle_wmlM2 _4th_Ventricle_wmlM2 Right_Accumbens_Area_wmlM2 Left_Accumbens_Area_wmlM2 Right_Amygdala_wmlM2 Left_Amygdala_wmlM2 Brain_Stem_wmlM2 Right_Caudate_wmlM2 Left_Caudate_wmlM2 Right_Cerebellum_Exterior_wmlM2 Left_Cerebellum_Exterior_wmlM2 Right_Cerebellum_White_Matter_wm Left_Cerebellum_White_Matter_wml Right_Hippocampus_wmlM2 Left_Hippocampus_wmlM2 Right_Inf_Lat_Vent_wmlM2 Left_Inf_Lat_Vent_wmlM2 Right_Lateral_Ventricle_wmlM2 Left_Lateral_Ventricle_wmlM2 Right_Pallidum_wmlM2 Left_Pallidum_wmlM2 Right_Putamen_wmlM2 Left_Putamen_wmlM2 Right_Thalamus_Proper_wmlM2 Left_Thalamus_Proper_wmlM2 Right_Ventral_DC_wmlM2 Left_Ventral_DC_wmlM2 Cerebellar_Vermal_Lobules_I_V_wm Cerebellar_Vermal_Lobules_VI_VI0 Cerebellar_Vermal_Lobules_VIII_0 Left_Basal_Forebrain_wmlM2 Right_Basal_Forebrain_wmlM2 frontal_lobe_WM_right_wmlM2 frontal_lobe_WM_left_wmlM2 occipital_lobe_WM_right_wmlM2 occipital_lobe_WM_left_wmlM2 parietal_lobe_WM_right_wmlM2 parietal_lobe_WM_left_wmlM2 temporal_lobe_WM_right_wmlM2 temporal_lobe_WM_left_wmlM2 fornix_right_wmlM2 fornix_left_wmlM2 anterior_limb_of_internal_capsu1 anterior_limb_of_internal_capsu2 posterior_limb_of_internal_caps1 posterior_limb_of_internal_caps2 corpus_callosum_wmlM2 Right_ACgG_anterior_cingulate_g0 Left_ACgG_anterior_cingulate_gy0 Right_AIns_anterior_insula_wmlM2 Left_AIns_anterior_insula_wmlM2 Right_AOrG_anterior_orbital_gyr0 Left_AOrG_anterior_orbital_gyru0 Right_AnG_angular_gyrus_wmlM2 Left_AnG_angular_gyrus_wmlM2 Right_Calc_calcarine_cortex_wmlM Left_Calc_calcarine_cortex_wmlM2 Right_CO_central_operculum_wmlM2 Left_CO_central_operculum_wmlM2 Right_Cun_cuneus_wmlM2 Left_Cun_cuneus_wmlM2 Right_Ent_entorhinal_area_wmlM2 Left_Ent_entorhinal_area_wmlM2 Right_FO_frontal_operculum_wmlM2 Left_FO_frontal_operculum_wmlM2 Right_FRP_frontal_pole_wmlM2 Left_FRP_frontal_pole_wmlM2 Right_FuG_fusiform_gyrus_wmlM2 Left_FuG_fusiform_gyrus_wmlM2 Right_GRe_gyrus_rectus_wmlM2 Left_GRe_gyrus_rectus_wmlM2 Right_IOG_inferior_occipital_gy0 Left_IOG_inferior_occipital_gyr0 Right_ITG_inferior_temporal_gyr0 Left_ITG_inferior_temporal_gyru0 Right_LiG_lingual_gyrus_wmlM2 Left_LiG_lingual_gyrus_wmlM2 Right_LOrG_lateral_orbital_gyru0 Left_LOrG_lateral_orbital_gyrus0 Right_MCgG_middle_cingulate_gyr0 Left_MCgG_middle_cingulate_gyru0 Right_MFC_medial_frontal_cortex0 Left_MFC_medial_frontal_cortex_w Right_MFG_middle_frontal_gyrus_w Left_MFG_middle_frontal_gyrus_wm Right_MOG_middle_occipital_gyru0 Left_MOG_middle_occipital_gyrus0 Right_MOrG_medial_orbital_gyrus0 Left_MOrG_medial_orbital_gyrus_w Right_MPoG_postcentral_gyrus_me0 Left_MPoG_postcentral_gyrus_med0 Right_MPrG_precentral_gyrus_med0 Left_MPrG_precentral_gyrus_medi0 Right_MSFG_superior_frontal_gyr0 Left_MSFG_superior_frontal_gyru0 Right_MTG_middle_temporal_gyrus0 Left_MTG_middle_temporal_gyrus_w Right_OCP_occipital_pole_wmlM2 Left_OCP_occipital_pole_wmlM2 Right_OFuG_occipital_fusiform_g0 Left_OFuG_occipital_fusiform_gy0 Right_OpIFG_opercular_part_of_t0 Left_OpIFG_opercular_part_of_th0 Right_OrIFG_orbital_part_of_the0 Left_OrIFG_orbital_part_of_the_0 Right_PCgG_posterior_cingulate_0 Left_PCgG_posterior_cingulate_g0 Right_PCu_precuneus_wmlM2 Left_PCu_precuneus_wmlM2 Right_PHG_parahippocampal_gyrus0 Left_PHG_parahippocampal_gyrus_w Right_PIns_posterior_insula_wmlM Left_PIns_posterior_insula_wmlM2 Right_PO_parietal_operculum_wmlM Left_PO_parietal_operculum_wmlM2 Right_PoG_postcentral_gyrus_wmlM Left_PoG_postcentral_gyrus_wmlM2 Right_POrG_posterior_orbital_gy0 Left_POrG_posterior_orbital_gyr0 Right_PP_planum_polare_wmlM2 Left_PP_planum_polare_wmlM2 Right_PrG_precentral_gyrus_wmlM2 Left_PrG_precentral_gyrus_wmlM2 Right_PT_planum_temporale_wmlM2 Left_PT_planum_temporale_wmlM2 Right_SCA_subcallosal_area_wmlM2 Left_SCA_subcallosal_area_wmlM2 Right_SFG_superior_frontal_gyru0 Left_SFG_superior_frontal_gyrus0 Right_SMC_supplementary_motor_c0 Left_SMC_supplementary_motor_co0 Right_SMG_supramarginal_gyrus_wm Left_SMG_supramarginal_gyrus_wml Right_SOG_superior_occipital_gy0 Left_SOG_superior_occipital_gyr0 Right_SPL_superior_parietal_lob0 Left_SPL_superior_parietal_lobu0 Right_STG_superior_temporal_gyr0 Left_STG_superior_temporal_gyru0 Right_TMP_temporal_pole_wmlM2 Left_TMP_temporal_pole_wmlM2 Right_TrIFG_triangular_part_of_0 Left_TrIFG_triangular_part_of_t0 Right_TTG_transverse_temporal_g0 Left_TTG_transverse_temporal_gy0

*ANALYSIS A: GLOBAL BRAIN VOLUMES: GM/WM**

capture drop TOTALBRAIN
gen TOTALBRAIN=TOTALBrain_volM2

capture drop GM
gen GM=GM_volM2

capture drop WM
gen WM=WM_volM2

save, replace


capture drop zTOTALBRAIN zGM zWM

foreach var of varlist TOTALBRAIN GM WM  {
  egen z`var' = std(`var') if sample1==1
}


save, replace

*ANALYSIS Aprime: LARGE REGION, GM/WM*

capture drop FRONTAL_GM_L OCCIPITAL_GM_L PARIETAL_GM_L TEMPORAL_GM_L FRONTAL_WM_L OCCIPITAL_WM_L PARIETAL_WM_L TEMPORAL_WM_L ///
FRONTAL_GM_R OCCIPITAL_GM_R PARIETAL_GM_R TEMPORAL_GM_R FRONTAL_WM_R OCCIPITAL_WM_R PARIETAL_WM_R TEMPORAL_WM_R

gen FRONTAL_GM_L=FRONTAL_GM_L_volM2
gen OCCIPITAL_GM_L=OCCIPITAL_GM_L_volM2
gen PARIETAL_GM_L=PARIETAL_GM_L_volM2
gen TEMPORAL_GM_L=TEMPORAL_GM_L_volM2
gen FRONTAL_WM_L=FRONTAL_GM_L_volM2
gen OCCIPITAL_WM_L=OCCIPITAL_GM_L_volM2
gen PARIETAL_WM_L=PARIETAL_GM_L_volM2
gen TEMPORAL_WM_L=TEMPORAL_WM_R_volM2
gen FRONTAL_GM_R=FRONTAL_GM_R_volM2
gen OCCIPITAL_GM_R=OCCIPITAL_GM_R_volM2
gen PARIETAL_GM_R=PARIETAL_GM_L_volM2
gen TEMPORAL_GM_R=TEMPORAL_GM_L_volM2
gen FRONTAL_WM_R=FRONTAL_WM_R_volM2
gen OCCIPITAL_WM_R=OCCIPITAL_WM_R_volM2
gen PARIETAL_WM_R=PARIETAL_WM_R_volM2
gen TEMPORAL_WM_R=TEMPORAL_WM_R_volM2

save, replace


capture drop zFRONTAL_GM_L zOCCIPITAL_GM_L zPARIETAL_GM_L zTEMPORAL_GM_L zFRONTAL_WM_L zOCCIPITAL_WM_L zPARIETAL_WM_L zTEMPORAL_WM_L
capture drop zFRONTAL_GM_R zOCCIPITAL_GM_R zPARIETAL_GM_R zTEMPORAL_GM_R zFRONTAL_WM_R zOCCIPITAL_WM_R zPARIETAL_WM_R zTEMPORAL_WM_R


foreach var of varlist FRONTAL_GM_L OCCIPITAL_GM_L PARIETAL_GM_L TEMPORAL_GM_L FRONTAL_WM_L OCCIPITAL_WM_L PARIETAL_WM_L TEMPORAL_WM_L ///
FRONTAL_GM_R OCCIPITAL_GM_R PARIETAL_GM_R TEMPORAL_GM_R FRONTAL_WM_R OCCIPITAL_WM_R PARIETAL_WM_R TEMPORAL_WM_R  {
  egen z`var' = std(`var') if sample1==1
}

*ANALYSIS B: Hippocampus*

capture drop Right_Hippocampus
capture drop Left_Hippocampus 
capture drop Right_Hippocampuspct
capture drop Left_Hippocampuspct

gen Right_Hippocampus=Right_Hippocampus_volM2
gen Left_Hippocampus=Left_Hippocampus_volM2
gen Right_Hippocampuspct=Right_Hippocampus*100/ICV_volM2
gen Left_Hippocampuspct=Left_Hippocampus_volM2*100/ICV_volM2

save, replace

capture drop zRight_Hippocampus
egen zRight_Hippocampus=std(Right_Hippocampus) if sample1==1


capture drop zLeft_Hippocampus
egen zLeft_Hippocampus=std(Left_Hippocampus) if sample1==1

capture drop zRight_Hippocampuspct
egen zRight_Hippocampuspct=std(Right_Hippocampuspct) if sample1==1


capture drop zLeft_Hippocampuspct
egen zLeft_Hippocampuspct=std(Left_Hippocampuspct) if sample1==1


save, replace


*ANALYSIS C: WHITE MATTER LESIONS**

capture drop Lesion_Volume Lesion_Volumepct

gen Lesion_Volume=TOTALBrain_wmlM2
replace Lesion_Volume=0.00000001 if TOTALBrain_wmlM2==0
gen Lesion_Volumepct=Lesion_Volume*100/ICV_volM2

capture drop zLesion_Volume
egen zLesion_Volume=std(Lesion_Volume) if sample1==1


capture drop LnLesion_Volume
gen LnLesion_Volume=ln(Lesion_Volume)


capture drop LnLesion_Volumepct
gen LnLesion_Volumepct=ln((Lesion_Volume*100/ICV_volM2))


save, replace


****************************************************COVARIATE DATA MANAGEMENT**************************************

/////////////////////EDUCATION//////////////

**EDUCATION**

capture drop edubr
gen edubr=.
replace edubr=1 if Education>=1 & Education<=8
replace edubr=2 if Education>=9 & Education<=12
replace edubr=3 if Education>=13 & Education~=.

tab edubr if HNDwave==1
tab edubr Education

save, replace


**Current drug use**

tab1 MarijCurr CokeCurr OpiateCurr

capture drop currdrugs
gen currdrugs=.
replace currdrugs=1 if MarijCurr==1 | CokeCurr==1 | OpiateCurr==1
replace currdrugs=0 if currdrugs~=1 & MarijCurr~=. & CokeCurr~=. & OpiateCurr~=.
replace currdrugs=9 if currdrugs==.

tab currdrugs

tab currdrugs MarijCurr
tab currdrugs CokeCurr
tab currdrugs OpiateCurr


save, replace

**Current smoking status**

tab  CigaretteCurr
su CigaretteCurr

capture drop smoke
gen smoke=.
replace smoke=1 if CigaretteCurr==1 
replace smoke=0 if CigaretteCurr~=1 & CigaretteCurr~=.
replace smoke=9 if smoke==.

tab smoke CigaretteCurr

capture drop smoke1 smoke9
gen smoke1=1 if smoke==1
replace smoke1=0 if smoke~=1

gen smoke9=1 if smoke==9
replace smoke9=0 if smoke~=9

sort HNDID

save, replace


**Self-rated health**

capture drop SRH
gen SRH=.
replace SRH=1 if SF01==1 | SF01==2
replace SRH=2 if SF01==3
replace SRH=3 if SF01==4 | SF01==5


tab SRH

save, replace

**BMI and WAIST SIZE**

su BMI

su WaistSize

su WHR

save, replace

**Hypertension and Diabetes**

su dxHTN dxDiabetes

tab1 dxHTN dxDiabetes

destring dxHTN, replace

recode dxHTN 2=1 1=0
label define yesno 0 "no" 1 "yes"

label val dxHTN yesno

recode dxDiabetes 1=0 2=1 3=2
capture label drop diabetes
label define diabetes 0 "no" 1 "pre-diabetes" 2 "diabetes"
label val dxDiabetes diabetes

**Dyslipidemia**
su CVhighChol rxLipid
tab1 CVhighChol rxLipid

capture drop HypercholesterolemiarxCV
gen HypercholesterolemiarxCV=.
replace HypercholesterolemiarxCV=1 if CVhighChol==2 | rxLipid==2
replace HypercholesterolemiarxCV=0 if HypercholesterolemiarxCV~=1 & CVhighChol~=. & rxLipid~=.

tab HypercholesterolemiarxCV 

tab HypercholesterolemiarxCV  CVhighChol
tab HypercholesterolemiarxCV rxLipid

save, replace


**CVD**

tab1  CVaFib CVangina CVcad CVchf CVmi

capture drop cvdbr
gen cvdbr=.
replace cvdbr=1 if CVaFib==2 | CVangina==2 | CVcad==2 | CVchf==2 | CVmi==2
replace cvdbr=0 if cvdbr~=1 & CVaFib~=. & CVangina~=. & CVcad~=. & CVchf~=. & CVmi~=.


tab cvdbr

**Inflammatory conditions***

tab1 MultipleSclerosis Lupus Gout RheumatoidArthritis Psoriasis ThyroidHyper ThyroidHypo Crohns

capture drop inflamcond
gen inflamcond=.
replace inflamcond=1 if MultipleSclerosis==2 | Lupus==2 | Gout==2 | RheumatoidArthritis==2 | Psoriasis==2 | ThyroidHyper==2 | ThyroidHypo==2 | Crohns==2
replace inflamcond=0 if inflamcond~=1 & MultipleSclerosis~=. & Lupus~=. & Gout~=. & RheumatoidArthritis~=. & Psoriasis~=. & ThyroidHyper~=. & ThyroidHypo~=. & Crohns~=.

tab inflamcond

save, replace

**NSAIDs**

tab rxNSAID

**HEI**

su hei2010_total_score
histogram hei2010_total_score


save, replace



*****DASH******

su DASH_TotalScore
histogram DASH_TotalScore

save, replace

******MAR*******

capture drop mar
egen mar=rmean(NutrientAdequacyRatio*)

capture drop mar100
gen mar100=mar*100
histogram mar100

save, replace


///////////////////////BIOCHEMICAL AND HEMATOLOGICAL INDICES////////////////////////

su AlbGlob Albumin ALP ALT Amylase AST B12 BasoPct BiliDir BiliTot BUN Calcium Chol Cl CO2 Creatinine EosinPct ESR Fe Ferritin FeSat Folate GGT Globulin Glucose HBV Hct HCV HDL Hgb HgbA1C Insulin K LDH LDLcalc LymphPct MCH MCHC MCV Mg MonoPct Na NeutPct Phosphate Platelets Protein RBC RDW RPR T3uptake T4Feee T4tot TIBC Triglyc TSH UpH UricAcid USpecGrav WBC HIV ChoHDLRat CRP Eosinophils Lipase Lymphocytes Monoocytes Neutrophils UCreatinine UMicroAlb VLDL TotalD




save, replace




****REMOVE "9" AND CHANGE TO MISSING FROM EDUBR, SMOKE, CURRDRUGS AND MARRIED****

recode edubr 9=.
recode smoke 9=.
recode currdrugs 9=.


save, replace

/////////////////////////////////FULL LIST OF WAVE 1 COVARIATES, ADDSTUB////////////////////////

addstub Age edubr currdrugs smoke SRH BMI WaistSize WHR dxHTN dxDiabetes HypercholesterolemiarxCV cvdbr inflamcond rxNSAID, stub(w1)
addstub hei2010_total_score DASH_TotalScore mar100 AlbGlob Albumin ALP ALT Amylase AST B12 BasoPct BiliDir BiliTot BUN Calcium Chol Cl CO2 Creatinine EosinPct, stub(w1) //
addstub ESR Fe Ferritin FeSat Folate GGT Globulin Glucose HBV Hct HCV HDL Hgb HgbA1C Insulin K LDH LDLcalc LymphPct MCH MCHC MCV Mg MonoPct Na NeutPct, stub(w1)
addstub Phosphate Platelets Protein RBC RDW RPR T3uptake T4Feee T4tot TIBC Triglyc TSH UpH UricAcid USpecGrav WBC HIV ChoHDLRat CRP Eosinophils Lipase, stub(w1)
addstub Lymphocytes Monoocytes Neutrophils UCreatinine UMicroAlb VLDL TotalD, stub(w1)

save, replace 



///////////////////////////////////EMPIRICAL BAYES ESTIMATORS/////////////////////////////////////

**NFL_w1w3
use NFL_w1w3,clear
sort HNDID
save, replace


**Take Loge of NFL**
capture drop LnNFL
gen LnNFL=ln(NFL)

histogram LnNFL

save, replace


**COVARIATESW1: Extract Age Sex Race and PovStat and merge with NFLw1w3

use COVARIATESw1,clear
sort HNDID
capture drop _merge
save, replace

keep HNDID HNDwave Age Sex Race PovStat
keep if HNDwave==1
save DEMOw1, replace

addstub Age, stub(w1)
sort HNDID
save, replace

**Merge with w3Age and LnNFLw3**
use COVARIATESW1W3, clear
keep HNDID HNDwave Age
keep if HNDwave==3
capture drop HNDwave
sort HNDID
save w3Age, replace
addstub Age, stub(w3)
save, replace

use NFL_w1w3, clear
sort HNDID
save, replace
keep HNDID HNDwave NFL
keep if HNDwave==1
capture rename NFL NFLw1
sort HNDID
save NFLw1, replace


use NFL_w1w3, clear
sort HNDID
save, replace
keep HNDID HNDwave NFL
keep if HNDwave==3
capture rename NFL NFLw3
sort HNDID
save NFLw3, replace

use NFL_w1w3,clear
merge HNDID using DEMOw1
tab _merge
capture drop _merge
sort HNDID
save NFL_MIXEDMODEL_LONG, replace
merge HNDID using w3Age
tab _merge
capture drop _merge
sort HNDID
save, replace
merge HNDID using NFLw1
tab _merge
capture drop _merge
sort HNDID
save, replace
merge HNDID using NFLw3
tab _merge
capture drop _merge
sort HNDID
save, replace

capture drop LnNFLw3
gen LnNFLw3=ln(NFLw3)
capture drop LnNFLw1
gen LnNFLw1=ln(NFLw1)

save, replace

**Generate observed delaNFL, Loge transformed**


capture drop deltaLnNFL
gen deltaLnNFL=(LnNFLw3-LnNFLw1)/(w3Age-w1Age)
su deltaLnNFL
histogram deltaLnNFL

**Restrict to non-missing longitudinal NFL**

keep if NFL~=.

save, replace

**Re-generate Age longitudinal using w1Age, w3Age and HNDwave**
capture drop Agelong
gen Agelong=.
replace Agelong=w1Age if HNDwave==1
replace Agelong=w3Age if HNDwave==3

save, replace

**Generate timew1w3**

capture drop timew1w3
gen timew1w3=Agelong-w1Age

save, replace

**Run mixed models and extract bayes1NFL**

**MEAN TIME +/- SD BETWEEN VISITS 1 AND 2**

su timew1w3 if LnNFL~=. & HNDwave==3

capture drop w1Agecenter50
gen w1Agecenter50=w1Age-50

*************LnNFL mixed model**************


xi:xtmixed LnNFL c.timew1w3##c.w1Agecenter50 c.timew1w3##Sex c.timew1w3##Race c.timew1w3##PovStat || HNDID: timew1w3 ,  cov(un) variance mle


**Performing EM optimization: 
**
**Performing gradient-based optimization: 
**
**Iteration 0:   log likelihood = -743.74158  
**Iteration 1:   log likelihood = -743.69449  
**Iteration 2:   log likelihood = -743.69446  
**
**Computing standard errors:
**
**Mixed-effects ML regression                     Number of obs     =      1,403
**Group variable: HNDID                           Number of groups  =        729
**                                                Obs per group:
**                                                              min =          1
**                                                              avg =        1.9
**                                                              max =          2
**                                                Wald chi2(9)      =     452.69
**Log likelihood = -743.69446                     Prob > chi2       =     0.0000
**
**--------------------------------------------------------------------------------------------
**                     LnNFL | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
**---------------------------+----------------------------------------------------------------
**                  timew1w3 |   .0367367   .0069329     5.30   0.000     .0231485    .0503249
**             w1Agecenter50 |   .0260676   .0018511    14.08   0.000     .0224396    .0296956
**                           |
**c.timew1w3#c.w1Agecenter50 |   .0005684   .0004078     1.39   0.163    -.0002309    .0013676
**                           |
**                  timew1w3 |          0  (omitted)
**                           |
**                       Sex |
**                      Men  |   .0944587   .0334616     2.82   0.005     .0288753    .1600422
**                           |
**            Sex#c.timew1w3 |
**                      Men  |   .0080145    .007381     1.09   0.278     -.006452     .022481
**                           |
**                  timew1w3 |          0  (omitted)
**                           |
**                      Race |
**                    AfrAm  |  -.1836598   .0333263    -5.51   0.000    -.2489781   -.1183414
**                           |
**           Race#c.timew1w3 |
**                    AfrAm  |   .0075765   .0074284     1.02   0.308     -.006983    .0221359
**                           |
**                  timew1w3 |          0  (omitted)
**                           |
**                   PovStat |
**                    Below  |   .0373991    .036728     1.02   0.309    -.0345865    .1093846
**                           |
**        PovStat#c.timew1w3 |
**                    Below  |   .0044728   .0080016     0.56   0.576      -.01121    .0201555
**                           |
**                     _cons |   2.087789   .0301934    69.15   0.000     2.028611    2.146967
**--------------------------------------------------------------------------------------------
**
**------------------------------------------------------------------------------
**  Random-effects parameters  |   Estimate   Std. err.     [95% conf. interval]
**-----------------------------+------------------------------------------------
**HNDID: Unstructured          |
*               var(timew1w3) |   .0062101   .0009051       .004667    .0082635
**                  var(_cons) |   .1672103   .0124954      .1444288    .1935851
**         cov(timew1w3,_cons) |  -.0083442   .0023229      -.012897   -.0037914
**-----------------------------+------------------------------------------------
**               var(Residual) |   .0249631   .0072747      .0141008     .044193
**------------------------------------------------------------------------------
**LR test vs. linear model: chi2(3) = 345.12                Prob > chi2 = 0.0000
**
**Note: LR test is conservative and provided only for reference.

capture drop e_TIMELnNFL e_consLnNFL 
predict  e_TIMELnNFL e_consLnNFL, reffects level(HNDID)

estat ic

***bayes1: empirical bayes estimator of level-1 coefficient**

capture drop bayes1LnNFL
gen bayes1LnNFL=     .0083177  +    .0005684*w1Agecenter50 +    .0080145*Sex +   .0075765*Race +      .0044728*PovStat+ e_TIMELnNFL



save   , replace


**Histogram of empirical bayes etimator for LnNFL**************************************

histogram bayes1LnNFL

**Scatterplot of LnNFLw1 LnNFLw3 and bayes1LnNFL**

graph matrix LnNFLw1 LnNFLw3 bayes1LnNFL

save, replace


keep HNDID HNDwave bayes1LnNFL

keep if HNDwave==1

save bayes1LnNFL, replace




**Take mean of bayes1NFL across waves. 

collapse (mean) bayes1LnNFL, by(HNDID)

save bayes1LnNFL_collapse, replace
capture drop HNDwave
sort HNDID
save, replace




**Merge resulting file with with final dataset. 

use HANDLS_paper51_NFLBRAINSCANFINALIZED, clear
sort HNDID
capture drop _merge
save, replace

merge HNDID using bayes1LnNFL_collapse
tab _merge
capture drop _merge
save, replace

**Generate annual rate of change as percentage of baseline LnNFL**************************************

capture drop LnNFLw3
gen LnNFLw3=ln(NFLw3)
capture drop LnNFLw1
gen LnNFLw1=ln(NFLw1)


capture drop bayes1LnNFLpct
gen bayes1LnNFLpct=bayes1LnNFL*100/LnNFLw1

save, replace


capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\IMPUTATION.smcl", replace



//////////////////////////////////////IMPUTED DATA//////////////////////

use HANDLS_paper51_NFLBRAINSCANFINALIZED, clear
sort HNDID
save, replace


merge HNDID using w3Age
save, replace

capture rename dxDiabetes w1dxDiabetes
capture rename dxHTN w1dxHTN
capture rename BMI w1BMI
capture rename WaistSize w1WaistSize
capture rename WHR w1WHR
capture rename hei2010_total_score w1hei2010_total_score
capture rename DASH_TotalScore w1DASH_TotalScore

save, replace

capture drop w1CES
capture rename CES w1CES

save, replace

**Impute data for wave 1 covariates. 

su w1edubr w1currdrugs w1smoke w1SRH w1BMI w1WaistSize w1WHR  w1dxDiabetes w1dxHTN  w1HypercholesterolemiarxCV w1cvdbr w1inflamcond w1rxNSAID w1hei2010_total_score w1mar100 w1DASH_TotalScore w1CES w1AlbGlob w1Albumin w1ALP w1ALT w1Amylase w1AST w1B12 w1BasoPct w1BiliDir w1BiliTot w1BUN w1Calcium w1Chol w1Cl w1CO2 w1Creatinine w1EosinPct w1ESR w1Fe w1Ferritin w1FeSat w1Folate w1GGT w1Globulin w1Glucose w1HBV w1Hct w1HCV w1HDL w1Hgb w1HgbA1C w1Insulin w1K w1LDH w1LDLcalc w1LymphPct w1MCH w1MCHC w1MCV w1Mg w1MonoPct w1Na w1NeutPct w1Phosphate w1Platelets w1RBC w1RDW w1RPR w1T3uptake w1T4Feee w1T4tot w1TIBC w1Triglyc w1TSH w1UpH w1UricAcid w1USpecGrav w1WBC w1HIV w1ChoHDLRat w1CRP w1Eosinophils w1Lipase w1Lymphocytes w1Monoocytes w1Neutrophils w1UCreatinine w1UMicroAlb w1VLDL w1TotalD LnNFLw1 LnNFLw3 bayes1LnNFL, det





save finaldata_imputed, replace


set matsize 11000

mi set flong

save, replace

mi unregister HNDID HNDwave w1edubr w1currdrugs w1smoke w1SRH w1BMI w1WaistSize w1WHR  w1dxDiabetes w1dxHTN  w1HypercholesterolemiarxCV w1cvdbr w1inflamcond w1rxNSAID w1hei2010_total_score w1mar100 w1DASH_TotalScore w1CES w1AlbGlob w1Albumin w1ALP w1ALT w1Amylase w1AST w1B12 w1BasoPct w1BiliDir w1BiliTot w1BUN w1Calcium w1Chol w1Cl w1CO2 w1Creatinine w1EosinPct w1ESR w1Fe w1Ferritin w1FeSat w1Folate w1GGT w1Globulin w1Glucose w1HBV w1Hct w1HCV w1HDL w1Hgb w1HgbA1C w1Insulin w1K w1LDH w1LDLcalc w1LymphPct w1MCH w1MCHC w1MCV w1Mg w1MonoPct w1Na w1NeutPct w1Phosphate w1Platelets w1RBC w1RDW w1RPR w1T3uptake w1T4Feee w1T4tot w1TIBC w1Triglyc w1TSH w1UpH w1UricAcid w1USpecGrav w1WBC w1HIV w1ChoHDLRat w1CRP w1Eosinophils w1Lipase w1Lymphocytes w1Monoocytes w1Neutrophils w1UCreatinine w1UMicroAlb w1VLDL w1TotalD LnNFLw1 LnNFLw3 bayes1LnNFL 


mi register imputed  w1edubr w1currdrugs w1smoke w1SRH w1BMI w1WaistSize w1WHR  w1dxDiabetes w1dxHTN  w1HypercholesterolemiarxCV w1cvdbr w1inflamcond w1rxNSAID w1hei2010_total_score w1mar100 w1DASH_TotalScore w1CES w1AlbGlob w1Albumin w1ALP w1ALT w1Amylase w1AST w1B12 w1BasoPct w1BiliDir w1BiliTot w1BUN w1Calcium w1Chol w1Cl w1CO2 w1Creatinine w1EosinPct w1ESR w1Fe w1Ferritin w1FeSat w1Folate w1GGT w1Globulin w1Glucose w1HBV w1Hct w1HCV w1HDL w1Hgb w1HgbA1C w1Insulin w1K w1LDH w1LDLcalc w1LymphPct w1MCH w1MCHC w1MCV w1Mg w1MonoPct w1Na w1NeutPct w1Phosphate w1Platelets w1RBC w1RDW w1RPR w1T3uptake w1T4Feee w1T4tot w1TIBC w1Triglyc w1TSH w1UpH w1UricAcid w1USpecGrav w1WBC w1HIV w1ChoHDLRat w1CRP w1Eosinophils w1Lipase w1Lymphocytes w1Monoocytes w1Neutrophils w1UCreatinine w1UMicroAlb w1VLDL w1TotalD


replace HNDwave=1 if HNDwave==.
save, replace


mi register passive HNDID HNDwave LnNFLw1 LnNFLw3 bayes1LnNFL 


mi impute chained (mlogit) w1edubr w1currdrugs w1smoke w1SRH w1dxDiabetes w1dxHTN  w1HypercholesterolemiarxCV w1cvdbr w1inflamcond w1rxNSAID (regress) w1BMI w1WaistSize w1WHR w1hei2010_total_score w1mar100 w1DASH_TotalScore w1CES w1AlbGlob w1Albumin w1ALP w1ALT w1Amylase w1AST w1B12 w1BasoPct w1BiliDir w1BiliTot w1BUN w1Calcium w1Chol w1Cl w1CO2 w1Creatinine w1EosinPct w1ESR w1Fe w1Ferritin w1FeSat w1Folate w1GGT w1Globulin w1Glucose w1HBV w1Hct w1HCV w1HDL w1Hgb w1HgbA1C w1Insulin w1K w1LDH w1LDLcalc w1LymphPct w1MCH w1MCHC w1MCV w1Mg w1MonoPct w1Na w1NeutPct w1Phosphate w1Platelets w1RBC w1RDW w1RPR w1T3uptake w1T4Feee w1T4tot w1TIBC w1Triglyc w1TSH w1UpH w1UricAcid w1USpecGrav w1WBC w1HIV w1ChoHDLRat w1CRP w1Eosinophils w1Lipase w1Lymphocytes w1Monoocytes w1Neutrophils w1UCreatinine w1UMicroAlb w1VLDL w1TotalD=Sex PovStat Race w1Age if w1Age~=., force augment noisily  add(5) rseed(1234) savetrace(tracefile, replace) 


save finaldata_imputed, replace


**Extract first imputation**

mi extract 1

save finaldata_imputed_one, replace





capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\LASSO.smcl", replace

/////////////////////////////////////LASSO//////////////////////////////

**LASSO COMMANDS FOR LnNFLw1**

use finaldata_imputed_one, clear


capture drop sample_LASSO
sort HNDID
save, replace
splitsample , generate(sample_LASSO) nsplit(2) rseed(1234)
tab sample_LASSO

save, replace



******************************************************LASSO MODELS*************************************************


lasso linear LnNFLw1 (w1Age Sex PovStat Race) i.w1edubr w1currdrugs w1smoke i.w1SRH i.w1dxDiabetes w1dxHTN  w1HypercholesterolemiarxCV w1cvdbr w1inflamcond w1rxNSAID w1BMI w1WaistSize w1WHR w1hei2010_total_score w1mar100 w1DASH_TotalScore w1CES w1AlbGlob w1Albumin w1ALP w1ALT w1Amylase w1AST w1B12 w1BasoPct w1BiliDir w1BiliTot w1BUN w1Calcium w1Chol w1Cl w1CO2 w1Creatinine w1EosinPct w1ESR w1Fe w1Ferritin w1FeSat w1Folate w1GGT w1Globulin w1Glucose w1HBV w1Hct w1HCV w1HDL w1Hgb w1HgbA1C w1Insulin w1K w1LDH w1LDLcalc w1LymphPct w1MCH w1MCHC w1MCV w1Mg w1MonoPct w1Na w1NeutPct w1Phosphate w1Platelets w1RBC w1RDW w1RPR w1T3uptake w1T4Feee w1T4tot w1TIBC w1Triglyc w1TSH w1UpH w1UricAcid w1USpecGrav w1WBC w1HIV w1ChoHDLRat w1CRP w1Eosinophils w1Lipase w1Lymphocytes w1Monoocytes w1Neutrophils w1UCreatinine w1UMicroAlb w1VLDL w1TotalD if sample_LASSO==1, rseed(1234)


cvplot 
estimates store cvLnNFLw1

lassoknots, display(nonzero osr2 bic)

lassoselect id=1

cvplot

estimates store minBICLnNFLw1

lasso linear LnNFLw1 (w1Age Sex PovStat Race)  i.w1edubr w1currdrugs w1smoke i.w1SRH i.w1dxDiabetes w1dxHTN  w1HypercholesterolemiarxCV w1cvdbr w1inflamcond w1rxNSAID w1BMI w1WaistSize w1WHR w1hei2010_total_score w1mar100 w1DASH_TotalScore w1CES w1AlbGlob w1Albumin w1ALP w1ALT w1Amylase w1AST w1B12 w1BasoPct w1BiliDir w1BiliTot w1BUN w1Calcium w1Chol w1Cl w1CO2 w1Creatinine w1EosinPct w1ESR w1Fe w1Ferritin w1FeSat w1Folate w1GGT w1Globulin w1Glucose w1HBV w1Hct w1HCV w1HDL w1Hgb w1HgbA1C w1Insulin w1K w1LDH w1LDLcalc w1LymphPct w1MCH w1MCHC w1MCV w1Mg w1MonoPct w1Na w1NeutPct w1Phosphate w1Platelets w1RBC w1RDW w1RPR w1T3uptake w1T4Feee w1T4tot w1TIBC w1Triglyc w1TSH w1UpH w1UricAcid w1USpecGrav w1WBC w1HIV w1ChoHDLRat w1CRP w1Eosinophils w1Lipase w1Lymphocytes w1Monoocytes w1Neutrophils w1UCreatinine w1UMicroAlb w1VLDL w1TotalD, selection(adaptive) rseed(1234)


estimates store adaptiveLnNFLw1

lassocoef cvLnNFLw1 minBICLnNFLw1 adaptiveLnNFLw1, sort(coef, standardized) nofvlabel

**----------------------------------------------------------
**               | cvLnNFLw1 minBICLnNFLw1  adaptiveLnNFLw1 
**---------------+------------------------------------------
**         w1Age |     x           x               x        
**          Race |     x           x               x        
**         w1BMI |     x                           x        
**   w1currdrugs |     x                           x        
**     w1Insulin |     x    
**   w1Phosphate |     x                           x        
**          w1Cl |     x    
**         w1CO2 |     x                           x        
**       PovStat |     x           x               x        
**   w1USpecGrav |     x                           x        
**  w1Creatinine |     x                           x        
**           Sex |     x           x               x        
**           w1K |     x    
**      w1Lipase |     x                           x        
**   w1Platelets |     x                           x        
**     3.w1edubr |     x                           x        
**               |
**         w1SRH |
**            3  |     x    
**               |
**   w1ChoHDLRat |     x                           x        
**        w1TIBC |     x                           x        
**         _cons |     x           x               x        
**     w1Albumin |                                 x        
**      w1HgbA1C |                                 x        
**               |
**         w1SRH |
**            1  |                                 x        
**               |
**      w1TotalD |                                 x        
**         w1HDL |                                 x        
**     w1Glucose |                                 x        
**     w1LDLcalc |                                 x        
**         w1ALT |                                 x        
**         w1LDH |                                 x        
**     w1BasoPct |                                 x        
**         w1BUN |                                 x        
** 0.w1dxDiabetes |                                 x        
**     w1Amylase |                                 x        
**         w1TSH |                                 x        
**       w1dxHTN |                                 x        
**     w1AlbGlob |                                 x        
**         w1CES |                                 x        
** w1Neutrophils |                                 x        
**         w1MCV |                                 x        
**   w1UMicroAlb |                                 x        
**----------------------------------------------------------
** Legend:
**  b - base level
**  e - empty cell
**  o - omitted
**  x - estimated





lassogof cvLnNFLw1 minBICLnNFLw1 adaptiveLnNFLw1, over(sample_LASSO) postselection


**Postselection coefficients
**-------------------------------------------------------------
**Name         sample_L~O |         MSE    R-squared        Obs
**------------------------+------------------------------------
**cvLnNFLw1               |
**                      1 |    .1359842       0.4197        348
**                      2 |     .186248       0.3296        346
**------------------------+------------------------------------
**minBICLnNFLw1           |
**                      1 |    .1750809       0.2529        348
**                      2 |    .2070778       0.2546        346
**------------------------+------------------------------------
**adaptiveLnNFLw1         |
**                      1 |    .1369333       0.4157        348
**                      2 |    .1313442       0.5272        346
**-------------------------------------------------------------



//////////////////RUN MODEL ON IMPUTED DATA: ADAPTIVE LASSO AS STARTING POINT///////////////////////////

use finaldata_imputed, clear


**Selected model: Adaptive**

mi estimate: reg LnNFLw1 w1Age Race w1BMI w1currdrugs w1Phosphate w1CO2 PovStat w1USpecGrav w1Creatinine Sex w1Lipase w1Platelets  i.w1edubr w1ChoHDLRat w1TIBC w1Albumin w1HgbA1C i.w1SRH w1TotalD w1HDL w1Glucose w1LDLcalc w1ALT w1LDH w1BasoPct w1BUN i.w1dxDiabetes w1Amylase w1TSH w1dxHTN w1AlbGlob w1CES w1Neutrophils w1MCV w1UMicroAlb

**Multiple-imputation estimates                   Imputations       =          5
**Linear regression                               Number of obs     =        694
**                                                Average RVI       =     0.0627
**                                                Largest FMI       =     0.5141
**                                                Complete DF       =        655
**DF adjustment:   Small sample                   DF:     min       =      17.67
**                                                        avg       =     459.07
**                                                        max       =     645.45
**Model F test:       Equal FMI                   F(  38,  642.1)   =      14.47
**Within VCE type:          OLS                   Prob > F          =     0.0000

**-------------------------------------------------------------------------------
**      LnNFLw1 | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
**--------------+----------------------------------------------------------------
**        w1Age |   .0206254   .0019414    10.62   0.000     .0168129    .0244379
**         Race |  -.2571507   .0418943    -6.14   0.000    -.3394497   -.1748517
**        w1BMI |  -.0149823   .0025984    -5.77   0.000    -.0200869   -.0098776
**  w1currdrugs |    .095642   .0443862     2.15   0.032     .0082717    .1830123
**  w1Phosphate |    .041439   .0300657     1.38   0.169    -.0176001    .1004782
**        w1CO2 |   .0085833   .0060529     1.42   0.157    -.0033036    .0204701
**      PovStat |   .0228562   .0350715     0.65   0.515    -.0460611    .0917735
**  w1USpecGrav |  -7.783302   2.421102    -3.21   0.001    -12.53784   -3.028768
** w1Creatinine |   .2193194   .0562754     3.90   0.000     .1060456    .3325932
**          Sex |   .0299499   .0389516     0.77   0.442    -.0465869    .1064867
**     w1Lipase |   .0007097   .0006674     1.06   0.289    -.0006051    .0020244
**  w1Platelets |  -.0004899   .0002472    -1.98   0.048    -.0009755   -4.32e-06
**              |
**      w1edubr |
**           2  |  -.0092962    .066896    -0.14   0.890    -.1407016    .1221091
**           3  |   .0440571   .0690753     0.64   0.524    -.0916295    .1797436
**              |
**  w1ChoHDLRat |  -.0114831   .0251946    -0.46   0.649    -.0612619    .0382957
**       w1TIBC |   -.000288   .0003153    -0.91   0.361    -.0009073    .0003313
**    w1Albumin |   -.253668   .0740989    -3.42   0.001    -.3991967   -.1081394
**     w1HgbA1C |   .0542668   .0311854     1.74   0.082    -.0069911    .1155246
**              |
**        w1SRH |
**           2  |  -.0745658   .0428724    -1.74   0.082    -.1587553    .0096237
**           3  |  -.1158557   .0456543    -2.54   0.011     -.205505   -.0262064
**              |
**     w1TotalD |   .0034205   .0016848     2.03   0.043      .000106     .006735
**        w1HDL |   .0017982   .0016258     1.11   0.270    -.0014055    .0050019
**    w1Glucose |   .0015101   .0010527     1.43   0.152    -.0005572    .0035775
**    w1LDLcalc |  -.0007784   .0006192    -1.26   0.210    -.0019989    .0004421
**        w1ALT |  -.0020316   .0009104    -2.23   0.026    -.0038196   -.0002435
**        w1LDH |   .0009866   .0005865     1.68   0.095    -.0001725    .0021457
**    w1BasoPct |   .1182973   .0545222     2.17   0.031     .0111628    .2254317
**        w1BUN |   .0065282   .0043948     1.49   0.138    -.0021159    .0151724
**              |
** w1dxDiabetes |
**pre-diabetes  |  -.0688421   .0424794    -1.62   0.106    -.1523137    .0146295
**    diabetes  |  -.1476909   .0700361    -2.11   0.035    -.2852196   -.0101621
**              |
**    w1Amylase |     .00118   .0008345     1.41   0.158    -.0004605    .0028205
**        w1TSH |    .006038   .0050152     1.20   0.229    -.0038133    .0158892
**      w1dxHTN |   .0489709   .0362335     1.35   0.178    -.0224181    .1203598
**    w1AlbGlob |  -.0799191   .0828021    -0.97   0.335    -.2425187    .0826806
**        w1CES |   .0014843   .0014136     1.05   0.294    -.0012918    .0042604
**w1Neutrophils |   9.20e-06   7.32e-06     1.26   0.210    -5.20e-06    .0000236
**        w1MCV |   .0024464    .002631     0.93   0.353    -.0027201    .0076128
**  w1UMicroAlb |   .0016303   .0015968     1.02   0.321     -.001729    .0049895
**        _cons |   9.549438   2.499104     3.82   0.000      4.64162    14.45725
**-------------------------------------------------------------------------------




**Reduced: Remove w1Phosphate w1CO2 w1Lipase i.w1edubr w1ChoHDLRat w1TIBC w1HgbA1C w1HDL w1Glucose w1LDLcalc w1LDH w1BUN w1Amylase w1TSH w1dxHTN  w1AlbGlob w1CES w1Neutrophils w1MCV w1UMicroAlb

 
mi estimate: reg LnNFLw1 w1Age Race w1BMI w1currdrugs PovStat w1USpecGrav w1Creatinine Sex  w1Platelets  w1Albumin  i.w1SRH w1TotalD  w1ALT w1BasoPct i.w1dxDiabetes 

**Multiple-imputation estimates                   Imputations       =          5
**Linear regression                               Number of obs     =        694
**                                                Average RVI       =     0.0503
**                                                Largest FMI       =     0.2731
**                                                Complete DF       =        676
**DF adjustment:   Small sample                   DF:     min       =      56.76
**                                                        avg       =     496.96
**                                                        max       =     669.61
**Model F test:       Equal FMI                   F(  17,  655.1)   =      27.07
**Within VCE type:          OLS                   Prob > F          =     0.0000
**
**-------------------------------------------------------------------------------
**      LnNFLw1 | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
**--------------+----------------------------------------------------------------
**        w1Age |   .0246135   .0017682    13.92   0.000     .0211413    .0280856
**         Race |  -.1929861   .0345389    -5.59   0.000    -.2608143   -.1251579
**        w1BMI |  -.0147849   .0023223    -6.37   0.000    -.0193451   -.0102247
** w1currdrugs |   .1144469   .0448129     2.55   0.011     .0260818    .2028121
**      PovStat |   .0121685   .0347124     0.35   0.726    -.0560017    .0803387
**  w1USpecGrav |  -6.332564   2.303314    -2.75   0.006    -10.85541   -1.809716
** w1Creatinine |   .3121058   .0466882     6.68   0.000     .2186055    .4056062
**          Sex |   -.009861   .0368254    -0.27   0.789     -.082238     .062516
**  w1Platelets |  -.0005506   .0002402    -2.29   0.022    -.0010222   -.0000789
**    w1Albumin |  -.2527239   .0598923    -4.22   0.000    -.3704225   -.1350254
**              |
**        w1SRH |
**           2  |  -.1037221   .0423314    -2.45   0.015    -.1868402   -.0206039
**           3  |   -.143625   .0438826    -3.27   0.001    -.2297927   -.0574574
**              |
**     w1TotalD |   .0036827   .0017536     2.10   0.037     .0002212    .0071441
**        w1ALT |  -.0010145   .0008714    -1.16   0.245    -.0027256    .0006967
**    w1BasoPct |   .0889644   .0557595     1.60   0.111    -.0206391    .1985678
**              |
* w1dxDiabetes |
**pre-diabetes  |  -.0524902   .0419486    -1.25   0.212    -.1350196    .0300393
**    diabetes  |   .0409333   .0532342     0.77   0.442    -.0636483    .1455148
**             |
**        _cons |    8.95895    2.38137     3.76   0.000      4.28246    13.63544
**-------------------------------------------------------------------------------


save, replace

***REMOVE w1ALT w1BasoPct and w1dxDiabetes**

 
mi estimate: reg LnNFLw1 w1Age Race w1BMI w1currdrugs PovStat w1USpecGrav w1Creatinine Sex  w1Platelets  w1Albumin  i.w1SRH w1TotalD   



**Multiple-imputation estimates                   Imputations       =          5
**Linear regression                               Number of obs     =        694
**                                                Average RVI       =     0.0490
**                                                Largest FMI       =     0.2576
**                                                Complete DF       =        680
**DF adjustment:   Small sample                   DF:     min       =      62.79
**                                                        avg       =     523.48
**                                                        max       =     671.77
**Model F test:       Equal FMI                   F(  13,  653.0)   =      34.69
**Within VCE type:          OLS                   Prob > F          =     0.0000

**------------------------------------------------------------------------------
**     LnNFLw1 | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
**-------------+----------------------------------------------------------------
**       w1Age |   .0247122   .0017217    14.35   0.000     .0213316    .0280928
**        Race |  -.1863256   .0342279    -5.44   0.000     -.253535   -.1191161
**       w1BMI |  -.0156944   .0022247    -7.05   0.000    -.0200627   -.0113261
** w1currdrugs |   .1111015   .0444687     2.50   0.013     .0235101    .1986929
**     PovStat |   .0134033   .0347012     0.39   0.699    -.0547456    .0815522
** w1USpecGrav |  -6.328351   2.287697    -2.77   0.006     -10.8203   -1.836404
** w1Creatinine |   .3141079   .0463553     6.78   0.000     .2214682    .4067475
**         Sex |  -.0222659   .0358619    -0.62   0.535    -.0927372    .0482053
** w1Platelets |  -.0004652    .000237    -1.96   0.050    -.0009305    1.29e-07
**   w1Albumin |  -.2674781   .0593553    -4.51   0.000    -.3840966   -.1508595
**             |
**       w1SRH |
**          2  |  -.1052666   .0416435    -2.53   0.012    -.1870337   -.0234995
**          3  |  -.1448167   .0426358    -3.40   0.001    -.2285343   -.0610991
**             |
**    w1TotalD |   .0037931   .0017544     2.16   0.032     .0003298    .0072564
**       _cons |   9.032853   2.358898     3.83   0.000     4.400942    13.66476
**------------------------------------------------------------------------------


**Remove w1Platelets**

mi estimate: reg LnNFLw1 w1Age Race w1BMI w1currdrugs PovStat w1USpecGrav w1Creatinine Sex   w1Albumin  i.w1SRH w1TotalD   

**Multiple-imputation estimates                   Imputations       =          5
**Linear regression                               Number of obs     =        694
**                                                Average RVI       =     0.0511
**                                                Largest FMI       =     0.2558
**                                                Complete DF       =        681
**DF adjustment:   Small sample                   DF:     min       =      63.57
**                                                        avg       =     515.58
**                                                        max       =     673.38
**Model F test:       Equal FMI                   F(  12,  649.1)   =      37.02
**Within VCE type:          OLS                   Prob > F          =     0.0000

**------------------------------------------------------------------------------
**     LnNFLw1 | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
**-------------+----------------------------------------------------------------
**       w1Age |   .0251809   .0017083    14.74   0.000     .0218265    .0285352
**        Race |  -.1890119   .0342617    -5.52   0.000    -.2562872   -.1217367
**       w1BMI |  -.0158338   .0022279    -7.11   0.000    -.0202083   -.0114593
** w1currdrugs |   .1110841   .0446752     2.49   0.014     .0230595    .1991086
**     PovStat |   .0166845   .0347086     0.48   0.631    -.0514771    .0848461
** w1USpecGrav |  -5.967447   2.285962    -2.61   0.009    -10.45599     -1.4789
**w1Creatinine |   .3180144   .0463473     6.86   0.000     .2254131    .4106156
**         Sex |  -.0068933   .0350987    -0.20   0.844    -.0758714    .0620848
**   w1Albumin |  -.2644972   .0594266    -4.45   0.000    -.3812516   -.1477427
**             |
**       w1SRH |
**          2  |   -.107733   .0417067    -2.58   0.010    -.1896239   -.0258421
**          3  |   -.143665   .0427294    -3.36   0.001    -.2275663   -.0597637
**             |
**    w1TotalD |   .0037885   .0017478     2.17   0.031     .0003411     .007236
**       _cons |   8.487892   2.347961     3.62   0.000     3.877453    13.09833
**------------------------------------------------------------------------------


******************************************************LASSO MODELS FOR WAVE 3 NFL*************************************************


**LASSO COMMANDS FOR bayes1LnNFLpct**

use finaldata_imputed_one, clear


capture drop sample_LASSO
sort HNDID
save, replace
splitsample , generate(sample_LASSO) nsplit(2) rseed(1234)
tab sample_LASSO

save, replace



******************************************************LASSO MODELS*************************************************


lasso linear LnNFLw3 (w1Age Sex PovStat Race) i.w1edubr w1currdrugs w1smoke i.w1SRH i.w1dxDiabetes w1dxHTN  w1HypercholesterolemiarxCV w1cvdbr w1inflamcond w1rxNSAID w1BMI w1WaistSize w1WHR w1hei2010_total_score w1mar100 w1DASH_TotalScore w1CES w1AlbGlob w1Albumin w1ALP w1ALT w1Amylase w1AST w1B12 w1BasoPct w1BiliDir w1BiliTot w1BUN w1Calcium w1Chol w1Cl w1CO2 w1Creatinine w1EosinPct w1ESR w1Fe w1Ferritin w1FeSat w1Folate w1GGT w1Globulin w1Glucose w1HBV w1Hct w1HCV w1HDL w1Hgb w1HgbA1C w1Insulin w1K w1LDH w1LDLcalc w1LymphPct w1MCH w1MCHC w1MCV w1Mg w1MonoPct w1Na w1NeutPct w1Phosphate w1Platelets w1RBC w1RDW w1RPR w1T3uptake w1T4Feee w1T4tot w1TIBC w1Triglyc w1TSH w1UpH w1UricAcid w1USpecGrav w1WBC w1HIV w1ChoHDLRat w1CRP w1Eosinophils w1Lipase w1Lymphocytes w1Monoocytes w1Neutrophils w1UCreatinine w1UMicroAlb w1VLDL w1TotalD if sample_LASSO==1, rseed(1234)


cvplot 
estimates store cvLnNFLw3

lassoknots, display(nonzero osr2 bic)

lassoselect id=1

cvplot

estimates store minBICLnNFLw3

lasso linear LnNFLw3 (w1Age Sex PovStat Race)  i.w1edubr w1currdrugs w1smoke i.w1SRH i.w1dxDiabetes w1dxHTN  w1HypercholesterolemiarxCV w1cvdbr w1inflamcond w1rxNSAID w1BMI w1WaistSize w1WHR w1hei2010_total_score w1mar100 w1DASH_TotalScore w1CES w1AlbGlob w1Albumin w1ALP w1ALT w1Amylase w1AST w1B12 w1BasoPct w1BiliDir w1BiliTot w1BUN w1Calcium w1Chol w1Cl w1CO2 w1Creatinine w1EosinPct w1ESR w1Fe w1Ferritin w1FeSat w1Folate w1GGT w1Globulin w1Glucose w1HBV w1Hct w1HCV w1HDL w1Hgb w1HgbA1C w1Insulin w1K w1LDH w1LDLcalc w1LymphPct w1MCH w1MCHC w1MCV w1Mg w1MonoPct w1Na w1NeutPct w1Phosphate w1Platelets w1RBC w1RDW w1RPR w1T3uptake w1T4Feee w1T4tot w1TIBC w1Triglyc w1TSH w1UpH w1UricAcid w1USpecGrav w1WBC w1HIV w1ChoHDLRat w1CRP w1Eosinophils w1Lipase w1Lymphocytes w1Monoocytes w1Neutrophils w1UCreatinine w1UMicroAlb w1VLDL w1TotalD, selection(adaptive) rseed(1234)


estimates store adaptiveLnNFLw3

lassocoef cvLnNFLw3 minBICLnNFLw3 adaptiveLnNFLw3, sort(coef, standardized) nofvlabel



**----------------------------------------------------------
**               | cvLnNFLw3 minBICLnNFLw3  adaptiveLnNFLw3 
**---------------+------------------------------------------
**         w1Age |     x           x               x        
**          Race |     x           x               x        
**         w1GGT |     x                           x        
**         w1BMI |     x                           x        
**  w1Creatinine |     x                           x        
**          w1Cl |     x    
**     w1Insulin |     x                           x        
**      w1Lipase |     x                           x        
**     w1Amylase |     x                           x        
**      w1HgbA1C |     x                           x        
**       PovStat |     x           x               x        
**           Sex |     x           x               x        
**               |
**         w1SRH |
**            3  |     x    
**               |
**     w1rxNSAID |     x    
**         w1ALP |     x                           x        
**         w1WHR |     x    
**         w1WBC |     x    
**     w1NeutPct |     x    
**   w1USpecGrav |     x                           x        
**         w1BUN |     x                           x        
**       w1smoke |     x    
**         w1B12 |     x    
**               |
**         w1SRH |
**            1  |     x                           x        
**               |
**         w1RPR |     x    
**         w1CRP |     x    
**         w1HIV |     x    
**       w1dxHTN |     x    
**         _cons |     x           x               x        
**     w1Glucose |                                 x        
**    w1UricAcid |                                 x        
**    w1EosinPct |                                 x        
**         w1HDL |                                 x        
**0.w1dxDiabetes |                                 x        
**   w1currdrugs |                                 x        
**         w1CES |                                 x        
**         w1ESR |                                 x        
**         w1LDH |                                 x        
**----------------------------------------------------------
**Legend:
**  b - base level
**  e - empty cell
**  o - omitted
**  x - estimated


lassogof cvLnNFLw3 minBICLnNFLw3 adaptiveLnNFLw3, over(sample_LASSO) postselection




**Postselection coefficients
**-------------------------------------------------------------
**Name         sample_L~O |         MSE    R-squared        Obs
**------------------------+------------------------------------
**cvLnNFLw3               |
**                      1 |    .1551697       0.4401        354
**                      2 |    .4142864      -0.1388        355
**------------------------+------------------------------------
**minBICLnNFLw3           |
**                      1 |    .2081272       0.2490        354
**                      2 |    .2836157       0.2204        355
**------------------------+------------------------------------
**adaptiveLnNFLw3         |
**                      1 |    .1756309       0.3663        354
**                      2 |    .2022459       0.4440        355
**-------------------------------------------------------------





//////////////RUN MODEL ON IMPUTED DATA//////////////////

***ADAPTIVE LASSO AS STARTING POINT****

use finaldata_imputed, clear


mi estimate: reg LnNFLw3 w1Age Race w1GGT w1BMI w1Creatinine w1Insulin w1Lipase w1Amylase w1HgbA1C PovStat Sex w1ALP w1USpecGrav w1BUN i.w1SRH w1Glucose w1UricAcid w1EosinPct w1HDL i.w1dxDiabetes w1currdrugs w1CES w1ESR w1LDH

**Multiple-imputation estimates                   Imputations       =          5
**Linear regression                               Number of obs     =        709
**                                                Average RVI       =     0.0837
**                                                Largest FMI       =     0.3864
**                                                Complete DF       =        682
**DF adjustment:   Small sample                   DF:     min       =      30.52
**                                                        avg       =     400.26
**                                                        max       =     672.80
**Model F test:       Equal FMI                   F(  26,  650.7)   =      16.64
**Within VCE type:          OLS                   Prob > F          =     0.0000
**
**-------------------------------------------------------------------------------
**      LnNFLw3 | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
**--------------+----------------------------------------------------------------
**        w1Age |   .0230277    .002152    10.70   0.000     .0188006    .0272548
**         Race |  -.1486843   .0411694    -3.61   0.000    -.2295257   -.0678429
**        w1GGT |    .000457   .0003765     1.21   0.226    -.0002852    .0011993
**        w1BMI |  -.0101147    .002831    -3.57   0.000     -.015676   -.0045534
** w1Creatinine |    .224075   .0625206     3.58   0.001     .0964817    .3516682
**    w1Insulin |   -.003059   .0016202    -1.89   0.059    -.0062402    .0001222
**     w1Lipase |   .0011063   .0007795     1.42   0.157    -.0004302    .0026429
**    w1Amylase |   .0014437   .0009809     1.47   0.142    -.0004884    .0033759
**     w1HgbA1C |    .034876   .0336058     1.04   0.300    -.0311154    .1008673
**      PovStat |   .0476348    .038738     1.23   0.219    -.0284309    .1237004
**          Sex |   -.011214   .0481885    -0.23   0.816    -.1060721    .0836442
**        w1ALP |    .002384   .0008367     2.85   0.005     .0007332    .0040347
**  w1USpecGrav |   -9.88734   2.851523    -3.47   0.001    -15.49335   -4.281326
**        w1BUN |   .0121905    .005151     2.37   0.019     .0020141    .0223669
**              |
**        w1SRH |
**           2  |  -.0445514   .0495764    -0.90   0.369    -.1419333    .0528305
**           3  |  -.0696356   .0536921    -1.30   0.196    -.1752781    .0360068
**              |
**    w1Glucose |   .0035985   .0011794     3.05   0.002     .0012824    .0059145
**   w1UricAcid |   .0407004   .0152045     2.68   0.008     .0107068     .070694
**   w1EosinPct |   .0274093   .0093342     2.94   0.003     .0090677    .0457508
**        w1HDL |   .0013418   .0013506     0.99   0.322    -.0013247    .0040083
**              |
** w1dxDiabetes |
**pre-diabetes  |  -.0803385   .0487673    -1.65   0.100    -.1761682    .0154913
**    diabetes  |  -.2222999    .077587    -2.87   0.004    -.3747355   -.0698643
**              |
**  w1currdrugs |   .0780181    .049201     1.59   0.113    -.0186405    .1746768
**        w1CES |    .002077   .0016708     1.24   0.214    -.0012037    .0053578
**        w1ESR |   .0014405   .0014504     0.99   0.321    -.0014102    .0042912
**        w1LDH |   .0006464   .0006804     0.95   0.345    -.0007063    .0019991
**        _cons |   9.993284   2.884871     3.46   0.001     4.322987    15.66358
**-------------------------------------------------------------------------------



**keep  w1Age Sex Race PovStat w1BMI w1ALP w1UMicroAlb w1HgbA1C w1USpecGrav w1UricAcid w1Creatinine w1EosinPct w1currdrugs**


mi estimate: reg LnNFLw3 w1Age Sex Race PovStat w1BMI w1Creatinine  w1ALP w1USpecGrav w1BUN w1Glucose w1UricAcid w1EosinPct i.w1dxDiabetes

**Multiple-imputation estimates                   Imputations       =          5
**Linear regression                               Number of obs     =        709
**                                                Average RVI       =     0.1021
**                                                Largest FMI       =     0.4265
**                                                Complete DF       =        694
**DF adjustment:   Small sample                   DF:     min       =      25.40
**                                                        avg       =     421.08
**                                                        max       =     688.84
**Model F test:       Equal FMI                   F(  14,  609.4)   =      27.33
**Within VCE type:          OLS                   Prob > F          =     0.0000
**
**-------------------------------------------------------------------------------
**      LnNFLw3 | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
**--------------+----------------------------------------------------------------
**        w1Age |   .0234049    .002118    11.05   0.000     .0192454    .0275644
**          Sex |  -.0204517   .0432522    -0.47   0.637    -.1055021    .0645987
**         Race |  -.0987301   .0361961    -2.73   0.007    -.1697981   -.0276621
**      PovStat |    .060181   .0383666     1.57   0.117    -.0151496    .1355116
**        w1BMI |  -.0124843   .0026263    -4.75   0.000    -.0176426   -.0073259
** w1Creatinine |   .2528599   .0619328     4.08   0.000     .1254091    .3803107
**        w1ALP |   .0029093   .0008022     3.63   0.000     .0013264    .0044923
**  w1USpecGrav |  -11.02134   2.842592    -3.88   0.000    -16.60622   -5.436465
**        w1BUN |   .0129297   .0052072     2.48   0.015     .0026101    .0232494
**    w1Glucose |   .0042431   .0007989     5.31   0.000     .0026729    .0058133
**   w1UricAcid |   .0421852   .0150355     2.81   0.006     .0124621    .0719083
**   w1EosinPct |   .0263639   .0093716     2.81   0.005     .0079557    .0447722
**              |
** w1dxDiabetes |
**pre-diabetes  |  -.0954143   .0492707    -1.94   0.054    -.1923514    .0015227
**    diabetes  |   -.200183   .0757791    -2.64   0.008    -.3490059   -.0513601
**              |
*        _cons |   11.46363   2.868967     4.00   0.000     5.827406    17.09985
**-------------------------------------------------------------------------------



////////////////////////////COMMON SET OF COVARIATES/////////////////////////////

**LnNFLw3: w1Age Sex Race PovStat w1BMI w1Creatinine  w1ALP w1USpecGrav w1BUN w1Glucose w1UricAcid w1EosinPct i.w1dxDiabetes
**LnNFLw1:  w1Age Race w1BMI w1currdrugs PovStat w1USpecGrav w1Creatinine Sex   w1Albumin  i.w1SRH w1TotalD

**Additional covariates to basic covariates: 
**w1BMI w1Creatinine  w1ALP w1USpecGrav w1BUN w1Glucose w1UricAcid w1EosinPct i.w1dxDiabetes  w1currdrugs  w1USpecGrav  w1Albumin  i.w1SRH w1TotalD


**Model 2: Adiposity: BMI**
**w1BMI

**Model 3: BMI+ Diabetes/blood glucose **
**w1BMI+w1dxDiabetes+w1Glucose 

**Modle 4: BMI+ liver/kidney disease **
**w1BMI+w1Creatinine+w1USpecGrav+w1BUN+w1ALP+w1UricAcid

**Model 5: BMI+ inflammation**
**w1BMI+TotalD w1Albumin w1EosinPct


**Model 6: BMI+ lifestyle/health-related factors**
**w1BMI+w1currdrugs+w1SRH




********************************************MAIN ANALYSIS***************************************************

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLE1PRELIM.smcl", replace





//TABLE 1: DESCRIPTIVES//


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

tab sample_final


capture drop w1Agebr
gen w1Agebr=.
replace w1Agebr=0 if w1Age<=50
replace w1Agebr=1 if w1Age>50 & w1Age~=.

save, replace

capture drop deltaLnNFL
gen deltaLnNFL=(LnNFLw3-LnNFLw1)/(w3Age-w1Age)
su deltaLnNFL
histogram deltaLnNFL

save, replace

**Correlation between exposures**

corr bayes1LnNFLpct bayes1LnNFL LnNFLw1 LnNFLw3 deltaLnNFL if sample_final==1
pwcorr bayes1LnNFLpct bayes1LnNFL LnNFLw1 LnNFLw3 deltaLnNFL if sample_final==1, sig

graph matrix  LnNFLw1 LnNFLw3 deltaLnNFL bayes1LnNFL if sample_final==1

graph save scatterplotNFL.gph, replace


save, replace


tab1 Sex w1Agebr Race PovStat if sample_final==1

**Selection bias**

**Total population**

ttest w1Age, by(sample_final)
tab Sex sample_final, row col chi
tab Race sample_final, row col chi
tab PovStat sample_final, row col chi

///CREATE BINARY TRACKING NFL above or below 8 pg/mL variable VARIABLE////

capture drop NFLw1w3trackhigh
gen NFLw1w3trackhigh=.
replace NFLw1w3trackhigh=1 if NFLw1>8 & NFLw3>8 & NFLw1~=. & NFLw3~=.
replace NFLw1w3trackhigh=0 if NFLw1w3trackhigh~=1 & NFLw1~=. & NFLw3~=.


tab NFLw1w3trackhigh

reg Left_Hippocampus NFLw1w3trackhigh w1Age Sex Race PovStat ICV_volM2 if sample_final==1
reg Right_Hippocampus NFLw1w3trackhigh w1Age Sex Race PovStat ICV_volM2 if sample_final==1
reg LnLesion_Volume NFLw1w3trackhigh w1Age Sex Race PovStat ICV_volM2 if sample_final==1


reg Left_Hippocampus NFLw1w3trackhigh w1Age Sex Race PovStat ICV_volM2 if sample_final==1 & Sex==1
reg Right_Hippocampus NFLw1w3trackhigh w1Age Sex Race PovStat ICV_volM2 if sample_final==1 & Sex==1
reg LnLesion_Volume NFLw1w3trackhigh w1Age Sex Race PovStat ICV_volM2 if sample_final==1 & Sex==1


reg Left_Hippocampus NFLw1w3trackhigh w1Age Sex Race PovStat ICV_volM2 if sample_final==1 & Sex==2
reg Right_Hippocampus NFLw1w3trackhigh w1Age Sex Race PovStat ICV_volM2 if sample_final==1 & Sex==2
reg LnLesion_Volume NFLw1w3trackhigh w1Age Sex Race PovStat ICV_volM2 if sample_final==1 & Sex==2



capture drop NFLw1w3tracklow
gen NFLw1w3tracklow=.
replace NFLw1w3tracklow=1 if NFLw1<=8 & NFLw3<=8 & NFLw1~=. & NFLw3~=.
replace NFLw1w3tracklow=0 if NFLw1w3tracklow~=1 & NFLw1~=. & NFLw3~=.


tab NFLw1w3tracklow

reg Left_Hippocampus NFLw1w3tracklow w1Age Sex Race PovStat ICV_volM2 if sample_final==1
reg Right_Hippocampus NFLw1w3tracklow w1Age Sex Race PovStat ICV_volM2 if sample_final==1
reg LnLesion_Volume NFLw1w3tracklow w1Age Sex Race PovStat ICV_volM2 if sample_final==1


reg Left_Hippocampus NFLw1w3tracklow w1Age Sex Race PovStat ICV_volM2 if sample_final==1 & Sex==1
reg Right_Hippocampus NFLw1w3tracklow w1Age Sex Race PovStat ICV_volM2 if sample_final==1  & Sex==1
reg LnLesion_Volume NFLw1w3tracklow w1Age Sex Race PovStat ICV_volM2 if sample_final==1  & Sex==1


reg Left_Hippocampus NFLw1w3tracklow w1Age Sex Race PovStat ICV_volM2 if sample_final==1 & Sex==2
reg Right_Hippocampus NFLw1w3tracklow w1Age Sex Race PovStat ICV_volM2 if sample_final==1  & Sex==2
reg LnLesion_Volume NFLw1w3tracklow w1Age Sex Race PovStat ICV_volM2 if sample_final==1  & Sex==2


save,replace


capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLE1.smcl", replace




**********************************************TABLE 1: VARIABLES BY SEX******************************


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear
sort HNDID
capture drop _merge
save, replace

**Overall**

tab Sex if sample_final==1
su w1Age if sample_final==1
tab w1Agebr if sample_final==1
tab Race if sample_final==1
tab PovStat if sample_final==1

su TIME_V1SCAN if sample_final==1, det
su TIME_V2SCAN if sample_final==1, det

save, replace

****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear


save, replace



****w1BMI w1dxDiabetes w1Glucose w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid w1Albumin w1EosinPct w1TotalD w1currdrugs w1SRH

mi estimate: mean w1BMI if sample_final==1
mi estimate: prop w1dxDiabetes if sample_final==1
mi estimate: mean w1Glucose if sample_final==1
mi estimate: mean w1Creatinine if sample_final==1
mi estimate: mean w1USpecGrav if sample_final==1
mi estimate: mean w1BUN if sample_final==1
mi estimate: mean w1ALP if sample_final==1
mi estimate: mean w1UricAcid if sample_final==1
mi estimate: mean w1Albumin if sample_final==1
mi estimate: mean w1EosinPct if sample_final==1
mi estimate: mean w1TotalD if sample_final==1
mi estimate: prop w1currdrugs if sample_final==1
mi estimate: prop w1SRH if sample_final==1

save, replace

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear



su LnNFLw1 LnNFLw3 bayes1LnNFL deltaLnNFL if sample_final==1, det


su TOTALBRAIN if sample_final==1
su GM if sample_final==1
su WM if sample_final==1


su FRONTAL_GM_L_volM2  if sample_final==1
su FRONTAL_WM_L_volM2  if sample_final==1
su TEMPORAL_GM_L_volM2  if sample_final==1
su TEMPORAL_WM_L_volM2  if sample_final==1
su PARIETAL_GM_L_volM2  if sample_final==1
su PARIETAL_WM_L_volM2  if sample_final==1
su OCCIPITAL_GM_L_volM2  if sample_final==1
su OCCIPITAL_WM_L_volM2  if sample_final==1

su FRONTAL_GM_R_volM2  if sample_final==1
su FRONTAL_WM_R_volM2  if sample_final==1
su TEMPORAL_GM_R_volM2  if sample_final==1
su TEMPORAL_WM_R_volM2  if sample_final==1
su PARIETAL_GM_R_volM2  if sample_final==1
su PARIETAL_WM_R_volM2  if sample_final==1
su OCCIPITAL_GM_R_volM2  if sample_final==1
su OCCIPITAL_WM_R_volM2  if sample_final==1

su Left_Hippocampus if sample_final==1
su Right_Hippocampus if sample_final==1


su LnLesion_Volume if sample_final==1


su Left_Hippocampuspct if sample_final==1
su Right_Hippocampuspct if sample_final==1


su LnLesion_Volumepct if sample_final==1

su ICV_volM2 if sample_final==1


save, replace

****************************T-tests by Sex**********************

**by Sex group**
ttest w1Age if sample_final==1, by(Sex)
tab Race Sex if sample_final==1, row col chi
tab PovStat Sex if sample_final==1, row col chi
ttest TIME_V1SCAN if sample_final==1, by(Sex)
ttest TIME_V2SCAN if sample_final==1, by(Sex)




***Follow-up time by other important covariates**
ttest TIME_V1SCAN if sample_final==1, by(w1Agebr)
ttest TIME_V1SCAN if sample_final==1, by(PovStat)
ttest TIME_V1SCAN if sample_final==1, by(Race)

ttest TIME_V2SCAN if sample_final==1, by(w1Agebr)
ttest TIME_V2SCAN if sample_final==1, by(PovStat)
ttest TIME_V2SCAN if sample_final==1, by(Race)

****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear
save, replace

*****w1Age at w1, categorical**
capture drop w1Agebr
mi passive: gen w1Agebr=.
mi passive: replace w1Agebr=0 if w1Age<=50
mi passive: replace w1Agebr=1 if w1Age>50 & w1Age~=.

save, replace

mi estimate: reg w1BMI Sex if sample_final==1
mi estimate: mlogit w1dxDiabetes Sex if sample_final==1,baseoutcome(0)
mi estimate: reg w1Glucose Sex if sample_final==1
mi estimate: reg w1Creatinine Sex if sample_final==1
mi estimate: reg w1USpecGrav Sex if sample_final==1
mi estimate: reg w1BUN Sex if sample_final==1
mi estimate: reg w1ALP Sex if sample_final==1
mi estimate: reg w1UricAcid Sex if sample_final==1
mi estimate: reg w1Albumin Sex if sample_final==1
mi estimate: reg w1EosinPct Sex if sample_final==1
mi estimate: reg w1TotalD Sex if sample_final==1
mi estimate: reg w1currdrugs Sex if sample_final==1
mi estimate: mlogit w1SRH Sex if sample_final==1, baseoutcome(1)



mi estimate: reg w1BMI Sex w1Age Race PovStat if sample_final==1
mi estimate: mlogit w1dxDiabetes Sex w1Age Race PovStat if sample_final==1,baseoutcome(0)
mi estimate: reg w1Glucose Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1Creatinine Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1USpecGrav Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1BUN Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1ALP Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1UricAcid Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1Albumin Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1EosinPct Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1TotalD Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1currdrugs Sex w1Age Race PovStat if sample_final==1
mi estimate: mlogit w1SRH Sex w1Age Race PovStat if sample_final==1, baseoutcome(1)

**Further adjusted for ICV_volM2**


mi estimate: reg w1BMI Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: mlogit w1dxDiabetes Sex w1Age Race PovStat ICV_volM2 if sample_final==1,baseoutcome(0)
mi estimate: reg w1Glucose Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1Creatinine Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1USpecGrav Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1BUN Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1ALP Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1UricAcid Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1Albumin Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1EosinPct Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1TotalD Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1currdrugs Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: mlogit w1SRH Sex w1Age Race PovStat ICV_volM2 if sample_final==1, baseoutcome(1)


****w1BMI w1dxDiabetes w1Glucose w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid w1Albumin w1EosinPct w1TotalD w1currdrugs w1SRH

**Females**

mi estimate: mean w1BMI if sample_final==1 & Sex==1
mi estimate: prop w1dxDiabetes if sample_final==1 & Sex==1
mi estimate: mean w1Glucose if sample_final==1 & Sex==1
mi estimate: mean w1Creatinine if sample_final==1 & Sex==1
mi estimate: mean w1USpecGrav if sample_final==1 & Sex==1
mi estimate: mean w1BUN if sample_final==1 & Sex==1
mi estimate: mean w1ALP if sample_final==1 & Sex==1
mi estimate: mean w1UricAcid if sample_final==1 & Sex==1
mi estimate: mean w1Albumin if sample_final==1  & Sex==1
mi estimate: mean w1EosinPct if sample_final==1 & Sex==1
mi estimate: mean w1TotalD if sample_final==1 & Sex==1
mi estimate: prop w1currdrugs if sample_final==1 & Sex==1
mi estimate: prop w1SRH if sample_final==1 & Sex==1



**Males**

mi estimate: mean w1BMI if sample_final==1 & Sex==2
mi estimate: prop w1dxDiabetes if sample_final==1 & Sex==2
mi estimate: mean w1Glucose if sample_final==1 & Sex==2
mi estimate: mean w1Creatinine if sample_final==1 & Sex==2
mi estimate: mean w1USpecGrav if sample_final==1 & Sex==2
mi estimate: mean w1BUN if sample_final==1 & Sex==2
mi estimate: mean w1ALP if sample_final==1 & Sex==2
mi estimate: mean w1UricAcid if sample_final==1 & Sex==2
mi estimate: mean w1Albumin if sample_final==1  & Sex==2
mi estimate: mean w1EosinPct if sample_final==1 & Sex==2
mi estimate: mean w1TotalD if sample_final==1 & Sex==2
mi estimate: prop w1currdrugs if sample_final==1 & Sex==2
mi estimate: prop w1SRH if sample_final==1 & Sex==2




use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


ttest LnNFLw1 if sample_final==1, by(Sex)
ttest LnNFLw3 if sample_final==1, by(Sex)
ttest bayes1LnNFL if sample_final==1, by(Sex)
ttest deltaLnNFL if sample_final==1, by(Sex)

tab NFLw1w3trackhigh Sex if sample_final==1, row col chi
tab NFLw1w3tracklow Sex if sample_final==1, row col chi

mlogit NFLw1w3trackhigh Sex w1Age Race PovStat ICV_volM2 if sample_final==1, baseoutcome(0)
mlogit NFLw1w3tracklow Sex w1Age Race PovStat ICV_volM2 if sample_final==1, baseoutcome(0)


**Median/IQR by Sex**

su LnNFLw1 if sample_final==1 & Sex==1, det
su LnNFLw1 if sample_final==1 & Sex==2, det

su LnNFLw3 if sample_final==1 & Sex==1, det
su LnNFLw3 if sample_final==1 & Sex==2, det
 

su bayes1LnNFL if sample_final==1 & Sex==1, det
su bayes1LnNFL if sample_final==1 & Sex==2, det
 
 
su deltaLnNFL if sample_final==1 & Sex==1, det
su deltaLnNFL if sample_final==1 & Sex==2, det

tab1 NFLw1w3trackhigh NFLw1w3tracklow if sample_final==1 
tab1 NFLw1w3trackhigh NFLw1w3tracklow if sample_final==1 & Sex==1
tab1 NFLw1w3trackhigh NFLw1w3tracklow if sample_final==1 & Sex==2


ttest TOTALBRAIN  if sample_final==1, by(Sex)
ttest GM  if sample_final==1, by(Sex)
ttest WM  if sample_final==1, by(Sex)


ttest FRONTAL_GM_L_volM2   if sample_final==1, by(Sex)
ttest FRONTAL_WM_L_volM2   if sample_final==1, by(Sex)
ttest TEMPORAL_GM_L_volM2   if sample_final==1, by(Sex)
ttest TEMPORAL_WM_L_volM2   if sample_final==1, by(Sex)
ttest PARIETAL_GM_L_volM2   if sample_final==1, by(Sex)
ttest PARIETAL_WM_L_volM2   if sample_final==1, by(Sex)
ttest OCCIPITAL_GM_L_volM2   if sample_final==1, by(Sex)
ttest OCCIPITAL_WM_L_volM2   if sample_final==1, by(Sex)
ttest FRONTAL_GM_R_volM2   if sample_final==1, by(Sex)
ttest FRONTAL_WM_R_volM2   if sample_final==1, by(Sex)
ttest TEMPORAL_GM_R_volM2   if sample_final==1, by(Sex)
ttest TEMPORAL_WM_R_volM2   if sample_final==1, by(Sex)
ttest PARIETAL_GM_R_volM2   if sample_final==1, by(Sex)
ttest PARIETAL_WM_R_volM2   if sample_final==1, by(Sex)
ttest OCCIPITAL_GM_R_volM2   if sample_final==1, by(Sex)
ttest OCCIPITAL_WM_R_volM2   if sample_final==1, by(Sex)

ttest Left_Hippocampus  if sample_final==1, by(Sex)
ttest Right_Hippocampus  if sample_final==1, by(Sex)

ttest LnLesion_Volume  if sample_final==1, by(Sex)


ttest Left_Hippocampuspct  if sample_final==1, by(Sex)
ttest Right_Hippocampuspct  if sample_final==1, by(Sex)

ttest LnLesion_Volumepct  if sample_final==1, by(Sex)

ttest ICV_volM2  if sample_final==1, by(Sex)


save, replace

**************Adjusted model for other variables: Race, w1Age and PovStat***********************

reg TIME_V1SCAN Sex w1Age Race PovStat  if sample_final==1
reg w1Age Sex Race PovStat if sample_final==1
mlogit Race w1Age Sex PovStat if sample_final==1, baseoutcome(1)
mlogit PovStat w1Age Sex Race if sample_final==1, baseoutcome(1)

reg LnNFLw1 Sex w1Age Race PovStat if sample_final==1
reg LnNFLw3 Sex w1Age Race PovStat if sample_final==1
reg bayes1LnNFL Sex w1Age Race PovStat if sample_final==1


reg TOTALBRAIN Sex w1Age Race PovStat if sample_final==1
reg GM Sex w1Age Race PovStat if sample_final==1
reg WM Sex w1Age Race PovStat if sample_final==1


reg FRONTAL_GM_L_volM2  Sex w1Age Race PovStat if sample_final==1
reg FRONTAL_WM_L_volM2  Sex w1Age Race PovStat if sample_final==1
reg TEMPORAL_GM_L_volM2  Sex w1Age Race PovStat if sample_final==1
reg TEMPORAL_WM_L_volM2  Sex w1Age Race PovStat if sample_final==1
reg PARIETAL_GM_L_volM2  Sex w1Age Race PovStat if sample_final==1
reg PARIETAL_WM_L_volM2  Sex w1Age Race PovStat if sample_final==1
reg OCCIPITAL_GM_L_volM2  Sex w1Age Race PovStat if sample_final==1
reg OCCIPITAL_WM_L_volM2  Sex w1Age Race PovStat if sample_final==1


reg FRONTAL_GM_R_volM2  Sex w1Age Race PovStat if sample_final==1
reg FRONTAL_WM_R_volM2  Sex w1Age Race PovStat if sample_final==1
reg TEMPORAL_GM_R_volM2  Sex w1Age Race PovStat if sample_final==1
reg TEMPORAL_WM_R_volM2  Sex w1Age Race PovStat if sample_final==1
reg PARIETAL_GM_R_volM2  Sex w1Age Race PovStat if sample_final==1
reg PARIETAL_WM_R_volM2  Sex w1Age Race PovStat if sample_final==1
reg OCCIPITAL_GM_R_volM2  Sex w1Age Race PovStat if sample_final==1
reg OCCIPITAL_WM_R_volM2  Sex w1Age Race PovStat if sample_final==1


reg Left_Hippocampus Sex w1Age Race PovStat if sample_final==1
reg Right_Hippocampus Sex w1Age Race PovStat if sample_final==1

reg LnLesion_Volume Sex w1Age Race PovStat if sample_final==1



**************Adjusted model for other variables: Race, w1Age, PovStat and ICV***********************
reg TIME_V1SCAN Sex w1Age Race PovStat ICV_volM2  if sample_final==1
reg w1Age Sex Race PovStat ICV_volM2 if sample_final==1
mlogit Race w1Age Sex PovStat ICV_volM2 if sample_final==1, baseoutcome(1)
mlogit PovStat w1Age Sex Race ICV_volM2 if sample_final==1, baseoutcome(1)

reg LnNFLw1 Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg LnNFLw3 Sex w1Age Race PovStat ICV_volM2 if sample_final==1

reg TOTALBRAIN Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg GM Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg WM Sex w1Age Race PovStat ICV_volM2 if sample_final==1


reg FRONTAL_GM_L_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg FRONTAL_WM_L_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg TEMPORAL_GM_L_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg TEMPORAL_WM_L_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg PARIETAL_GM_L_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg PARIETAL_WM_L_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg OCCIPITAL_GM_L_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg OCCIPITAL_WM_L_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1


reg FRONTAL_GM_R_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg FRONTAL_WM_R_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg TEMPORAL_GM_R_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg TEMPORAL_WM_R_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg PARIETAL_GM_R_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg PARIETAL_WM_R_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg OCCIPITAL_GM_R_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg OCCIPITAL_WM_R_volM2  Sex w1Age Race PovStat ICV_volM2 if sample_final==1


reg Left_Hippocampus Sex w1Age Race PovStat ICV_volM2 if sample_final==1
reg Right_Hippocampus Sex w1Age Race PovStat ICV_volM2 if sample_final==1

reg LnLesion_Volume Sex w1Age Race PovStat ICV_volM2 if sample_final==1



save, replace


reg Left_Hippocampuspct Sex w1Age Race PovStat   if sample_final==1
reg Right_Hippocampuspct Sex w1Age Race PovStat if sample_final==1

reg LnLesion_Volumepct Sex w1Age Race PovStat if sample_final==1



save, replace




*******************************TABLE S1. LIST OF ROIs for WML ANALYSIS, ANALYSIS C'*********
capture log close


*******************************TABLE S2. VARIABLES BY TERTILE OF LNNFL AT VISITS 1 AND 2***************



capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLES2.smcl", replace



***********************TABLE S2: BY LnNFL visits 1 and 2 TERTILES********************************* 




/////////////////////////////////BY LnNFL at visit 1, tertiles/////////////////////


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear
sort HNDID
capture drop _merge
save, replace

tab sample_final



capture drop LnNFLw1tert

xtile LnNFLw1tert=LnNFLw1 if sample_final==1, nq(3)
bysort LnNFLw1tert: su LnNFLw1
bysort LnNFLw1tert: su LnNFLw3
bysort LnNFLw1tert: su bayes1LnNFL
bysort LnNFLw1tert: su deltaLnNFL


save, replace


**************LnNFL, v1, T1******************************

tab Sex if sample_final==1 & LnNFLw1tert==1
su w1Age if sample_final==1 & LnNFLw1tert==1
tab w1Agebr if sample_final==1 & LnNFLw1tert==1
tab Race if sample_final==1 & LnNFLw1tert==1
tab PovStat if sample_final==1 & LnNFLw1tert==1

su TIME_V1SCAN if sample_final==1 & LnNFLw1tert==1
su TIME_V2SCAN if sample_final==1 & LnNFLw1tert==1



****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear


capture drop LnNFLw1tert

xtile LnNFLw1tert=LnNFLw1 if sample_final==1, nq(3)
bysort LnNFLw1tert: su LnNFLw1

save, replace


****w1BMI w1dxDiabetes w1Glucose w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid w1Albumin w1EosinPct w1TotalD w1currdrugs w1SRH

mi estimate: mean w1BMI if sample_final==1 & LnNFLw1tert==1
mi estimate: prop w1dxDiabetes if sample_final==1 & LnNFLw1tert==1
mi estimate: mean w1Glucose if sample_final==1 & LnNFLw1tert==1
mi estimate: mean w1Creatinine if sample_final==1 & LnNFLw1tert==1
mi estimate: mean w1USpecGrav if sample_final==1 & LnNFLw1tert==1
mi estimate: mean w1BUN if sample_final==1 & LnNFLw1tert==1
mi estimate: mean w1ALP if sample_final==1 & LnNFLw1tert==1
mi estimate: mean w1UricAcid if sample_final==1 & LnNFLw1tert==1
mi estimate: mean w1Albumin if sample_final==1 & LnNFLw1tert==1
mi estimate: mean w1EosinPct if sample_final==1 & LnNFLw1tert==1
mi estimate: mean w1TotalD if sample_final==1 & LnNFLw1tert==1
mi estimate: prop w1currdrugs if sample_final==1 & LnNFLw1tert==1
mi estimate: prop w1SRH if sample_final==1 & LnNFLw1tert==1




use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear




su LnNFLw1 LnNFLw3 if sample_final==1 & LnNFLw1tert==1, det
su LnNFLw3 if sample_final==1 & LnNFLw1tert==1, det
su bayes1LnNFL if sample_final==1 & LnNFLw1tert==1, det
su deltaLnNFL if sample_final==1 & LnNFLw1tert==1, det

su ICV_volM2 if sample_final==1 & LnNFLw1tert==1

su TOTALBRAIN if sample_final==1 & LnNFLw1tert==1
su GM if sample_final==1 & LnNFLw1tert==1
su WM if sample_final==1 & LnNFLw1tert==1


su FRONTAL_GM_L_volM2  if sample_final==1 & LnNFLw1tert==1 
su FRONTAL_WM_L_volM2  if sample_final==1 & LnNFLw1tert==1
su TEMPORAL_GM_L_volM2  if sample_final==1 & LnNFLw1tert==1
su TEMPORAL_WM_L_volM2  if sample_final==1 & LnNFLw1tert==1
su PARIETAL_GM_L_volM2  if sample_final==1 & LnNFLw1tert==1
su PARIETAL_WM_L_volM2  if sample_final==1 & LnNFLw1tert==1
su OCCIPITAL_GM_L_volM2  if sample_final==1 & LnNFLw1tert==1
su OCCIPITAL_WM_L_volM2  if sample_final==1 & LnNFLw1tert==1


su FRONTAL_GM_R_volM2  if sample_final==1 & LnNFLw1tert==1
su FRONTAL_WM_R_volM2  if sample_final==1 & LnNFLw1tert==1
su TEMPORAL_GM_R_volM2  if sample_final==1 & LnNFLw1tert==1
su TEMPORAL_WM_R_volM2  if sample_final==1 & LnNFLw1tert==1
su PARIETAL_GM_R_volM2  if sample_final==1 & LnNFLw1tert==1
su PARIETAL_WM_R_volM2  if sample_final==1 & LnNFLw1tert==1
su OCCIPITAL_GM_R_volM2 if sample_final==1 & LnNFLw1tert==1
su OCCIPITAL_WM_R_volM2  if sample_final==1 & LnNFLw1tert==1


su Left_Hippocampus if sample_final==1 & LnNFLw1tert==1
su Right_Hippocampus if sample_final==1 & LnNFLw1tert==1

su LnLesion_Volume if sample_final==1 & LnNFLw1tert==1

su Left_Hippocampuspct if sample_final==1 & LnNFLw1tert==1
su Right_Hippocampuspct if sample_final==1 & LnNFLw1tert==1

su LnLesion_Volumepct if sample_final==1 & LnNFLw1tert==1

*************************LnNFL, V1, second tertile****************

**************LnNFL, v1, T2******************************
use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

tab Sex if sample_final==1 & LnNFLw1tert==2
su w1Age if sample_final==1 & LnNFLw1tert==2
tab w1Agebr if sample_final==1 & LnNFLw1tert==2
tab Race if sample_final==1 & LnNFLw1tert==2
tab PovStat if sample_final==1 & LnNFLw1tert==2

su TIME_V1SCAN if sample_final==1 & LnNFLw1tert==2
su TIME_V2SCAN if sample_final==1 & LnNFLw1tert==2



****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear


capture drop LnNFLw1tert

xtile LnNFLw1tert=LnNFLw1 if sample_final==1, nq(3)
bysort LnNFLw1tert: su LnNFLw1

save, replace


****w1BMI w1dxDiabetes w1Glucose w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid w1Albumin w1EosinPct w1TotalD w1currdrugs w1SRH

mi estimate: mean w1BMI if sample_final==1 & LnNFLw1tert==2
mi estimate: prop w1dxDiabetes if sample_final==1 & LnNFLw1tert==2
mi estimate: mean w1Glucose if sample_final==1 & LnNFLw1tert==2
mi estimate: mean w1Creatinine if sample_final==1 & LnNFLw1tert==2
mi estimate: mean w1USpecGrav if sample_final==1 & LnNFLw1tert==2
mi estimate: mean w1BUN if sample_final==1 & LnNFLw1tert==2
mi estimate: mean w1ALP if sample_final==1 & LnNFLw1tert==2
mi estimate: mean w1UricAcid if sample_final==1 & LnNFLw1tert==2
mi estimate: mean w1Albumin if sample_final==1 & LnNFLw1tert==2
mi estimate: mean w1EosinPct if sample_final==1 & LnNFLw1tert==2
mi estimate: mean w1TotalD if sample_final==1 & LnNFLw1tert==2
mi estimate: prop w1currdrugs if sample_final==1 & LnNFLw1tert==2
mi estimate: prop w1SRH if sample_final==1 & LnNFLw1tert==2




use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear




su LnNFLw1 LnNFLw3 if sample_final==1 & LnNFLw1tert==2, det
su LnNFLw3 if sample_final==1 & LnNFLw1tert==2, det
su bayes1LnNFL if sample_final==1 & LnNFLw1tert==2, det
su deltaLnNFL if sample_final==1 & LnNFLw1tert==2, det

su ICV_volM2 if sample_final==1 & LnNFLw1tert==2

su TOTALBRAIN if sample_final==1 & LnNFLw1tert==2
su GM if sample_final==1 & LnNFLw1tert==2
su WM if sample_final==1 & LnNFLw1tert==2


su FRONTAL_GM_L_volM2 if sample_final==1 & LnNFLw1tert==2 
su FRONTAL_WM_L_volM2 if sample_final==1 & LnNFLw1tert==2
su TEMPORAL_GM_L_volM2 if sample_final==1 & LnNFLw1tert==2
su TEMPORAL_WM_L_volM2 if sample_final==1 & LnNFLw1tert==2
su PARIETAL_GM_L_volM2 if sample_final==1 & LnNFLw1tert==2
su PARIETAL_WM_L_volM2 if sample_final==1 & LnNFLw1tert==2
su OCCIPITAL_GM_L_volM2 if sample_final==1 & LnNFLw1tert==2
su OCCIPITAL_WM_L_volM2 if sample_final==1 & LnNFLw1tert==2


su FRONTAL_GM_R_volM2 if sample_final==1 & LnNFLw1tert==2
su FRONTAL_WM_R_volM2 if sample_final==1 & LnNFLw1tert==2
su TEMPORAL_GM_R_volM2 if sample_final==1 & LnNFLw1tert==2
su TEMPORAL_WM_R_volM2 if sample_final==1 & LnNFLw1tert==2
su PARIETAL_GM_R_volM2 if sample_final==1 & LnNFLw1tert==2
su PARIETAL_WM_R_volM2 if sample_final==1 & LnNFLw1tert==2
su OCCIPITAL_GM_R_volM2 if sample_final==1 & LnNFLw1tert==2
su OCCIPITAL_WM_R if sample_final==1 & LnNFLw1tert==2


su Left_Hippocampus if sample_final==1 & LnNFLw1tert==2
su Right_Hippocampus if sample_final==1 & LnNFLw1tert==2

su LnLesion_Volume if sample_final==1 & LnNFLw1tert==2

su Left_Hippocampuspct if sample_final==1 & LnNFLw1tert==2
su Right_Hippocampuspct if sample_final==1 & LnNFLw1tert==2

su LnLesion_Volumepct if sample_final==1 & LnNFLw1tert==2


**************LnNFL, v1, T3******************************
use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

tab Sex if sample_final==1 & LnNFLw1tert==3
su w1Age if sample_final==1 & LnNFLw1tert==3
tab w1Agebr if sample_final==1 & LnNFLw1tert==3
tab Race if sample_final==1 & LnNFLw1tert==3
tab PovStat if sample_final==1 & LnNFLw1tert==3

su TIME_V1SCAN if sample_final==1 & LnNFLw1tert==3
su TIME_V2SCAN if sample_final==1 & LnNFLw1tert==3



****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear


capture drop LnNFLw1tert

xtile LnNFLw1tert=LnNFLw1 if sample_final==1, nq(3)

save, replace


****w1BMI w1dxDiabetes w1Glucose w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid w1Albumin w1EosinPct w1TotalD w1currdrugs w1SRH

mi estimate: mean w1BMI if sample_final==1 & LnNFLw1tert==3
mi estimate: prop w1dxDiabetes if sample_final==1 & LnNFLw1tert==3
mi estimate: mean w1Glucose if sample_final==1 & LnNFLw1tert==3
mi estimate: mean w1Creatinine if sample_final==1 & LnNFLw1tert==3
mi estimate: mean w1USpecGrav if sample_final==1 & LnNFLw1tert==3
mi estimate: mean w1BUN if sample_final==1 & LnNFLw1tert==3
mi estimate: mean w1ALP if sample_final==1 & LnNFLw1tert==3
mi estimate: mean w1UricAcid if sample_final==1 & LnNFLw1tert==3
mi estimate: mean w1Albumin if sample_final==1 & LnNFLw1tert==3
mi estimate: mean w1EosinPct if sample_final==1 & LnNFLw1tert==3
mi estimate: mean w1TotalD if sample_final==1 & LnNFLw1tert==3
mi estimate: prop w1currdrugs if sample_final==1 & LnNFLw1tert==3
mi estimate: prop w1SRH if sample_final==1 & LnNFLw1tert==3




use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear




su LnNFLw1 LnNFLw3 if sample_final==1 & LnNFLw1tert==3, det
su LnNFLw3 if sample_final==1 & LnNFLw1tert==3, det
su bayes1LnNFL if sample_final==1 & LnNFLw1tert==3, det
su deltaLnNFL if sample_final==1 & LnNFLw1tert==3, det


su ICV_volM2 if sample_final==1 & LnNFLw1tert==3


su TOTALBRAIN if sample_final==1 & LnNFLw1tert==3
su GM if sample_final==1 & LnNFLw1tert==3
su WM if sample_final==1 & LnNFLw1tert==3


su FRONTAL_GM_L_volM2 if sample_final==1 & LnNFLw1tert==3 
su FRONTAL_WM_L_volM2 if sample_final==1 & LnNFLw1tert==3
su TEMPORAL_GM_L_volM2 if sample_final==1 & LnNFLw1tert==3
su TEMPORAL_WM_L_volM2 if sample_final==1 & LnNFLw1tert==3
su PARIETAL_GM_L_volM2 if sample_final==1 & LnNFLw1tert==3
su PARIETAL_WM_L_volM2 if sample_final==1 & LnNFLw1tert==3
su OCCIPITAL_GM_L_volM2 if sample_final==1 & LnNFLw1tert==3
su OCCIPITAL_WM_L_volM2 if sample_final==1 & LnNFLw1tert==3


su FRONTAL_GM_R_volM2 if sample_final==1 & LnNFLw1tert==3
su FRONTAL_WM_R if sample_final==1 & LnNFLw1tert==3
su TEMPORAL_GM_R_volM2 if sample_final==1 & LnNFLw1tert==3
su TEMPORAL_WM_R if sample_final==1 & LnNFLw1tert==3
su PARIETAL_GM_R_volM2 if sample_final==1 & LnNFLw1tert==3
su PARIETAL_WM_R if sample_final==1 & LnNFLw1tert==3
su OCCIPITAL_GM_R_volM2 if sample_final==1 & LnNFLw1tert==3
su OCCIPITAL_WM_R if sample_final==1 & LnNFLw1tert==3


su Left_Hippocampus if sample_final==1 & LnNFLw1tert==3
su Right_Hippocampus if sample_final==1 & LnNFLw1tert==3

su LnLesion_Volume if sample_final==1 & LnNFLw1tert==3

su Left_Hippocampuspct if sample_final==1 & LnNFLw1tert==3
su Right_Hippocampuspct if sample_final==1 & LnNFLw1tert==3

su LnLesion_Volumepct if sample_final==1 & LnNFLw1tert==3

****************BY LnNFL tertile, V1*************************
use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

**Overall**

tab Sex LnNFLw1tert if sample_final==1, row col chi 
reg w1Age LnNFLw1tert  if sample_final==1
tab w1Agebr LnNFLw1tert if sample_final==1, row col chi
tab Race LnNFLw1tert if sample_final==1, row col chi 
tab PovStat LnNFLw1tert if sample_final==1, row col chi 

reg TIME_V1SCAN LnNFLw1tert  if sample_final==1
reg TIME_V2SCAN LnNFLw1tert  if sample_final==1


****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear


mi estimate: reg w1BMI LnNFLw1tert if sample_final==1
mi estimate: mlogit w1dxDiabetes LnNFLw1tert if sample_final==1,baseoutcome(0)
mi estimate: reg w1Glucose LnNFLw1tert if sample_final==1
mi estimate: reg w1Creatinine LnNFLw1tert if sample_final==1
mi estimate: reg w1USpecGrav LnNFLw1tert if sample_final==1
mi estimate: reg w1BUN LnNFLw1tert if sample_final==1
mi estimate: reg w1ALP LnNFLw1tert if sample_final==1
mi estimate: reg w1UricAcid LnNFLw1tert if sample_final==1
mi estimate: reg w1Albumin LnNFLw1tert if sample_final==1
mi estimate: reg w1EosinPct LnNFLw1tert if sample_final==1
mi estimate: reg w1TotalD LnNFLw1tert if sample_final==1
mi estimate: reg w1currdrugs LnNFLw1tert if sample_final==1
mi estimate: mlogit w1SRH LnNFLw1tert if sample_final==1, baseoutcome(1)



mi estimate: reg w1BMI LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: mlogit w1dxDiabetes LnNFLw1tert Sex w1Age Race PovStat if sample_final==1,baseoutcome(0)
mi estimate: reg w1Glucose LnNFLw1tert Sex  w1Age Race PovStat if sample_final==1
mi estimate: reg w1Creatinine LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1USpecGrav LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1BUN LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1ALP LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1UricAcid LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1Albumin LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1EosinPct LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1TotalD LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1currdrugs LnNFLw1tert Sex w1Age Race PovStat if sample_final==1
mi estimate: mlogit w1SRH LnNFLw1tert Sex w1Age Race PovStat if sample_final==1, baseoutcome(1)

**Further adjusted for ICV_volM2**


mi estimate: reg w1BMI LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: mlogit w1dxDiabetes LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1,baseoutcome(0)
mi estimate: reg w1Glucose LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1Creatinine LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1USpecGrav LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1BUN LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1ALP LnNFLw1tert  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1UricAcid LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1Albumin LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1EosinPct LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1TotalD LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1currdrugs LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: mlogit w1SRH LnNFLw1tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1, baseoutcome(1)




use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

reg LnNFLw1 LnNFLw1tert  if sample_final==1

reg LnNFLw3 LnNFLw1tert if sample_final==1

reg bayes1LnNFL LnNFLw1tert  if sample_final==1
 
reg deltaLnNFL LnNFLw1tert  if sample_final==1

save, replace

reg ICV_volM2 LnNFLw1tert  if sample_final==1


reg TOTALBRAIN LnNFLw1tert  if sample_final==1
reg GM LnNFLw1tert if sample_final==1
reg WM LnNFLw1tert  if sample_final==1


reg FRONTAL_GM_L_volM2 LnNFLw1tert  if sample_final==1 
reg FRONTAL_WM_L_volM2 LnNFLw1tert  if sample_final==1
reg TEMPORAL_GM_L_volM2 LnNFLw1tert  if sample_final==1
reg TEMPORAL_WM_L_volM2 LnNFLw1tert  if sample_final==1
reg PARIETAL_GM_L_volM2 LnNFLw1tert  if sample_final==1
reg PARIETAL_WM_L_volM2 LnNFLw1tert  if sample_final==1
reg OCCIPITAL_GM_L_volM2 LnNFLw1tert  if sample_final==1
reg OCCIPITAL_WM_L_volM2 LnNFLw1tert  if sample_final==1


reg FRONTAL_GM_R_volM2 LnNFLw1tert  if sample_final==1
reg FRONTAL_WM_R_volM2 LnNFLw1tert  if sample_final==1
reg TEMPORAL_GM_R_volM2 LnNFLw1tert  if sample_final==1
reg TEMPORAL_WM_R_volM2 LnNFLw1tert  if sample_final==1
reg PARIETAL_GM_R_volM2 LnNFLw1tert   if sample_final==1
reg PARIETAL_WM_R_volM2 LnNFLw1tert   if sample_final==1
reg OCCIPITAL_GM_R_volM2 LnNFLw1tert  if sample_final==1
reg OCCIPITAL_WM_R_volM2 LnNFLw1tert  if sample_final==1


reg Left_Hippocampus LnNFLw1tert   if sample_final==1
reg Right_Hippocampus LnNFLw1tert  if sample_final==1

reg Left_Hippocampuspct LnNFLw1tert   if sample_final==1
reg Right_Hippocampuspct LnNFLw1tert  if sample_final==1


reg LnLesion_Volume LnNFLw1tert  if sample_final==1
reg LnLesion_Volumepct LnNFLw1tert  if sample_final==1


tab NFLw1w3trackhigh LnNFLw1tert if sample_final==1, row col chi
tab NFLw1w3tracklow LnNFLw1tert if sample_final==1, row col chi

****

reg LnNFLw1 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1

reg LnNFLw3 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1

reg bayes1LnNFL LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
 
reg deltaLnNFL LnNFLw1tert w1Age Sex Race PovStat if sample_final==1

save, replace

reg ICV_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1

reg TOTALBRAIN LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg GM LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg WM LnNFLw1tert w1Age Sex Race PovStat if sample_final==1


reg FRONTAL_GM_L_volM2 LnNFLw1tert w1Age Sex Race PovStat  if sample_final==1 
reg FRONTAL_WM_L_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg TEMPORAL_GM_L_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg TEMPORAL_WM_L_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg PARIETAL_GM_L_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg PARIETAL_WM_L_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg OCCIPITAL_GM_L_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg OCCIPITAL_WM_L_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1


reg FRONTAL_GM_R_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg FRONTAL_WM_R_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg TEMPORAL_GM_R_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg TEMPORAL_WM_R_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg PARIETAL_GM_R_volM2 LnNFLw1tert w1Age Sex Race PovStat  if sample_final==1
reg PARIETAL_WM_R_volM2 LnNFLw1tert w1Age Sex Race PovStat  if sample_final==1
reg OCCIPITAL_GM_R_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg OCCIPITAL_WM_R_volM2 LnNFLw1tert w1Age Sex Race PovStat if sample_final==1


reg Left_Hippocampus LnNFLw1tert w1Age Sex Race PovStat  if sample_final==1
reg Right_Hippocampus LnNFLw1tert w1Age Sex Race PovStat if sample_final==1

reg Left_Hippocampuspct LnNFLw1tert  w1Age Sex Race PovStat if sample_final==1
reg Right_Hippocampuspct LnNFLw1tert w1Age Sex Race PovStat if sample_final==1


reg LnLesion_Volume LnNFLw1tert w1Age Sex Race PovStat if sample_final==1
reg LnLesion_Volumepct LnNFLw1tert w1Age Sex Race PovStat if sample_final==1



mlogit NFLw1w3trackhigh LnNFLw1tert w1Age Sex Race PovStat  if sample_final==1, baseoutcome(0)
mlogit NFLw1w3tracklow LnNFLw1tert w1Age Sex Race PovStat  if sample_final==1, baseoutcome(0)


mlogit NFLw1w3trackhigh LnNFLw1tert w1Age Sex Race PovStat ICV_volM2 if sample_final==1, baseoutcome(0)
mlogit NFLw1w3tracklow LnNFLw1tert w1Age Sex Race PovStat ICV_volM2 if sample_final==1, baseoutcome(0)


save, replace




/////////////////////////////////BY LnNFL at visit 2, tertiles/////////////////////


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear
sort HNDID
capture drop _merge
save, replace

tab sample_final



capture drop LnNFLw3tert

xtile LnNFLw3tert=LnNFLw3 if sample_final==1, nq(3)
bysort LnNFLw3tert: su LnNFLw1
bysort LnNFLw3tert: su LnNFLw3
bysort LnNFLw3tert: su bayes1LnNFL
bysort LnNFLw3tert: su deltaLnNFL


save, replace


**************LnNFL, v2, T1******************************

tab Sex if sample_final==1 & LnNFLw3tert==1
su w1Age if sample_final==1 & LnNFLw3tert==1
tab w1Agebr if sample_final==1 & LnNFLw3tert==1
tab Race if sample_final==1 & LnNFLw3tert==1
tab PovStat if sample_final==1 & LnNFLw3tert==1

su TIME_V1SCAN if sample_final==1 & LnNFLw3tert==1
su TIME_V2SCAN if sample_final==1 & LnNFLw3tert==1



****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear


capture drop LnNFLw3tert

xtile LnNFLw3tert=LnNFLw3 if sample_final==1, nq(3)
bysort LnNFLw3tert: su LnNFLw1

save, replace


****w1BMI w1dxDiabetes w1Glucose w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid w1Albumin w1EosinPct w1TotalD w1currdrugs w1SRH

mi estimate: mean w1BMI if sample_final==1 & LnNFLw3tert==1
mi estimate: prop w1dxDiabetes if sample_final==1 & LnNFLw3tert==1
mi estimate: mean w1Glucose if sample_final==1 & LnNFLw3tert==1
mi estimate: mean w1Creatinine if sample_final==1 & LnNFLw3tert==1
mi estimate: mean w1USpecGrav if sample_final==1 & LnNFLw3tert==1
mi estimate: mean w1BUN if sample_final==1 & LnNFLw3tert==1
mi estimate: mean w1ALP if sample_final==1 & LnNFLw3tert==1
mi estimate: mean w1UricAcid if sample_final==1 & LnNFLw3tert==1
mi estimate: mean w1Albumin if sample_final==1 & LnNFLw3tert==1
mi estimate: mean w1EosinPct if sample_final==1 & LnNFLw3tert==1
mi estimate: mean w1TotalD if sample_final==1 & LnNFLw3tert==1
mi estimate: prop w1currdrugs if sample_final==1 & LnNFLw3tert==1
mi estimate: prop w1SRH if sample_final==1 & LnNFLw3tert==1




use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear




su LnNFLw1 LnNFLw3 if sample_final==1 & LnNFLw3tert==1, det
su LnNFLw3 if sample_final==1 & LnNFLw3tert==1, det
su bayes1LnNFL if sample_final==1 & LnNFLw3tert==1, det
su deltaLnNFL if sample_final==1 & LnNFLw3tert==1, det

su ICV_volM2 if sample_final==1 & LnNFLw3tert==1

su TOTALBRAIN if sample_final==1 & LnNFLw3tert==1
su GM if sample_final==1 & LnNFLw3tert==1
su WM if sample_final==1 & LnNFLw3tert==1


su FRONTAL_GM_L_volM2  if sample_final==1 & LnNFLw3tert==1 
su FRONTAL_WM_L_volM2  if sample_final==1 & LnNFLw3tert==1
su TEMPORAL_GM_L_volM2  if sample_final==1 & LnNFLw3tert==1
su TEMPORAL_WM_L_volM2  if sample_final==1 & LnNFLw3tert==1
su PARIETAL_GM_L_volM2  if sample_final==1 & LnNFLw3tert==1
su PARIETAL_WM_L_volM2  if sample_final==1 & LnNFLw3tert==1
su OCCIPITAL_GM_L_volM2  if sample_final==1 & LnNFLw3tert==1
su OCCIPITAL_WM_L_volM2  if sample_final==1 & LnNFLw3tert==1


su FRONTAL_GM_R_volM2  if sample_final==1 & LnNFLw3tert==1
su FRONTAL_WM_R_volM2  if sample_final==1 & LnNFLw3tert==1
su TEMPORAL_GM_R_volM2  if sample_final==1 & LnNFLw3tert==1
su TEMPORAL_WM_R_volM2  if sample_final==1 & LnNFLw3tert==1
su PARIETAL_GM_R_volM2  if sample_final==1 & LnNFLw3tert==1
su PARIETAL_WM_R_volM2  if sample_final==1 & LnNFLw3tert==1
su OCCIPITAL_GM_R_volM2 if sample_final==1 & LnNFLw3tert==1
su OCCIPITAL_WM_R_volM2  if sample_final==1 & LnNFLw3tert==1


su Left_Hippocampus if sample_final==1 & LnNFLw3tert==1
su Right_Hippocampus if sample_final==1 & LnNFLw3tert==1

su LnLesion_Volume if sample_final==1 & LnNFLw3tert==1

su Left_Hippocampuspct if sample_final==1 & LnNFLw3tert==1
su Right_Hippocampuspct if sample_final==1 & LnNFLw3tert==1

su LnLesion_Volumepct if sample_final==1 & LnNFLw3tert==1

*************************LnNFL, V2, second tertile****************

**************LnNFL, v2, T2******************************
use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

tab Sex if sample_final==1 & LnNFLw3tert==2
su w1Age if sample_final==1 & LnNFLw3tert==2
tab w1Agebr if sample_final==1 & LnNFLw3tert==2
tab Race if sample_final==1 & LnNFLw3tert==2
tab PovStat if sample_final==1 & LnNFLw3tert==2

su TIME_V1SCAN if sample_final==1 & LnNFLw3tert==2
su TIME_V2SCAN if sample_final==1 & LnNFLw3tert==2



****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear


capture drop LnNFLw3tert

xtile LnNFLw3tert=LnNFLw3 if sample_final==1, nq(3)
bysort LnNFLw3tert: su LnNFLw1

save, replace


****w1BMI w1dxDiabetes w1Glucose w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid w1Albumin w1EosinPct w1TotalD w1currdrugs w1SRH

mi estimate: mean w1BMI if sample_final==1 & LnNFLw3tert==2
mi estimate: prop w1dxDiabetes if sample_final==1 & LnNFLw3tert==2
mi estimate: mean w1Glucose if sample_final==1 & LnNFLw3tert==2
mi estimate: mean w1Creatinine if sample_final==1 & LnNFLw3tert==2
mi estimate: mean w1USpecGrav if sample_final==1 & LnNFLw3tert==2
mi estimate: mean w1BUN if sample_final==1 & LnNFLw3tert==2
mi estimate: mean w1ALP if sample_final==1 & LnNFLw3tert==2
mi estimate: mean w1UricAcid if sample_final==1 & LnNFLw3tert==2
mi estimate: mean w1Albumin if sample_final==1 & LnNFLw3tert==2
mi estimate: mean w1EosinPct if sample_final==1 & LnNFLw3tert==2
mi estimate: mean w1TotalD if sample_final==1 & LnNFLw3tert==2
mi estimate: prop w1currdrugs if sample_final==1 & LnNFLw3tert==2
mi estimate: prop w1SRH if sample_final==1 & LnNFLw3tert==2




use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear




su LnNFLw1 LnNFLw3 if sample_final==1 & LnNFLw3tert==2, det
su LnNFLw3 if sample_final==1 & LnNFLw3tert==2, det
su bayes1LnNFL if sample_final==1 & LnNFLw3tert==2, det
su deltaLnNFL if sample_final==1 & LnNFLw3tert==2, det

su ICV_volM2 if sample_final==1 & LnNFLw3tert==2

su TOTALBRAIN if sample_final==1 & LnNFLw3tert==2
su GM if sample_final==1 & LnNFLw3tert==2
su WM if sample_final==1 & LnNFLw3tert==2


su FRONTAL_GM_L_volM2 if sample_final==1 & LnNFLw3tert==2 
su FRONTAL_WM_L_volM2 if sample_final==1 & LnNFLw3tert==2
su TEMPORAL_GM_L_volM2 if sample_final==1 & LnNFLw3tert==2
su TEMPORAL_WM_L_volM2 if sample_final==1 & LnNFLw3tert==2
su PARIETAL_GM_L_volM2 if sample_final==1 & LnNFLw3tert==2
su PARIETAL_WM_L_volM2 if sample_final==1 & LnNFLw3tert==2
su OCCIPITAL_GM_L_volM2 if sample_final==1 & LnNFLw3tert==2
su OCCIPITAL_WM_L_volM2 if sample_final==1 & LnNFLw3tert==2


su FRONTAL_GM_R_volM2 if sample_final==1 & LnNFLw3tert==2
su FRONTAL_WM_R_volM2 if sample_final==1 & LnNFLw3tert==2
su TEMPORAL_GM_R_volM2 if sample_final==1 & LnNFLw3tert==2
su TEMPORAL_WM_R_volM2 if sample_final==1 & LnNFLw3tert==2
su PARIETAL_GM_R_volM2 if sample_final==1 & LnNFLw3tert==2
su PARIETAL_WM_R_volM2 if sample_final==1 & LnNFLw3tert==2
su OCCIPITAL_GM_R_volM2 if sample_final==1 & LnNFLw3tert==2
su OCCIPITAL_WM_R if sample_final==1 & LnNFLw3tert==2


su Left_Hippocampus if sample_final==1 & LnNFLw3tert==2
su Right_Hippocampus if sample_final==1 & LnNFLw3tert==2

su LnLesion_Volume if sample_final==1 & LnNFLw3tert==2

su Left_Hippocampuspct if sample_final==1 & LnNFLw3tert==2
su Right_Hippocampuspct if sample_final==1 & LnNFLw3tert==2

su LnLesion_Volumepct if sample_final==1 & LnNFLw3tert==2


**************LnNFL, v2, T3******************************
use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

tab Sex if sample_final==1 & LnNFLw3tert==3
su w1Age if sample_final==1 & LnNFLw3tert==3
tab w1Agebr if sample_final==1 & LnNFLw3tert==3
tab Race if sample_final==1 & LnNFLw3tert==3
tab PovStat if sample_final==1 & LnNFLw3tert==3

su TIME_V1SCAN if sample_final==1 & LnNFLw3tert==3
su TIME_V2SCAN if sample_final==1 & LnNFLw3tert==3



****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear


capture drop LnNFLw3tert

xtile LnNFLw3tert=LnNFLw3 if sample_final==1, nq(3)

save, replace


****w1BMI w1dxDiabetes w1Glucose w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid w1Albumin w1EosinPct w1TotalD w1currdrugs w1SRH

mi estimate: mean w1BMI if sample_final==1 & LnNFLw3tert==3
mi estimate: prop w1dxDiabetes if sample_final==1 & LnNFLw3tert==3
mi estimate: mean w1Glucose if sample_final==1 & LnNFLw3tert==3
mi estimate: mean w1Creatinine if sample_final==1 & LnNFLw3tert==3
mi estimate: mean w1USpecGrav if sample_final==1 & LnNFLw3tert==3
mi estimate: mean w1BUN if sample_final==1 & LnNFLw3tert==3
mi estimate: mean w1ALP if sample_final==1 & LnNFLw3tert==3
mi estimate: mean w1UricAcid if sample_final==1 & LnNFLw3tert==3
mi estimate: mean w1Albumin if sample_final==1 & LnNFLw3tert==3
mi estimate: mean w1EosinPct if sample_final==1 & LnNFLw3tert==3
mi estimate: mean w1TotalD if sample_final==1 & LnNFLw3tert==3
mi estimate: prop w1currdrugs if sample_final==1 & LnNFLw3tert==3
mi estimate: prop w1SRH if sample_final==1 & LnNFLw3tert==3




use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear




su LnNFLw1 LnNFLw3 if sample_final==1 & LnNFLw3tert==3, det
su LnNFLw3 if sample_final==1 & LnNFLw3tert==3, det
su bayes1LnNFL if sample_final==1 & LnNFLw3tert==3, det
su deltaLnNFL if sample_final==1 & LnNFLw3tert==3, det


su ICV_volM2 if sample_final==1 & LnNFLw3tert==3


su TOTALBRAIN if sample_final==1 & LnNFLw3tert==3
su GM if sample_final==1 & LnNFLw3tert==3
su WM if sample_final==1 & LnNFLw3tert==3


su FRONTAL_GM_L_volM2 if sample_final==1 & LnNFLw3tert==3 
su FRONTAL_WM_L_volM2 if sample_final==1 & LnNFLw3tert==3
su TEMPORAL_GM_L_volM2 if sample_final==1 & LnNFLw3tert==3
su TEMPORAL_WM_L_volM2 if sample_final==1 & LnNFLw3tert==3
su PARIETAL_GM_L_volM2 if sample_final==1 & LnNFLw3tert==3
su PARIETAL_WM_L_volM2 if sample_final==1 & LnNFLw3tert==3
su OCCIPITAL_GM_L_volM2 if sample_final==1 & LnNFLw3tert==3
su OCCIPITAL_WM_L_volM2 if sample_final==1 & LnNFLw3tert==3


su FRONTAL_GM_R_volM2 if sample_final==1 & LnNFLw3tert==3
su FRONTAL_WM_R if sample_final==1 & LnNFLw3tert==3
su TEMPORAL_GM_R_volM2 if sample_final==1 & LnNFLw3tert==3
su TEMPORAL_WM_R if sample_final==1 & LnNFLw3tert==3
su PARIETAL_GM_R_volM2 if sample_final==1 & LnNFLw3tert==3
su PARIETAL_WM_R if sample_final==1 & LnNFLw3tert==3
su OCCIPITAL_GM_R_volM2 if sample_final==1 & LnNFLw3tert==3
su OCCIPITAL_WM_R if sample_final==1 & LnNFLw3tert==3


su Left_Hippocampus if sample_final==1 & LnNFLw3tert==3
su Right_Hippocampus if sample_final==1 & LnNFLw3tert==3

su LnLesion_Volume if sample_final==1 & LnNFLw3tert==3

su Left_Hippocampuspct if sample_final==1 & LnNFLw3tert==3
su Right_Hippocampuspct if sample_final==1 & LnNFLw3tert==3

su LnLesion_Volumepct if sample_final==1 & LnNFLw3tert==3

****************BY LnNFL tertile, V2*************************
use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

**Overall**

tab Sex LnNFLw3tert if sample_final==1, row col chi 
reg w1Age LnNFLw3tert  if sample_final==1
tab w1Agebr LnNFLw3tert if sample_final==1, row col chi
tab Race LnNFLw3tert if sample_final==1, row col chi 
tab PovStat LnNFLw3tert if sample_final==1, row col chi 

reg TIME_V1SCAN LnNFLw3tert  if sample_final==1
reg TIME_V2SCAN LnNFLw3tert  if sample_final==1


****IMPUTED DATA COVARIATES*****
use finaldata_imputed, clear


mi estimate: reg w1BMI LnNFLw3tert if sample_final==1
mi estimate: mlogit w1dxDiabetes LnNFLw3tert if sample_final==1,baseoutcome(0)
mi estimate: reg w1Glucose LnNFLw3tert if sample_final==1
mi estimate: reg w1Creatinine LnNFLw3tert if sample_final==1
mi estimate: reg w1USpecGrav LnNFLw3tert if sample_final==1
mi estimate: reg w1BUN LnNFLw3tert if sample_final==1
mi estimate: reg w1ALP LnNFLw3tert if sample_final==1
mi estimate: reg w1UricAcid LnNFLw3tert if sample_final==1
mi estimate: reg w1Albumin LnNFLw3tert if sample_final==1
mi estimate: reg w1EosinPct LnNFLw3tert if sample_final==1
mi estimate: reg w1TotalD LnNFLw3tert if sample_final==1
mi estimate: reg w1currdrugs LnNFLw3tert if sample_final==1
mi estimate: mlogit w1SRH LnNFLw3tert if sample_final==1, baseoutcome(1)



mi estimate: reg w1BMI LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: mlogit w1dxDiabetes LnNFLw3tert Sex w1Age Race PovStat if sample_final==1,baseoutcome(0)
mi estimate: reg w1Glucose LnNFLw3tert Sex  w1Age Race PovStat if sample_final==1
mi estimate: reg w1Creatinine LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1USpecGrav LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1BUN LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1ALP LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1UricAcid LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1Albumin LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1EosinPct LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1TotalD LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: reg w1currdrugs LnNFLw3tert Sex w1Age Race PovStat if sample_final==1
mi estimate: mlogit w1SRH LnNFLw3tert Sex w1Age Race PovStat if sample_final==1, baseoutcome(1)

**Further adjusted for ICV_volM2**


mi estimate: reg w1BMI LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: mlogit w1dxDiabetes LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1,baseoutcome(0)
mi estimate: reg w1Glucose LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1Creatinine LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1USpecGrav LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1BUN LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1ALP LnNFLw1tert  Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1UricAcid LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1Albumin LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1EosinPct LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1TotalD LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: reg w1currdrugs LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1
mi estimate: mlogit w1SRH LnNFLw3tert Sex w1Age Race PovStat ICV_volM2 if sample_final==1, baseoutcome(1)




use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

reg LnNFLw1 LnNFLw3tert  if sample_final==1

reg LnNFLw3 LnNFLw3tert if sample_final==1

reg bayes1LnNFL LnNFLw3tert  if sample_final==1
 
reg deltaLnNFL LnNFLw3tert  if sample_final==1

save, replace

reg ICV_volM2 LnNFLw3tert  if sample_final==1


reg TOTALBRAIN LnNFLw3tert  if sample_final==1
reg GM LnNFLw3tert if sample_final==1
reg WM LnNFLw3tert  if sample_final==1


reg FRONTAL_GM_L_volM2 LnNFLw3tert  if sample_final==1 
reg FRONTAL_WM_L_volM2 LnNFLw3tert  if sample_final==1
reg TEMPORAL_GM_L_volM2 LnNFLw3tert  if sample_final==1
reg TEMPORAL_WM_L_volM2 LnNFLw3tert  if sample_final==1
reg PARIETAL_GM_L_volM2 LnNFLw3tert  if sample_final==1
reg PARIETAL_WM_L_volM2 LnNFLw3tert  if sample_final==1
reg OCCIPITAL_GM_L_volM2 LnNFLw3tert  if sample_final==1
reg OCCIPITAL_WM_L_volM2 LnNFLw3tert  if sample_final==1


reg FRONTAL_GM_R_volM2 LnNFLw3tert  if sample_final==1
reg FRONTAL_WM_R_volM2 LnNFLw3tert  if sample_final==1
reg TEMPORAL_GM_R_volM2 LnNFLw3tert  if sample_final==1
reg TEMPORAL_WM_R_volM2 LnNFLw3tert  if sample_final==1
reg PARIETAL_GM_R_volM2 LnNFLw3tert   if sample_final==1
reg PARIETAL_WM_R_volM2 LnNFLw3tert   if sample_final==1
reg OCCIPITAL_GM_R_volM2 LnNFLw3tert  if sample_final==1
reg OCCIPITAL_WM_R_volM2 LnNFLw3tert  if sample_final==1


reg Left_Hippocampus LnNFLw3tert   if sample_final==1
reg Right_Hippocampus LnNFLw3tert  if sample_final==1

reg Left_Hippocampuspct LnNFLw3tert   if sample_final==1
reg Right_Hippocampuspct LnNFLw3tert  if sample_final==1


reg LnLesion_Volume LnNFLw3tert  if sample_final==1
reg LnLesion_Volumepct LnNFLw3tert  if sample_final==1


tab NFLw1w3trackhigh LnNFLw3tert if sample_final==1, row col chi
tab NFLw1w3tracklow LnNFLw3tert if sample_final==1, row col chi

****

reg LnNFLw1 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1

reg LnNFLw3 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1

reg bayes1LnNFL LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
 
reg deltaLnNFL LnNFLw3tert w1Age Sex Race PovStat if sample_final==1

save, replace

reg ICV_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1

reg TOTALBRAIN LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg GM LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg WM LnNFLw3tert w1Age Sex Race PovStat if sample_final==1


reg FRONTAL_GM_L_volM2 LnNFLw3tert w1Age Sex Race PovStat  if sample_final==1 
reg FRONTAL_WM_L_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg TEMPORAL_GM_L_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg TEMPORAL_WM_L_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg PARIETAL_GM_L_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg PARIETAL_WM_L_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg OCCIPITAL_GM_L_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg OCCIPITAL_WM_L_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1


reg FRONTAL_GM_R_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg FRONTAL_WM_R_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg TEMPORAL_GM_R_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg TEMPORAL_WM_R_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg PARIETAL_GM_R_volM2 LnNFLw3tert w1Age Sex Race PovStat  if sample_final==1
reg PARIETAL_WM_R_volM2 LnNFLw3tert w1Age Sex Race PovStat  if sample_final==1
reg OCCIPITAL_GM_R_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg OCCIPITAL_WM_R_volM2 LnNFLw3tert w1Age Sex Race PovStat if sample_final==1


reg Left_Hippocampus LnNFLw3tert w1Age Sex Race PovStat  if sample_final==1
reg Right_Hippocampus LnNFLw3tert w1Age Sex Race PovStat if sample_final==1

reg Left_Hippocampuspct LnNFLw3tert  w1Age Sex Race PovStat if sample_final==1
reg Right_Hippocampuspct LnNFLw3tert w1Age Sex Race PovStat if sample_final==1


reg LnLesion_Volume LnNFLw3tert w1Age Sex Race PovStat if sample_final==1
reg LnLesion_Volumepct LnNFLw3tert w1Age Sex Race PovStat if sample_final==1



mlogit NFLw1w3trackhigh LnNFLw3tert w1Age Sex Race PovStat  if sample_final==1, baseoutcome(0)
mlogit NFLw1w3tracklow LnNFLw3tert w1Age Sex Race PovStat  if sample_final==1, baseoutcome(0)


mlogit NFLw1w3trackhigh LnNFLw3tert w1Age Sex Race PovStat ICV_volM2 if sample_final==1, baseoutcome(0)
mlogit NFLw1w3tracklow LnNFLw3tert w1Age Sex Race PovStat ICV_volM2 if sample_final==1, baseoutcome(0)


save, replace






***************************************TABLE 2-3: MODELS WITH HIPPOCAMPUS AND LESION VOLUME ADJUSTED FOR ICV*******************************************************

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLE2_3.smcl", replace


/////////////////////////////////////////////OVERALL//////////////////////////////////////////////////////////////////

/////////////////////////////////////ANALYSIS A: Larger brain volumes/////////////////////////////

************************
**Create Output LOG for Table 2, larger volumes**


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear



save, replace

**ANALYSIS A**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

global tflist ""
global modseq=0
foreach X of var LnNFLw1 LnNFLw3 {
foreach Y of var TOTALBRAIN GM WM {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X'  Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1, beta", label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save Outputdata_overall_A, replace

**ANALYSIS Aprime**


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

global tflist ""
global modseq=0
foreach X of var LnNFLw1 LnNFLw3 {
foreach Y of var FRONTAL_GM_L_volM2 FRONTAL_WM_L_volM2 OCCIPITAL_GM_L_volM2 OCCIPITAL_WM_L_volM2 PARIETAL_GM_L_volM2 PARIETAL_WM_L_volM2 TEMPORAL_GM_L_volM2 TEMPORAL_WM_L_volM2 ///
 FRONTAL_GM_R_volM2 FRONTAL_WM_R_volM2 OCCIPITAL_GM_R_volM2 OCCIPITAL_WM_R_volM2 PARIETAL_GM_R_volM2 PARIETAL_WM_R_volM2 TEMPORAL_GM_R_volM2 TEMPORAL_WM_R_volM2 {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X'  Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1, beta", label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save Outputdata_overall_Aprime, replace



**ANALYSIS B**


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

global tflist ""
global modseq=0
foreach X of var LnNFLw1 LnNFLw3 {
foreach Y of var Left_Hippocampus Right_Hippocampus {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X'  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1, beta", label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save Outputdata_overall_B, replace




**ANALYSIS C**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


global tflist ""
global modseq=0
foreach X of var LnNFLw1 LnNFLw3 {
foreach Y of var LnLesion_Volume {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X'  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1, beta", label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save Outputdata_overall_C, replace


**ANALYSIS D**
use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

capture rename anterior_limb_of_internal_capsul  ALIC_rightvl
capture rename anterior_limb_of_internal_capsu0 ALIC_leftvl
capture rename posterior_limb_of_internal_capsu PLIC_rightvl
capture rename posterior_limb_of_internal_caps0 PLIC_leftvl

save HANDLS_paper51_NFLBRAINSCANFINALIZED_fin, replace


foreach X of varlist _3rd_Ventricle_volM2 _4th_Ventricle_volM2 Right_Accumbens_Area_volM2 Left_Accumbens_Area_volM2 Right_Amygdala_volM2 Left_Amygdala_volM2 Brain_Stem_volM2 Right_Caudate_volM2 Left_Caudate_volM2 Right_Cerebellum_Exterior_volM2 Left_Cerebellum_Exterior_volM2 Right_Cerebellum_White_Matter_vo Left_Cerebellum_White_Matter_vol Right_Hippocampus_volM2 Left_Hippocampus_volM2 Right_Inf_Lat_Vent_volM2 Left_Inf_Lat_Vent_volM2 Right_Lateral_Ventricle_volM2 Left_Lateral_Ventricle_volM2 Right_Pallidum_volM2 Left_Pallidum_volM2 Right_Putamen_volM2 Left_Putamen_volM2 Right_Thalamus_Proper_volM2 Left_Thalamus_Proper_volM2 Right_Ventral_DC_volM2 Left_Ventral_DC_volM2 Cerebellar_Vermal_Lobules_I_V_vo Cerebellar_Vermal_Lobules_VI_VII Cerebellar_Vermal_Lobules_VIII_X Left_Basal_Forebrain_volM2 Right_Basal_Forebrain_volM2 frontal_lobe_WM_right_volM2 frontal_lobe_WM_left_volM2 occipital_lobe_WM_right_volM2 occipital_lobe_WM_left_volM2 parietal_lobe_WM_right_volM2 parietal_lobe_WM_left_volM2 temporal_lobe_WM_right_volM2 temporal_lobe_WM_left_volM2 fornix_right_volM2 fornix_left_volM2 ALIC_rightvl ALIC_leftvl PLIC_rightvl PLIC_leftvl corpus_callosum_volM2 Right_ACgG_anterior_cingulate_gy Left_ACgG_anterior_cingulate_gyr Right_AIns_anterior_insula_volM2 Left_AIns_anterior_insula_volM2 Right_AOrG_anterior_orbital_gyru Left_AOrG_anterior_orbital_gyrus Right_AnG_angular_gyrus_volM2 Left_AnG_angular_gyrus_volM2 Right_Calc_calcarine_cortex_volM Left_Calc_calcarine_cortex_volM2 Right_CO_central_operculum_volM2 Left_CO_central_operculum_volM2 Right_Cun_cuneus_volM2 Left_Cun_cuneus_volM2 Right_Ent_entorhinal_area_volM2 Left_Ent_entorhinal_area_volM2 Right_FO_frontal_operculum_volM2 Left_FO_frontal_operculum_volM2 Right_FRP_frontal_pole_volM2 Left_FRP_frontal_pole_volM2 Right_FuG_fusiform_gyrus_volM2 Left_FuG_fusiform_gyrus_volM2 Right_GRe_gyrus_rectus_volM2 Left_GRe_gyrus_rectus_volM2 Right_IOG_inferior_occipital_gyr Left_IOG_inferior_occipital_gyru Right_ITG_inferior_temporal_gyru Left_ITG_inferior_temporal_gyrus Right_LiG_lingual_gyrus_volM2 Left_LiG_lingual_gyrus_volM2 Right_LOrG_lateral_orbital_gyrus Left_LOrG_lateral_orbital_gyrus_ Right_MCgG_middle_cingulate_gyru Left_MCgG_middle_cingulate_gyrus Right_MFC_medial_frontal_cortex_ Left_MFC_medial_frontal_cortex_v Right_MFG_middle_frontal_gyrus_v Left_MFG_middle_frontal_gyrus_vo Right_MOG_middle_occipital_gyrus Left_MOG_middle_occipital_gyrus_ Right_MOrG_medial_orbital_gyrus_ Left_MOrG_medial_orbital_gyrus_v Right_MPoG_postcentral_gyrus_med Left_MPoG_postcentral_gyrus_medi Right_MPrG_precentral_gyrus_medi Left_MPrG_precentral_gyrus_media Right_MSFG_superior_frontal_gyru Left_MSFG_superior_frontal_gyrus Right_MTG_middle_temporal_gyrus_ Left_MTG_middle_temporal_gyrus_v Right_OCP_occipital_pole_volM2 Left_OCP_occipital_pole_volM2 Right_OFuG_occipital_fusiform_gy Left_OFuG_occipital_fusiform_gyr Right_OpIFG_opercular_part_of_th Left_OpIFG_opercular_part_of_the Right_OrIFG_orbital_part_of_the_ Left_OrIFG_orbital_part_of_the_i Right_PCgG_posterior_cingulate_g Left_PCgG_posterior_cingulate_gy Right_PCu_precuneus_volM2 Left_PCu_precuneus_volM2 Right_PHG_parahippocampal_gyrus_ Left_PHG_parahippocampal_gyrus_v Right_PIns_posterior_insula_volM Left_PIns_posterior_insula_volM2 Right_PO_parietal_operculum_volM Left_PO_parietal_operculum_volM2 Right_PoG_postcentral_gyrus_volM Left_PoG_postcentral_gyrus_volM2 Right_POrG_posterior_orbital_gyr Left_POrG_posterior_orbital_gyru Right_PP_planum_polare_volM2 Left_PP_planum_polare_volM2 Right_PrG_precentral_gyrus_volM2 Left_PrG_precentral_gyrus_volM2 Right_PT_planum_temporale_volM2 Left_PT_planum_temporale_volM2 Right_SCA_subcallosal_area_volM2 Left_SCA_subcallosal_area_volM2 Right_SFG_superior_frontal_gyrus Left_SFG_superior_frontal_gyrus_ Right_SMC_supplementary_motor_co Left_SMC_supplementary_motor_cor Right_SMG_supramarginal_gyrus_vo Left_SMG_supramarginal_gyrus_vol Right_SOG_superior_occipital_gyr Left_SOG_superior_occipital_gyru Right_SPL_superior_parietal_lobu Left_SPL_superior_parietal_lobul Right_STG_superior_temporal_gyru Left_STG_superior_temporal_gyrus Right_TMP_temporal_pole_volM2 Left_TMP_temporal_pole_volM2 Right_TrIFG_triangular_part_of_t Left_TrIFG_triangular_part_of_th Right_TTG_transverse_temporal_gy Left_TTG_transverse_temporal_gyr {
renvars `X', postdrop(1) 
}


foreach X2 of varlist _3rd_Ventricle_volM _4th_Ventricle_volM Right_Accumbens_Area_volM Left_Accumbens_Area_volM Right_Amygdala_volM Left_Amygdala_volM Brain_Stem_volM Right_Caudate_volM Left_Caudate_volM Right_Cerebellum_Exterior_volM Left_Cerebellum_Exterior_volM Right_Cerebellum_White_Matter_v Left_Cerebellum_White_Matter_vo Right_Hippocampus_volM Left_Hippocampus_volM Right_Inf_Lat_Vent_volM Left_Inf_Lat_Vent_volM Right_Lateral_Ventricle_volM Left_Lateral_Ventricle_volM Right_Pallidum_volM Left_Pallidum_volM Right_Putamen_volM Left_Putamen_volM Right_Thalamus_Proper_volM Left_Thalamus_Proper_volM Right_Ventral_DC_volM Left_Ventral_DC_volM Cerebellar_Vermal_Lobules_I_V_v Cerebellar_Vermal_Lobules_VI_VI Cerebellar_Vermal_Lobules_VIII_ Left_Basal_Forebrain_volM Right_Basal_Forebrain_volM frontal_lobe_WM_right_volM frontal_lobe_WM_left_volM occipital_lobe_WM_right_volM occipital_lobe_WM_left_volM parietal_lobe_WM_right_volM parietal_lobe_WM_left_volM temporal_lobe_WM_right_volM temporal_lobe_WM_left_volM fornix_right_volM fornix_left_volM ALIC_rightv ALIC_leftv PLIC_rightv PLIC_leftv corpus_callosum_volM Right_ACgG_anterior_cingulate_g Left_ACgG_anterior_cingulate_gy Right_AIns_anterior_insula_volM Left_AIns_anterior_insula_volM Right_AOrG_anterior_orbital_gyr Left_AOrG_anterior_orbital_gyru Right_AnG_angular_gyrus_volM Left_AnG_angular_gyrus_volM Right_Calc_calcarine_cortex_vol Left_Calc_calcarine_cortex_volM Right_CO_central_operculum_volM Left_CO_central_operculum_volM Right_Cun_cuneus_volM Left_Cun_cuneus_volM Right_Ent_entorhinal_area_volM Left_Ent_entorhinal_area_volM Right_FO_frontal_operculum_volM Left_FO_frontal_operculum_volM Right_FRP_frontal_pole_volM Left_FRP_frontal_pole_volM Right_FuG_fusiform_gyrus_volM Left_FuG_fusiform_gyrus_volM Right_GRe_gyrus_rectus_volM Left_GRe_gyrus_rectus_volM Right_IOG_inferior_occipital_gy Left_IOG_inferior_occipital_gyr Right_ITG_inferior_temporal_gyr Left_ITG_inferior_temporal_gyru Right_LiG_lingual_gyrus_volM Left_LiG_lingual_gyrus_volM Right_LOrG_lateral_orbital_gyru Left_LOrG_lateral_orbital_gyrus Right_MCgG_middle_cingulate_gyr Left_MCgG_middle_cingulate_gyru Right_MFC_medial_frontal_cortex Left_MFC_medial_frontal_cortex_ Right_MFG_middle_frontal_gyrus_ Left_MFG_middle_frontal_gyrus_v Right_MOG_middle_occipital_gyru Left_MOG_middle_occipital_gyrus Right_MOrG_medial_orbital_gyrus Left_MOrG_medial_orbital_gyrus_ Right_MPoG_postcentral_gyrus_me Left_MPoG_postcentral_gyrus_med Right_MPrG_precentral_gyrus_med Left_MPrG_precentral_gyrus_medi Right_MSFG_superior_frontal_gyr Left_MSFG_superior_frontal_gyru Right_MTG_middle_temporal_gyrus Left_MTG_middle_temporal_gyrus_ Right_OCP_occipital_pole_volM Left_OCP_occipital_pole_volM Right_OFuG_occipital_fusiform_g Left_OFuG_occipital_fusiform_gy Right_OpIFG_opercular_part_of_t Left_OpIFG_opercular_part_of_th Right_OrIFG_orbital_part_of_the Left_OrIFG_orbital_part_of_the_ Right_PCgG_posterior_cingulate_ Left_PCgG_posterior_cingulate_g Right_PCu_precuneus_volM Left_PCu_precuneus_volM Right_PHG_parahippocampal_gyrus Left_PHG_parahippocampal_gyrus_ Right_PIns_posterior_insula_vol Left_PIns_posterior_insula_volM Right_PO_parietal_operculum_vol Left_PO_parietal_operculum_volM Right_PoG_postcentral_gyrus_vol Left_PoG_postcentral_gyrus_volM Right_POrG_posterior_orbital_gy Left_POrG_posterior_orbital_gyr Right_PP_planum_polare_volM Left_PP_planum_polare_volM Right_PrG_precentral_gyrus_volM Left_PrG_precentral_gyrus_volM Right_PT_planum_temporale_volM Left_PT_planum_temporale_volM Right_SCA_subcallosal_area_volM Left_SCA_subcallosal_area_volM Right_SFG_superior_frontal_gyru Left_SFG_superior_frontal_gyrus Right_SMC_supplementary_motor_c Left_SMC_supplementary_motor_co Right_SMG_supramarginal_gyrus_v Left_SMG_supramarginal_gyrus_vo Right_SOG_superior_occipital_gy Left_SOG_superior_occipital_gyr Right_SPL_superior_parietal_lob Left_SPL_superior_parietal_lobu Right_STG_superior_temporal_gyr Left_STG_superior_temporal_gyru Right_TMP_temporal_pole_volM Left_TMP_temporal_pole_volM Right_TrIFG_triangular_part_of_ Left_TrIFG_triangular_part_of_t Right_TTG_transverse_temporal_g Left_TTG_transverse_temporal_gy {
capture drop z`X2' 
	egen z`X2'=std(`X2') if sample_final==1 
}

save, replace

capture drop zLnNFLw1
egen zLnNFLw1=std(LnNFLw1) if sample_final==1

capture drop zLnNFLw3
egen zLnNFLw3=std(LnNFLw3) if sample_final==1

save, replace

global tflist ""
global modseq=0
foreach X of var zLnNFLw1 zLnNFLw3 {
foreach Y of var z_3rd_Ventricle_volM- zLeft_TTG_transverse_temporal_gy {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X'  Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1, beta", label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save zOutputdata_overall_D, replace


**ANALYSIS E: Small ROIs, WML**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

capture rename anterior_limb_of_internal_capsu1  ALIC_rightwml
capture rename anterior_limb_of_internal_capsu2 ALIC_leftwml
capture rename posterior_limb_of_internal_caps1 PLIC_rightwml
capture rename posterior_limb_of_internal_caps2 PLIC_leftwml

save HANDLS_paper51_NFLBRAINSCANFINALIZED_fin, replace


foreach X of varlist _3rd_Ventricle_wmlM2 _4th_Ventricle_wmlM2 Right_Accumbens_Area_wmlM2 Left_Accumbens_Area_wmlM2 Right_Amygdala_wmlM2 Left_Amygdala_wmlM2 Brain_Stem_wmlM2 Right_Caudate_wmlM2 Left_Caudate_wmlM2 Right_Cerebellum_Exterior_wmlM2 Left_Cerebellum_Exterior_wmlM2 Right_Cerebellum_White_Matter_wm Left_Cerebellum_White_Matter_wml Right_Hippocampus_wmlM2 Left_Hippocampus_wmlM2 Right_Inf_Lat_Vent_wmlM2 Left_Inf_Lat_Vent_wmlM2 Right_Lateral_Ventricle_wmlM2 Left_Lateral_Ventricle_wmlM2 Right_Pallidum_wmlM2 Left_Pallidum_wmlM2 Right_Putamen_wmlM2 Left_Putamen_wmlM2 Right_Thalamus_Proper_wmlM2 Left_Thalamus_Proper_wmlM2 Right_Ventral_DC_wmlM2 Left_Ventral_DC_wmlM2 Cerebellar_Vermal_Lobules_I_V_wm Cerebellar_Vermal_Lobules_VI_VI0 Cerebellar_Vermal_Lobules_VIII_0 Left_Basal_Forebrain_wmlM2 Right_Basal_Forebrain_wmlM2 frontal_lobe_WM_right_wmlM2 frontal_lobe_WM_left_wmlM2 occipital_lobe_WM_right_wmlM2 occipital_lobe_WM_left_wmlM2 parietal_lobe_WM_right_wmlM2 parietal_lobe_WM_left_wmlM2 temporal_lobe_WM_right_wmlM2 temporal_lobe_WM_left_wmlM2 fornix_right_wmlM2 fornix_left_wmlM2 ALIC_rightwml ALIC_leftwml PLIC_rightwml PLIC_leftwml corpus_callosum_wmlM2 Right_ACgG_anterior_cingulate_g0 Left_ACgG_anterior_cingulate_gy0 Right_AIns_anterior_insula_wmlM2 Left_AIns_anterior_insula_wmlM2 Right_AOrG_anterior_orbital_gyr0 Left_AOrG_anterior_orbital_gyru0 Right_AnG_angular_gyrus_wmlM2 Left_AnG_angular_gyrus_wmlM2 Right_Calc_calcarine_cortex_wmlM Left_Calc_calcarine_cortex_wmlM2 Right_CO_central_operculum_wmlM2 Left_CO_central_operculum_wmlM2 Right_Cun_cuneus_wmlM2 Left_Cun_cuneus_wmlM2 Right_Ent_entorhinal_area_wmlM2 Left_Ent_entorhinal_area_wmlM2 Right_FO_frontal_operculum_wmlM2 Left_FO_frontal_operculum_wmlM2 Right_FRP_frontal_pole_wmlM2 Left_FRP_frontal_pole_wmlM2 Right_FuG_fusiform_gyrus_wmlM2 Left_FuG_fusiform_gyrus_wmlM2 Right_GRe_gyrus_rectus_wmlM2 Left_GRe_gyrus_rectus_wmlM2 Right_IOG_inferior_occipital_gy0 Left_IOG_inferior_occipital_gyr0 Right_ITG_inferior_temporal_gyr0 Left_ITG_inferior_temporal_gyru0 Right_LiG_lingual_gyrus_wmlM2 Left_LiG_lingual_gyrus_wmlM2 Right_LOrG_lateral_orbital_gyru0 Left_LOrG_lateral_orbital_gyrus0 Right_MCgG_middle_cingulate_gyr0 Left_MCgG_middle_cingulate_gyru0 Right_MFC_medial_frontal_cortex0 Left_MFC_medial_frontal_cortex_w Right_MFG_middle_frontal_gyrus_w Left_MFG_middle_frontal_gyrus_wm Right_MOG_middle_occipital_gyru0 Left_MOG_middle_occipital_gyrus0 Right_MOrG_medial_orbital_gyrus0 Left_MOrG_medial_orbital_gyrus_w Right_MPoG_postcentral_gyrus_me0 Left_MPoG_postcentral_gyrus_med0 Right_MPrG_precentral_gyrus_med0 Left_MPrG_precentral_gyrus_medi0 Right_MSFG_superior_frontal_gyr0 Left_MSFG_superior_frontal_gyru0 Right_MTG_middle_temporal_gyrus0 Left_MTG_middle_temporal_gyrus_w Right_OCP_occipital_pole_wmlM2 Left_OCP_occipital_pole_wmlM2 Right_OFuG_occipital_fusiform_g0 Left_OFuG_occipital_fusiform_gy0 Right_OpIFG_opercular_part_of_t0 Left_OpIFG_opercular_part_of_th0 Right_OrIFG_orbital_part_of_the0 Left_OrIFG_orbital_part_of_the_0 Right_PCgG_posterior_cingulate_0 Left_PCgG_posterior_cingulate_g0 Right_PCu_precuneus_wmlM2 Left_PCu_precuneus_wmlM2 Right_PHG_parahippocampal_gyrus0 Left_PHG_parahippocampal_gyrus_w Right_PIns_posterior_insula_wmlM Left_PIns_posterior_insula_wmlM2 Right_PO_parietal_operculum_wmlM Left_PO_parietal_operculum_wmlM2 Right_PoG_postcentral_gyrus_wmlM Left_PoG_postcentral_gyrus_wmlM2 Right_POrG_posterior_orbital_gy0 Left_POrG_posterior_orbital_gyr0 Right_PP_planum_polare_wmlM2 Left_PP_planum_polare_wmlM2 Right_PrG_precentral_gyrus_wmlM2 Left_PrG_precentral_gyrus_wmlM2 Right_PT_planum_temporale_wmlM2 Left_PT_planum_temporale_wmlM2 Right_SCA_subcallosal_area_wmlM2 Left_SCA_subcallosal_area_wmlM2 Right_SFG_superior_frontal_gyru0 Left_SFG_superior_frontal_gyrus0 Right_SMC_supplementary_motor_c0 Left_SMC_supplementary_motor_co0 Right_SMG_supramarginal_gyrus_wm Left_SMG_supramarginal_gyrus_wml Right_SOG_superior_occipital_gy0 Left_SOG_superior_occipital_gyr0 Right_SPL_superior_parietal_lob0 Left_SPL_superior_parietal_lobu0 Right_STG_superior_temporal_gyr0 Left_STG_superior_temporal_gyru0 Right_TMP_temporal_pole_wmlM2 Left_TMP_temporal_pole_wmlM2 Right_TrIFG_triangular_part_of_0 Left_TrIFG_triangular_part_of_t0 Right_TTG_transverse_temporal_g0 Left_TTG_transverse_temporal_gy0 {
renvars `X', postdrop(2) 
}

save HANDLS_paper51_NFLBRAINSCANFINALIZED_fin, replace



foreach X1 of varlist _3rd_Ventricle_wml _4th_Ventricle_wml Right_Accumbens_Area_wml Left_Accumbens_Area_wml Right_Amygdala_wml Left_Amygdala_wml Brain_Stem_wml Right_Caudate_wml Left_Caudate_wml Right_Cerebellum_Exterior_wml Left_Cerebellum_Exterior_wml Right_Cerebellum_White_Matter_ Left_Cerebellum_White_Matter_w Right_Hippocampus_wml Left_Hippocampus_wml Right_Inf_Lat_Vent_wml Left_Inf_Lat_Vent_wml Right_Lateral_Ventricle_wml Left_Lateral_Ventricle_wml Right_Pallidum_wml Left_Pallidum_wml Right_Putamen_wml Left_Putamen_wml Right_Thalamus_Proper_wml Left_Thalamus_Proper_wml Right_Ventral_DC_wml Left_Ventral_DC_wml Cerebellar_Vermal_Lobules_I_V_ Cerebellar_Vermal_Lobules_VI_V Cerebellar_Vermal_Lobules_VIII Left_Basal_Forebrain_wml Right_Basal_Forebrain_wml frontal_lobe_WM_right_wml frontal_lobe_WM_left_wml occipital_lobe_WM_right_wml occipital_lobe_WM_left_wml parietal_lobe_WM_right_wml parietal_lobe_WM_left_wml temporal_lobe_WM_right_wml temporal_lobe_WM_left_wml fornix_right_wml fornix_left_wml ALIC_rightw ALIC_leftw PLIC_rightw PLIC_leftw corpus_callosum_wml Right_ACgG_anterior_cingulate_ Left_ACgG_anterior_cingulate_g Right_AIns_anterior_insula_wml Left_AIns_anterior_insula_wml Right_AOrG_anterior_orbital_gy Left_AOrG_anterior_orbital_gyr Right_AnG_angular_gyrus_wml Left_AnG_angular_gyrus_wml Right_Calc_calcarine_cortex_wm Left_Calc_calcarine_cortex_wml Right_CO_central_operculum_wml Left_CO_central_operculum_wml Right_Cun_cuneus_wml Left_Cun_cuneus_wml Right_Ent_entorhinal_area_wml Left_Ent_entorhinal_area_wml Right_FO_frontal_operculum_wml Left_FO_frontal_operculum_wml Right_FRP_frontal_pole_wml Left_FRP_frontal_pole_wml Right_FuG_fusiform_gyrus_wml Left_FuG_fusiform_gyrus_wml Right_GRe_gyrus_rectus_wml Left_GRe_gyrus_rectus_wml Right_IOG_inferior_occipital_g Left_IOG_inferior_occipital_gy Right_ITG_inferior_temporal_gy Left_ITG_inferior_temporal_gyr Right_LiG_lingual_gyrus_wml Left_LiG_lingual_gyrus_wml Right_LOrG_lateral_orbital_gyr Left_LOrG_lateral_orbital_gyru Right_MCgG_middle_cingulate_gy Left_MCgG_middle_cingulate_gyr Right_MFC_medial_frontal_corte Left_MFC_medial_frontal_cortex Right_MFG_middle_frontal_gyrus Left_MFG_middle_frontal_gyrus_ Right_MOG_middle_occipital_gyr Left_MOG_middle_occipital_gyru Right_MOrG_medial_orbital_gyru Left_MOrG_medial_orbital_gyrus Right_MPoG_postcentral_gyrus_m Left_MPoG_postcentral_gyrus_me Right_MPrG_precentral_gyrus_me Left_MPrG_precentral_gyrus_med Right_MSFG_superior_frontal_gy Left_MSFG_superior_frontal_gyr Right_MTG_middle_temporal_gyru Left_MTG_middle_temporal_gyrus Right_OCP_occipital_pole_wml Left_OCP_occipital_pole_wml Right_OFuG_occipital_fusiform_ Left_OFuG_occipital_fusiform_g Right_OpIFG_opercular_part_of_ Left_OpIFG_opercular_part_of_t Right_OrIFG_orbital_part_of_th Left_OrIFG_orbital_part_of_the Right_PCgG_posterior_cingulate Left_PCgG_posterior_cingulate_ Right_PCu_precuneus_wml Left_PCu_precuneus_wml Right_PHG_parahippocampal_gyru Left_PHG_parahippocampal_gyrus Right_PIns_posterior_insula_wm Left_PIns_posterior_insula_wml Right_PO_parietal_operculum_wm Left_PO_parietal_operculum_wml Right_PoG_postcentral_gyrus_wm Left_PoG_postcentral_gyrus_wml Right_POrG_posterior_orbital_g Left_POrG_posterior_orbital_gy Right_PP_planum_polare_wml Left_PP_planum_polare_wml Right_PrG_precentral_gyrus_wml Left_PrG_precentral_gyrus_wml Right_PT_planum_temporale_wml Left_PT_planum_temporale_wml Right_SCA_subcallosal_area_wml Left_SCA_subcallosal_area_wml Right_SFG_superior_frontal_gyr Left_SFG_superior_frontal_gyru Right_SMC_supplementary_motor_ Left_SMC_supplementary_motor_c Right_SMG_supramarginal_gyrus_ Left_SMG_supramarginal_gyrus_w Right_SOG_superior_occipital_g Left_SOG_superior_occipital_gy Right_SPL_superior_parietal_lo Left_SPL_superior_parietal_lob Right_STG_superior_temporal_gy Left_STG_superior_temporal_gyr Right_TMP_temporal_pole_wml Left_TMP_temporal_pole_wml Right_TrIFG_triangular_part_of Left_TrIFG_triangular_part_of_ Right_TTG_transverse_temporal_ Left_TTG_transverse_temporal_g {
capture drop z`X1' 
	gen L`X1'=(`X1')^(1/3) if sample_final==1
	
}

save HANDLS_paper51_NFLBRAINSCANFINALIZED_fin, replace

**Select non-missing from _3rd_Ventricle_wmlM through Left_TTG_transverse_temporal_gy: N=60**


tab1 _3rd_Ventricle_wml _4th_Ventricle_wml Right_Accumbens_Area_wml Left_Accumbens_Area_wml Brain_Stem_wml Right_Caudate_wml Left_Caudate_wml Right_Cerebellum_Exterior_wml Right_Hippocampus_wml Left_Hippocampus_wml Right_Inf_Lat_Vent_wml Left_Inf_Lat_Vent_wml Right_Lateral_Ventricle_wml Left_Lateral_Ventricle_wml Right_Putamen_wml Left_Putamen_wml Right_Thalamus_Proper_wml Left_Thalamus_Proper_wml Right_Ventral_DC_wml Left_Ventral_DC_wml  Left_Basal_Forebrain_wml frontal_lobe_WM_right_wml frontal_lobe_WM_left_wml occipital_lobe_WM_right_wml occipital_lobe_WM_left_wml parietal_lobe_WM_right_wml parietal_lobe_WM_left_wml temporal_lobe_WM_right_wml temporal_lobe_WM_left_wml fornix_right_wml fornix_left_wml ALIC_rightw ALIC_leftw PLIC_rightw PLIC_leftw corpus_callosum_wml Right_ACgG_anterior_cingulate_ Right_AIns_anterior_insula_wml Right_AOrG_anterior_orbital_gy Right_AnG_angular_gyrus_wml Left_AnG_angular_gyrus_wml Right_Calc_calcarine_cortex_wm Left_Calc_calcarine_cortex_wml Right_Cun_cuneus_wml Right_LiG_lingual_gyrus_wml Left_LiG_lingual_gyrus_wml Right_MFC_medial_frontal_corte Right_MOG_middle_occipital_gyr Right_MOrG_medial_orbital_gyru Right_MSFG_superior_frontal_gy Right_MTG_middle_temporal_gyru Right_PCu_precuneus_wml Left_PCu_precuneus_wml Left_PIns_posterior_insula_wml Right_PO_parietal_operculum_wm Right_PoG_postcentral_gyrus_wm Right_PrG_precentral_gyrus_wml Left_SFG_superior_frontal_gyru Right_SMC_supplementary_motor_ Right_SMG_supramarginal_gyrus_ Right_SOG_superior_occipital_g if sample_final==1


su L_3rd_Ventricle_wml L_4th_Ventricle_wml LRight_Accumbens_Area_wml LLeft_Accumbens_Area_wml LBrain_Stem_wml LRight_Caudate_wml LLeft_Caudate_wml LRight_Cerebellum_Exterior_wml LRight_Hippocampus_wml LLeft_Hippocampus_wml LRight_Inf_Lat_Vent_wml LLeft_Inf_Lat_Vent_wml LRight_Lateral_Ventricle_wml LLeft_Lateral_Ventricle_wml LRight_Putamen_wml LLeft_Putamen_wml LRight_Thalamus_Proper_wml LLeft_Thalamus_Proper_wml LRight_Ventral_DC_wml LLeft_Ventral_DC_wml  LLeft_Basal_Forebrain_wml Lfrontal_lobe_WM_right_wml Lfrontal_lobe_WM_left_wml Loccipital_lobe_WM_right_wml Loccipital_lobe_WM_left_wml Lparietal_lobe_WM_right_wml Lparietal_lobe_WM_left_wml Ltemporal_lobe_WM_right_wml Ltemporal_lobe_WM_left_wml Lfornix_right_wml Lfornix_left_wml LALIC_rightw LALIC_leftw LPLIC_rightw LPLIC_leftw Lcorpus_callosum_wml LRight_ACgG_anterior_cingulate_ LRight_AIns_anterior_insula_wml LRight_AOrG_anterior_orbital_gy LRight_AnG_angular_gyrus_wml LLeft_AnG_angular_gyrus_wml LRight_Calc_calcarine_cortex_wm LLeft_Calc_calcarine_cortex_wml LRight_Cun_cuneus_wml LRight_LiG_lingual_gyrus_wml LLeft_LiG_lingual_gyrus_wml LRight_MFC_medial_frontal_corte LRight_MOG_middle_occipital_gyr LRight_MOrG_medial_orbital_gyru LRight_MSFG_superior_frontal_gy LRight_MTG_middle_temporal_gyru LRight_PCu_precuneus_wml LLeft_PCu_precuneus_wml LLeft_PIns_posterior_insula_wml LRight_PO_parietal_operculum_wm LRight_PoG_postcentral_gyrus_wm LRight_PrG_precentral_gyrus_wml LLeft_SFG_superior_frontal_gyru LRight_SMC_supplementary_motor_ LRight_SMG_supramarginal_gyrus_ LRight_SOG_superior_occipital_g if sample_final==1, det

**Create all z-scores**
foreach X2 of varlist L_3rd_Ventricle_wml L_4th_Ventricle_wml LRight_Accumbens_Area_wml LLeft_Accumbens_Area_wml LBrain_Stem_wml LRight_Caudate_wml LLeft_Caudate_wml LRight_Cerebellum_Exterior_wml LRight_Hippocampus_wml LLeft_Hippocampus_wml LRight_Inf_Lat_Vent_wml LLeft_Inf_Lat_Vent_wml LRight_Lateral_Ventricle_wml LLeft_Lateral_Ventricle_wml LRight_Putamen_wml LLeft_Putamen_wml LRight_Thalamus_Proper_wml LLeft_Thalamus_Proper_wml LRight_Ventral_DC_wml LLeft_Ventral_DC_wml  LLeft_Basal_Forebrain_wml Lfrontal_lobe_WM_right_wml Lfrontal_lobe_WM_left_wml Loccipital_lobe_WM_right_wml Loccipital_lobe_WM_left_wml Lparietal_lobe_WM_right_wml Lparietal_lobe_WM_left_wml Ltemporal_lobe_WM_right_wml Ltemporal_lobe_WM_left_wml Lfornix_right_wml Lfornix_left_wml LALIC_rightw LALIC_leftw LPLIC_rightw LPLIC_leftw Lcorpus_callosum_wml LRight_ACgG_anterior_cingulate_ LRight_AIns_anterior_insula_wml LRight_AOrG_anterior_orbital_gy LRight_AnG_angular_gyrus_wml LLeft_AnG_angular_gyrus_wml LRight_Calc_calcarine_cortex_wm LLeft_Calc_calcarine_cortex_wml LRight_Cun_cuneus_wml LRight_LiG_lingual_gyrus_wml LLeft_LiG_lingual_gyrus_wml LRight_MFC_medial_frontal_corte LRight_MOG_middle_occipital_gyr LRight_MOrG_medial_orbital_gyru LRight_MSFG_superior_frontal_gy LRight_MTG_middle_temporal_gyru LRight_PCu_precuneus_wml LLeft_PCu_precuneus_wml LLeft_PIns_posterior_insula_wml LRight_PO_parietal_operculum_wm LRight_PoG_postcentral_gyrus_wm LRight_PrG_precentral_gyrus_wml LLeft_SFG_superior_frontal_gyru LRight_SMC_supplementary_motor_ LRight_SMG_supramarginal_gyrus_ LRight_SOG_superior_occipital_g  { 
capture drop z`X2' 
	egen z`X2'=std(`X2') if sample_final==1
	
}


**Exclude volumes with <5% have non-zero volumes, 24 volumes left**

foreach X3 of varlist LRight_Caudate_wml  LLeft_Caudate_wml LRight_Hippocampus_wml LLeft_Hippocampus_wml  LRight_Lateral_Ventricle_wml LLeft_Lateral_Ventricle_wml  LRight_Thalamus_Proper_wml LLeft_Thalamus_Proper_wml Lfrontal_lobe_WM_right_wml Lfrontal_lobe_WM_left_wml Loccipital_lobe_WM_right_wml Loccipital_lobe_WM_left_wml Lparietal_lobe_WM_right_wml Lparietal_lobe_WM_left_wml Ltemporal_lobe_WM_right_wml Ltemporal_lobe_WM_left_wml  Lfornix_left_wml LALIC_rightw LALIC_leftw LPLIC_rightw LPLIC_leftw Lcorpus_callosum_wml   LRight_Calc_calcarine_cortex_wm LLeft_Calc_calcarine_cortex_wml   { 
capture drop z`X3' 
	egen z`X3'=std(`X3') if sample_final==1 
	
}




save HANDLS_paper51_NFLBRAINSCANFINALIZED_fin, replace


**N=24 ROIs from previous step, Exclude volumes that have gray matter: Exclude xx, xx, xx, xx and xx**
**LRight_Lateral_Ventricle_wml LLeft_Lateral_Ventricle_wml 


use HANDLS_paper51_NFLBRAINSCANFINALIZED_fin,clear


capture drop zLnNFLw1
egen zLnNFLw1=std(LnNFLw1) if sample_final==1

capture drop zLnNFLw3
egen zLnNFLw3=std(LnNFLw3) if sample_final==1

save, replace




global tflist ""
global modseq=0
foreach X of var zLnNFLw1 zLnNFLw3 {
foreach Y of var zLRight_Caudate_wml  zLLeft_Caudate_wml zLRight_Hippocampus_wml zLLeft_Hippocampus_wml  zLRight_Lateral_Ventricle_wml zLLeft_Lateral_Ventricle_wml  zLRight_Thalamus_Proper_wml zLLeft_Thalamus_Proper_wml zLfrontal_lobe_WM_right_wml zLfrontal_lobe_WM_left_wml zLoccipital_lobe_WM_right_wml zLoccipital_lobe_WM_left_wml zLparietal_lobe_WM_right_wml zLparietal_lobe_WM_left_wml zLtemporal_lobe_WM_right_wml zLtemporal_lobe_WM_left_wml  zLfornix_left_wml zLALIC_rightw zLALIC_leftw zLPLIC_rightw zLPLIC_leftw zLcorpus_callosum_wml   zLRight_Calc_calcarine_cortex_wm zLLeft_Calc_calcarine_cortex_wml {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X'  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2  if sample_final==1, beta", label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }

capture drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save zOutputdata_overall_E, replace

keep if parmseq==1
save zOutputdata_overall_Eparm1,replace


////////////////////////////////////////////////////////////////////////BY Sex ////////////////////////////////////////////////////////////////

**ANALYSIS A**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


global tflist ""
global modseq=0
foreach X of var LnNFLw1 LnNFLw3 {
foreach Y of var TOTALBRAIN GM WM  {
foreach Z of var Sex {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X' Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1, beta", by(`Z') label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 }
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save Outputdata_bysociodem_A, replace


**ANALYSIS Aprime**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


global tflist ""
global modseq=0
foreach X of var LnNFLw1 LnNFLw3 {
foreach Y of var FRONTAL_GM_L_volM2 FRONTAL_WM_L_volM2 OCCIPITAL_GM_L_volM2 OCCIPITAL_WM_L_volM2 PARIETAL_GM_L_volM2 PARIETAL_WM_L_volM2 TEMPORAL_GM_L_volM2 TEMPORAL_WM_L_volM2 ///
 FRONTAL_GM_R_volM2 FRONTAL_WM_R_volM2 OCCIPITAL_GM_R_volM2 OCCIPITAL_WM_R_volM2 PARIETAL_GM_R_volM2 PARIETAL_WM_R_volM2 TEMPORAL_GM_R_volM2 TEMPORAL_WM_R_volM2  {
foreach Z of var Sex {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X' Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1, beta", by(`Z') label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 }
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save Outputdata_bysociodem_Aprime, replace

**ANALYSIS B**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


global tflist ""
global modseq=0
foreach X of var LnNFLw1 LnNFLw3 {
foreach Y of var Left_Hippocampus Right_Hippocampus  {
foreach Z of var Sex {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X' Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1, beta", by(`Z') label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 }
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save Outputdata_bysociodem_B, replace

**ANALYSIS C**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


global tflist ""
global modseq=0
foreach X of var LnNFLw1 LnNFLw3 {
foreach Y of var LnLesion_Volume  {
foreach Z of var Sex {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X' Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1, beta", by(`Z') label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 }
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save Outputdata_bysociodem_C, replace


**ANALYSIS D**
use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

capture rename anterior_limb_of_internal_capsul  ALIC_rightvl
capture rename anterior_limb_of_internal_capsu0 ALIC_leftvl
capture rename posterior_limb_of_internal_capsu PLIC_rightvl
capture rename posterior_limb_of_internal_caps0 PLIC_leftvl

save HANDLS_paper51_NFLBRAINSCANFINALIZED_fin, replace


foreach X of varlist _3rd_Ventricle_volM2 _4th_Ventricle_volM2 Right_Accumbens_Area_volM2 Left_Accumbens_Area_volM2 Right_Amygdala_volM2 Left_Amygdala_volM2 Brain_Stem_volM2 Right_Caudate_volM2 Left_Caudate_volM2 Right_Cerebellum_Exterior_volM2 Left_Cerebellum_Exterior_volM2 Right_Cerebellum_White_Matter_vo Left_Cerebellum_White_Matter_vol Right_Hippocampus_volM2 Left_Hippocampus_volM2 Right_Inf_Lat_Vent_volM2 Left_Inf_Lat_Vent_volM2 Right_Lateral_Ventricle_volM2 Left_Lateral_Ventricle_volM2 Right_Pallidum_volM2 Left_Pallidum_volM2 Right_Putamen_volM2 Left_Putamen_volM2 Right_Thalamus_Proper_volM2 Left_Thalamus_Proper_volM2 Right_Ventral_DC_volM2 Left_Ventral_DC_volM2 Cerebellar_Vermal_Lobules_I_V_vo Cerebellar_Vermal_Lobules_VI_VII Cerebellar_Vermal_Lobules_VIII_X Left_Basal_Forebrain_volM2 Right_Basal_Forebrain_volM2 frontal_lobe_WM_right_volM2 frontal_lobe_WM_left_volM2 occipital_lobe_WM_right_volM2 occipital_lobe_WM_left_volM2 parietal_lobe_WM_right_volM2 parietal_lobe_WM_left_volM2 temporal_lobe_WM_right_volM2 temporal_lobe_WM_left_volM2 fornix_right_volM2 fornix_left_volM2 ALIC_rightvl ALIC_leftvl PLIC_rightvl PLIC_leftvl corpus_callosum_volM2 Right_ACgG_anterior_cingulate_gy Left_ACgG_anterior_cingulate_gyr Right_AIns_anterior_insula_volM2 Left_AIns_anterior_insula_volM2 Right_AOrG_anterior_orbital_gyru Left_AOrG_anterior_orbital_gyrus Right_AnG_angular_gyrus_volM2 Left_AnG_angular_gyrus_volM2 Right_Calc_calcarine_cortex_volM Left_Calc_calcarine_cortex_volM2 Right_CO_central_operculum_volM2 Left_CO_central_operculum_volM2 Right_Cun_cuneus_volM2 Left_Cun_cuneus_volM2 Right_Ent_entorhinal_area_volM2 Left_Ent_entorhinal_area_volM2 Right_FO_frontal_operculum_volM2 Left_FO_frontal_operculum_volM2 Right_FRP_frontal_pole_volM2 Left_FRP_frontal_pole_volM2 Right_FuG_fusiform_gyrus_volM2 Left_FuG_fusiform_gyrus_volM2 Right_GRe_gyrus_rectus_volM2 Left_GRe_gyrus_rectus_volM2 Right_IOG_inferior_occipital_gyr Left_IOG_inferior_occipital_gyru Right_ITG_inferior_temporal_gyru Left_ITG_inferior_temporal_gyrus Right_LiG_lingual_gyrus_volM2 Left_LiG_lingual_gyrus_volM2 Right_LOrG_lateral_orbital_gyrus Left_LOrG_lateral_orbital_gyrus_ Right_MCgG_middle_cingulate_gyru Left_MCgG_middle_cingulate_gyrus Right_MFC_medial_frontal_cortex_ Left_MFC_medial_frontal_cortex_v Right_MFG_middle_frontal_gyrus_v Left_MFG_middle_frontal_gyrus_vo Right_MOG_middle_occipital_gyrus Left_MOG_middle_occipital_gyrus_ Right_MOrG_medial_orbital_gyrus_ Left_MOrG_medial_orbital_gyrus_v Right_MPoG_postcentral_gyrus_med Left_MPoG_postcentral_gyrus_medi Right_MPrG_precentral_gyrus_medi Left_MPrG_precentral_gyrus_media Right_MSFG_superior_frontal_gyru Left_MSFG_superior_frontal_gyrus Right_MTG_middle_temporal_gyrus_ Left_MTG_middle_temporal_gyrus_v Right_OCP_occipital_pole_volM2 Left_OCP_occipital_pole_volM2 Right_OFuG_occipital_fusiform_gy Left_OFuG_occipital_fusiform_gyr Right_OpIFG_opercular_part_of_th Left_OpIFG_opercular_part_of_the Right_OrIFG_orbital_part_of_the_ Left_OrIFG_orbital_part_of_the_i Right_PCgG_posterior_cingulate_g Left_PCgG_posterior_cingulate_gy Right_PCu_precuneus_volM2 Left_PCu_precuneus_volM2 Right_PHG_parahippocampal_gyrus_ Left_PHG_parahippocampal_gyrus_v Right_PIns_posterior_insula_volM Left_PIns_posterior_insula_volM2 Right_PO_parietal_operculum_volM Left_PO_parietal_operculum_volM2 Right_PoG_postcentral_gyrus_volM Left_PoG_postcentral_gyrus_volM2 Right_POrG_posterior_orbital_gyr Left_POrG_posterior_orbital_gyru Right_PP_planum_polare_volM2 Left_PP_planum_polare_volM2 Right_PrG_precentral_gyrus_volM2 Left_PrG_precentral_gyrus_volM2 Right_PT_planum_temporale_volM2 Left_PT_planum_temporale_volM2 Right_SCA_subcallosal_area_volM2 Left_SCA_subcallosal_area_volM2 Right_SFG_superior_frontal_gyrus Left_SFG_superior_frontal_gyrus_ Right_SMC_supplementary_motor_co Left_SMC_supplementary_motor_cor Right_SMG_supramarginal_gyrus_vo Left_SMG_supramarginal_gyrus_vol Right_SOG_superior_occipital_gyr Left_SOG_superior_occipital_gyru Right_SPL_superior_parietal_lobu Left_SPL_superior_parietal_lobul Right_STG_superior_temporal_gyru Left_STG_superior_temporal_gyrus Right_TMP_temporal_pole_volM2 Left_TMP_temporal_pole_volM2 Right_TrIFG_triangular_part_of_t Left_TrIFG_triangular_part_of_th Right_TTG_transverse_temporal_gy Left_TTG_transverse_temporal_gyr {
renvars `X', postdrop(1) 
}


foreach X2 of varlist _3rd_Ventricle_volM _4th_Ventricle_volM Right_Accumbens_Area_volM Left_Accumbens_Area_volM Right_Amygdala_volM Left_Amygdala_volM Brain_Stem_volM Right_Caudate_volM Left_Caudate_volM Right_Cerebellum_Exterior_volM Left_Cerebellum_Exterior_volM Right_Cerebellum_White_Matter_v Left_Cerebellum_White_Matter_vo Right_Hippocampus_volM Left_Hippocampus_volM Right_Inf_Lat_Vent_volM Left_Inf_Lat_Vent_volM Right_Lateral_Ventricle_volM Left_Lateral_Ventricle_volM Right_Pallidum_volM Left_Pallidum_volM Right_Putamen_volM Left_Putamen_volM Right_Thalamus_Proper_volM Left_Thalamus_Proper_volM Right_Ventral_DC_volM Left_Ventral_DC_volM Cerebellar_Vermal_Lobules_I_V_v Cerebellar_Vermal_Lobules_VI_VI Cerebellar_Vermal_Lobules_VIII_ Left_Basal_Forebrain_volM Right_Basal_Forebrain_volM frontal_lobe_WM_right_volM frontal_lobe_WM_left_volM occipital_lobe_WM_right_volM occipital_lobe_WM_left_volM parietal_lobe_WM_right_volM parietal_lobe_WM_left_volM temporal_lobe_WM_right_volM temporal_lobe_WM_left_volM fornix_right_volM fornix_left_volM ALIC_rightv ALIC_leftv PLIC_rightv PLIC_leftv corpus_callosum_volM Right_ACgG_anterior_cingulate_g Left_ACgG_anterior_cingulate_gy Right_AIns_anterior_insula_volM Left_AIns_anterior_insula_volM Right_AOrG_anterior_orbital_gyr Left_AOrG_anterior_orbital_gyru Right_AnG_angular_gyrus_volM Left_AnG_angular_gyrus_volM Right_Calc_calcarine_cortex_vol Left_Calc_calcarine_cortex_volM Right_CO_central_operculum_volM Left_CO_central_operculum_volM Right_Cun_cuneus_volM Left_Cun_cuneus_volM Right_Ent_entorhinal_area_volM Left_Ent_entorhinal_area_volM Right_FO_frontal_operculum_volM Left_FO_frontal_operculum_volM Right_FRP_frontal_pole_volM Left_FRP_frontal_pole_volM Right_FuG_fusiform_gyrus_volM Left_FuG_fusiform_gyrus_volM Right_GRe_gyrus_rectus_volM Left_GRe_gyrus_rectus_volM Right_IOG_inferior_occipital_gy Left_IOG_inferior_occipital_gyr Right_ITG_inferior_temporal_gyr Left_ITG_inferior_temporal_gyru Right_LiG_lingual_gyrus_volM Left_LiG_lingual_gyrus_volM Right_LOrG_lateral_orbital_gyru Left_LOrG_lateral_orbital_gyrus Right_MCgG_middle_cingulate_gyr Left_MCgG_middle_cingulate_gyru Right_MFC_medial_frontal_cortex Left_MFC_medial_frontal_cortex_ Right_MFG_middle_frontal_gyrus_ Left_MFG_middle_frontal_gyrus_v Right_MOG_middle_occipital_gyru Left_MOG_middle_occipital_gyrus Right_MOrG_medial_orbital_gyrus Left_MOrG_medial_orbital_gyrus_ Right_MPoG_postcentral_gyrus_me Left_MPoG_postcentral_gyrus_med Right_MPrG_precentral_gyrus_med Left_MPrG_precentral_gyrus_medi Right_MSFG_superior_frontal_gyr Left_MSFG_superior_frontal_gyru Right_MTG_middle_temporal_gyrus Left_MTG_middle_temporal_gyrus_ Right_OCP_occipital_pole_volM Left_OCP_occipital_pole_volM Right_OFuG_occipital_fusiform_g Left_OFuG_occipital_fusiform_gy Right_OpIFG_opercular_part_of_t Left_OpIFG_opercular_part_of_th Right_OrIFG_orbital_part_of_the Left_OrIFG_orbital_part_of_the_ Right_PCgG_posterior_cingulate_ Left_PCgG_posterior_cingulate_g Right_PCu_precuneus_volM Left_PCu_precuneus_volM Right_PHG_parahippocampal_gyrus Left_PHG_parahippocampal_gyrus_ Right_PIns_posterior_insula_vol Left_PIns_posterior_insula_volM Right_PO_parietal_operculum_vol Left_PO_parietal_operculum_volM Right_PoG_postcentral_gyrus_vol Left_PoG_postcentral_gyrus_volM Right_POrG_posterior_orbital_gy Left_POrG_posterior_orbital_gyr Right_PP_planum_polare_volM Left_PP_planum_polare_volM Right_PrG_precentral_gyrus_volM Left_PrG_precentral_gyrus_volM Right_PT_planum_temporale_volM Left_PT_planum_temporale_volM Right_SCA_subcallosal_area_volM Left_SCA_subcallosal_area_volM Right_SFG_superior_frontal_gyru Left_SFG_superior_frontal_gyrus Right_SMC_supplementary_motor_c Left_SMC_supplementary_motor_co Right_SMG_supramarginal_gyrus_v Left_SMG_supramarginal_gyrus_vo Right_SOG_superior_occipital_gy Left_SOG_superior_occipital_gyr Right_SPL_superior_parietal_lob Left_SPL_superior_parietal_lobu Right_STG_superior_temporal_gyr Left_STG_superior_temporal_gyru Right_TMP_temporal_pole_volM Left_TMP_temporal_pole_volM Right_TrIFG_triangular_part_of_ Left_TrIFG_triangular_part_of_t Right_TTG_transverse_temporal_g Left_TTG_transverse_temporal_gy {
capture drop z`X2' 
	egen z`X2'=std(`X2') if sample_final==1, by(Sex) 
}    

save, replace

capture drop zLnNFLw1
egen zLnNFLw1=std(LnNFLw1) if sample_final==1, by(Sex)

capture drop zLnNFLw3
egen zLnNFLw3=std(LnNFLw3) if sample_final==1,by(Sex)

save, replace

global tflist ""
global modseq=0
foreach X of var zLnNFLw1 zLnNFLw3 {
foreach Y of var z_3rd_Ventricle_volM- zLeft_TTG_transverse_temporal_gy {
foreach Z of var Sex {	
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X'  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1, beta", by(`Z')  label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
 }
 
drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save zOutputdata_overall_Dbysociodem, replace



**ANALYSIS E: Small ROIs, WML**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

capture rename anterior_limb_of_internal_capsu1  ALIC_rightwml
capture rename anterior_limb_of_internal_capsu2 ALIC_leftwml
capture rename posterior_limb_of_internal_caps1 PLIC_rightwml
capture rename posterior_limb_of_internal_caps2 PLIC_leftwml

save HANDLS_paper51_NFLBRAINSCANFINALIZED_finBYSEX, replace


foreach X of varlist _3rd_Ventricle_wmlM2 _4th_Ventricle_wmlM2 Right_Accumbens_Area_wmlM2 Left_Accumbens_Area_wmlM2 Right_Amygdala_wmlM2 Left_Amygdala_wmlM2 Brain_Stem_wmlM2 Right_Caudate_wmlM2 Left_Caudate_wmlM2 Right_Cerebellum_Exterior_wmlM2 Left_Cerebellum_Exterior_wmlM2 Right_Cerebellum_White_Matter_wm Left_Cerebellum_White_Matter_wml Right_Hippocampus_wmlM2 Left_Hippocampus_wmlM2 Right_Inf_Lat_Vent_wmlM2 Left_Inf_Lat_Vent_wmlM2 Right_Lateral_Ventricle_wmlM2 Left_Lateral_Ventricle_wmlM2 Right_Pallidum_wmlM2 Left_Pallidum_wmlM2 Right_Putamen_wmlM2 Left_Putamen_wmlM2 Right_Thalamus_Proper_wmlM2 Left_Thalamus_Proper_wmlM2 Right_Ventral_DC_wmlM2 Left_Ventral_DC_wmlM2 Cerebellar_Vermal_Lobules_I_V_wm Cerebellar_Vermal_Lobules_VI_VI0 Cerebellar_Vermal_Lobules_VIII_0 Left_Basal_Forebrain_wmlM2 Right_Basal_Forebrain_wmlM2 frontal_lobe_WM_right_wmlM2 frontal_lobe_WM_left_wmlM2 occipital_lobe_WM_right_wmlM2 occipital_lobe_WM_left_wmlM2 parietal_lobe_WM_right_wmlM2 parietal_lobe_WM_left_wmlM2 temporal_lobe_WM_right_wmlM2 temporal_lobe_WM_left_wmlM2 fornix_right_wmlM2 fornix_left_wmlM2 ALIC_rightwml ALIC_leftwml PLIC_rightwml PLIC_leftwml corpus_callosum_wmlM2 Right_ACgG_anterior_cingulate_g0 Left_ACgG_anterior_cingulate_gy0 Right_AIns_anterior_insula_wmlM2 Left_AIns_anterior_insula_wmlM2 Right_AOrG_anterior_orbital_gyr0 Left_AOrG_anterior_orbital_gyru0 Right_AnG_angular_gyrus_wmlM2 Left_AnG_angular_gyrus_wmlM2 Right_Calc_calcarine_cortex_wmlM Left_Calc_calcarine_cortex_wmlM2 Right_CO_central_operculum_wmlM2 Left_CO_central_operculum_wmlM2 Right_Cun_cuneus_wmlM2 Left_Cun_cuneus_wmlM2 Right_Ent_entorhinal_area_wmlM2 Left_Ent_entorhinal_area_wmlM2 Right_FO_frontal_operculum_wmlM2 Left_FO_frontal_operculum_wmlM2 Right_FRP_frontal_pole_wmlM2 Left_FRP_frontal_pole_wmlM2 Right_FuG_fusiform_gyrus_wmlM2 Left_FuG_fusiform_gyrus_wmlM2 Right_GRe_gyrus_rectus_wmlM2 Left_GRe_gyrus_rectus_wmlM2 Right_IOG_inferior_occipital_gy0 Left_IOG_inferior_occipital_gyr0 Right_ITG_inferior_temporal_gyr0 Left_ITG_inferior_temporal_gyru0 Right_LiG_lingual_gyrus_wmlM2 Left_LiG_lingual_gyrus_wmlM2 Right_LOrG_lateral_orbital_gyru0 Left_LOrG_lateral_orbital_gyrus0 Right_MCgG_middle_cingulate_gyr0 Left_MCgG_middle_cingulate_gyru0 Right_MFC_medial_frontal_cortex0 Left_MFC_medial_frontal_cortex_w Right_MFG_middle_frontal_gyrus_w Left_MFG_middle_frontal_gyrus_wm Right_MOG_middle_occipital_gyru0 Left_MOG_middle_occipital_gyrus0 Right_MOrG_medial_orbital_gyrus0 Left_MOrG_medial_orbital_gyrus_w Right_MPoG_postcentral_gyrus_me0 Left_MPoG_postcentral_gyrus_med0 Right_MPrG_precentral_gyrus_med0 Left_MPrG_precentral_gyrus_medi0 Right_MSFG_superior_frontal_gyr0 Left_MSFG_superior_frontal_gyru0 Right_MTG_middle_temporal_gyrus0 Left_MTG_middle_temporal_gyrus_w Right_OCP_occipital_pole_wmlM2 Left_OCP_occipital_pole_wmlM2 Right_OFuG_occipital_fusiform_g0 Left_OFuG_occipital_fusiform_gy0 Right_OpIFG_opercular_part_of_t0 Left_OpIFG_opercular_part_of_th0 Right_OrIFG_orbital_part_of_the0 Left_OrIFG_orbital_part_of_the_0 Right_PCgG_posterior_cingulate_0 Left_PCgG_posterior_cingulate_g0 Right_PCu_precuneus_wmlM2 Left_PCu_precuneus_wmlM2 Right_PHG_parahippocampal_gyrus0 Left_PHG_parahippocampal_gyrus_w Right_PIns_posterior_insula_wmlM Left_PIns_posterior_insula_wmlM2 Right_PO_parietal_operculum_wmlM Left_PO_parietal_operculum_wmlM2 Right_PoG_postcentral_gyrus_wmlM Left_PoG_postcentral_gyrus_wmlM2 Right_POrG_posterior_orbital_gy0 Left_POrG_posterior_orbital_gyr0 Right_PP_planum_polare_wmlM2 Left_PP_planum_polare_wmlM2 Right_PrG_precentral_gyrus_wmlM2 Left_PrG_precentral_gyrus_wmlM2 Right_PT_planum_temporale_wmlM2 Left_PT_planum_temporale_wmlM2 Right_SCA_subcallosal_area_wmlM2 Left_SCA_subcallosal_area_wmlM2 Right_SFG_superior_frontal_gyru0 Left_SFG_superior_frontal_gyrus0 Right_SMC_supplementary_motor_c0 Left_SMC_supplementary_motor_co0 Right_SMG_supramarginal_gyrus_wm Left_SMG_supramarginal_gyrus_wml Right_SOG_superior_occipital_gy0 Left_SOG_superior_occipital_gyr0 Right_SPL_superior_parietal_lob0 Left_SPL_superior_parietal_lobu0 Right_STG_superior_temporal_gyr0 Left_STG_superior_temporal_gyru0 Right_TMP_temporal_pole_wmlM2 Left_TMP_temporal_pole_wmlM2 Right_TrIFG_triangular_part_of_0 Left_TrIFG_triangular_part_of_t0 Right_TTG_transverse_temporal_g0 Left_TTG_transverse_temporal_gy0 {
renvars `X', postdrop(2) 
}




foreach X1 of varlist _3rd_Ventricle_wml _4th_Ventricle_wml Right_Accumbens_Area_wml Left_Accumbens_Area_wml Right_Amygdala_wml Left_Amygdala_wml Brain_Stem_wml Right_Caudate_wml Left_Caudate_wml Right_Cerebellum_Exterior_wml Left_Cerebellum_Exterior_wml Right_Cerebellum_White_Matter_ Left_Cerebellum_White_Matter_w Right_Hippocampus_wml Left_Hippocampus_wml Right_Inf_Lat_Vent_wml Left_Inf_Lat_Vent_wml Right_Lateral_Ventricle_wml Left_Lateral_Ventricle_wml Right_Pallidum_wml Left_Pallidum_wml Right_Putamen_wml Left_Putamen_wml Right_Thalamus_Proper_wml Left_Thalamus_Proper_wml Right_Ventral_DC_wml Left_Ventral_DC_wml Cerebellar_Vermal_Lobules_I_V_ Cerebellar_Vermal_Lobules_VI_V Cerebellar_Vermal_Lobules_VIII Left_Basal_Forebrain_wml Right_Basal_Forebrain_wml frontal_lobe_WM_right_wml frontal_lobe_WM_left_wml occipital_lobe_WM_right_wml occipital_lobe_WM_left_wml parietal_lobe_WM_right_wml parietal_lobe_WM_left_wml temporal_lobe_WM_right_wml temporal_lobe_WM_left_wml fornix_right_wml fornix_left_wml ALIC_rightw ALIC_leftw PLIC_rightw PLIC_leftw corpus_callosum_wml Right_ACgG_anterior_cingulate_ Left_ACgG_anterior_cingulate_g Right_AIns_anterior_insula_wml Left_AIns_anterior_insula_wml Right_AOrG_anterior_orbital_gy Left_AOrG_anterior_orbital_gyr Right_AnG_angular_gyrus_wml Left_AnG_angular_gyrus_wml Right_Calc_calcarine_cortex_wm Left_Calc_calcarine_cortex_wml Right_CO_central_operculum_wml Left_CO_central_operculum_wml Right_Cun_cuneus_wml Left_Cun_cuneus_wml Right_Ent_entorhinal_area_wml Left_Ent_entorhinal_area_wml Right_FO_frontal_operculum_wml Left_FO_frontal_operculum_wml Right_FRP_frontal_pole_wml Left_FRP_frontal_pole_wml Right_FuG_fusiform_gyrus_wml Left_FuG_fusiform_gyrus_wml Right_GRe_gyrus_rectus_wml Left_GRe_gyrus_rectus_wml Right_IOG_inferior_occipital_g Left_IOG_inferior_occipital_gy Right_ITG_inferior_temporal_gy Left_ITG_inferior_temporal_gyr Right_LiG_lingual_gyrus_wml Left_LiG_lingual_gyrus_wml Right_LOrG_lateral_orbital_gyr Left_LOrG_lateral_orbital_gyru Right_MCgG_middle_cingulate_gy Left_MCgG_middle_cingulate_gyr Right_MFC_medial_frontal_corte Left_MFC_medial_frontal_cortex Right_MFG_middle_frontal_gyrus Left_MFG_middle_frontal_gyrus_ Right_MOG_middle_occipital_gyr Left_MOG_middle_occipital_gyru Right_MOrG_medial_orbital_gyru Left_MOrG_medial_orbital_gyrus Right_MPoG_postcentral_gyrus_m Left_MPoG_postcentral_gyrus_me Right_MPrG_precentral_gyrus_me Left_MPrG_precentral_gyrus_med Right_MSFG_superior_frontal_gy Left_MSFG_superior_frontal_gyr Right_MTG_middle_temporal_gyru Left_MTG_middle_temporal_gyrus Right_OCP_occipital_pole_wml Left_OCP_occipital_pole_wml Right_OFuG_occipital_fusiform_ Left_OFuG_occipital_fusiform_g Right_OpIFG_opercular_part_of_ Left_OpIFG_opercular_part_of_t Right_OrIFG_orbital_part_of_th Left_OrIFG_orbital_part_of_the Right_PCgG_posterior_cingulate Left_PCgG_posterior_cingulate_ Right_PCu_precuneus_wml Left_PCu_precuneus_wml Right_PHG_parahippocampal_gyru Left_PHG_parahippocampal_gyrus Right_PIns_posterior_insula_wm Left_PIns_posterior_insula_wml Right_PO_parietal_operculum_wm Left_PO_parietal_operculum_wml Right_PoG_postcentral_gyrus_wm Left_PoG_postcentral_gyrus_wml Right_POrG_posterior_orbital_g Left_POrG_posterior_orbital_gy Right_PP_planum_polare_wml Left_PP_planum_polare_wml Right_PrG_precentral_gyrus_wml Left_PrG_precentral_gyrus_wml Right_PT_planum_temporale_wml Left_PT_planum_temporale_wml Right_SCA_subcallosal_area_wml Left_SCA_subcallosal_area_wml Right_SFG_superior_frontal_gyr Left_SFG_superior_frontal_gyru Right_SMC_supplementary_motor_ Left_SMC_supplementary_motor_c Right_SMG_supramarginal_gyrus_ Left_SMG_supramarginal_gyrus_w Right_SOG_superior_occipital_g Left_SOG_superior_occipital_gy Right_SPL_superior_parietal_lo Left_SPL_superior_parietal_lob Right_STG_superior_temporal_gy Left_STG_superior_temporal_gyr Right_TMP_temporal_pole_wml Left_TMP_temporal_pole_wml Right_TrIFG_triangular_part_of Left_TrIFG_triangular_part_of_ Right_TTG_transverse_temporal_ Left_TTG_transverse_temporal_g {
capture drop L`X1' 
	gen L`X1'=ln(`X1')^1/3 if sample_final==1
	
}



**Select non-missing from z_3rd_Ventricle_wmlM through zLeft_TTG_transverse_temporal_gy**


foreach X2 of varlist LRight_Caudate_wml  LLeft_Caudate_wml LRight_Hippocampus_wml LLeft_Hippocampus_wml  LRight_Lateral_Ventricle_wml LLeft_Lateral_Ventricle_wml  LRight_Thalamus_Proper_wml LLeft_Thalamus_Proper_wml Lfrontal_lobe_WM_right_wml Lfrontal_lobe_WM_left_wml Loccipital_lobe_WM_right_wml Loccipital_lobe_WM_left_wml Lparietal_lobe_WM_right_wml Lparietal_lobe_WM_left_wml Ltemporal_lobe_WM_right_wml Ltemporal_lobe_WM_left_wml  Lfornix_left_wml LALIC_rightw LALIC_leftw LPLIC_rightw LPLIC_leftw Lcorpus_callosum_wml   LRight_Calc_calcarine_cortex_wm LLeft_Calc_calcarine_cortex_wml  { 
capture drop z`X2' 
	egen z`X2'=std(`X2') if sample_final==1, by(Sex)
	
}

save HANDLS_paper51_NFLBRAINSCANFINALIZED_finBYSEX, replace

capture drop zLnNFLw1
egen zLnNFLw1=std(LnNFLw1) if sample_final==1

capture drop zLnNFLw3
egen zLnNFLw3=std(LnNFLw3) if sample_final==1

save, replace



global tflist ""
global modseq=0
foreach X of var zLnNFLw1 zLnNFLw3 {
foreach Y of var zLRight_Caudate_wml  zLLeft_Caudate_wml zLRight_Hippocampus_wml zLLeft_Hippocampus_wml  zLRight_Lateral_Ventricle_wml zLLeft_Lateral_Ventricle_wml  zLRight_Thalamus_Proper_wml zLLeft_Thalamus_Proper_wml zLfrontal_lobe_WM_right_wml zLfrontal_lobe_WM_left_wml zLoccipital_lobe_WM_right_wml zLoccipital_lobe_WM_left_wml zLparietal_lobe_WM_right_wml zLparietal_lobe_WM_left_wml zLtemporal_lobe_WM_right_wml zLtemporal_lobe_WM_left_wml  zLfornix_left_wml zLALIC_rightw zLALIC_leftw zLPLIC_rightw zLPLIC_leftw zLcorpus_callosum_wml   zLRight_Calc_calcarine_cortex_wm zLLeft_Calc_calcarine_cortex_wml {
foreach Z of var Sex {
global modseq=$modseq+1
tempfile tfcur
parmby "regr `Y' `X'  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2  if sample_final==1, beta", by(`Z')  label command format(estimate min95 max95 %8.2f p %8.1e) idn($modseq) saving(`"`tfcur'"', replace) flist(tflist)
 }
 }
}
capture drop _all
append using $tflist
sort idnum command `Y' parmseq
describe
by idnum command:list `Y' parm label estimate min95 max95 p, noobs

save zOutputdata_overall_Ebysociodem, replace



//Manhattan plots@@//

**    manhattan chromosome base-pair pvalue [if] [, options]

**    options               Description
    **---------------------------------------------------------------------------------------------------------------------------------------------------------**-----------------------------------
**    Plot options
**      title(string)       display a title
**      caption(string)     display a caption
**      xlabel(string)      set x label; default is xlabel(Chromosome)
**      width(#)            set width of plot; default is width(15)
**      height(#)           set height of plot; default is height(5)

**    Chromosome options
**      x(#)                specify chromosome number to be labeled as X; default is x(23)
**      y(#)                specify chromosome number to be labeled as Y; default is y(24)
**      mito(#)             specify chromosome number to be labeled as M; default is mito(25)

**    Graph options
**      bonferroni(h|v|n)   draw a line at Bonferroni significance level; label line with horizontal (h), vertical (v), or no (n) labels
**      mlabel(var)         set a variable to use for labeling markers
**      mthreshold(#|b)     set a -log(p-value) above which markers will be labeled, or use b to set your threshold at the Bonferroni significance level
**      yline(#)            set log(p-value) at which to draw a line
**      labelyline(h|v)     label line specified with yline() by using horizontal (h) or vertical (v) labels
**     addmargin           add a margin to the left and right of the plot, leaving room for labels

**    Style options
**      color1(color)       set first color of markers
**      color2(color)       set second color of markers
**      linecolor(color)    set the color of Bonferroni line and label or y line and label
**    **---------------------------------------------------------------------------------------------------------------------------------------------------------**-----------------------------------


capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLE2_3OUTPUTS.smcl", replace


********************************ANALYSIS A*******************************************


//Total population//

use Outputdata_overall_A, clear

keep if parmseq==1


manhattan idnum parmseq  p , bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.26) labelyline(h) 

graph save manhattan_overall_A.gph,replace

**Smile plot**
use Outputdata_overall_A, clear

keep if parmseq==1


sort parm
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp) ptl(command)  t1(" ")
list if signif, nodisp



//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)
sort myqvallargervolumes

list parm myqvallargervolumes command p  estimate stderr if p<0.05


////////////////////////By Sex/////////////////////////////
	
use Outputdata_bysociodem_A, clear

keep if dof>=20 & parmseq==1

manhattan idnum parmseq  p  ,  bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.56)  labelyline(h)  

graph save manhattanbysociodem_A.gph,replace


save, replace

**Smile plot**

//By Sex//
use Outputdata_bysociodem_A, clear

keep if parmseq==1

sort Sex
capture drop uncp corp signif
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp)by(Sex) ptl(command)  t1(" ")
list if signif, nodisp




//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command estimate stderr p Sex if p<0.05


save, replace

**--> No significant findings or q-values**



***********************************ANALYSIS Aprime*******************************


//Total population//

use Outputdata_overall_Aprime, clear

keep if parmseq==1


manhattan idnum parmseq  p , bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.26) labelyline(h) 

graph save manhattan_Aprime.gph, replace

**Smile plot**
use Outputdata_overall_Aprime, clear

keep if parmseq==1


sort parm
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp) ptl(command)  t1(" ")
list if signif, nodisp


//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command p  estimate stderr if p<0.05


//By Sex//
	
use Outputdata_bysociodem_Aprime, clear

keep if dof>=20 & parmseq==1

manhattan idnum parmseq  p  ,  bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.56)  labelyline(h)  

graph save manhattan_bysociodem_Aprime.gph, replace


save, replace

**Smile plot**

//By Sex//
use Outputdata_bysociodem_Aprime, clear

keep if parmseq==1

sort Sex
capture drop uncp corp signif
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp)by(Sex) ptl(command)  t1(" ")
list if signif, nodisp




//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command estimate stderr p Sex if p<0.05


save, replace

**--> No significant findings or q-values**



**************************************ANALYSIS B**********************************



//Total population//

use Outputdata_overall_B, clear

keep if parmseq==1


manhattan idnum parmseq  p , bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.96) labelyline(h) 

graph save manhattan_overall_B.gph,replace

**Smile plot**
use Outputdata_overall_B, clear

keep if parmseq==1


sort parm
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp) ptl(command)  t1(" ")
list if signif, nodisp


**     +--------------------------------------------------------------------------------------------
**> ------------------------------------------------------------------------------------------------*
**> ----------------------------------+
**     | parmseq   command                                                                          
**>            idnum      parm   label   estimate      stderr   dof            t         p     min95
**>     max95   uncp    corp   signif |
**     |--------------------------------------------------------------------------------------------
**> ------------------------------------------------------------------------------------------------
**> ----------------------------------|
**  3. |       1   regr Left_Hippocampus LnNFLw3  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sa
**> mple_f..       3   LnNFLw3            -124.98   42.412936   171   -2.9467604   3.7e-03   -208.70
**>    -41.26    .05   .0125        1 |
**     +--------------------------------------------------------------------------------------------
**> ------------------------------------------------------------------------------------------------
**> ----------------------------------+

. 



//q-value: FDR//

use Outputdata_overall_B, clear

keep if parmseq==1


capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)

sort myqvallargervolumes
list parm myqvallargervolumes command estimate stderr p  if p<0.05


**     +-------------------------------------------------------------------------------------------+
**  1. |                     parm                   |                  myqvall~s                   |
**     |                  LnNFLw3                   |                  .01463801                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr Left_Hippocampus LnNFLw3  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_f.. |
**     |-------------------------------------------------------------------------------------------|
**     |          estimate           |              stderr           |                 p           |
**     |           -124.98           |           42.412936           |           3.7e-03           |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  2. |                     parm                   |                  myqvall~s                   |
**     |                  LnNFLw3                   |                  .04868519                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr Right_Hippocampus LnNFLw3  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_.. |
**     |-------------------------------------------------------------------------------------------|
**     |          estimate           |              stderr           |                 p           |
**     |           -100.45           |           44.214994           |           2.4e-02           |
**     +---------------------------------------------------------------------------------------


//By Sex//  ///
	
use Outputdata_bysociodem_B, clear

keep if dof>=20 & parmseq==1

manhattan idnum parmseq  p  ,  bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.56)  labelyline(h)  


graph save manhattan_bysociodem_B.gph,replace


save, replace

**Smile plot**

//By Sex//
use Outputdata_bysociodem_B, clear

keep if parmseq==1

sort Sex
capture drop uncp corp signif
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp)by(Sex) ptl(command)  t1(" ")
list if signif, nodisp


//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command estimate stderr p Sex if p<0.05



**     +-------------------------------------------------------------------------------------------+
**  1. |                     parm                   |                  myqvall~s                   |
**     |                  LnNFLw3                   |                  .19341837                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_fi.. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate        |          stderr        |             p        |        Sex        |
**     |        -145.36        |       63.140141        |       2.4e-02        |        Men        |
**     +-------------------------------------------------------------------------------------------+



****************************ANALYSIS C********************************


//Total population//

use Outputdata_overall_C, clear

keep if parmseq==1


manhattan idnum parmseq  p , bonferroni (v) xlabel (idnum/parmseq) mlabel (command) labelyline(h) 

graph save manhattan_overall_C.gph, replace

**Smile plot**
use Outputdata_overall_C, clear

keep if parmseq==1


sort parm
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp) ptl(command)  t1(" ")
list if signif, nodisp


//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command p  estimate stderr if p<0.05


**     +-------------------------------------------------------------------------------------------+
**  1. |                     parm                   |                  myqvall~s                   |
**     |                  LnNFLw1                   |                  .00299568                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr LnLesion_Volume LnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_fi.. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          1.5e-03           |               2.13           |           .66041491           |
**     +-------------------------------------------------------------------------------------------+



//By Sex//
	
use Outputdata_bysociodem_C, clear

keep if dof>=20 & parmseq==1

manhattan idnum parmseq  p  ,  bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.56)  labelyline(h)  

graph save manhattan_bysociodem_C.gph, replace

save, replace

**Smile plot**

//By Sex//
use Outputdata_bysociodem_C, clear

keep if parmseq==1

sort Sex
capture drop uncp corp signif
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp)by(Sex) ptl(command)  t1(" ")
list if signif, nodisp


//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command estimate stderr p Sex if p<0.05


**     +-------------------------------------------------------------------------------------------+
**  1. |                     parm                   |                  myqvall~s                   |
**     |                  LnNFLw1                   |                  .05278914                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_fin.. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate       |          stderr        |             p        |         Sex        |
**     |           2.90       |       1.1459241        |       1.3e-02        |       Women        |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  2. |                     parm                   |                  myqvall~s                   |
**     |                  LnNFLw1                   |                   .0804528                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_fin.. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate       |          stderr        |             p        |         Sex        |
**     |           1.48       |       .70685723        |       4.0e-02        |         Men        |
**     +-----------------------------------------------------------------------------------------





save, replace


****************************ANALYSIS D********************************


//Total population//

use zOutputdata_overall_D, clear

keep if parmseq==1


manhattan idnum parmseq  p , bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.26) labelyline(h) 

graph save manhattan_overall_D.gph, replace

**Smile plot**
use zOutputdata_overall_D, clear

keep if parmseq==1


sort parm
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp) ptl(command)  t1(" ")
list if signif, nodisp


//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command p  estimate stderr if p<0.05

**no significant findings**

//By Sex//
	
use zOutputdata_overall_Dbysociodem, clear

keep if dof>=20 & parmseq==1

manhattan idnum parmseq  p  ,  bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.96)  labelyline(h)  

graph save manhattan_Dbysociodem.gph, replace

save, replace

**Smile plot**

//By Sex//
use zOutputdata_overall_Dbysociodem, clear

keep if parmseq==1

sort Sex
capture drop uncp corp signif
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp)by(Sex) ptl(command)  t1(" ")
list if signif, nodisp


//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command estimate stderr p Sex if p<0.05

**no significant findings**


save, replace



****************************ANALYSIS E********************************


//Total population//

use zOutputdata_overall_E, clear

keep if parmseq==1


manhattan idnum parmseq  p , bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.26) labelyline(h) 

graph save manhattan_overall_E.gph, replace

**Smile plot**
use zOutputdata_overall_E, clear

keep if parmseq==1


sort parm
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp) ptl(command)  t1(" ")
list if signif, nodisp


//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command p  estimate stderr if p<0.05


**     +-------------------------------------------------------------------------------------------+
**  1. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                   .0014513                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLPLIC_rightw zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2  if sample_fi.. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          3.0e-05           |               0.38           |           .08937904           |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  2. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                  .06088352                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLfrontal_lobe_WM_right_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2.. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          2.5e-03           |               0.26           |           .08622859           |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  3. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                  .09920288                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLfrontal_lobe_WM_left_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 .. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          6.4e-03           |               0.24           |           .08636262           |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  4. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                  .09920288                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLRight_Lateral_Ventricle_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_vol.. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          8.3e-03           |               0.24           |           .09063563           |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  5. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                   .1016253                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLparietal_lobe_WM_right_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM.. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          1.3e-02           |               0.23           |            .0905931           |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  6. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                   .1016253                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLLeft_Lateral_Ventricle_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM.. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          1.1e-02           |               0.23           |           .09085247           |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  7. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                  .11051753                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLRight_Caudate_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2  if sam.. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          1.6e-02           |               0.22           |           .09204836           |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  8. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                  .16535894                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLcorpus_callosum_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2  if s.. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          2.8e-02           |               0.21           |           .09293203           |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  9. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                  .19580818                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLALIC_rightw zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2  if sample_fi.. |
**     |-------------------------------------------------------------------------------------------|
**     |                p           |           estimate           |              stderr           |
**     |          3.7e-02           |               0.19           |           .09169732           |
**     +-------------------------------------------------------------------------------------------+

**



//By Sex//
	
use zOutputdata_overall_Ebysociodem, clear

keep if dof>=20 & parmseq==1

manhattan idnum parmseq  p  ,  bonferroni (v) xlabel (idnum/parmseq) mlabel (command) mthreshold(1.96)  labelyline(h)  

graph save manhattan_Ebysociodem.gph, replace

save, replace

**Smile plot**

//By Sex//
use zOutputdata_overall_Ebysociodem, clear

keep if parmseq==1

sort Sex
capture drop uncp corp signif
multproc, pval(p) meth(bonferroni) gpunc(uncp) gpcor(corp) rej(signif)
smileplot7, est(estimate) pval(p) punc(uncp) pcor(corp)by(Sex) ptl(command)  t1(" ")
list if signif, nodisp


//q-value: FDR//

capture drop myqvallargervolumes

qqvalue p, method(simes) qvalue(myqvallargervolumes)


sort myqvallargervolumes
list parm myqvallargervolumes command estimate stderr p Sex if p<0.05

**     +-------------------------------------------------------------------------------------------+
**  1. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw3                  |                   .3365128                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLRight_Lateral_Ventricle_wml zLnNFLw3  Sex w1Age Race PovStat TIME_V1SCAN ICV_vol.. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate       |          stderr        |             p        |         Sex        |
**     |           0.37       |       .16857698        |       2.9e-02        |       Women        |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  2. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                   .3365128                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLfrontal_lobe_WM_left_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 .. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate       |          stderr        |             p        |         Sex        |
**     |           0.35       |       .14066642        |       1.5e-02        |         Men        |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  3. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                   .3365128                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLoccipital_lobe_WM_left_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM.. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate       |          stderr        |             p        |         Sex        |
**     |          -0.35       |       .14739173        |       2.0e-02        |         Men        |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  4. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw3                  |                   .3365128                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLfrontal_lobe_WM_right_wml zLnNFLw3  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2.. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate       |          stderr        |             p        |         Sex        |
**     |          -0.35       |       .14758572        |       2.1e-02        |         Men        |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  5. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw3                  |                   .3365128                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLLeft_Lateral_Ventricle_wml zLnNFLw3  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM.. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate       |          stderr        |             p        |         Sex        |
**     |          -0.25       |       .10899528        |       2.8e-02        |         Men        |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  6. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw1                  |                  .37496042                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLRight_Caudate_wml zLnNFLw1  Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2  if sam.. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate       |          stderr        |             p        |         Sex        |
**     |           0.28       |       .13555283        |       4.5e-02        |         Men        |
**     +-------------------------------------------------------------------------------------------+
**
**     +-------------------------------------------------------------------------------------------+
**  7. |                      parm                  |                  myqvall~s                   |
**     |                  zLnNFLw3                  |                  .37496042                   |
**     |-------------------------------------------------------------------------------------------|
**     | command                                                                                   |
**     | regr zLoccipital_lobe_WM_right_wml zLnNFLw3  Sex w1Age Race PovStat TIME_V1SCAN ICV_vol.. |
**     |-------------------------------------------------------------------------------------------|
**     |       estimate       |          stderr        |             p        |         Sex        |
**     |           0.50       |       .24179921        |       4.5e-02        |       Women        |
**     +-------------------------------------------------------------------------------------------+



save, replace


**************************************ALL FINDINGS***********************************************

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLE2_S3_ALLFINDINGS.smcl",replace


////////////////////////////////////////TABLES 2 AND S3/////////////////////////////////////////////


***********************TABLE 3: LnNFLw1, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta



//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1, beta 



**Model 2: BMI-Adjusted**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2 

save, replace





*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta



//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1 

save, replace



//INTERACTION BY Sex//
use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace


*********************TABLE S3: LnNFLw1, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace


//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==2

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1

save, replace


//INTERACTION BY Sex//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace




***********MODEL 4: MODEL 2+liver/kidney disease******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//

mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 

save, replace

//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1 

save, replace

**INTERACTION BY Sex**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  if sample_final==1 

save, replace

***********MODEL 5: MODEL 2+oxidative stress******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1 

save, replace


**********INTERACTION BY Sex***************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==2



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1 

save, replace

******************INTERACTION BY Sex***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace



capture log close

capture log using "D:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLE3_S4_ALLFINDINGS.smcl", replace



////////////////////////////////////////TABLES 3 AND S4/////////////////////////////////////////////


***********************TABLE 3: LnNFLw3, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1, beta 



**Model 2: BMI-Adjusted**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2 

save, replace





*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1 

save, replace



//INTERACTION BY Sex//
use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace


*********************TABLE S3: LnNFLw3, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace


//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1

save, replace


//INTERACTION BY Sex//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace




***********MODEL 4: MODEL 2+liver/kidney disease******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//

mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 

save, replace

//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1 

save, replace

**INTERACTION BY Sex**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  if sample_final==1 

save, replace

***********MODEL 5: MODEL 2+oxidative stress******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1 

save, replace


**********INTERACTION BY Sex***************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1 

save, replace

******************INTERACTION BY Sex***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////TABLES 4 AND S4-S5/////////////////////////////////////////////

***********************TABLE 3: LnNFLw3, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2 

save, replace


*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1 

save, replace



//INTERACTION BY Sex//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace


capture log close


capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLE4.smcl",replace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

***********************TABLE 4: bayes1LnNFL, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg GM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg WM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta

//ANALYSIS B//
reg Left_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta
reg Right_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta

//ANALYSIS C//
reg LnLesion_Volume bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg GM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final==1
mi estimate: reg WM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg GM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg WM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta

//ANALYSIS B//
reg Left_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta
reg Right_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2
mi estimate: reg GM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==2
mi estimate: reg WM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2

//ANALYSIS B//
mi estimate: reg Left_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2 

save, replace


*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg GM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg WM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta


//ANALYSIS B//
reg Left_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta
reg Right_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1
mi estimate: reg GM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==1
mi estimate: reg WM bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume bayes1LnNFL LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1 

save, replace



//INTERACTION BY Sex//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.bayes1LnNFL##Sex LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg GM c.bayes1LnNFL##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg WM c.bayes1LnNFL##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.bayes1LnNFL##Sex  LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.bayes1LnNFL##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.bayes1LnNFL##Sex LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace




capture log close



capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLES5.smcl",replace


//////////////////////////////////////////////////TABLE S5: NFL at v1 AND v3 vs. Hippocampus and LnLesion volume adjusted for ICV_volM2,  by race///////////////////////////////////////////////////////////////////////////////////
	  
	  ////////////////////////////////////////LnNFLw1 EXPOSURE/////////////////////////////////////////////


***********************LnNFLw1, MODELS 1 AND 2*********************

*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==2,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==2 

save, replace





*********************WHITES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==1,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==1,beta 


**Model 2**

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==1 

save, replace



//INTERACTION BY Race//


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final==1 

save, replace


*********************LnNFLw1, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==2 

save, replace



//Whites//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==1 

save, replace


//INTERACTION by Race//



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==1 

save, replace

**INTERACTION by Race**


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 

save, replace


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1  

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1  

save, replace

***********MODEL 5: MODEL2+OXIDATIVE STRESS******


//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==2 

save, replace



//WHITE//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==1 

save, replace


**********INTERACTION by Race***************



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace



***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  



//AFRICAN-AMERICAN//


use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==1 

save, replace

******************INTERACTION by Race***********************


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace


////////////////////////////////////////TABLE S5. LnNFLw3 exposure/////////////////////////////////////////////




***********************LnNFLw3, MODELS 1 AND 2*********************


*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==2,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==2 

save, replace





*********************WHITE******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final==1 & Race==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Race==1 

save, replace



//INTERACTION by Race//


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final==1 

save, replace


*********************LnNFLw3, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==2 

save, replace



//WHITE//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Race==1 

save, replace


//INTERACTION by Race//



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Race==1 

save, replace

**INTERACTION by Race**



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 

save, replace


***********MODEL 5: MODEL2+OXIDATIVE STRESS******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Race==1 

save, replace


**********INTERACTION by Race***************


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  

//AFRICAN-AMERICAN//

use finaldata_imputed,clear




//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Race==1 

save, replace

******************INTERACTION by Race***********************



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace



  

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLES6.smcl",replace


/////////////////////////////////////////////////TABLE S6: TOTAL AND GM/WM VOLUMES VS. NFL exposures by Race***********************
////////////////////////////////////////LnNFLw1 EXPOSURE/////////////////////////////////////////////


***********************LnNFLw1, MODELS 1 AND 2*********************

*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==2,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==2,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==2,beta


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==2

save, replace

*********************WHITES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==1,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==1,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==1,beta


**Model 2**

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==1

save, replace



//INTERACTION BY Race//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1

save, replace


*********************LnNFLw1, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==2

save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==1


save, replace


//INTERACTION BY Race//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==2

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==1

save, replace

**INTERACTION BY Race**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1


save, replace


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1  

save, replace

***********MODEL 5: MODEL2+OXIDATIVE STRESS******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

save, replace


//AFRICAN-AMERICAN//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==2

mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==2


save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==1

mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==1

save, replace


**********INTERACTION BY Race***************



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

save, replace



***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1

save, replace


//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==2


save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==1


save, replace

******************INTERACTION BY Race***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1

save, replace


////////////////////////////////////////LnNFLw3 exposure/////////////////////////////////////////////


***********************LnNFLw3, MODELS 1 AND 2*********************

*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==2,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==2,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==2,beta


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==2

save, replace

*********************WHITES******************************************************

**Model 1**


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final==1 & Race==1,beta


**Model 2**


use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Race==1

save, replace



//INTERACTION BY Race//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1

save, replace


*********************LnNFLw3, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//


use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==2

save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Race==1


save, replace


//INTERACTION BY Race//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==2

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Race==1

save, replace

**INTERACTION BY Race**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1


save, replace


***********MODEL 5: MODEL2+OXIDATIVE STRESS******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==2

mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==2


save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==1

mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Race==1

save, replace


**********INTERACTION BY Race***************



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

save, replace



***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==2


save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Race==1


save, replace

******************INTERACTION BY Race***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1

save, replace

capture log close



////////////////////////////////////CREATE A CSV FILE FOR GURAY FOR ANALYSIS E////////////////////////////

cd "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\OUTPUTS_PARMBY"

use zOutputdata_overall_E,clear

keep if parmseq==1

save zOutputdata_overall_E_firstparm, replace

sort parm p

save, replace

**--> Change to CSV file******
**export delimited command parm estimate stderr dof p min95 max95 using  **"D:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\MayBaydoun_folder\HANDLS_PAPER51_NFLSMRI\OUTPUT\OUTPUTS_PARMBY\GURAY_PARMBY_FILEANALYSISE_OVERALL.csv",novarnames replace

**--> Manipulate the file and create two CSV files to run for v1 and v2 NFL using R code for volcano plot**

**-->ADD VARIABLES FOR MEAN AND SD FOR KEY WM AND WML VOLUME VARIABLES***


cd "E:\HANDLS_PAPER51_NFLSMRI\DATA"

use HANDLS_paper51_NFLBRAINSCANFINALIZED_fin,clear

su WM if sample_final==1
su Lesion_Volume if sample_final==1

su PLIC_rightv frontal_lobe_WM_right_wml frontal_lobe_WM_left_wml Right_Lateral_Ventricle_wml Left_Lateral_Ventricle_wml parietal_lobe_WM_right_wml Right_Caudate_wml corpus_callosum_wml ALIC_rightv if sample_final==1



**--> Send files to Guray to create the brain image plots**




////////////////////////////////////////////////////////////////SENSITIVITY ANALYSES///////////////////////////////////////////////////////




/////////////////////////////SENSITIVITY ANALYSIS #1: TRACKING HIGH VS. TRACKING LOW:TABLE S7/////////////////////////////////////////////////////////////////////////////

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLES7_ALLFINDINGS.smcl",replace


////////////////////////////////////////TABLES S7/////////////////////////////////////////////

cd "E:\HANDLS_PAPER51_NFLSMRI\DATA"

***********************TABLE 3: NFLw1w3trackhigh, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta



//ANALYSIS B//
reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta
reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta

//ANALYSIS C//
reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1, beta 



**Model 2: BMI-Adjusted**

use finaldata_imputed,clear

capture drop NFLw1w3trackhigh
gen NFLw1w3trackhigh=.
replace NFLw1w3trackhigh=1 if NFLw1>8 & NFLw3>8 & NFLw1~=. & NFLw3~=.
replace NFLw1w3trackhigh=0 if NFLw1w3trackhigh~=1 & NFLw1~=. & NFLw3~=.


capture drop NFLw1w3tracklow
gen NFLw1w3tracklow=.
replace NFLw1w3tracklow=1 if NFLw1<=8 & NFLw3<=8 & NFLw1~=. & NFLw3~=.
replace NFLw1w3tracklow=0 if NFLw1w3tracklow~=1 & NFLw1~=. & NFLw3~=.

save, replace

//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta


//ANALYSIS B//
reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta
reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2 

save, replace





*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta



//ANALYSIS B//
reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta
reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1 

save, replace



//INTERACTION BY Sex//
use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg GM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg WM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace


*********************TABLE S3: NFLw1w3trackhigh, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace


//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==2

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1

save, replace


//INTERACTION BY Sex//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1
mi estimate: reg GM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace




***********MODEL 4: MODEL 2+liver/kidney disease******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//

mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 

save, replace

//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1 

save, replace

**INTERACTION BY Sex**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1
mi estimate: reg GM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  if sample_final==1 

save, replace

***********MODEL 5: MODEL 2+oxidative stress******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1 

save, replace


**********INTERACTION BY Sex***************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg GM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==2



//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3trackhigh Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1 

save, replace

******************INTERACTION BY Sex***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1
mi estimate: reg GM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3trackhigh##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace




***********************TABLE S7: NFLw1w3tracklow, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta

//ANALYSIS B//
reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta
reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta

//ANALYSIS C//
reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1, beta 



**Model 2: BMI-Adjusted**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta


//ANALYSIS B//
reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta
reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2 

save, replace





*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta


//ANALYSIS B//
reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta
reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1 

save, replace



//INTERACTION BY Sex//
use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg GM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg WM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace


*********************TABLE S3: NFLw1w3tracklow, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace


//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 & Sex==1

save, replace


//INTERACTION BY Sex//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1
mi estimate: reg GM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final==1
mi estimate: reg WM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final==1 

save, replace




***********MODEL 4: MODEL 2+liver/kidney disease******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//

mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 

save, replace

//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1 & Sex==1 

save, replace

**INTERACTION BY Sex**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1
mi estimate: reg GM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final==1
mi estimate: reg WM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  if sample_final==1 

save, replace

***********MODEL 5: MODEL 2+oxidative stress******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==2



//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 & Sex==1 

save, replace


**********INTERACTION BY Sex***************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg GM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final==1 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 & Sex==1 

save, replace

******************INTERACTION BY Sex***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1
mi estimate: reg GM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final==1
mi estimate: reg WM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final==1 

save, replace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////TABLES 4 AND S4-S5/////////////////////////////////////////////

***********************TABLE 3: NFLw1w3tracklow, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta
reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1,beta

//ANALYSIS B//
reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta
reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta

//ANALYSIS C//
reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta
reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==2,beta

//ANALYSIS B//
reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta
reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==2
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==2

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==2 

save, replace


*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta
reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN if sample_final==1 & Sex==1,beta


//ANALYSIS B//
reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta
reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1
mi estimate: reg GM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1 & Sex==1
mi estimate: reg WM NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1
mi estimate: reg Right_Hippocampus NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume NFLw1w3tracklow Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 & Sex==1 

save, replace



//INTERACTION BY Sex//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1
mi estimate: reg GM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final==1
mi estimate: reg WM c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1
mi estimate: reg Right_Hippocampus c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.NFLw1w3tracklow##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final==1 

save, replace




capture log close

 
 
///////////////////////////////////////////////////////SENSITIVITY ANALYSIS # 2: LARGEST SAMPLE SIZE//////////////////////////////////////////////////

 
**************************************ALL FINDINGS***********************************************

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLE2_S3_ALLFINDINGS_SENS2.smcl",replace


////////////////////////////////////////TABLES 2 AND S3/////////////////////////////////////////////


***********************TABLE 3: LnNFLw1, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ,beta



//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 ,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 ,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 , beta 



**Model 2: BMI-Adjusted**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI 
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if Sex==2,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if Sex==2,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if Sex==2,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==2,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==2 

save, replace





*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if Sex==1,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if Sex==1,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if Sex==1,beta



//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==1,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==1 

save, replace



//INTERACTION BY Sex//
use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

save, replace


*********************TABLE S3: LnNFLw1, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose 
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

save, replace


//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if Sex==2

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==1

save, replace


//INTERACTION BY Sex//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose 
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

save, replace




***********MODEL 4: MODEL 2+liver/kidney disease******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid 
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

//ANALYSIS C//

mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

save, replace

//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==1 

save, replace

**INTERACTION BY Sex**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid 
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  

save, replace

***********MODEL 5: MODEL 2+oxidative stress******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct 
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct 
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct 

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==1 

save, replace


**********INTERACTION BY Sex***************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH 
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH 

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if Sex==2



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==1 

save, replace

******************INTERACTION BY Sex***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH 
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

save, replace



capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLE3_S4_ALLFINDINGS_SENS2.smcl", replace



////////////////////////////////////////TABLES 3 AND S4/////////////////////////////////////////////


***********************TABLE 3: LnNFLw3, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ,beta

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 ,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 ,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 , beta 



**Model 2: BMI-Adjusted**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI 
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI 

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==2,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==2,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==2,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==2,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==2 

save, replace





*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==1,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==1 

save, replace



//INTERACTION BY Sex//
use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

save, replace


*********************TABLE S3: LnNFLw3, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose 
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

save, replace


//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Sex==1

save, replace


//INTERACTION BY Sex//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose 
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

save, replace




***********MODEL 4: MODEL 2+liver/kidney disease******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid 
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid 

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

//ANALYSIS C//

mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

save, replace

//Males//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Sex==1 

save, replace

**INTERACTION BY Sex**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid 
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  

save, replace

***********MODEL 5: MODEL 2+oxidative stress******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct 
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct 
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct 


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==2



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Sex==1 

save, replace


**********INTERACTION BY Sex***************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH 
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH 

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

save, replace


//Males//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==2 

save, replace



//Females//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Sex==1 

save, replace

******************INTERACTION BY Sex***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH 
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH 

//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

save, replace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////TABLES 4 AND S4-S5/////////////////////////////////////////////

***********************TABLE 3: LnNFLw3, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ,beta

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 ,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 ,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 ,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==2,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==2,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==2,beta

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==2,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==2

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==2 

save, replace


*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if Sex==1,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if Sex==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Sex==1 

save, replace



//INTERACTION BY Sex//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 

save, replace




capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLES5_SENS2.smcl",replace


//////////////////////////////////////////////////TABLE S5: NFL at v1 AND v3 vs. Hippocampus and LnLesion volume adjusted for ICV_volM2,  by race///////////////////////////////////////////////////////////////////////////////////
	  
	  ////////////////////////////////////////LnNFLw1 EXPOSURE/////////////////////////////////////////////


***********************LnNFLw1, MODELS 1 AND 2*********************

*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==2,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==2 

save, replace





*********************WHITES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==1,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==1,beta 


**Model 2**

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==1 

save, replace



//INTERACTION BY Race//


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 

save, replace


*********************LnNFLw1, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==2 

save, replace



//Whites//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==1 

save, replace


//INTERACTION by Race//



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==1 

save, replace

**INTERACTION by Race**


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

save, replace


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  

save, replace

***********MODEL 5: MODEL2+OXIDATIVE STRESS******


//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==2 

save, replace



//WHITE//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==1 

save, replace


**********INTERACTION by Race***************



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

save, replace



***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  



//AFRICAN-AMERICAN//


use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==1 

save, replace

******************INTERACTION by Race***********************


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

save, replace


////////////////////////////////////////TABLE S5. LnNFLw3 exposure/////////////////////////////////////////////




***********************LnNFLw3, MODELS 1 AND 2*********************


*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==2,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==2,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==2 

save, replace





*********************WHITE******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if Race==1,beta 


**Model 2**

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if Race==1 

save, replace



//INTERACTION by Race//


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 

save, replace


*********************LnNFLw3, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==2 

save, replace



//WHITE//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if Race==1 

save, replace


//INTERACTION by Race//



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if Race==1 

save, replace

**INTERACTION by Race**



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 

save, replace


***********MODEL 5: MODEL2+OXIDATIVE STRESS******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if Race==1 

save, replace


**********INTERACTION by Race***************


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  

//AFRICAN-AMERICAN//

use finaldata_imputed,clear




//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==2 

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if Race==1 

save, replace

******************INTERACTION by Race***********************



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 

save, replace



  

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT\TABLES6_SENS2.smcl",replace


/////////////////////////////////////////////////TABLE S6: TOTAL AND GM/WM VOLUMES VS. NFL exposures by Race***********************
////////////////////////////////////////LnNFLw1 EXPOSURE/////////////////////////////////////////////


***********************LnNFLw1, MODELS 1 AND 2*********************

*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if Race==2,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if Race==2,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if Race==2,beta


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==2

save, replace

*********************WHITES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if Race==1,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if Race==1,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if Race==1,beta


**Model 2**

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==1

save, replace



//INTERACTION BY Race//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   

save, replace


*********************LnNFLw1, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==2

save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==1


save, replace


//INTERACTION BY Race//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==2

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==1

save, replace

**INTERACTION BY Race**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  


save, replace


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid   

save, replace

***********MODEL 5: MODEL2+OXIDATIVE STRESS******

//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

save, replace


//AFRICAN-AMERICAN//


use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==2

mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==2


save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==1

mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==1

save, replace


**********INTERACTION BY Race***************



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

save, replace



***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  

save, replace


//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==2


save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==1


save, replace

******************INTERACTION BY Race***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  

save, replace


////////////////////////////////////////LnNFLw3 exposure/////////////////////////////////////////////


***********************LnNFLw3, MODELS 1 AND 2*********************

*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if Race==2,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if Race==2,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if Race==2,beta


**Model 2**

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==2

save, replace

*********************WHITES******************************************************

**Model 1**


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if Race==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if Race==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if Race==1,beta


**Model 2**


use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if Race==1

save, replace



//INTERACTION BY Race//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   

save, replace


*********************LnNFLw3, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//


use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==2

save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if Race==1


save, replace


//INTERACTION BY Race//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==2

save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if Race==1

save, replace

**INTERACTION BY Race**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  


save, replace


***********MODEL 5: MODEL2+OXIDATIVE STRESS******

//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==2

mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==2


save, replace



//WHITES//

use finaldata_imputed,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==1

mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  Race==1

save, replace


**********INTERACTION BY Race***************



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final==1

save, replace



***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//AFRICAN-AMERICAN//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==2


save, replace



//WHITES//

use finaldata_imputed,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if Race==1


save, replace

******************INTERACTION BY Race***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  

save, replace



capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT2\TBI_MMSE_DATAMANGEMENT.smcl",replace

**********************************SENSITIVITY ANALYSIS FOR TBI AND MMSE SCORE****************************

cd "E:\HANDLS_PAPER51_NFLSMRI\DATA"


use finaldata_imputed,clear
sort HNDID
capture drop _merge
save finaldata_imputed_final, replace

use 2022-02-02_w1mms,clear
capture rename HNDid HNDID
sort HNDID


save, replace


use finaldata_imputed_final,clear
merge HNDID using 2022-02-02_w1mms
save, replace


mi estimate: mean  MMStot if sample_final==1 & tagMMS==1

mi estimate: prop  MMStot if sample_final==1 & tagMMS==1

mi estimate: prop HeadInjury if sample_final==1

capture drop exclusion
gen exclusion=.
replace exclusion=1 if MMStot<23 & MMStot~=. & tagMMS==1 | HeadInjury==2
replace exclusion=0 if exclusion~=1



capture drop sample_final2
gen sample_final2=.
replace sample_final2=1 if sample_final==1 & exclusion==0
replace sample_final2=0 if sample_final2~=1

tab sample_final2 if _mi_m==1

save, replace

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear
capture drop _merge
sort HNDID
save, replace

merge HNDID using 2022-02-02_w1mms
save, replace

capture drop exclusion
gen exclusion=.
replace exclusion=1 if MMStot<23 & MMStot~=. & tagMMS==1 | HeadInjury==2
replace exclusion=0 if exclusion~=1


capture drop sample_final2
gen sample_final2=.
replace sample_final2=1 if sample_final==1 & exclusion==0
replace sample_final2=0 if sample_final2~=1

tab sample_final2 

save, replace

**************************************ALL FINDINGS***********************************************

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT2\TABLE2_S3_ALLFINDINGS.smcl",replace


////////////////////////////////////////TABLES 2 AND S3/////////////////////////////////////////////


***********************TABLE 3: LnNFLw1, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1,beta



//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1, beta 



**Model 2: BMI-Adjusted**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==2,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==2,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==2,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==2,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==2 

save, replace





*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==1,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==1,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==1,beta



//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==1,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==1 

save, replace



//INTERACTION BY Sex//
use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final2==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 

save, replace


*********************TABLE S3: LnNFLw1, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 

save, replace


//Males//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1 & Sex==2

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==1

save, replace


//INTERACTION BY Sex//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 

save, replace




***********MODEL 4: MODEL 2+liver/kidney disease******

//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1

//ANALYSIS C//

mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 

save, replace

//Males//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==1 

save, replace

**INTERACTION BY Sex**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  if sample_final2==1 

save, replace

***********MODEL 5: MODEL 2+oxidative stress******

//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final2==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final2==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final2==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 

save, replace


//Males//


use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==1 

save, replace


**********INTERACTION BY Sex***************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 

save, replace


//Males//


use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1 & Sex==2



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==1 

save, replace

******************INTERACTION BY Sex***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1
mi estimate: reg GM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg WM c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 

save, replace



capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT2\TABLE3_S4_ALLFINDINGS.smcl", replace



////////////////////////////////////////TABLES 3 AND S4/////////////////////////////////////////////


***********************TABLE 3: LnNFLw3, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1,beta

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1, beta 



**Model 2: BMI-Adjusted**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==2,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==2,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==2,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==2,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==2 

save, replace





*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==1,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==1 

save, replace



//INTERACTION BY Sex//
use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final2==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 

save, replace


*********************TABLE S3: LnNFLw3, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 

save, replace


//Males//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Sex==1

save, replace


//INTERACTION BY Sex//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 

save, replace




***********MODEL 4: MODEL 2+liver/kidney disease******

//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1

//ANALYSIS C//

mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 

save, replace

//Males//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Sex==1 

save, replace

**INTERACTION BY Sex**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2  if sample_final2==1 

save, replace

***********MODEL 5: MODEL 2+oxidative stress******

//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final2==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final2==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 

save, replace


//Males//


use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==2



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Sex==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Sex==1 

save, replace


**********INTERACTION BY Sex***************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 

save, replace


//Males//


use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1 & Sex==2


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==2 

save, replace



//Females//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Sex==1 

save, replace

******************INTERACTION BY Sex***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH if sample_final2==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 

save, replace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////TABLES 4 AND S4-S5/////////////////////////////////////////////

***********************TABLE 3: LnNFLw3, MODELS 1 AND 2*********************

**ANALYSES A-C, TOTAL POPULATION**

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1,beta

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final2==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 

save, replace

*********************MALES*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==2,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==2,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==2,beta

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==2,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==2,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Sex==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==2

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==2 

save, replace


*********************FEMALES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN if sample_final2==1 & Sex==1,beta


//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN ICV_volM2 if sample_final2==1 & Sex==1,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Sex==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI if sample_final2==1 & Sex==1

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Sex==1 

save, replace



//INTERACTION BY Sex//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1
mi estimate: reg GM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final2==1
mi estimate: reg WM c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Sex Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 

save, replace




capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT2\TABLES5.smcl",replace


//////////////////////////////////////////////////TABLE S5: NFL at v1 AND v3 vs. Hippocampus and LnLesion volume adjusted for ICV_volM2,  by race///////////////////////////////////////////////////////////////////////////////////
	  
	  ////////////////////////////////////////LnNFLw1 EXPOSURE/////////////////////////////////////////////


***********************LnNFLw1, MODELS 1 AND 2*********************

*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==2,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==2,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==2 

save, replace





*********************WHITES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS B//
reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==1,beta
reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==1,beta 


**Model 2**

use finaldata_imputed_final,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==1 

save, replace



//INTERACTION BY Race//


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final2==1 

save, replace


*********************LnNFLw1, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==2 

save, replace



//Whites//

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==1 

save, replace


//INTERACTION by Race//



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed_final,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==1 

save, replace

**INTERACTION by Race**


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 

save, replace


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1  

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1  

save, replace

***********MODEL 5: MODEL2+OXIDATIVE STRESS******


//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==2 

save, replace



//WHITE//

use finaldata_imputed_final,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==1 

save, replace


**********INTERACTION by Race***************



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 

save, replace



***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  



//AFRICAN-AMERICAN//


use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==1 

save, replace

******************INTERACTION by Race***********************


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 

save, replace


////////////////////////////////////////TABLE S5. LnNFLw3 exposure/////////////////////////////////////////////




***********************LnNFLw3, MODELS 1 AND 2*********************


*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==2,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==2,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==2,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==2 

save, replace





*********************WHITE******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear

//ANALYSIS B//
reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==1,beta
reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==1,beta

//ANALYSIS C//
reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN    ICV_volM2 if sample_final2==1 & Race==1,beta 


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI ICV_volM2 if sample_final2==1 & Race==1 

save, replace



//INTERACTION by Race//


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  ICV_volM2 if sample_final2==1 

save, replace


*********************LnNFLw3, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==2 

save, replace



//WHITE//

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 & Race==1 

save, replace


//INTERACTION by Race//



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose ICV_volM2 if sample_final2==1 

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear

//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 & Race==1 

save, replace

**INTERACTION by Race**



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid ICV_volM2 if sample_final2==1 

save, replace


***********MODEL 5: MODEL2+OXIDATIVE STRESS******

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear


//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed_final,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 & Race==1 

save, replace


**********INTERACTION by Race***************


//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct ICV_volM2 if sample_final2==1 

save, replace

***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear




//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==2
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==2

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==2 

save, replace



//WHITES//

use finaldata_imputed_final,clear



//ANALYSIS B//
mi estimate: reg Left_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==1
mi estimate: reg Right_Hippocampus LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 & Race==1 

save, replace

******************INTERACTION by Race***********************



//ANALYSIS B//
mi estimate: reg Left_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1
mi estimate: reg Right_Hippocampus c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1

//ANALYSIS C//
mi estimate: reg LnLesion_Volume c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH ICV_volM2 if sample_final2==1 

save, replace



  

capture log close

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT2\TABLES6.smcl",replace


/////////////////////////////////////////////////TABLE S6: TOTAL AND GM/WM VOLUMES VS. NFL exposures by Race***********************
////////////////////////////////////////LnNFLw1 EXPOSURE/////////////////////////////////////////////


***********************LnNFLw1, MODELS 1 AND 2*********************

*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==2,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==2,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==2,beta


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==2

save, replace

*********************WHITES******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==1,beta
reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==1,beta
reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==1,beta


**Model 2**

use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==1

save, replace



//INTERACTION BY Race//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final2==1
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final2==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final2==1

save, replace


*********************LnNFLw1, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==2

save, replace



//WHITES//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==1


save, replace


//INTERACTION BY Race//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==2

save, replace



//WHITES//

use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==1

save, replace

**INTERACTION BY Race**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1


save, replace


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1  

save, replace

***********MODEL 5: MODEL2+OXIDATIVE STRESS******

//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1

save, replace


//AFRICAN-AMERICAN//


use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==2

mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==2


save, replace



//WHITES//

use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==1

mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==1

save, replace


**********INTERACTION BY Race***************



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1

mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1

save, replace



***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//Overall//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1

save, replace


//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==2
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==2


save, replace



//WHITES//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==1
mi estimate: reg GM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw1 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==1


save, replace

******************INTERACTION BY Race***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg GM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg WM c.LnNFLw1##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1

save, replace


////////////////////////////////////////LnNFLw3 exposure/////////////////////////////////////////////


***********************LnNFLw3, MODELS 1 AND 2*********************

*********************AFRICAN-AMERICAN*******************************************************

**Model 1**

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==2,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==2,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==2,beta


**Model 2**

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==2

save, replace

*********************WHITES******************************************************

**Model 1**


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear


//ANALYSIS A//
reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==1,beta
reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==1,beta
reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN     if sample_final2==1 & Race==1,beta


**Model 2**


use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  if sample_final2==1 & Race==1

save, replace



//INTERACTION BY Race//


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final2==1
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final2==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI   if sample_final2==1

save, replace


*********************LnNFLw3, MODELS 3-6***************************

***********MODEL 3: MODEL 2+w1dxDiabetes w1Glucose******

//AFRICAN-AMERICAN//


use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==2

save, replace



//WHITES//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1 & Race==1


save, replace


//INTERACTION BY Race//



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1dxDiabetes w1Glucose  if sample_final2==1

save, replace


***********MODEL 4: MODEL 2+liver/kidney disease******

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==2

save, replace



//WHITES//

use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1 & Race==1

save, replace

**INTERACTION BY Race**


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1Creatinine w1USpecGrav w1BUN w1ALP w1UricAcid  if sample_final2==1


save, replace


***********MODEL 5: MODEL2+OXIDATIVE STRESS******

//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==2

mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==2


save, replace



//WHITES//

use finaldata_imputed_final,clear



//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==1

mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1 & Race==1

save, replace


**********INTERACTION BY Race***************



//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1

mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI  w1TotalD w1Albumin w1EosinPct if  sample_final2==1

save, replace



***********MODEL 6: MODEL 2+lifestyle/health-related factors*******  


//AFRICAN-AMERICAN//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==2
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==2
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==2


save, replace



//WHITES//

use finaldata_imputed_final,clear


//ANALYSIS A//
mi estimate: reg TOTALBRAIN LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==1
mi estimate: reg GM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==1
mi estimate: reg WM LnNFLw3 Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1 & Race==1


save, replace

******************INTERACTION BY Race***********************


//ANALYSIS A//
mi estimate: reg TOTALBRAIN c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg GM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1
mi estimate: reg WM c.LnNFLw3##Race Sex w1Age Race PovStat TIME_V1SCAN w1BMI w1currdrugs w1SRH  if sample_final2==1

save, replace

capture log close


****************************************************REVISION OF MANUSCRIPT: ADDITIONAL ANALYSES***************************************************

capture log using "E:\HANDLS_PAPER51_NFLSMRI\OUTPUT3\SUPPFIG_ABETA_TAU.smcl",replace


use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear
sort HNDID
save, replace

use 2020-12-22_abeta,clear
capture rename HNDid HNDID
save, replace
keep HNDID HNDwave tau b40 b42 Ratio4240
save abeta_tau_long, replace

keep if HNDwave==1
save abeta_tau_w1, replace
capture rename tau w1tau
capture rename b40 w1b40
capture rename b42 w1b42
capture rename Ratio4240 w1Ratio4240
save abeta_tau_w1, replace
sort HNDID
save, replace

use abeta_tau_long,clear
keep if HNDwave==3
save abeta_tau_w3, replace
capture rename tau w3tau
capture rename b40 w3b40
capture rename b42 w3b42
capture rename Ratio4240 w3Ratio4240
save abeta_tau_w3, replace
sort HNDID
save, replace

use HANDLS_paper51_NFLBRAINSCANFINALIZED,clear
capture drop _merge
save, replace
merge HNDID using abeta_tau_w1
tab _merge
capture drop _merge
sort HNDID
save, replace
merge HNDID using abeta_tau_w3
tab _merge
capture drop _merge
sort HNDID
save, replace

capture drop LNw1tau LNw1b40 LNw1b42 LNw1Ratio4240 LNw3tau LNw3b40 LNw3b42 LNw3Ratio4240
foreach x of varlist w1tau w1b40 w1b42 w1Ratio4240 w3tau w3b40 w3b42 w3Ratio4240 {
	gen LN`x'=ln(`x')
}
 

su LnNFLw1 LnNFLw3 bayes1LnNFL LNw1tau LNw1b40 LNw1b42 LNw1Ratio4240 LNw3tau LNw3b40 LNw3b42 LNw3Ratio4240 if sample_final==1
graph matrix LnNFLw1 LnNFLw3 bayes1LnNFL LNw1tau LNw1b40 LNw1b42 LNw1Ratio4240 LNw3tau LNw3b40 LNw3b42 LNw3Ratio4240 if sample_final==1

graph save SUPPFIG_ABETA_TAU.gph,replace

pwcorr LnNFLw1 LnNFLw3 bayes1LnNFL LNw1tau LNw1b40 LNw1b42 LNw1Ratio4240 LNw3tau LNw3b40 LNw3b42 LNw3Ratio4240 if sample_final==1,sig

save, replace

capture log close





