#################################################################
##  List of functionals from libxc available in the interface  ##
##  and map for generic functional names to libxc definitions  ##
#################################################################

xc2libxc = {
## LDAs ##
'gl'        :   'LDA_X,LDA_C_GL',
'hl'        :   'LDA_X,LDA_C_HL',
'pw92'      :   'LDA_X,LDA_C_PW',
'pz'        :   'LDA_X,LDA_C_PZ',
'svwn'      :   'LDA_X,LDA_C_VWN',
'teter93'   :   'LDA_XC_TETER93',
'zlp'       :   'LDA_XC_ZLP',
## GGAs ##
'b96pbe'    :   'GGA_X_B86_MGC,GGA_C_PBE',
'beefvdw'   :   'GGA_XC_BEEFVDW',
'blyp'      :   'GGA_X_B88,GGA_C_LYP',
'bp86'      :   'GGA_X_B88,GGA_C_P86',
'gam'       :   'GGA_X_GAM,GGA_C_GAM',
'pbe'       :   'GGA_X_PBE,GGA_C_PBE',
'pbeint'    :   'GGA_X_PBEINT,GGA_C_PBEINT',
'pw86pbe'   :   'GGA_X_PW86,GGA_C_PBE',
'pw91'      :   'GGA_X_PW91,GGA_C_PW91',
'revpbe'    :   'GGA_X_PBE_R,GGA_C_PBE',
'rpbe'      :   'GGA_X_RPBE,GGA_C_PBE',
'sogga'     :   'GGA_X_SOGGA,GGA_C_PBE',
'sogga11'   :   'GGA_X_SOGGA11,GGA_C_SOGGA11',
'vv10'      :   'GGA_XC_VV10',
'xlyp'      :   'GGA_XC_XLYP',
## Hybrid GGAs ##
#'b1lyp'     :   'HYB_GGA_XC_B1LYP',
#'b1pw91'    :   'HYB_GGA_XC_B1PW91',
#'b3lyp'     :   'HYB_GGA_XC_B3LYP',
#'b3p86'     :   'HYB_GGA_XC_B3P86',
#'b3pw91'    :   'HYB_GGA_XC_B3PW91',
#'b97'       :   'HYB_GGA_XC_B97',
#'cam-b3lyp' :   'HYB_GGA_XC_CAM_B3LYP',
#'hjs-b88'   :   'HYB_GGA_XC_HJS_B88',
#'hjs-b97x'  :   'HYB_GGA_XC_HJS_B97X',
#'hjs-pbe'   :   'HYB_GGA_XC_HJS_PBE',
#'hjs-pbesol':   'HYB_GGA_XC_HJS_PBE_SOL',
#'hse03'     :   'HYB_GGA_XC_HSE03',
#'hse06'     :   'HYB_GGA_XC_HSE06',
#'lc-vv10'   :   'HYB_GGA_XC_LC_VV10',
#'o3lyp'     :   'HYB_GGA_XC_O3LYP',
#'pbe0'      :   'HYB_GGA_XC_PBEH',
#'pbeh'      :   'HYB_GGA_XC_PBEH',
#'rev-b3lyp' :   'HYB_GGA_XC_REVB3LYP',
#'wb97x'     :   'HYB_GGA_XC_WB97X',
#'x3lyp'     :   'HYB_GGA_XC_X3LYP',
}

for key, val in xc2libxc.items():
    if not key.islower(): xc2libxc[key.lower()] = xc2libxc.pop(key)

available_xcs = [xc for xc in xc2libxc.keys()]


#--EOF--#
