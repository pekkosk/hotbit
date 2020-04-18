from .pylibxc_functionals import available_xcs, xc2libxc

###################
## error messages #
###################
_ERR_NOLIBXC  = "Could not load python module of libxc. "
_ERR_NOLIBXC += "Only PW92 functional is supported natively. \n"
_ERR_NOLIBXC += "Please choose PW92 or add install libxc and its "
_ERR_NOLIBXC += "python module and add it to your PYTHONPATH "
_ERR_NOLIBXC += "environment variable! See README for further info."

_ERR_NOSTRXC  = "Definition of exchange-correlation functional "
_ERR_NOSTRXC += "has to be of type string. Aborting."

_ERR_UNKNOWNXC  = "The specified xc-functional '{0:s}' is not known"
_ERR_UNKNOWNXC += " or not supported yet.\nAvailable functionals:\n"
_ERR_UNKNOWNXC += ", ".join(["'"+f.upper()+"'" for f in available_xcs])

_ERR_NOGRID  = "GGA and Hybrid calculations require to set the grid "
_ERR_NOGRID += "via set_grid(grid)! I didn't find self.grid. Aborting."

_ERR_NORNL  = "Hybrid calculations require to the the radial wavefunctions "
_ERR_NORNL += "via set_Rnl(Rnl)! I didn't find self.Rnl. Aborting."


####################
# helper functions #
####################

## convert generic xc name to libxc definition
## TODO: allow for direct input of libxc functional definition?
def translate_xc(xc):
    if not (type(xc) == str): raise ValueError(_ERR_NOSTRXC)
    xc_l = xc.lower()
    if xc_l in available_xcs: return xc2libxc[xc_l]
    raise ValueError(_ERR_UNKNOWNXC.format(xc))
    



#--EOF--#
