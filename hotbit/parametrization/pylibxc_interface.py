from numpy import clip, where
from .pylibxc_helper import _ERR_NOLIBXC, _ERR_NOGRID, _ERR_NORNL, translate_xc
try:
    import pylibxc
except ImportError:
    raise ImportError(_ERR_NOLIBXC)


_ggaIDs = [pylibxc.flags.XC_FAMILY_GGA, pylibxc.flags.XC_FAMILY_HYB_GGA]
_hybIDs = [pylibxc.flags.XC_FAMILY_HYB_GGA, pylibxc.flags.XC_FAMILY_HYB_MGGA]

##########################################################################
##  Interface to python module of libxc (www.tddft.org/programs/libxc)  ##
##########################################################################
#                                                                        #
#  . libxc assumes Exc = \int fxc dr, fxc and xc energy density thus     #
#    refer to the exchange-correlation energy density per unit volume    #
#    as usually defined for GGA functionals and not the energy density   #
#    per unit particle as formulated in, e.g. the LDA.                   #
#  . this implementation currently supports a select set of LDA and GGA  #
#    functionals (see pylibxc_functionals).                              #
#  . Further support of meta-GGA functionals and Hybrid functionals is   #
#    planned for future versions.                                        #
#                                                                        #
##########################################################################



class libXCFunctional:
    """ libxc wrapper. """
    def __init__(self, xc_name_in):
        self.mini = 1e-90
        self._xGGA, self._cGGA, self._xcGGA = False, False, False
        self._EXXfrac = 0.
        
        xc_name = translate_xc(xc_name_in)
        self._sep_xc = (',' in xc_name)
        if self._sep_xc:
            [x_name, c_name] = xc_name.split(',')
            self._pylibxc_x = pylibxc.LibXCFunctional(x_name, "unpolarized")
            self._pylibxc_c = pylibxc.LibXCFunctional(c_name, "unpolarized")
            self._xGGA = (self._pylibxc_x.get_family() in _ggaIDs)
            self._cGGA = (self._pylibxc_c.get_family() in _ggaIDs)
            self._doEXX = (self._pylibxc_x.get_family() in _hybIDs)
            if self._doEXX: self._EXXfrac = self._pylibxc_x.get_hyb_exx_coef()
        else:
            self._pylibxc_xc = pylibxc.LibXCFunctional(xc_name, "unpolarized")
            self._xcGGA = (self._pylibxc_xc.get_family() in _ggaIDs)
            self._doEXX = (self._pylibxc_xc.get_family() in _hybIDs)
            if self._doEXX: self._EXXfrac = self._pylibxc_xc.get_hyb_exx_coef()
        
        self._anyGGA = any([self._xGGA, self._cGGA, self._xcGGA])
        if self._doEXX:
            raise NotImplementedError("Exact exchange not implemented yet.")
        
    
    def set_grid(self, grid):
        """ Communicates grid to xc instance (Needed for GGAs). """
        self.grid = grid
        
    
    def set_Rnl(self, Rnl):
        """
        Communicates radial wave functions to xc instance (Needed for hybrids).
        """
        self.Rnl = Rnl
        
    
    def exc(self, rho):
        """ Returns exchange-correlation energy for density rho. """
        inp = {'rho':rho}
        if self._anyGGA:
            if not hasattr(self, 'grid'): raise ValueError(_ERR_NOGRID)
            drho_dr = where( rho<self.mini, 0.0, self.grid.derivative(rho))
            inp['sigma'] = drho_dr * drho_dr
        
        if self._sep_xc:
            # calculate exchange energy (density)
            res = self._pylibxc_x.compute(inp, do_exc=True, do_vxc=False)
            e_x = (1. - self._EXXfrac) * res['zk'][0]
            # calculate correlation energy (density)
            res = self._pylibxc_c.compute(inp, do_exc=True, do_vxc=False)
            e_c = res['zk'][0]
            e_xc = e_x + e_c
        else:
            res = self._pylibxc_xc.compute(inp, do_exc=True, do_vxc=False)
            e_xc = res['zk'][0]
        
        if self._doEXX:
            if not hasattr(self, 'grid'): raise ValueError(_ERR_NOGRID)
            if not hasattr(self, 'Rnl'): raise ValueError(_ERR_NORNL)
            exc += self._EXXfrac * self.get_exx()
        
        return e_xc
        
    
    def vxc(self, rho):
        """ Returns exchange-correlation potential for density rho. """
        inp = {'rho':rho}
        if self._anyGGA:
            if not hasattr(self, 'grid'): raise ValueError(_ERR_NOGRID)
            drho_dr = where( rho<self.mini, 0.0, self.grid.derivative(rho))
            sigma = drho_dr * drho_dr
            inp['sigma'] = sigma

        if self._sep_xc:
            # calculate exchange potential
            res = self._pylibxc_x.compute(inp, do_exc=False, do_vxc=True)
            v_x = res['vrho'][0]
            # GGA contribution to vx
            if self._xGGA: v_x -= self.v_gga(res['vsigma'][0], drho_dr)
            v_x *= (1. - self._EXXfrac)
            
            # calculate correlation potential
            res = self._pylibxc_c.compute(inp, do_exc=False, do_vxc=True)
            v_c = res['vrho'][0]
            # GGA contribution to vc
            if self._cGGA: v_c -= self.v_gga(res['vsigma'][0], drho_dr)
            
            v_xc = v_x + v_c
        else:
            # calculate exchange-correlation potential
            res = self._pylibxc_xc.compute(inp, do_exc=False, do_vxc=True)
            v_xc = res['vrho'][0]
            # GGA contribution to vxc
            if self._xcGGA: v_xc -= self.v_gga(res['vsigma'][0], drho_dr)
        
        if self._doEXX:
            if not hasattr(self, 'grid'): raise ValueError(_ERR_NOGRID)
            if not hasattr(self, 'Rnl'): raise ValueError(_ERR_NORNL)
            vxc += self._EXXfrac * self.get_vxx()
        
        return v_xc
        
    
    def v_gga(self, de_dsigma, grad):
        """
        Returns GGA contribution to radial vx, vc, or vxc, as given by:
            2/r * dfxc/dgrad + d/dr(dfxc/dgrad),
        where
            grad = drho/dr
            dfxc/dgrad = 2 * dfxc/dsigma * grad
            sigma = grad * grad
        
        Arguments:
        ==========
            . de_dsigma: [nd-array] partial derivative of the xc energy density
                                    w.r.t. the square of the density gradient
            . grad:      [nd-array] partial derivative of the density w.r.t. the
                                    position in the radial grid, r
        """
        de_dgrad = 2. * de_dsigma * grad
        d2e_drdgrad = self.grid.derivative(de_dgrad)
        r = clip(self.grid.get_grid(), self.mini, None)
        TWOde_rdgrad = 2. * de_dgrad / r
        return TWOde_rdgrad + d2e_drdgrad
        
    

#########################################################################
##  Hybrids: libxc provides the density-based parts of Exc, vxc, etc.  ##
#########################################################################
#                                                                       #
#  . EXX energy contribution:                                           #
#      exx = -1/2 \sum_i,j \int n_HF(r,r') / |r-r'| dr' / rho           #
#    with n_HF(r,r') = Rnl_i(r) Rnl_i(r') Rnl_j(r) Rnl_j(r').           #
#  . EXX potential contribution:                                        #
#      vxx = dexx/dn = dexx/dRnl dRnl/dveff dveff/drho = ??             #
#                                                                       #
#########################################################################
    
    def get_exx(self):
        """
        Returns exact exchange contribution to xc energy density:
            exx = -1/2 \sum_i,j \int n_HF(r,r') / |r-r'| dr' / rho
        with n_HF(r,r') = Rnl_i(r) Rnl_i(r') Rnl_j(r) Rnl_j(r')
        where Rnl_i are the radial wave functions.
        """
        raise NotImplementedError("Exact exchange not implemented yet.")
        
    
    def get_vxx(self):
        """
        Returns exact exchange contribution to radial vxc or vx:
            vxx = ???
        """
        raise NotImplementedError("Exact exchange not implemented yet.")
        
    

#--EOF--#
