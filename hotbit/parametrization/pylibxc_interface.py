from numpy import clip, where
from .pylibxc_helper import _ERR_NOLIBXC, _ERR_NOGRID, translate_xc
try:
    import pylibxc
except ImportError:
    raise ImportError(_ERR_NOLIBXC)


_ggaIDs = [pylibxc.flags.XC_FAMILY_GGA, pylibxc.flags.XC_FAMILY_HYB_GGA]

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
        
        xc_name = translate_xc(xc_name_in)
        self._sep_xc = (',' in xc_name)
        if self._sep_xc:
            [x_name, c_name] = xc_name.split(',')
            self._pylibxc_x = pylibxc.LibXCFunctional(x_name, "unpolarized")
            self._pylibxc_c = pylibxc.LibXCFunctional(c_name, "unpolarized")
            self._xGGA = (self._pylibxc_x.get_family() in _ggaIDs)
            self._cGGA = (self._pylibxc_c.get_family() in _ggaIDs)
        else:
            self._pylibxc_xc = pylibxc.LibXCFunctional(xc_name, "unpolarized")
            self._xcGGA = (self._pylibxc_xc.get_family() in _ggaIDs)
        
        self._anyGGA = any([self._xGGA, self._cGGA, self._xcGGA])
        
    
    def set_grid(self, grid):
        """ Communicates grid to xc instance (Needed for GGAs). """
        self.grid = grid
        
    
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
            # calculate correlation energy (density)
            res = self._pylibxc_c.compute(inp, do_exc=True, do_vxc=False)
            e_c = res['zk'][0]
            e_xc = e_x + e_c
        else:
            res = self._pylibxc_xc.compute(inp, do_exc=True, do_vxc=False)
            e_xc = res['zk'][0]
        
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


#--EOF--#
