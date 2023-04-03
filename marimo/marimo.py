import numpy as np
import astropy.constants as cst
from scipy.special import voigt_profile
import matplotlib.pyplot as plt
import sys
import threading
sys.path.append('./')
import core
import time

import joblib
from scipy.interpolate import LinearNDInterpolator
from scipy import ndimage

from tqdm import tqdm

from galario.double import sampleImage, threads

c = cst.c.cgs.value
k_B = cst.k_B.cgs.value
h = cst.h.cgs.value
m_p = cst.m_p.cgs.value
au = cst.au.cgs.value
G = cst.G.cgs.value
GHz = 1e9
mu_gas = 2.37

threads(num=1)

class functions:

    def gaussian(x, x0, A, FWHM):
        sigma = FWHM/(2*np.sqrt(2*np.log(2)))
        return A*np.exp(-(x-x0)**2/2/sigma)
    
    def voigt(v, sigma, gamma):
        return voigt_profile( v, sigma, gamma )/voigt_profile( 0, sigma, gamma )

    def C(T):
    
        '''
        # HITRAN database CO 3-2
        gamma_H2 = 0.0703
        r_H2 = 0.66
        gamma_He = 0.0467
        r_He = 0.51
        Tref = 296
        '''

        # HITRAN database CO 2-1
        gamma_H2 = 0.0711
        r_H2 = 0.66
        gamma_He = 0.0472
        r_He = 0.51
        Tref = 296
    
        atm = 1013250 # Ba (cgsの圧力)
    
        return ( gamma_H2 * (Tref/T)**r_H2 + 0.19 * gamma_He * (Tref/T)**r_He ) * c*k_B*T/atm

    def B(nu, T):
        return 2*h*nu**3/c**2 / ( np.exp(h*nu/k_B/T) - 1 )

    def line_slab(v, v0, nu, cs, tau0, T, I0):
        tau = tau0 * np.exp(-(v-v0)**2/2/cs**2)
        return I0 * np.exp(-tau) + B(nu, T) * ( 1 - np.exp(-tau) )

    def cont_slab(v, v0, nu, cs, tau, T, I0):
        return I0 * np.exp(-tau) + B(nu, T) * ( 1 - np.exp(-tau) )

    def I2Tb(I, nu):
        return c**2/(2*k_B*nu**2) * I

    def Af( line_property, T ):
        # absorption cross section / line profile
        nu0 = line_property['nu0']
        Jl = line_property['Jl']
        Aul = line_property['Aul']
        m = line_property['molec_weight'] * m_p

        A1 = c**3 / ( 8*np.pi*np.sqrt(2*np.pi)*nu0**3 ) * (Jl+3/2)/(Jl+1)
        A2 = Aul * h*nu0/(k_B*T) * np.sqrt(m/(k_B*T))
        A3 = (1 - np.exp(-h*nu0/k_B/T))
        A4 = np.exp(-h*nu0/2/k_B/T * (Jl+1/(3*(Jl+1))))

        return A1*A2*A3*A4
    

class readdata:

    def molecdata(molec_df, molec_lev):

        data = open(molec_df, 'r')

        lines = data.readlines()

        N_enelev = int(lines[5])

        energy = float((lines[(6+molec_lev)].split())[1])
        g = float((lines[(6+molec_lev)].split())[2])
        Jl = float((lines[(6+molec_lev)].split())[3])
        Ju = Jl + 1

        Aul = float((lines[(6+N_enelev+3 + molec_lev)].split())[3])
        freq = float((lines[(6+N_enelev+3 + molec_lev)].split())[4])

        #print(energy)

        line_property = { 'datafile'      : molec_df,
                          'name'          : lines[1].strip(),
                          'molec_weight'  : float(lines[3]),
                          'energy'        : energy,
                          'g'             : g,
                          'Jl'            : Jl,
                          'Ju'            : Ju,
                          'Aul'           : Aul,
                          'nu0'          : freq * 1e9
                          }

        #line_property = 1
        return line_property

class disk_config:

    # set physical parameters as a function of (r, z)
    def __init__(self, r0, rin, rout, rexp, Nco_0, Nco_g, Tmid0, Tatm0, pmid, patm, Xco, M_star, rcav = None, dcav = None,  vturb = 0.0, d2g = 0.01, chi_d = 0.5, k_d = 1.0,):
        """
        Configure your disk model
    
        Parameters
        ----------
        r0 : float
            The normalized radius.
        rin : float
            Disk inner radius
        rout : float
            Disk outer radius
        rexp : float
            the exponential cutoff
        rcav : Optional[float]
            cavity radius if needed
        dcav : Optional[float]
            Depth of the cavity if needed
        d2g : Optional[float]
            dust-to-gas mass ratio
            default: 0.01
        chi_d : Optional[float]
            the intensity reduction factor as in Zhu+19
            default: 0.5
        k_d : Optional[float]
            dust opacity
            default: 1.0
        Nco_0 : float
            CO surface density at r0
        Nco_g : float
            power law index of the CO surface density
        Tmid0 : float
            The midplane temperature at r0
        Tatm0 : float
            The atmosphere temperature at r0  
        pmid : float
            The power-law index of the midplane temperature
        patm : float
            The power-law index of the atmosphere temperature
        vturb : float
            the turbulence velosity normalized by the local sound speed
            default: 0.0
        Mstar : float
            the stellar mass
        """

        # all quantities are given in cgs

        # normalize radius
        self.r0 = r0
        # inner radius
        self.rin = rin
        # outer radius
        self.rout = rout
        # expornential cutoff
        self.rexp = rexp
        # cavity radius

        if rcav != None:
            self.rcav = rcav
            # cavity depth
            self.dcav = dcav
        else:
            self.rcav = rin
            self.dcav = 1.0

        # gas surface Density at r0
        self.Sigma_g0 = Nco_0 / Xco * mu_gas * m_p
        # Gas surface density power law index
        self.gamma_g = Nco_g
        # dust/gas
        self.d2g = d2g
        # dust intensity reduction
        self.chi_d = chi_d
        # dust single scattering albedo
        self.k_d = k_d

        # dust and gas midplane temperature at r0
        self.Tmid0 = Tmid0
        # gas surface temperature at r0
        self.Tatm0 = Tatm0
        # midplane temperature power law index
        self.pmid = pmid
        # surface temperature power law index
        self.patm = patm
        # CO/H2 ratio
        self.Xco = Xco

        self.vturb = vturb

        # stellar mass
        self.M_star = M_star

        # calc. region in maximum scale height
        #self.NHg = NHg

        return None

    def tau_d(self, r):
        # dust total optical depth
        return self.chi_d * self.Sigma_d0 * ( r / self.r0 ) ** (-self.gamma_d) * ( r >= self.rcav  )
    
class calc_config:

    def __init__(self, line_files, line_levels, incl, pa, D, V, xmax, Nx, Rmin = 1e-4*au, subpix_max = 1, subpix_pow = 0.0, x0 = 0.0, y0 = 0.0, nz = 100, z2r = 0.5, bmaj = 0.1, bmin = 0.1, bpa = 0.0, voigt = True, kind = 1, dust = True, contsub = False, nthreads = 16, output = 'intensity'):
        """
        Configure calculations
    
        Parameters
        ----------
        line_files : array_like
            list of names of LAMDA database files (e.g., co.dat)).
            If voigt = 1, corresponding pressure broadening parameter files are needed (e.g. co_pb.dat).
        line_levels : array_like
            list of transition levels to be calculated.
        incl : float
            inclination angle of the disk in degrees.
        pa : float
            position angle of the disk in degrees.
        D : float
            distance to the source
        V : 1d array
            velosity grid in cm/s
        Rmin : Optional[float]
            inner calculation radius
        xmax : float
            outer calculation x coordinate
        Nx : int
            number of x-grids
        subpix_max : int
            number of sub pixels at Rmin
        subpix_pow : float
            power law index of the subpix
            -> subpix = int( subpix_max * ( R/Rmin )**(-subpix_pow) ) + 1 
        x0, y0 : Optional[float]
            beam center
        nz : Optional[int]
            number of the calculation grid along the line of sight.
        z2r : Optional[float]
            calculation boundary in z/r.
        bmaj : Optional[float]
            beam major axis in arcsec
        bmin : Optional[float]
            beam minor axis in arcsec
        bpa : Optional[float]
            beam position angle in degree
        voigt : Optional[bool]
            If ``True'', the Voigt function is used as a line profile.
        kind : Optional[int]
            1 : cpf12 by Humlicek (1979)
            2 : w4 by Humlicek (1982)
            3 : zpf16 by Humlicek (1979, 1982)
        contsub : Optional[bool]
        nthreads : Optional[int]
            nthreads for OpenMp
        dust : Optional[bool]
            If ``True'', dust continuum emission is included.
        output : Optional[str]
            'intensity'(default) : output intensity cube
            'tau' : output total optical depth
            'surf' : output height of tau=1 surface
        """

        self.line_files = line_files
        self.line_levels = line_levels
        self.incl = np.deg2rad(incl)
        self.pa = np.deg2rad(pa)
        self.V = V

        self.D = D

        #self.R_pos = np.logspace(np.log10(Rmin), np.log10(Rmax), NR)
        #self.T_pos = np.linspace(0.0, 2*np.pi, NT)
        self.X = np.linspace( -xmax, xmax, Nx )
        self.Y = self.X


        self.dX = self.X[1] - self.X[0]

        self.Rmin = Rmin
        #if Y is None:
        #    self.Y = X
        #else:
        #    self.Y = Y

        self.x0 = x0
        self.y0 = y0

        #RR_pos, TT_pos = np.meshgrid( self.R_pos, self.T_pos, indexing='ij' )

        self.XX, self.YY = np.meshgrid( self.X, self.Y, indexing='xy' )


        self.subpix_max = subpix_max
        self.subpix_pow = subpix_pow

        self.nz = nz

        if z2r >= 1/np.tan(self.incl):
            print('Error: z2r shoud be smaller than 1/tan(incl)', file=sys.stderr)
            sys.exit(1)
        else:
            self.z2r = z2r

        self.bmaj = bmaj
        self.bmin = bmin
        self.bpa = np.deg2rad(bpa)
        
        if voigt:
            self.voigt = 1
        else:
            self.voigt = 0

        self.kind = kind

        if contsub:
            self.contsub = 1
        else:
            self.contsub = 0

        self.dust = int(dust)

        self.nthreads = nthreads
        self.output = output



        # make a central mask

        xx = (self.XX/au/self.D - self.x0)*np.cos(self.bpa) - (self.YY/au/self.D - self.y0)*np.sin(self.bpa)
        yy = (self.XX/au/self.D - self.x0)*np.sin(self.bpa) + (self.YY/au/self.D - self.y0)*np.cos(self.bpa)
        Beam = np.exp( -(xx/self.fwhm2sig(self.bmaj)/2)**2 -(yy/self.fwhm2sig(self.bmin)/2)**2 )
        mask_beam = np.tile( Beam, (len(self.V), 1, 1) )

        #_dr = np.gradient(np.deg2rad(RR_pos/au/self.D/3600))
        #_dt = self.T_pos[1] - self.T_pos[0]
        #area_tmp =  np.deg2rad(RR_pos/au/self.D/3600) * _dr * _dt

        #print(np.shape(area_tmp))
        #area = np.tile( area_tmp, (len(self.V), 1, 1) ) 

        self.spectrum_mask = mask_beam #*area
    
    def fwhm2sig(self, fwhm):
        return fwhm / ( 2*np.sqrt(2*np.log(2)) )

class run:

    # set imagecube
    def __init__(self, disk, calc, progress = False):
        """
        Calculating the radiative transfer

        Parameters
        ----------
        disk : disk_config instance
        calc : calc_config instance
        """
        

        self.d = disk
        self.c = calc

        self.progress = progress

        self.res_dict = {}

        id = 0
        for line_file, line_level in zip(self.c.line_files, self.c.line_levels):
            
            self.lp = readdata.molecdata(line_file, line_level)

            self.res_dict[id] = {'transition':'{} J={}-{}'.format(self.lp['name'], int(self.lp['Ju']), int(self.lp['Jl'])),
                                 self.c.output: self.integrate()}

        return None

    def background(self):
        "background query/monitoring thread"
        time.sleep(0.2)  # wait a bit for foreground code to start
        
        total = core.core.max_prog
        #pbar = tqdm(total = total)

        while self.bg_runn:

            now = core.core.progress 

            if now >= total: 
                print("\r10/10 {}".format("🏃‍♂️"*int(10)), end="")
                break
            
            print("\r{}/10 {}".format(int(now/total*10), '🏃‍♂️'*int(now/total*10) + '  '*int(10-now/total*10+1)), end="")
            time.sleep(0.5)

        return  



    def integrate(self):
        '''
        calclate a image cube
        '''
   

        self.gamma_h2 = 0.0711
        self.r_h2 = 0.66
        self.gamma_he = 0.0472
        self.r_he = 0.51
        self.tref = 296

        

        if self.progress:
            self.bg_runn = True
            thrd = threading.Thread(target=self.background)
            thrd.start()

        #print(self.lp['molec_weight'], self.lp['nu0'], self.lp['Jl'], self.lp['Aul'])

        I = core.core.integrate( self.c.V, self.c.X, self.c.Y, self.c.dX, self.c.Rmin, self.c.subpix_max, self.c.subpix_pow, self.c.pa, self.c.incl, self.c.z2r, self.lp['molec_weight'], self.lp['nu0'], self.lp['Jl'], self.lp['Aul'], self.c.nz, self.d.Sigma_g0, self.d.gamma_g, self.d.r0, self.d.rexp, self.d.rin, self.d.rout, self.d.rcav, self.d.dcav,self.d.d2g, self.d.k_d, self.d.Tmid0, self.d.pmid, self.d.Tatm0, self.d.patm, mu_gas, self.d.M_star, self.d.Xco, self.d.vturb, self.gamma_h2, self.r_h2, self.gamma_he, self.r_he, self.tref, self.c.dust, self.d.chi_d, self.c.contsub, self.c.voigt, self.c.kind, self.c.D, self.c.nthreads)
        
        if self.progress:
            self.bg_runn = False
            thrd.join()

        return I.T


    def vis_cube(self, id, u, v):

        # vis = joblib.Parallel(n_jobs=-1, verbose=50)(
        #    joblib.delayed(self.vis_sample)(self.res_dict[id][self.c.output][nch], self.c.dX/au/self.c.D, u, v) for nch in range(len(self.c.V))
        #    )

        vis = np.zeros((len(self.c.V), len(u)), dtype='complex')
        dxy = np.deg2rad(self.c.dX/au/self.c.D/3600)
        fac = dxy**2 * 1e23

        #print(self.res_dict[0][self.c.output][10])

        for nch in range(len(self.c.V)):
            vis[nch] = sampleImage(self.res_dict[id][self.c.output][nch] * fac, dxy, u, v, origin='upper')

        return vis

    def show_channelmap(self, id, nch, clim=None, unit = 'cgs'):

        if unit == 'cgs':
            fac = 1.0
        elif unit == 'jy/pix':
            dxy = np.deg2rad(self.c.dX/au/self.c.D/3600)
            fac = dxy**2 * 1e23
        
        plt.pcolormesh(self.c.XX/au, self.c.YY/au, np.flipud(np.fliplr(self.res_dict[id][self.c.output][nch]*fac)), shading='auto')

        if clim != None:
             plt.clim(clim)
        plt.gca().set_aspect('equal')

        plt.title('v = {:.2f} km/s'.format(self.c.V[nch]/1e5 ))
        plt.xlabel('au')
        plt.ylabel('au')
        plt.colorbar(label=unit)
        plt.show()

    def itp_cube(self, _img,  points, Xnew, Ynew, sigma_smo, ch):

        _img_itp = LinearNDInterpolator( points, np.ravel(_img), fill_value=0.0)
        _img_new = _img_itp(Xnew, Ynew)

        _img_new_smo = ndimage.gaussian_filter(_img_new, sigma=sigma_smo)

        return _img_new_smo

    
    def channel(self, id, nch):
        #dxy = np.deg2rad(self.c.dX/au/self.c.D/3600)
        #fac = dxy**2 * 1e23
        _img = self.res_dict[id][self.c.output][nch] #* fac
        return self.c.XX/au/self.c.D, self.c.YY/au/self.c.D, _img

    def cube(self, id):
        #dxy = np.deg2rad(self.c.dX/au/self.c.D/3600)
        #fac = dxy**2 * 1e23
        _cube = self.res_dict[id][self.c.output] #* fac
        return self.c.XX/au/self.c.D, self.c.YY/au/self.c.D, _cube

    def show_spectrum(self, id, scale='linear', mode='F', pos = None ):
        dxy = np.deg2rad(self.c.dX/au/self.c.D/3600)
        fac = dxy**2 * 1e23

        if mode == 'F':
            F = np.sum( self.res_dict[id][self.c.output] * fac , axis=(1,2))  # to Jy

            plt.plot(self.c.V/1e5, F)
            plt.xlabel('v km/s')
            plt.ylabel('Jy')
        
        elif mode == 'I':
            F = self.res_dict[id][self.c.output][:, pos[0], pos[1]]
            plt.plot(self.c.V/1e5, F)
            plt.xlabel('v km/s')
            plt.ylabel('intensity')

        plt.yscale(scale)

        return self.c.V/1e5, F
    
    

    '''
    def Fsca(self, tau_d, omega_d, mu):
        eps = np.sqrt(1-omega_d)
        t1 = -np.sqrt(3) * eps * tau_d

        return 1/( np.exp(t1)*(eps-1) - (eps+1) ) * (  (1-np.exp(-(np.sqrt(3)*eps+1/mu)*tau_d) )/(np.sqrt(3)*eps*mu+1) + (np.exp(-tau_d/mu)-np.exp(t1))/(np.sqrt(3)*eps*mu-1)  )

    def I_dust(self, nu, T, tau_d, omega_d, i ):
        mu = np.cos(i)
        return functions.B(nu, T) * ( (1-np.exp(-tau_d/mu)) + omega_d * self.Fsca(tau_d, omega_d, mu) )
    '''
