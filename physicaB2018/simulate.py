import numpy as np
from datetime import datetime
import bornagain as ba
from bornagain import deg, nm

# SLDs
sld_Si = 2.074e-6
sld_Si_im = 0.0     # -2.3819e-11
sld_D2O = 6.356e-6
sld_D2O_im = -1.1295e-13

# microgel particle parameters
b = 9.0e+07
xi = 9.0
xiz = 2.3

suf = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
sim_datafile = "sim{}.txt".format(suf)
exp_datafile = "data/ME3O5/I-SM-1181-KWS2.txt"

# particle dimensions
core_radius = 53.0*nm
core_height = 53.0*nm
radius = 135.0*nm
height = 97.0*nm

# KWS-2 detector parameters
npx, npy = 128, 128   # number of detector pixels
det_width, det_height = 675.0, 675.0  # mm, detector size
psize = det_width / npx  # pixel size, mm
sdd = 3710.0     # mm, sample-detector distance

# direct beam position
beam_xpos, beam_ypos = 64.0, 64.5       # pixel

# incident angle
ai = 0.73   # degree
wavelength = 5.0  # angstrom

# beam
beam_intensity = 500.0

# integration over the bin
mc_integration=False


class FormFactorMicrogel(ba.IFormFactorBorn):
    """
    A custom defined form factor for a MicroGel particle.
    """
    def __init__(self, b, xi, xiz):
        ba.IFormFactorBorn.__init__(self)
        # parameters describing the form factor
        self.b = b          # scaling parameter
        self.xi = xi        # correlation length in XY-plane
        self.xiz = xiz      # correlation length in Z direction

    def clone(self):
        """
        IMPORTANT NOTE:
        The clone method needs to call transferToCPP() on the cloned object
        to transfer the ownership of the clone to the cpp code
        """
        cloned_ff = FormFactorMicrogel(self.b, self.xi, self.xiz)
        cloned_ff.transferToCPP()
        return cloned_ff

    def evaluate_for_q(self, q):
        return np.sqrt(self.b/(1.0 + self.xi*self.xi*(q.y()*q.y() + q.x()*q.x()) + self.xiz*self.xiz*q.z()*q.z()))


def vol(r, h):
    """
    Calculate volume of the truncated sphere
    :param r: radius
    :param h: height
    :return: volume
    """
    v = 2.0/3.0 + (h-r)/r - ((h-r)/r)**3/3.0
    return v*np.pi*r*r*r


def get_sample():
    """
    Returns a sample
    """
    # defining materials
    m_si = ba.MaterialBySLD("Si", sld_Si, sld_Si_im)
    m_d2o = ba.MaterialBySLD("D2O", sld_D2O, sld_D2O_im)
    m_core = ba.MaterialBySLD("Me3O5:D2O2", 2.0*1.0e-06, 0.0)
    m_shell = ba.MaterialBySLD("Me3O5:D2O", 3.9*1.0e-06, 0.0)

    # layer with particles
    # calculate average SLD
    Vcore = vol(core_radius, core_height)
    Vshell = vol(radius, height) - Vcore
    f_d2o = 0.7
    f_core = (1.0 - f_d2o)/(1 + Vshell/Vcore)
    f_shell = (1.0 - f_d2o)/(1 + Vcore/Vshell)
    sld_mix = f_d2o*sld_D2O + f_shell*3.9*1.0e-06 + f_core*2.0*1.0e-06
    m_mix = ba.MaterialBySLD("mix", sld_mix, 0.0)

    # fluctuation component
    ff_microgel = FormFactorMicrogel(b, xi, xiz)
    microgel = ba.Particle(m_core, ff_microgel)
    microgel_layout = ba.ParticleLayout()
    microgel_layout.addParticle(microgel, 1.0)

    # collection of particles
    ff = ba.FormFactorTruncatedSphere(radius=radius, height=height)
    ff_core = ba.FormFactorTruncatedSphere(radius=core_radius, height=core_height)
    transform = ba.RotationY(180.0 * deg)
    shell_particle = ba.Particle(m_shell, ff)
    core_particle = ba.Particle(m_core, ff_core)
    core_position = ba.kvector_t(0.0, 0.0, 0.0)
    particle = ba.ParticleCoreShell(shell_particle, core_particle, core_position)
    particle.setPosition(ba.kvector_t(0.0, 0.0, 0.0))
    particle.setRotation(transform)

    nparticles = 2  # the larger is this number, the more slow will be the simulation. 10 is usually enough
    sigma = 0.2*radius

    gauss_distr = ba.DistributionGaussian(radius, sigma)

    sigma_factor = 2.0
    par_distr = ba.ParameterDistribution(
        "/ParticleCoreShell/Particle1/TruncatedSphere/Radius", gauss_distr, nparticles, sigma_factor,
        ba.RealLimits.lowerLimited(core_radius + 1.0))
    par_distr.linkParameter("/ParticleCoreShell/Particle1/TruncatedSphere/Height")
    par_distr.linkParameter("/ParticleCoreShell/Particle0/TruncatedSphere/Height")
    par_distr.linkParameter("/ParticleCoreShell/Particle0/TruncatedSphere/Radius")
    part_coll = ba.ParticleDistribution(particle, par_distr)

    microgel_layout.addParticle(part_coll, 1.2e-05)

    # interference can be neglected
    interference = ba.InterferenceFunctionNone()
    microgel_layout.setInterferenceFunction(interference)

    # describe layer roughness
    roughness = ba.LayerRoughness()
    roughness.setSigma(1.2 * ba.nm)
    roughness.setHurstParameter(0.8)
    roughness.setLatteralCorrLength(570.0 * ba.nm)

    # create layers
    d2o_layer = ba.Layer(m_d2o)
    mix_layer = ba.Layer(m_mix, 2.0*height)
    mix_layer.addLayout(microgel_layout)
    si_layer = ba.Layer(m_si)
    multi_layer = ba.MultiLayer()
    multi_layer.addLayer(si_layer)
    multi_layer.addLayer(mix_layer)
    multi_layer.addLayerWithTopRoughness(d2o_layer, roughness)

    return multi_layer


def create_detector():
    """
    Creates and returns KWS-2 detector
    """
    u0 = beam_xpos*psize  # in mm
    v0 = beam_ypos*psize  # in mm
    detector = ba.RectangularDetector(npx, det_width, npy, det_height)
    detector.setPerpendicularToDirectBeam(sdd, u0, v0)
    return detector


def get_simulation(wl=5.0, alpha_i=ai):
    """
    Returns a GISAXS simulation with beam and detector defined
    """
    simulation = ba.GISASSimulation()
    simulation.setBeamParameters(wl*ba.angstrom, alpha_i*ba.deg, 0.0*ba.deg)
    simulation.setDetector(create_detector())
    simulation.setBeamIntensity(beam_intensity)
    simulation.getOptions().setMonteCarloIntegration(mc_integration, 50)
    return simulation


def run_simulation():
    """
    Runs simulation and returns resulting intensity map.
    """
    sample = get_sample()
    simulation = get_simulation(wavelength, ai)
    simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(1.5*psize, 0.7*psize))
    simulation.setSample(sample)
    # set region of interest in mm
    simulation.setRegionOfInterest(150.0, 390.0, 520.0, 650.0)
    # uncomment lines below to add the beam divergence
    # beware: this will increase the simulation time a lot!
    # wavelength_distr = ba.DistributionLogNormal(wavelength*angstrom, 0.2*wavelength*angstrom)
    # alpha_distr = ba.DistributionGaussian(ai * deg, 0.05 * deg)
    # phi_distr = ba.DistributionGaussian(0.0 * deg, 0.15 * deg)
    # simulation.addParameterDistribution("*/Beam/Wavelength", wavelength_distr, 5)
    # simulation.addParameterDistribution("*/Beam/InclinationAngle", alpha_distr, 3)
    # simulation.addParameterDistribution("*/Beam/AzimuthalAngle", phi_distr, 5)
    # options
    # simulation.getOptions().setUseAvgMaterials(True)   # does not work for the microgel FF in 1.11
    simulation.getOptions().setIncludeSpecular(True)    # include specular peak
    simulation.getOptions().setNumberOfThreads(-1)      # custom FF can be calculaten only in a single thread
    simulation.setTerminalProgressMonitor()             # show progress

    simulation.runSimulation()                          # run simulation
    return simulation.result()                          # return result


if __name__ == '__main__':
    result = run_simulation()
    # to save the simulated matrix to a file, uncomment the line below
    # np.savetxt(sim_datafile, result.histogram2d(ba.AxesUnits.QSPACE).getArray())
    ba.plot_simulation_result(result, units=ba.AxesUnits.QSPACE, intensity_min=1.0e-05, intensity_max=0.2)
