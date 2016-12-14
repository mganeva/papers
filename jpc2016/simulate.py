import numpy as np
import bornagain as ba
from datetime import datetime
from plot import plot_data

# cylinders
radius = 3.0 * ba.nanometer
length = 1.0 * ba.micrometer

# hexagonal ordering
distance = 7.2 * ba.nanometer   # 5.95 for FCC lattice
hex_angle = 60.0 * ba.degree      # 64.55 for the FCC lattice

# mesocrystal
nbr_layers = 30             # defines height
meso_width = 72.0*distance  # defines width

# detector parameters
sdd = 2400.0  # in mm
pixel_size = 200.0/240.0  # in mm
npx, npy = 240, 240  # number of pixels
beam_center_y, beam_center_z = 120.0, 103.8  # pixels

# detector extension to see the side peaks in absence of the angular beam divergence
extend_det = 100    # pixel

# beam parameters
ai = 0.82*ba.degree
pi = 0.0*ba.degree
wavelength = 4.28*ba.angstrom
# beam intensity represents the total amount of neutrons/photons = flux*time*irradiated_area
# it acts as a constant intensity scaling factor
beam_intensity = 1.7590e+10     # 8.0090e+7 for the FCC lattice
# phi_div = 5.0*ba.degree       # not considered
# wl_div = 0.01*wavelength      # not considered

# rotation distribution parameters
# absolute value of the Euler angle Gamma |Gamma| must be below 90 degree
# |Gamma| is in the range from 0 to 90 with the step of  6.02 degree
rz_range = 84.3     # 84.6 for the FCC lattice
rz_num = 28         # 30 for the FCC lattice

# output file name
suf = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
outfile = "sim" + suf + ".dat"


# ------------
# Create sample with rotational distribution of domain orientations
# ------------
def get_sample_rot_distr():
    # create materials
    substrate_material = ba.HomogeneousMaterial("Substrate", 6.05878e-6, 6.94321e-11)
    particle_material = ba.HomogeneousMaterial("ParticleCore", -0.609e-6, 0.183e-6)
    layer_material = ba.HomogeneousMaterial("D2O", 1.85762e-5, 3.31309e-11)

    # create mesocrystal
    # lattice
    lat = create_hex_lattice(distance, length)
    # basis particle
    ff = ba.FormFactorCylinder(radius, length)
    particle = ba.Particle(particle_material, ff)
    # rotate cylinder to make it parallel to the substrate
    rotation = ba.RotationY(90.0*ba.degree)
    particle.setRotation(rotation)
    basis = ba.ParticleComposition(particle)

    # mesocrystal
    total_height = nbr_layers*distance*np.sin(hex_angle)
    # using a Gaussian form factor as the overall shape of where the cylinders are
    meso_ff = ba.FormFactorGauss(meso_width, total_height)
    # assemble them into mesocrystal
    npc = ba.Crystal(basis, lat)
    dw_factor = 0.2
    npc.setDWFactor(dw_factor)
    meso = ba.MesoCrystal(npc, meso_ff)
    # rotate mesocrystal
    rot_y = ba.RotationY(180.0 * ba.degree)
    meso.setRotation(rot_y)          # turn upside down
    rot_z = ba.RotationZ(0.1 * ba.degree)
    meso.applyRotation(rot_z)        # rotate around Z
    meso.setPosition(0.0, 0.0, 0.0)

    # add uniform distribution of the domain orientations
    rot_distr = ba.DistributionGate((90.0-rz_range)*ba.degree, (90.0+rz_range)*ba.degree)
    ang_distr = ba.ParameterDistribution("*MesoCrystal/EulerRotation/Gamma", rot_distr, rz_num)
    part_coll = ba.ParticleDistribution(meso, ang_distr)

    # Create multilayer
    multi_layer = ba.MultiLayer()
    d2o_layer = ba.Layer(layer_material)
    substrate_layer = ba.Layer(substrate_material)

    particle_layout = ba.ParticleLayout()
    particle_layout.addParticle(part_coll)
    d2o_layer.addLayout(particle_layout)

    # the sample is upside down
    multi_layer.addLayer(substrate_layer)
    multi_layer.addLayer(d2o_layer)

    return multi_layer


# -------------------------------------------------------------------------
# create lattice
# -------------------------------------------------------------------------
def create_hex_lattice(lat_length, cyl_length):
    ca = np.cos(0.0)
    sa = np.sin(0.0)
    ca3 = 0.5*np.cos(0.0)
    sa3 = 0.5*np.tan(0.0 + hex_angle)
    a1 = ba.kvector_t(cyl_length, 0.0, 0.0)
    a2 = ba.kvector_t(0.0, lat_length*ca, lat_length*sa)
    a3 = ba.kvector_t(0.0, lat_length*ca3, lat_length*sa3)
    result = ba.Lattice(a1, a2, a3)
    return result


def create_detector(large=False):
    """
    Creates and returns NREX detector
    """
    if large:       # if simulated detector larger as the real one
        u0 = (beam_center_y + extend_det) * pixel_size  # in mm
        v0 = beam_center_z * pixel_size  # in mm
        detector = ba.RectangularDetector(npx + 2*extend_det, (npx + 2*extend_det) * pixel_size, npy, npy * pixel_size)
    else:
        u0 = beam_center_y*pixel_size  # in mm
        v0 = beam_center_z*pixel_size  # in mm
        detector = ba.RectangularDetector(npx, npx*pixel_size, npy, npy*pixel_size)
    normal = ba.kvector_t(sdd * np.cos(ai), 0.0, 1.0 * sdd * np.sin(ai))
    detector.setPosition(normal, u0, v0)
    # detector.setPerpendicularToReflectedBeam(sdd, u0, v0)      # alternatively

    return detector


def create_simulation(large=False):
    """
    Creates and returns GISANS simulation with beam and detector defined
    """
    simulation = ba.GISASSimulation()
    simulation.setDetector(create_detector(large))
    # the broader detector resolution is taken instead of the angular beam divergence
    simulation.setDetectorResolutionFunction(ba.ResolutionFunction2DGaussian(36.0 * pixel_size, 3.0 * pixel_size))
    simulation.setBeamParameters(wavelength, ai, pi)
    simulation.setBeamIntensity(beam_intensity)
    # beam divergence is not considered
    # wl_distr = ba.DistributionGaussian(wavelength, wl_div / 2.355)
    # simulation.addParameterDistribution("*/Beam/Wavelength", wl_distr, 3)
    # phi_distr = ba.DistributionGaussian(0.0 * ba.degree, phi_div / 2.355)
    # simulation.addParameterDistribution("*/Beam/AzimuthalAngle", phi_distr, 300)
    return simulation


def run_simulation():
    sample = get_sample_rot_distr()
    sim = create_simulation(large=True)
    sim.setSample(sample)
    sim.runSimulation()
    sim.printParameters()
    result = sim.getIntensityData(ba.IDetector2D.NBINS)
    return result.crop(extend_det, 0, npx + extend_det - 1, npy)


def simulate():
    """
    run simulation, save and plot results
    """
    simul_data = run_simulation()
    np.savetxt(outfile, simul_data.getArray(), delimiter='\t', newline='\n')

    # plot data
    # plot_settings = {'title': "Hexagonal Lattice", 'filename': outfile, 'f': 1.0}
    plot_settings = {'title': r"FCC Lattice with $\beta=64.55^{\circ}$", 'filename': outfile, 'f': 1.0}
    plot_data(plot_settings)


if __name__ == '__main__':
    # run simulation, save and plot the result
    simulate()
