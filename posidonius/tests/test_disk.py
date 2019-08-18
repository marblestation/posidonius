import os
import inspect
import shutil
import filecmp
import posidonius
from posidonius.tests import common
from posidonius.tests.test_base import TestBase

class Disk(TestBase):

    def setUp(self):
        TestBase.setUp(self)
        self.current_filename, ignore = os.path.splitext(os.path.basename(__file__)) # Filename without extension
        self.current_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))

    def tearDown(self):
        TestBase.tearDown(self)

    def test_enabled_disk(self):
        current_function_name = inspect.stack()[0][3][5:] # Remove first 5 characters "test_"
        expected_json_filename, json_filename = common.setup(self.current_dirname, self.current_filename, current_function_name)

        initial_time, time_step, time_limit, historic_snapshot_period, recovery_snapshot_period = common.simulation_properties()
        consider_effects = posidonius.ConsiderEffects({
            "tides": True,
            "rotational_flattening": True,
            "general_relativity": True,
            "disk": True,
            "wind": True,
            "evolution": True,
        })
        general_relativity_implementation = "Kidder1995" # Assumes one central massive body
        #general_relativity_implementation = "Anderson1975" # Assumes one central massive body
        #general_relativity_implementation = "Newhall1983" # Considers all bodies
        universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

        star_mass = 1.0 # Solar masses
        star_rotation_period = 8. # hours
        star_dissipation_factor_scale = 1.0
        star_position = posidonius.Axes(0., 0., 0.)
        star_velocity = posidonius.Axes(0., 0., 0.)

        star_evolution = posidonius.Baraffe1998(star_mass) # M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
        #star_evolution = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
        #star_evolution = posidonius.Leconte2011(star_mass) # BrownDwarf (mass = 0.01 .. 0.08)
        #star_evolution = posidonius.BolmontMathis2016(star_mass) # SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
        #star_evolution = posidonius.GalletBolmont2017(star_mass) # SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
        #star_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
        #star_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
        #star_evolution = posidonius.NonEvolving()
        star = common.solar_like_with_disk(star_mass, star_dissipation_factor_scale, star_position, star_velocity, star_rotation_period, general_relativity_implementation, star_evolution, wind_k_factor=4.0e-18, wind_rotation_saturation=1.7592918860102842)
        universe.add_particle(star)
        common.basic_configuration(universe)

        ############################################################################
        whfast_alternative_coordinates="Jacobi"
        #whfast_alternative_coordinates="DemocraticHeliocentric"
        #whfast_alternative_coordinates="WHDS"
        #universe.write(json_filename, integrator="LeapFrog")
        #universe.write(json_filename, integrator="IAS15")
        universe.write(json_filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)

        if not os.path.exists(expected_json_filename):
            shutil.copyfile(json_filename, expected_json_filename)

        self.assertTrue(filecmp.cmp(json_filename, expected_json_filename, shallow=False), "Generated JSON case is not equal to the expected one: {} {}".format(json_filename, expected_json_filename))
        shutil.rmtree(os.path.dirname(json_filename))

    def test_disabled_disk(self):
        current_function_name = inspect.stack()[0][3][5:] # Remove first 5 characters "test_"
        expected_json_filename, json_filename = common.setup(self.current_dirname, self.current_filename, current_function_name)

        initial_time, time_step, time_limit, historic_snapshot_period, recovery_snapshot_period = common.simulation_properties()
        consider_effects = posidonius.ConsiderEffects({
            "tides": True,
            "rotational_flattening": True,
            "general_relativity": True,
            "disk": False,
            "wind": True,
            "evolution": True,
        })
        general_relativity_implementation = "Kidder1995" # Assumes one central massive body
        #general_relativity_implementation = "Anderson1975" # Assumes one central massive body
        #general_relativity_implementation = "Newhall1983" # Considers all bodies
        universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

        star_mass = 1.0 # Solar masses
        star_rotation_period = 8. # hours
        star_dissipation_factor_scale = 1.0
        star_position = posidonius.Axes(0., 0., 0.)
        star_velocity = posidonius.Axes(0., 0., 0.)

        star_evolution = posidonius.Baraffe1998(star_mass) # M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
        #star_evolution = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
        #star_evolution = posidonius.Leconte2011(star_mass) # BrownDwarf (mass = 0.01 .. 0.08)
        #star_evolution = posidonius.BolmontMathis2016(star_mass) # SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
        #star_evolution = posidonius.GalletBolmont2017(star_mass) # SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
        #star_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
        #star_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
        #star_evolution = posidonius.NonEvolving()
        star = common.solar_like(star_mass, star_dissipation_factor_scale, star_position, star_velocity, star_rotation_period, general_relativity_implementation, star_evolution, wind_k_factor=4.0e-18, wind_rotation_saturation=1.7592918860102842)
        universe.add_particle(star)
        common.basic_configuration(universe)

        ############################################################################
        whfast_alternative_coordinates="Jacobi"
        #whfast_alternative_coordinates="DemocraticHeliocentric"
        #whfast_alternative_coordinates="WHDS"
        #universe.write(json_filename, integrator="LeapFrog")
        #universe.write(json_filename, integrator="IAS15")
        universe.write(json_filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)

        if not os.path.exists(expected_json_filename):
            shutil.copyfile(json_filename, expected_json_filename)

        self.assertTrue(filecmp.cmp(json_filename, expected_json_filename, shallow=False), "Generated JSON case is not equal to the expected one: {} {}".format(json_filename, expected_json_filename))
        shutil.rmtree(os.path.dirname(json_filename))
