import os
import six
import datetime
from posidonius.particles.axes import Axes
from posidonius.integrator import WHFast, Ias15, LeapFrog
from posidonius.constants import *
from posidonius.effects.evolution import NonEvolving, Leconte2011, Baraffe2015, Baraffe1998, LeconteChabrier2013, BolmontMathis2016, GalletBolmont2017
from posidonius.tools import calculate_spin, mass_radius_relation, get_center_of_mass_of_pair
import posidonius.effects as effects
from posidonius.particles.particle import Particle, DummyParticle

class ConsiderEffects(object):
    def __init__(self, input_properties):
        self._data = {
            "tides": False,
            "rotational_flattening": False,
            "general_relativity": False,
            "disk": False,
            "wind": False,
            "evolution": False,
        }
        # Update default values, ignore non-recognised keys
        for key, value in six.iteritems(input_properties):
            if key in self._data:
                self._data[key] = value

    def get(self):
        if type(self._data) == str:
            return self._data
        else:
            return self._data.copy()


class Universe(object):
    def __init__(self, initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects):
        self._time_step = time_step
        self._recovery_snapshot_period = recovery_snapshot_period
        self._historic_snapshot_period = historic_snapshot_period
        self._data = {
            "consider_effects": consider_effects.get(),
            "general_relativity_implementation": "Disabled",
            "hosts": {
                "index": {
                    "most_massive": MAX_PARTICLES+1,
                    "tides": MAX_PARTICLES+1,
                    "rotational_flattening": MAX_PARTICLES+1,
                    "general_relativity": MAX_PARTICLES+1,
                    "disk": MAX_PARTICLES+1,
                },
                "most_massive": {
                    "all": False,
                    "tides": False,
                    "rotational_flattening": False,
                    "general_relativity": False,
                    "disk": False,
                }
            },
            "initial_time": float(initial_time),
            "n_particles": 0,
            "particles": [],
            "particles_evolvers": [],
            "star_planet_dependent_dissipation_factors": {},
            "temporary_copied_particle_positions": [],
            "temporary_copied_particles_masses": [],
            "temporary_copied_particles_radiuses": [],
            "time_limit": float(time_limit),
        }


    def add_dummy_particle(self):
        self.add_particle(DummyParticle())
        self._data['n_particles'] -= 1 # Compensate the addition from the previous add_particle call

    def add_particle(self, particle):
        if self._data['n_particles'] > MAX_PARTICLES:
            raise Exception("Maximum number of particles reached: {}".format(MAX_PARTICLES))

        if effects.tides.CentralBody in particle.effects():
            if self._data["consider_effects"]["tides"]:
                if self._data['hosts']['index']['tides'] == MAX_PARTICLES+1:
                    self._data['hosts']['index']['tides'] = self._data['n_particles']
                else:
                    raise Exception("There can only be one central body for tidal effects!")
            else:
                print("[WARNING {} UTC] Added a particle with tidal effect (central body) but the tidal effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if not self._data["consider_effects"]["tides"] and effects.tides.OrbitingBody in particle.effects():
            print("[WARNING {} UTC] Added a particle with tidal effect (orbiting body) but the tidal effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))

        if effects.rotational_flattening.CentralBody in particle.effects():
            if self._data["consider_effects"]["rotational_flattening"]:
                if self._data['hosts']['index']['rotational_flattening'] == MAX_PARTICLES+1:
                    self._data['hosts']['index']['rotational_flattening'] = self._data['n_particles']
                else:
                    raise Exception("There can only be one central body for rotational flattening effects!")
            else:
                print("[WARNING {} UTC] Added a particle with rotational flattening effect (central body) but the rotational flattening effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if not self._data["consider_effects"]["rotational_flattening"] and effects.rotational_flattening.OrbitingBody in particle.effects():
            print("[WARNING {} UTC] Added a particle with rotational flattening effect (orbiting body) but the rotational flattening effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))

        if effects.general_relativity.CentralBody in particle.effects():
            if self._data["consider_effects"]["general_relativity"]:
                if self._data['hosts']['index']['general_relativity'] == MAX_PARTICLES+1:
                    self._data['hosts']['index']['general_relativity'] = self._data['n_particles']
                else:
                    raise Exception("There can only be one central body for rotational flattening effects!")
                self._data['general_relativity_implementation'] = particle.general_relativity_implementation()
            else:
                print("[WARNING {} UTC] Added a particle with general relativity effect (central body) but the general relativity effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if not self._data["consider_effects"]["general_relativity"] and effects.general_relativity.OrbitingBody in particle.effects():
            print("[WARNING {} UTC] Added a particle with general relativity effect (orbiting body) but the general relativity effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))

        if effects.disk.CentralBody in particle.effects():
            if self._data["consider_effects"]["disk"]:
                if self._data['hosts']['index']['disk'] == MAX_PARTICLES+1:
                    self._data['hosts']['index']['disk'] = self._data['n_particles']
                else:
                    raise Exception("Only one body with a disk is allowed!")
            else:
                print("[WARNING {} UTC] Added a particle with disk effect (central body) but the disk effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if not self._data["consider_effects"]["disk"] and effects.disk.OrbitingBody in particle.effects():
            print("[WARNING {} UTC] Added a particle with disk effect (orbiting body) but the disk effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))

        if effects.wind.Interaction in particle.effects():
            if not self._data["consider_effects"]["wind"]:
                print("[WARNING {} UTC] Added a particle with wind effect but the wind effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))

        if particle.evolution() != effects.evolution.NonEvolving:
            if not self._data["consider_effects"]["evolution"]:
                print("[WARNING {} UTC] Added a particle with evolution effect but the evolution effect is disabled for this simulation".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))

        evolver = particle.get_evolver(self._data['initial_time'])
        if len(evolver['time']) > 0 and evolver['time'][0] > 0.:
            raise Exception("Your initial time ({} days | {:.2e} years) is smaller than the minimum allowed age of the star ({} days | {:.2e} years)".format(self._data['initial_time'], self._data['initial_time']/365.25, evolver['time'][0]+self._data['initial_time'], (evolver['time'][0]+self._data['initial_time'])/365.25));
        if len(evolver['time']) > 0 and evolver['time'][-1] < self._data['initial_time']:
            raise Exception("Your time limit ({} days | {:.2e} years) is greater than the maximum allowed age of the star ({} days | {:.2e} years)", self._data['initial_time'], self._data['initial_time']/365.25, evolver['time'][0], evolver['time'][0]/365.25)

        self._data['particles'].append(particle.get())
        self._data['particles_evolvers'].append(evolver)
        self._data['temporary_copied_particles_radiuses'].append(0.0)
        self._data['temporary_copied_particles_masses'].append(0.0)
        self._data['temporary_copied_particle_positions'].append(Axes(0.0, 0.0, 0.0).get())
        self._data['n_particles'] += 1


    def populate_inertial_frame(self):
        # Compute center of mass
        center_of_mass_position = Axes(0., 0., 0.)
        center_of_mass_velocity = Axes(0., 0., 0.)
        center_of_mass_mass = 0.;

        for particle in self._data['particles']:
            center_of_mass_mass = get_center_of_mass_of_pair(center_of_mass_position,
                                                                    center_of_mass_velocity,
                                                                    center_of_mass_mass,
                                                                    particle)

        for particle in self._data['particles']:
            particle['inertial_position']['x'] = particle['heliocentric_position']['x'] - center_of_mass_position.x()
            particle['inertial_position']['y'] = particle['heliocentric_position']['y'] - center_of_mass_position.y()
            particle['inertial_position']['z'] = particle['heliocentric_position']['z'] - center_of_mass_position.z()
            particle['inertial_velocity']['x'] = particle['heliocentric_velocity']['x'] - center_of_mass_velocity.x()
            particle['inertial_velocity']['y'] = particle['heliocentric_velocity']['y'] - center_of_mass_velocity.y()
            particle['inertial_velocity']['z'] = particle['heliocentric_velocity']['z'] - center_of_mass_velocity.z()



    def assign_id(self):
        for i, particle in enumerate(self._data['particles']):
            particle['id'] = i

    def find_most_massive_particle(self):
        max_mass_found = 0.
        for i, particle in enumerate(self._data['particles']):
            if particle['mass'] > max_mass_found:
                max_mass_found = particle['mass']
                self._data['hosts']['index']['most_massive'] = i
        if self._data['hosts']['index']['general_relativity'] != MAX_PARTICLES+1 and self._data['hosts']['index']['general_relativity'] != self._data['hosts']['index']['most_massive']:
                raise Exception("The most massive body should be the central body for general relativity effects!")

    def check_effects_vs_central_and_orbiting(self):
        found_tides_central_body = False
        found_tides_orbiting_body = False
        found_rotational_flattening_central_body = False
        found_rotational_flattening_orbiting_body = False
        found_general_relativity_central_body = False
        found_general_relativity_orbiting_body = False
        found_disk_central_body = False
        found_disk_orbiting_body = False
        found_wind = False
        found_evolution = False
        for i, particle in enumerate(self._data['particles']):
            if "CentralBody" in particle['tides']['effect']:
                found_tides_central_body = True
            elif "OrbitingBody" in particle['tides']['effect']:
                found_tides_orbiting_body = True
            elif "CentralBody" in particle['rotational_flattening']['effect']:
                found_rotational_flattening_central_body = True
            elif "OrbitingBody" in particle['rotational_flattening']['effect']:
                found_rotational_flattening_orbiting_body = True
            elif "CentralBody" in particle['general_relativity']['effect']:
                found_general_relativity_central_body = True
            elif "OrbitingBody" in particle['general_relativity']['effect']:
                found_general_relativity_orbiting_body = True
            elif "CentralBody" in particle['disk']['effect']:
                found_disk_central_body = True
            elif "OrbitingBody" in particle['disk']['effect']:
                found_disk_orbiting_body = True
            elif "Interaction" in particle['wind']['effect']:
                found_wind = True
            elif "NonEvolving" in particle['evolution']:
                found_evolution = True
        if self._data["consider_effects"]["tides"]:
            if not found_tides_central_body:
                print("[INFO {} UTC] No central body for tidal effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            if not found_tides_orbiting_body:
                print("[INFO {} UTC] No orbiting body for tidal effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if self._data["consider_effects"]["rotational_flattening"]:
            if not found_rotational_flattening_central_body:
                print("[INFO {} UTC] No central body for rotational flattening effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            if not found_rotational_flattening_orbiting_body:
                print("[INFO {} UTC] No orbiting body for rotational flattening effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if self._data["consider_effects"]["general_relativity"]:
            if not found_general_relativity_central_body:
                print("[INFO {} UTC] No central body for general relativity effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            if not found_general_relativity_orbiting_body:
                print("[INFO {} UTC] No orbiting body for general relativity effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if self._data["consider_effects"]["disk"]:
            if not found_disk_central_body:
                print("[INFO {} UTC] No central body for disk effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            if not found_disk_orbiting_body:
                print("[INFO {} UTC] No orbiting body for disk effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if self._data["consider_effects"]["wind"]:
            if not found_wind:
                print("[INFO {} UTC] No wind effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        if self._data["consider_effects"]["evolution"]:
            if not found_evolution:
                print("[INFO {} UTC] No evolution effects".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))

    def disable_unnecessary_effects(self):
        if self._data["consider_effects"]["tides"] and self._data['hosts']['index']['tides'] == MAX_PARTICLES+1:
            print("[INFO {} UTC] Disabled tides because no central host was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["tides"] = False
        if self._data["consider_effects"]["rotational_flattening"] and self._data['hosts']['index']['rotational_flattening'] == MAX_PARTICLES+1:
            print("[INFO {} UTC] Disabled rotational flattening because no central host was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["rotational_flattening"] = False
        if self._data["consider_effects"]["general_relativity"] and self._data['hosts']['index']['general_relativity'] == MAX_PARTICLES+1:
            print("[INFO {} UTC] Disabled general relativity because no central host was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["general_relativity"] = False
        if self._data["consider_effects"]["disk"] and self._data['hosts']['index']['disk'] == MAX_PARTICLES+1:
            print("[INFO {} UTC] Disabled disk because no central host was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["disk"] = False
        found_evolving_body = False
        found_wind = False
        for particle in self._data['particles']:
            if particle["evolution"] != "NonEvolving":
                found_evolving_body = True
            if particle["wind"]["effect"] != "Disabled":
                found_wind = True
        if self._data["consider_effects"]["wind"] and not found_wind:
            print("[INFO {} UTC] Disabled wind because no wind was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["wind"] = False
        if self._data["consider_effects"]["evolution"] and not found_evolving_body:
            print("[INFO {} UTC] Disabled evolution because no evolving body was included".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
            self._data["consider_effects"]["evolution"] = False


    def update_hosts(self):
        self._data["hosts"]["most_massive"] = {
            "all": False,
            "tides": self._data["consider_effects"]["tides"] and self._data['hosts']['index']['most_massive'] == self._data['hosts']['index']['tides'],
            "rotational_flattening": self._data["consider_effects"]["rotational_flattening"] and self._data['hosts']['index']['most_massive'] == self._data['hosts']['index']['rotational_flattening'],
            "general_relativity": self._data["consider_effects"]["general_relativity"] and self._data['hosts']['index']['most_massive'] == self._data['hosts']['index']['general_relativity'],
            "disk": self._data["consider_effects"]["disk"] and self._data['hosts']['index']['most_massive'] == self._data['hosts']['index']['disk'],
        }
        self._data["hosts"]["most_massive"]["all"] = ((self._data["consider_effects"]["tides"] and self._data["hosts"]["most_massive"]["tides"]) or not self._data["consider_effects"]["tides"]) \
            and ((self._data["consider_effects"]["rotational_flattening"] and self._data["hosts"]["most_massive"]["rotational_flattening"]) or not self._data["consider_effects"]["rotational_flattening"]) \
            and ((self._data["consider_effects"]["general_relativity"] and self._data["hosts"]["most_massive"]["general_relativity"]) or not self._data["consider_effects"]["general_relativity"]) \
            and ((self._data["consider_effects"]["disk"] and self._data["hosts"]["most_massive"]["disk"]) or not self._data["consider_effects"]["disk"])

    def compute_general_relativity_factor(self):
        most_massive_particle_index = self._data['hosts']['index']['most_massive']
        central_body_mass_g = self._data['particles'][most_massive_particle_index]['mass_g']
        for i, particle in enumerate(self._data['particles']):
            if not self._data["consider_effects"]["general_relativity"] \
                    or "CentralBody" in particle['general_relativity']['effect'] \
                    or particle['general_relativity']['effect'] == "Disabled" \
                    or particle['mass'] == 0:
                # Disabled GR or Enabled GR and central body
                particle['general_relativity']['parameters']['internal']['factor'] = 0.0
            else:
                particle['general_relativity']['parameters']['internal']['factor'] = central_body_mass_g*particle['mass_g'] / np.power(central_body_mass_g + particle['mass_g'], 2)

    def get(self):
        self.disable_unnecessary_effects()
        self.check_effects_vs_central_and_orbiting()
        self.populate_inertial_frame()
        self.assign_id()
        self.find_most_massive_particle()
        self.update_hosts()
        self.compute_general_relativity_factor()

        # Add dummy particles to fill the vector
        n_dummy_particles = MAX_PARTICLES - self._data['n_particles']
        for i in range(n_dummy_particles):
            self.add_dummy_particle()

        data = self._data.copy()

        # Reset indices according to enabled effects
        if data["consider_effects"]["rotational_flattening"] and not data["consider_effects"]["tides"]:
            # In practice, rotational_flattening will use tidal_host_particle_index
            # make sure it is the correct index, even if tides are disabled
            data['hosts']['index']['tides'] = data['hosts']['index']['rotational_flattening']
        elif not data["consider_effects"]["tides"]:
            data['hosts']['index']['tides'] = MAX_PARTICLES+1
        if not data["consider_effects"]["rotational_flattening"]:
            data['hosts']['index']['rotational_flattening'] = MAX_PARTICLES+1
        if not data["consider_effects"]["general_relativity"]:
            data['hosts']['index']['general_relativity'] = MAX_PARTICLES+1
        if not data["consider_effects"]["disk"]:
            data['hosts']['index']['disk'] = MAX_PARTICLES+1

        # Forget the dummy particles
        if n_dummy_particles > 0:
            self._data['particles'] = self._data['particles'][:-n_dummy_particles]
            self._data['particles_evolvers'] = self._data['particles_evolvers'][:-n_dummy_particles]
            #self._data['n_particles'] -= n_dummy_particles
        return data

    def write(self, filename, integrator="WHFast", whfast_alternative_coordinates="DemocraticHeliocentric"):
        if integrator.lower() == "whfast":
            universe_integrator = WHFast(whfast_alternative_coordinates, self._time_step, self._recovery_snapshot_period, self._historic_snapshot_period, self)
            universe_integrator.write(filename)
        elif integrator.lower() == "ias15":
            universe_integrator = Ias15(self._time_step, self._recovery_snapshot_period, self._historic_snapshot_period, self)
            universe_integrator.write(filename)
        elif integrator.lower() == "leapfrog":
            universe_integrator = LeapFrog(self._time_step, self._recovery_snapshot_period, self._historic_snapshot_period, self)
            universe_integrator.write(filename)
        else:
            raise Exception("Unknown integtrator '{}'".format(integrator))
        base_filename = os.path.splitext(filename)[0]
        print("[INFO {} UTC] Start the simulation with:".format(datetime.datetime.utcnow().strftime("%Y.%m.%d %H:%M:%S")))
        print("posidonius start {} {} {}".format(filename, base_filename+".bin", base_filename+"_history.bin"))


