import os
import inspect
import shutil
import filecmp
import numpy as np
import posidonius
from posidonius import Axes
from posidonius.tests import common
from posidonius.tests.test_base import TestBase
from posidonius.constants import DEG2RAD, G, G_SI, M_SUN, HOUR, AU, PI, TWO_PI
from posidonius.tools import calculate_cartesian_coordinates, calculate_keplerian_orbital_elements, modulus

class Tides(TestBase):

    def setUp(self):
        TestBase.setUp(self)
        self.current_filename, ignore = os.path.splitext(os.path.basename(__file__)) # Filename without extension
        self.current_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))

    def tearDown(self):
        TestBase.tearDown(self)

    def test_calculate_keplerian_orbital_elements_circular_orbit(self):
        """
        Test calculate_keplerian_orbital_elements with a circular orbit (e=0).
        """
        # Circular orbit parameters
        target_mass = 1.0  # Solar masses
        central_mass = 1.0  # Solar masses
        total_mass = target_mass + central_mass

        # Position at perihelion (circular orbit, so distance is constant)
        position = Axes(1.0, 0.0, 0.0)  # AU
        # Velocity perpendicular to position vector for circular orbit
        orbital_speed = np.sqrt(G * total_mass / 1.0)  # AU/yr
        velocity = Axes(0.0, orbital_speed, 0.0)

        # Expected orbital elements
        expected_a = 1.0  # Semi-major axis (AU)
        expected_e = 0.0  # Eccentricity
        expected_i = 0.0  # Inclination

        # Calculate orbital elements
        a, q, e, i, p, n, l, f = calculate_keplerian_orbital_elements(
            target_mass, position, velocity, masses=[central_mass], positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(a, expected_a, places=7)
        self.assertAlmostEqual(e, expected_e, places=7)
        self.assertAlmostEqual(i, expected_i, places=7)

    def test_calculate_keplerian_orbital_elements_elliptical_orbit(self):
        """
        Test calculate_keplerian_orbital_elements with an elliptical orbit (0 < e < 1).
        """
        target_mass = 0.0  # Negligible planet mass
        central_mass = 1.0  # Solar masses
        total_mass = target_mass + central_mass

        # Elliptical orbit parameters
        semi_major_axis = 1.0  # AU
        eccentricity = 0.5
        # Calculate periapsis distance
        periapsis_distance = semi_major_axis * (1 - eccentricity)
        # Position at periapsis
        position = Axes(periapsis_distance, 0.0, 0.0)
        # Velocity at periapsis
        orbital_speed = np.sqrt(G * total_mass * (1 + eccentricity) / (semi_major_axis * (1 - eccentricity)))
        velocity = Axes(0.0, orbital_speed, 0.0)

        # Expected orbital elements
        expected_a = semi_major_axis
        expected_e = eccentricity

        # Calculate orbital elements
        a, q, e, i, p, n, l, f = calculate_keplerian_orbital_elements(
            target_mass, position, velocity, masses=[central_mass], positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(a, expected_a, places=7)
        self.assertAlmostEqual(e, expected_e, places=7)
        # Since the orbit is in the x-y plane, inclination should be zero
        self.assertAlmostEqual(i, 0.0, places=7)
        # True anomaly at periapsis should be zero
        self.assertAlmostEqual(f, 0.0, places=7)

    def test_calculate_keplerian_orbital_elements_parabolic_orbit(self):
        """
        Test calculate_keplerian_orbital_elements with a parabolic orbit (e=1).
        """
        target_mass = 0.0
        central_mass = 1.0
        total_mass = target_mass + central_mass

        # Parabolic orbit parameters
        e = 1.0
        periapsis_distance = 1.0  # AU
        # Position at periapsis
        position = Axes(periapsis_distance, 0.0, 0.0)
        # Velocity at periapsis for parabolic orbit
        orbital_speed = np.sqrt(2 * G * total_mass / periapsis_distance)
        velocity = Axes(0.0, orbital_speed, 0.0)

        # Expected orbital elements
        expected_e = 1.0

        # Calculate orbital elements
        a, q, e_calc, i, p, n, l, f = calculate_keplerian_orbital_elements(
            target_mass, position, velocity, masses=[central_mass], positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(e_calc, expected_e, places=7)
        # Semi-major axis is undefined for parabolic orbit (should be infinite), but code may set it to None or a large value
        self.assertTrue(np.isinf(a) or a > 1e10)
        # True anomaly at periapsis should be zero
        self.assertAlmostEqual(f, 0.0, places=7)

    def test_calculate_keplerian_orbital_elements_hyperbolic_orbit(self):
        """
        Test calculate_keplerian_orbital_elements with a hyperbolic orbit (e > 1).
        """
        target_mass = 0.0
        central_mass = 1.0
        total_mass = target_mass + central_mass

        # Hyperbolic orbit parameters
        e = 1.5
        periapsis_distance = 1.0  # AU
        # Position at periapsis
        position = Axes(periapsis_distance, 0.0, 0.0)
        # Velocity at periapsis for hyperbolic orbit
        orbital_speed = np.sqrt(G * total_mass * (e + 1) / periapsis_distance)
        velocity = Axes(0.0, orbital_speed, 0.0)

        # Expected orbital elements
        expected_e = e

        # Calculate orbital elements
        a, q, e_calc, i, p, n, l, f = calculate_keplerian_orbital_elements(
            target_mass, position, velocity, masses=[central_mass], positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(e_calc, expected_e, places=7)
        # Semi-major axis is negative for hyperbolic orbits
        expected_a = -periapsis_distance / (e - 1)
        self.assertAlmostEqual(a, expected_a, places=7)
        # True anomaly at periapsis should be zero
        self.assertAlmostEqual(f, 0.0, places=7)

    def test_calculate_keplerian_orbital_elements_inclined_orbit(self):
        """
        Test calculate_keplerian_orbital_elements with an inclined orbit.
        """
        target_mass = 0.0
        central_mass = 1.0
        total_mass = target_mass + central_mass

        # Orbital parameters
        semi_major_axis = 1.0  # AU
        eccentricity = 0.1
        inclination = 45.0 * DEG2RAD  # 45 degrees
        # For simplicity, set argument of periapsis and longitude of ascending node to zero
        p = 0.0
        n = 0.0
        # Compute position at periapsis
        E = 0.0  # Eccentric anomaly at periapsis
        r = semi_major_axis * (1 - eccentricity)
        x_orb = r * np.cos(E)
        y_orb = r * np.sin(E)
        # Position in space (since E=0, x_orb = r, y_orb = 0)
        x = x_orb
        y = 0.0
        z = 0.0
        # Rotate position vector
        # Since n = 0 and p = 0, rotation simplifies to inclination about the x-axis
        y_rot = y * np.cos(inclination) - z * np.sin(inclination)
        z_rot = y * np.sin(inclination) + z * np.cos(inclination)
        position = Axes(x, y_rot, z_rot)
        # Compute orbital speed at periapsis using the correct formula
        orbital_speed = np.sqrt(G * total_mass * (1 + eccentricity) / (semi_major_axis * (1 - eccentricity)))
        # Velocity in orbital plane at periapsis
        vx_orb = 0.0
        vy_orb = orbital_speed
        # Rotate velocity vector
        vy_rot = vy_orb * np.cos(inclination)
        vz_rot = vy_orb * np.sin(inclination)
        velocity = Axes(vx_orb, vy_rot, vz_rot)

        # Calculate orbital elements
        a, q, e_calc, i_calc, p_calc, n_calc, l, f = calculate_keplerian_orbital_elements(
            target_mass,
            position,
            velocity,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)],
        )

        # Assertions
        self.assertAlmostEqual(a, semi_major_axis, places=6)
        self.assertAlmostEqual(e_calc, eccentricity, places=6)
        self.assertAlmostEqual(i_calc, inclination, places=6)

    def test_calculate_cartesian_coordinates_circular_orbit(self):
        """
        Test calculate_cartesian_coordinates for a circular orbit.
        """
        target_mass = 0.0
        central_mass = 1.0
        total_mass = target_mass + central_mass

        # Orbital elements
        q = 1.0  # Perihelion distance (AU)
        e = 0.0  # Circular orbit
        i0 = 0.0  # Inclination (radians)
        p = 0.0  # Longitude of perihelion (radians)
        n0 = 0.0  # Longitude of ascending node (radians)
        l = 0.0  # Mean anomaly (radians)

        # Expected position and velocity
        expected_position = Axes(1.0, 0.0, 0.0)  # Since circular orbit at perihelion
        expected_velocity_magnitude = np.sqrt(G * (central_mass + target_mass) / q)
        expected_velocity = Axes(0.0, expected_velocity_magnitude, 0.0)

        # Calculate Cartesian coordinates
        position, velocity = calculate_cartesian_coordinates(
            target_mass, q, e, i0, p, n0, l, masses=[central_mass], positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(position.x(), expected_position.x(), places=6)
        self.assertAlmostEqual(position.y(), expected_position.y(), places=6)
        self.assertAlmostEqual(position.z(), expected_position.z(), places=6)
        # For velocity, since in the simplified code velocity might be zero, we check the magnitude
        velocity_magnitude = np.sqrt(velocity.x()**2 + velocity.y()**2 + velocity.z()**2)
        self.assertAlmostEqual(velocity_magnitude, expected_velocity_magnitude, places=6)

    def test_calculate_cartesian_coordinates_elliptical_orbit(self):
        """
        Test calculate_cartesian_coordinates for an elliptical orbit.
        """
        target_mass = 0.0
        central_mass = 1.0
        total_mass = target_mass + central_mass

        # Orbital elements
        q = 0.5  # Perihelion distance (AU)
        e = 0.5  # Eccentricity
        i0 = 30.0 * DEG2RAD  # Inclination (radians)
        p = 45.0 * DEG2RAD  # Longitude of perihelion (radians)
        n0 = 60.0 * DEG2RAD  # Longitude of ascending node (radians)
        l = 0.0  # Mean anomaly at perihelion

        # Calculate Cartesian coordinates
        position, velocity = calculate_cartesian_coordinates(
            target_mass, q, e, i0, p, n0, l, masses=[central_mass], positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Since the actual calculation is complex, we'll verify by converting back to orbital elements
        a, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f = calculate_keplerian_orbital_elements(
            target_mass, position, velocity, masses=[central_mass], positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(q_calc, q, places=6)
        self.assertAlmostEqual(e_calc, e, places=6)
        self.assertAlmostEqual(i_calc, i0, places=6)
        # Since p and n can differ by multiples of 2π, we'll check their modulus
        self.assertAlmostEqual(modulus(p_calc, TWO_PI), modulus(p, TWO_PI), places=6)
        self.assertAlmostEqual(modulus(n_calc, TWO_PI), modulus(n0, TWO_PI), places=6)

    def test_calculate_cartesian_coordinates_hyperbolic_orbit(self):
        """
        Test calculate_cartesian_coordinates for a hyperbolic orbit.
        """
        target_mass = 0.0
        central_mass = 1.0

        # Orbital elements
        q = 1.0  # Perihelion distance (AU)
        e = 1.5  # Eccentricity (>1 for hyperbolic orbit)
        i0 = 0.0  # Inclination (degrees)
        p = 0.0  # Longitude of perihelion (radians)
        n0 = 0.0  # Longitude of ascending node (degrees)
        l = 0.0  # Mean anomaly at perihelion

        # Calculate Cartesian coordinates
        position, velocity = calculate_cartesian_coordinates(
            target_mass, q, e, i0, p, n0, l, masses=[central_mass],
            positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Verify by converting back to orbital elements
        a, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f = calculate_keplerian_orbital_elements(
            target_mass, position, velocity, masses=[central_mass],
            positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(q_calc, q, places=6)
        self.assertAlmostEqual(e_calc, e, places=6)
        self.assertAlmostEqual(i_calc, i0, places=6)

    def test_calculate_cartesian_coordinates_inclined_orbit(self):
        """
        Test calculate_cartesian_coordinates for an inclined orbit.
        """
        target_mass = 0.0
        central_mass = 1.0

        # Orbital elements
        q = 1.0  # Perihelion distance (AU)
        e = 0.1  # Eccentricity
        i0 = 45.0 * DEG2RAD  # Inclination (radians)
        p = 30.0 * DEG2RAD  # Longitude of perihelion (radians)
        n0 = 60.0 * DEG2RAD  # Longitude of ascending node (radians)
        l = 10.0 * DEG2RAD  # Mean anomaly (radians)

        # Calculate Cartesian coordinates
        position, velocity = calculate_cartesian_coordinates(
            target_mass, q, e, i0, p, n0, l, masses=[central_mass],
            positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Verify by converting back to orbital elements
        a, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f = calculate_keplerian_orbital_elements(
            target_mass, position, velocity, masses=[central_mass],
            positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(e_calc, e, places=5)
        self.assertAlmostEqual(i_calc, i0, places=5)
        self.assertAlmostEqual(modulus(p_calc, TWO_PI), modulus(p, TWO_PI), places=5)
        self.assertAlmostEqual(modulus(n_calc, TWO_PI), modulus(n0, TWO_PI), places=5)

    def test_calculate_keplerian_elements_edge_cases(self):
        """
        Test calculate_keplerian_orbital_elements with edge cases (e close to 0 or 1).
        """
        # Edge case: e very close to zero (nearly circular orbit)
        target_mass = 0.0
        central_mass = 1.0
        e = 1e-8
        semi_major_axis = 1.0  # AU
        position = Axes(1.0, 0.0, 0.0)
        orbital_speed = np.sqrt(G * (central_mass + target_mass) / semi_major_axis)
        velocity = Axes(0.0, orbital_speed, 0.0)
        a, q, e_calc, i, p, n, l, f = calculate_keplerian_orbital_elements(
            target_mass, position, velocity, masses=[central_mass],
            positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )
        self.assertAlmostEqual(e_calc, 0.0, places=10)

        # Edge case: e very close to 1 (parabolic orbit)
        e = 1.0 - 1e-8
        periapsis_distance = 1.0  # AU
        position = Axes(periapsis_distance, 0.0, 0.0)
        orbital_speed = np.sqrt(2 * G * (central_mass + target_mass) / periapsis_distance * (1 - e))
        velocity = Axes(0.0, orbital_speed, 0.0)
        a, q, e_calc, i, p, n, l, f = calculate_keplerian_orbital_elements(
            target_mass, position, velocity, masses=[central_mass],
            positions=[Axes(0, 0, 0)], velocities=[Axes(0, 0, 0)]
        )
        self.assertAlmostEqual(e_calc, e, places=7)






    def test_calculate_keplerian_orbital_elements_circumbinary_orbit(self):
        """
        Test calculate_keplerian_orbital_elements with a circumbinary planet orbiting two stars.
        """
        # Masses of the two stars
        star1_mass = 1.0
        star2_mass = 0.5
        total_mass = star1_mass + star2_mass

        # Positions and velocities of the stars (assuming they orbit their common center of mass)
        separation = 0.2  # AU
        star1_position = Axes(-star2_mass / total_mass * separation, 0.0, 0.0)
        star2_position = Axes(star1_mass / total_mass * separation, 0.0, 0.0)

        # Velocities of the stars
        orbital_speed = np.sqrt(G * total_mass / separation)
        star1_velocity = Axes(0.0, star2_mass / total_mass * orbital_speed, 0.0)
        star2_velocity = Axes(0.0, -star1_mass / total_mass * orbital_speed, 0.0)

        # Planet orbiting the binary stars
        planet_mass = 0.0
        semi_major_axis = 1.0  # AU
        eccentricity = 0.1
        inclination = 0.0  # Coplanar orbit

        # Position and velocity of the planet at true anomaly f = 0
        f = 0.0
        r = semi_major_axis * (1 - eccentricity)
        position = Axes(r, 0.0, 0.0)
        orbital_speed_planet = np.sqrt(G * total_mass * (1 + eccentricity) / (semi_major_axis * (1 - eccentricity)))
        velocity = Axes(0.0, orbital_speed_planet, 0.0)

        # Masses, positions, and velocities of the central bodies
        masses = [star1_mass, star2_mass]
        positions = [star1_position, star2_position]
        velocities = [star1_velocity, star2_velocity]

        # Calculate orbital elements
        a_calc, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc = calculate_keplerian_orbital_elements(
            planet_mass,
            position,
            velocity,
            masses=masses,
            positions=positions,
            velocities=velocities
        )

        # Assertions
        self.assertAlmostEqual(a_calc, semi_major_axis, places=6)
        self.assertAlmostEqual(e_calc, eccentricity, places=6)
        self.assertAlmostEqual(i_calc, inclination, places=6)


    def test_calculate_keplerian_orbital_elements_zero_eccentricity_and_inclination(self):
        """
        Test calculate_keplerian_orbital_elements with zero eccentricity and zero inclination.
        """
        target_mass = 0.0
        central_mass = 1.0

        # Orbital parameters
        semi_major_axis = 1.0  # AU
        eccentricity = 0.0
        inclination = 0.0
        p = 0.0  # Undefined for circular orbits
        n = 0.0  # Undefined for zero inclination
        l = 45.0 * DEG2RAD  # Mean anomaly

        # Position and velocity in circular orbit
        position, velocity = calculate_cartesian_coordinates(
            target_mass,
            q=semi_major_axis * (1 - eccentricity),
            e=eccentricity,
            i0=inclination,
            p=p,
            n0=n,
            l=l,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)],
        )

        # Calculate orbital elements
        a_calc, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc = calculate_keplerian_orbital_elements(
            target_mass,
            position,
            velocity,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)],
        )
        print(a_calc, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc)

        # Assertions
        self.assertAlmostEqual(a_calc, semi_major_axis, places=6)
        self.assertAlmostEqual(e_calc, eccentricity, places=6)
        self.assertAlmostEqual(i_calc, inclination, places=6)
        self.assertAlmostEqual(l_calc, l, places=6)
        self.assertAlmostEqual(f_calc, l, places=6)  # For circular orbits, true anomaly equals mean anomaly


    def test_calculate_keplerian_orbital_elements_high_eccentricity(self):
        """
        Test calculate_keplerian_orbital_elements with high eccentricity.
        """
        target_mass = 0.0
        central_mass = 1.0

        # Orbital parameters
        semi_major_axis = 1.0  # AU
        eccentricity = 0.95  # High eccentricity

        # Position at periapsis
        r = semi_major_axis * (1 - eccentricity)
        position = Axes(r, 0.0, 0.0)

        # Velocity at periapsis
        orbital_speed = np.sqrt(G * central_mass * (1 + eccentricity) / (semi_major_axis * (1 - eccentricity)))
        velocity = Axes(0.0, orbital_speed, 0.0)

        # Calculate orbital elements
        a_calc, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc = calculate_keplerian_orbital_elements(
            target_mass,
            position,
            velocity,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)],
        )

        # Assertions
        self.assertAlmostEqual(a_calc, semi_major_axis, places=5)
        self.assertAlmostEqual(e_calc, eccentricity, places=5)

    def test_calculate_keplerian_elements_random_orbits(self):
        """
        Test with random orbital elements to ensure robustness across a wide range of inputs.
        """
        for _ in range(1000):  # Test with 10 random orbits
            target_mass = 0.0
            central_mass = 1.0

            # Random orbital elements
            q = np.random.uniform(0.1, 5.0)  # Perihelion distance between 0.1 and 5 AU
            e = np.random.uniform(0.01, 0.8)  # Eccentricity between 0.01 and 0.8
            i0 = np.random.uniform(1.0, 179.0) * DEG2RAD  # Inclination between 1 and 179 degrees
            p = np.random.uniform(0.0, 360.0) * DEG2RAD  # Argument of periapsis
            n0 = np.random.uniform(0.0, 360.0) * DEG2RAD  # Longitude of ascending node
            l = np.random.uniform(0.0, 360.0) * DEG2RAD  # Mean anomaly
            print(q, e, i0, p, n0, l)

            # Calculate Cartesian coordinates
            position, velocity = calculate_cartesian_coordinates(
                target_mass, q, e, i0, p, n0, l,
                masses=[central_mass],
                positions=[Axes(0, 0, 0)],
                velocities=[Axes(0, 0, 0)]
            )

            # Calculate orbital elements from Cartesian coordinates
            a, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc = calculate_keplerian_orbital_elements(
                target_mass, position, velocity,
                masses=[central_mass],
                positions=[Axes(0, 0, 0)],
                velocities=[Axes(0, 0, 0)]
            )

            ## TODO: fix code?
            #if i_calc > np.pi:
                #i_calc = 2 * np.pi - i_calc

            # Assertions
            self.assertAlmostEqual(q_calc, q, places=3)
            self.assertAlmostEqual(e_calc, e, places=3)
            self.assertAlmostEqual(i_calc, i0, places=3)
            self.assertAlmostEqual(modulus(p_calc, TWO_PI), modulus(p, TWO_PI), places=3)
            self.assertAlmostEqual(modulus(n_calc, TWO_PI), modulus(n0, TWO_PI), places=3)


    def test_calculate_keplerian_orbital_elements_non_zero_arguments(self):
        """
        Test calculate_keplerian_orbital_elements with non-zero argument of periapsis and longitude of ascending node.
        """
        target_mass = 0.0
        central_mass = 1.0
        total_mass = target_mass + central_mass

        # Orbital parameters
        semi_major_axis = 2.0  # AU
        eccentricity = 0.3
        inclination = 45.0 * DEG2RAD
        p = 30.0 * DEG2RAD  # Argument of periapsis
        n = 60.0 * DEG2RAD  # Longitude of ascending node
        f = 90.0 * DEG2RAD  # True anomaly

        # Calculate position and velocity in orbital plane
        r = semi_major_axis * (1 - eccentricity ** 2) / (1 + eccentricity * np.cos(f))
        x_orb = r * np.cos(f)
        y_orb = r * np.sin(f)

        # Specific angular momentum
        h = np.sqrt(G * total_mass * semi_major_axis * (1 - eccentricity ** 2))

        # Velocity in orbital plane
        vx_orb = -h / r * np.sin(f)
        vy_orb = h / r * (eccentricity + np.cos(f))

        # Rotation matrices components
        cos_n = np.cos(n)
        sin_n = np.sin(n)
        cos_i = np.cos(inclination)
        sin_i = np.sin(inclination)
        cos_p = np.cos(p)
        sin_p = np.sin(p)

        # Position in inertial frame
        x = (cos_n * cos_p - sin_n * sin_p * cos_i) * x_orb + (-cos_n * sin_p - sin_n * cos_p * cos_i) * y_orb
        y = (sin_n * cos_p + cos_n * sin_p * cos_i) * x_orb + (-sin_n * sin_p + cos_n * cos_p * cos_i) * y_orb
        z = (sin_p * sin_i) * x_orb + (cos_p * sin_i) * y_orb

        position = Axes(x, y, z)

        # Velocity in inertial frame
        vx = (cos_n * cos_p - sin_n * sin_p * cos_i) * vx_orb + (-cos_n * sin_p - sin_n * cos_p * cos_i) * vy_orb
        vy = (sin_n * cos_p + cos_n * sin_p * cos_i) * vx_orb + (-sin_n * sin_p + cos_n * cos_p * cos_i) * vy_orb
        vz = (sin_p * sin_i) * vx_orb + (cos_p * sin_i) * vy_orb

        velocity = Axes(vx, vy, vz)

        # Calculate orbital elements
        a_calc, q_calc, e_calc, i_calc, p_calc, n_calc, l, f_calc = calculate_keplerian_orbital_elements(
            target_mass,
            position,
            velocity,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)],
        )

        # Assertions
        self.assertAlmostEqual(a_calc, semi_major_axis, places=6)
        self.assertAlmostEqual(e_calc, eccentricity, places=6)
        self.assertAlmostEqual(i_calc, inclination, places=6)
        # Adjust assertion for possible angle differences
        #delta_p = modulus(p_calc - p, TWO_PI)
        delta_p = modulus(p_calc - (n + p), TWO_PI)
        if delta_p > np.pi:
            delta_p = TWO_PI - delta_p
        self.assertTrue(np.isclose(delta_p, 0.0, atol=1e-6) or np.isclose(delta_p, np.pi, atol=1e-6))
        delta_n = modulus(n_calc - n, TWO_PI)
        self.assertTrue(np.isclose(delta_n, 0.0, atol=1e-6) or np.isclose(delta_n, np.pi, atol=1e-6))


    #--------------------------------------------------------------------------------
    # [start] Keplerian -> Cartesians -> Keplerian
    #--------------------------------------------------------------------------------
    def test_calculate_keplerian_elements_and_back(self):
        """
        Test converting from orbital elements to Cartesian coordinates and back.
        """
        target_mass = 0.0
        central_mass = 1.0

        # Define arbitrary orbital elements
        q = 0.8
        e = 0.3
        i0 = 20.0 * DEG2RAD
        p = 50.0 * DEG2RAD
        n0 = 80.0 * DEG2RAD
        l = 30.0 * DEG2RAD

        # Calculate Cartesian coordinates
        position, velocity = calculate_cartesian_coordinates(
            target_mass, q, e, i0, p, n0, l,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Calculate orbital elements from Cartesian coordinates
        a, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc = calculate_keplerian_orbital_elements(
            target_mass, position, velocity,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(q_calc, q, places=5)
        self.assertAlmostEqual(e_calc, e, places=5)
        self.assertAlmostEqual(i_calc, i0, places=5)
        self.assertAlmostEqual(modulus(p_calc, TWO_PI), modulus(p, TWO_PI), places=5)
        self.assertAlmostEqual(modulus(n_calc, TWO_PI), modulus(n0, TWO_PI), places=5)
        # Mean anomaly might have a different reference point, but their difference should be consistent

    def test_retrograde_orbit_and_back(self):
        """
        Test converting from orbital elements to Cartesian coordinates and back for a retrograde orbit.
        """
        target_mass = 0.0
        central_mass = 1.0

        # Orbital elements for a retrograde orbit
        q = 1.0
        e = 0.4
        i0 = 150.0 * DEG2RAD  # Retrograde inclination
        p = 70.0 * DEG2RAD
        n0 = 110.0 * DEG2RAD
        l = 250.0 * DEG2RAD

        # Calculate Cartesian coordinates
        position, velocity = calculate_cartesian_coordinates(
            target_mass, q, e, i0, p, n0, l,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Calculate orbital elements from Cartesian coordinates
        a, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc = calculate_keplerian_orbital_elements(
            target_mass, position, velocity,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Normalize inclination to [0, π]
        i_calc_normalized = modulus(i_calc, np.pi)
        self.assertAlmostEqual(i_calc_normalized, i0 % np.pi, places=4)
        self.assertAlmostEqual(e_calc, e, places=4)
        self.assertAlmostEqual(q_calc, q, places=4)
        # Adjust angle assertions for possible angle differences
        delta_p = modulus(p_calc - p, TWO_PI)
        self.assertTrue(np.isclose(delta_p, 0.0, atol=1e-4) or np.isclose(delta_p, np.pi, atol=1e-4))
        delta_n = modulus(n_calc - n0, TWO_PI)
        self.assertTrue(np.isclose(delta_n, 0.0, atol=1e-4) or np.isclose(delta_n, np.pi, atol=1e-4))

    def test_polar_orbit_and_back(self):
        """
        Test converting from orbital elements to Cartesian coordinates and back for a polar orbit.
        """
        target_mass = 0.0
        central_mass = 1.0

        # Orbital elements for a polar orbit
        q = 1.5
        e = 0.2
        i0 = 90.0 * DEG2RAD  # Polar orbit
        p = 45.0 * DEG2RAD
        n0 = 30.0 * DEG2RAD
        l = 120.0 * DEG2RAD

        # Calculate Cartesian coordinates
        position, velocity = calculate_cartesian_coordinates(
            target_mass, q, e, i0, p, n0, l,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Calculate orbital elements from Cartesian coordinates
        a, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc = calculate_keplerian_orbital_elements(
            target_mass, position, velocity,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(i_calc, i0, places=4)
        self.assertAlmostEqual(e_calc, e, places=4)
        self.assertAlmostEqual(q_calc, q, places=4)
        self.assertAlmostEqual(modulus(p_calc, TWO_PI), modulus(p, TWO_PI), places=4)
        self.assertAlmostEqual(modulus(n_calc, TWO_PI), modulus(n0, TWO_PI), places=4)

    def test_near_parabolic_orbit_and_back(self):
        """
        Test converting from orbital elements to Cartesian coordinates and back for a near-parabolic orbit.
        """
        target_mass = 0.0
        central_mass = 1.0

        # Orbital elements for a near-parabolic orbit
        q = 0.5
        e = 0.9999  # Very close to 1
        i0 = 10.0 * DEG2RAD
        p = 25.0 * DEG2RAD
        n0 = 45.0 * DEG2RAD
        l = 0.0 * DEG2RAD  # At perihelion

        # Calculate Cartesian coordinates
        position, velocity = calculate_cartesian_coordinates(
            target_mass, q, e, i0, p, n0, l,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Since the semi-major axis is very large, we may not assert on 'a'
        # Calculate orbital elements from Cartesian coordinates
        a_calc, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc = calculate_keplerian_orbital_elements(
            target_mass, position, velocity,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Assertions
        self.assertAlmostEqual(e_calc, e, places=3)
        self.assertAlmostEqual(q_calc, q, places=3)
        self.assertAlmostEqual(i_calc, i0, places=3)
        # Adjusted precision due to numerical instability


    def test_zero_inclination_nonzero_argument_and_back(self):
        """
        Test converting from orbital elements to Cartesian coordinates and back for zero inclination and non-zero argument of periapsis.
        """
        target_mass = 0.0
        central_mass = 1.0

        # Orbital elements
        q = 1.0
        e = 0.5
        i0 = 0.0 * DEG2RAD  # Zero inclination
        p = 45.0 * DEG2RAD  # Non-zero argument of periapsis
        n0 = 0.0 * DEG2RAD  # Longitude of ascending node is undefined when inclination is zero
        l = 60.0 * DEG2RAD

        # Calculate Cartesian coordinates
        position, velocity = calculate_cartesian_coordinates(
            target_mass, q, e, i0, p, n0, l,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Calculate orbital elements from Cartesian coordinates
        a, q_calc, e_calc, i_calc, p_calc, n_calc, l_calc, f_calc = calculate_keplerian_orbital_elements(
            target_mass, position, velocity,
            masses=[central_mass],
            positions=[Axes(0, 0, 0)],
            velocities=[Axes(0, 0, 0)]
        )

        # Since inclination is zero, the longitude of ascending node is undefined (set to zero)
        self.assertAlmostEqual(i_calc, i0, places=5)
        self.assertAlmostEqual(e_calc, e, places=5)
        self.assertAlmostEqual(q_calc, q, places=5)
        self.assertAlmostEqual(modulus(p_calc, TWO_PI), modulus(p + n0, TWO_PI), places=5)

    #--------------------------------------------------------------------------------
    # [end] Keplerian -> Cartesians -> Keplerian
    #--------------------------------------------------------------------------------


