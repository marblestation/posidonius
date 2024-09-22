extern crate time;
use std;
use time::{OffsetDateTime, format_description};
use super::constants::*;
use super::particles::Axes;
//use std::cmp::Ord;

pub fn calculate_inclination_orbital_equatorial_plane(position: Axes, velocity: Axes, spin: Axes) -> f64 {
    // Calculate the spin axis inclination described as
    // the inclination of the orbit plane with respect to the equatorial plane
    // --- Input
    let x = position.x;
    let y = position.y;
    let z = position.z;
    let u = velocity.x;
    let v = velocity.y;
    let w = velocity.z;
    let sx = spin.x;
    let sy = spin.y;
    let sz = spin.z;
    let s = (sx.powi(2) + sy.powi(2) + sz.powi(2)).sqrt();

    // --- Calculate the component of the orbital angular momentum
    let hx = y * w - z * v;
    let hy = z * u - x * w;
    let hz = x * v - y * u;
    let h = (hx.powi(2) + hy.powi(2) + hz.powi(2)).sqrt();

    let hx_rel = hx / h;
    let hy_rel = hy / h;
    let hz_rel = hz / h;
    let h_rel = (hx_rel.powi(2) + hy_rel.powi(2) + hz_rel.powi(2)).sqrt();

    let numerator = hx_rel * sx + hy_rel * sy + hz_rel * sz;
    let denominator = h_rel * s;
    //
    // println!("F1: {:?}, {:?}, {:?}", position, velocity, spin);
    // println!("F2: {:?}, {:?}, {:?}, {:?}", hx, hy, hz, h);
    // println!("F3: {:?}, {:?}", numerator, denominator);

    let cos_inclination = (numerator / denominator).clamp(-1., 1.);

    cos_inclination.acos()
}

pub fn calculate_eccentricity_vector(gm: f64, position: Axes, velocity: Axes) -> Axes {
    // --- Input --- //
    let x = position.x;
    let y = position.y;
    let z = position.z;
    let u = velocity.x;
    let v = velocity.y;
    let w = velocity.z;
    
    // --- Output --- //
    let mut eccentricity_vector = Axes{x:0., y:0., z:0. };
    
    // --- Local --- //
    // Angular momentum
    let hx = y * w  -  z * v;
    let hy = z * u  -  x * w;
    let hz = x * v  -  y * u;
    // v vectorial h
    let v_vect_h_x = v * hz - w * hy;
    let v_vect_h_y = w * hx - u * hz;
    let v_vect_h_z = u * hy - v * hx;
    // distance
    let r = (x*x + y*y + z*z).sqrt();
    // eccentricity componant
    eccentricity_vector.x = (v_vect_h_x / gm) - (x / r);
    eccentricity_vector.y = (v_vect_h_y / gm) - (y / r);
    eccentricity_vector.z = (v_vect_h_z / gm) - (z / r);

    eccentricity_vector
}

pub fn calculate_keplerian_orbital_elements(gm: f64, position: Axes, velocity: Axes) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    // ! Based on the implementation of Chambers in Mercury
    // Calculates Keplerian orbital elements given relative coordinates and
    // velocities, and GM = G times the sum of the masses.

    // --- Input
    let x = position.x;
    let y = position.y;
    let z = position.z;
    let u = velocity.x;
    let v = velocity.y;
    let w = velocity.z;
    // --- Output
    let a: f64; // semi-major axis (in AU)
    let q: f64; // perihelion distance
    let mut eccentricity: f64; // eccentricity
    let mut i: f64 = 0.; // inclination
    let mut _p: f64; // longitude of perihelion (NOT argument of perihelion!!)
    let mut _n: f64; // longitude of ascending node
    let mut _l: f64; // mean anomaly (or mean longitude if eccentricity < 1.e-8)
    let orbital_period: f64; // orbital_period (in days)
    // --- Local
    let eccentricity_vector = calculate_eccentricity_vector(gm, position, velocity);
    let ex = eccentricity_vector.x;
    let ey = eccentricity_vector.y;
    let ez = eccentricity_vector.z;
    let e_scal_r = ex*x + ey*y + ez*z;
    let hx = y * w  -  z * v;
    let hy = z * u  -  x * w;
    let hz = x * v  -  y * u;
    let h2 = hx.powf(2.) + hy.powf(2.) + hz.powf(2.);
    let v2 = u * u  +  v * v  +  w * w;
    let rv = x * u  +  y * v  +  z * w;
    let r = (x*x + y*y + z*z).sqrt();
    let h = h2.sqrt();
    let s = h2 / gm;

    // --- Semi-major axis
    a  = gm * r / (2.0 * gm  -  r * v2);

  
    // --- Inclination and longitude of ascending node
    let mut longitude_of_ascending_node:f64;
    let cos_i = (hz / h).clamp(-1.,1.);
    if cos_i.abs() < 1. {
        i = cos_i.acos();
        if hy > 0. {
            i = TWO_PI - i;
        }
        //longitude_of_ascending_node = (hx / -hy).atan();
        longitude_of_ascending_node = hx.atan2(-hy);
        if longitude_of_ascending_node < 0. {
            longitude_of_ascending_node = longitude_of_ascending_node + TWO_PI;
        }
        if i == 0. {
            longitude_of_ascending_node = 0.;
        }
    } else {
        if cos_i > 0. {
            i = 0.;
        }
        if cos_i < 0. { 
            i = PI;
        }
        longitude_of_ascending_node = 0.;
    }
    //
  
    // --- Eccentricity and perihelion distance
    let temp = 1.  +  s * (v2 / gm  -  2. / r);
    if temp <= 0. {
        eccentricity = 0.;
    } else {
        eccentricity = temp.sqrt();
    }
    q = s / (1. + eccentricity);
    if eccentricity < 3.0e-8 {
        eccentricity = 0.;
    }

    // --- Longitude of perihelion
    //let longitude_perihelion:f64;
    //if eccentricity == 0. {
        //longitude_perihelion= 0.;
    //} else {
        //longitude_perihelion = modulus(ey.atan2(ex) + TWO_PI + TWO_PI, TWO_PI);
    //}
    let argument_perihelion:f64;
    if eccentricity == 0. {
        argument_perihelion= 0.;
    } else {
        let ex_rot = ex*longitude_of_ascending_node.cos() +ey*longitude_of_ascending_node.sin();
        let ey_rot = -ex*longitude_of_ascending_node.sin() +ey*longitude_of_ascending_node.cos();
        let ez_rot = ez;
        let ex_rot_2 = ex_rot;
        let ey_rot_2 = ey_rot*i.cos() + ez_rot*i.sin();
        //let ez_rot_2 = -ey_rot*i.sin() + ez_rot*i.cos();

        // argument_perihelion = (ey_rot_2/ex_rot_2).atan();
        argument_perihelion = (ey_rot_2).atan2(ex_rot_2);
    }

    let longitude_perihelion = argument_perihelion + longitude_of_ascending_node;

    // --- True anomaly
    let mut true_anomaly:f64;
    let cos_f:f64;
    if eccentricity == 0. {
        let x_rot = x*longitude_of_ascending_node.cos() + y*longitude_of_ascending_node.sin();
        cos_f = (x_rot / r).clamp(-1.,1.);
        true_anomaly = cos_f.acos();
    } else {
        cos_f = (e_scal_r / (eccentricity * r)).clamp(-1., 1.);
        true_anomaly = cos_f.acos();
        if rv < 0. {
            true_anomaly = TWO_PI - true_anomaly;
        }
    }

    
    // --- Mean anomaly
    let mut mean_anomaly:f64;
    if eccentricity == 0. {
        mean_anomaly = true_anomaly;
    } else {
        let mut cos_bige = ((1./eccentricity)*( 1. - (r / a) )).clamp(-1.,1.);

        // Mean anomaly for ellipse
        if eccentricity < 1. {
            if cos_bige.abs() > 1. {
                cos_bige = cos_bige.signum();
            }
            let mut bige = cos_bige.acos();
            if rv < 0. {
                bige = TWO_PI - bige;
            }
            mean_anomaly = bige - eccentricity*(bige.sin());
        } else {
            // Mean anomaly for hyperbola
            if cos_bige < 1. {
                cos_bige = 1.;
            }
            let mut bige = (cos_bige + (cos_bige*cos_bige-1.).sqrt()).log(std::f64::consts::E);
            if rv < 0. {
                bige = - bige;
            }
            mean_anomaly = eccentricity * bige.sinh() - bige;
        }

    }
  
    if mean_anomaly < 0. {
        mean_anomaly = mean_anomaly + TWO_PI;
    }
    if mean_anomaly > TWO_PI {
        mean_anomaly = modulus(mean_anomaly, TWO_PI);
    }

    // Given relative coordinates and velocities (of the body 1 respect to body 2),
    // and GM = G times the sum of the masses (body 1 + body 2)
    orbital_period = (TWO_PI / gm.sqrt()) * a.powf(3./2.); // in days

    // a Semimajor axis
    // q perihelion distance
    // e eccentricity
    // i inclination
    // p longitude of perihelion (! not the argument of perihelion)
    // n longitude of ascending node
    // l mean anomaly
    // orbital period bah...

    return (a, q, eccentricity, i, longitude_perihelion, longitude_of_ascending_node, mean_anomaly, orbital_period)
}

pub fn calculate_perihelion_distance_and_eccentricity(gm: f64, position: Axes, velocity: Axes) -> (f64, f64) {
    // ! Based on the implementation of Chambers in Mercury
    // Calculates Keplerian orbital elements given relative coordinates and
    // velocities, and GM = G times the sum of the masses.

    // Input
    let x = position.x;
    let y = position.y;
    let z = position.z;
    let u = velocity.x;
    let v = velocity.y;
    let w = velocity.z;
    // Output
    let q: f64; // perihelion distance
    let eccentricity: f64; // eccentricity
    // Local
    let hx = y * w  -  z * v;
    let hy = z * u  -  x * w;
    let hz = x * v  -  y * u;
    let h2 = hx.powf(2.) + hy.powf(2.) + hz.powf(2.);
    let v2 = u * u  +  v * v  +  w * w;
    let r = (x*x + y*y + z*z).sqrt();
    let s = h2 / gm;

    // Eccentricity and perihelion distance
    let temp = 1.  +  s * (v2 / gm  -  2. / r);
    if temp <= 0. {
        eccentricity = 0.;
    } else {
        eccentricity = temp.sqrt();
    }
    q = s / (1. + eccentricity);

    return(q, eccentricity)
}



pub fn calculate_cartesian_coordinates(gm: f64, q: f64, e: f64, i0: f64, p: f64, n0: f64, l: f64)  -> (f64, f64, f64, f64, f64, f64) {
    // Calculates Cartesian coordinates and velocities given Keplerian orbital
    // elements (for elliptical, parabolic or hyperbolic orbits).
    // ! Based on the implementation of Chambers in Mercury, which is
    //   based on a routine from Levison and Duncan's SWIFT integrator.
    // WARNING: It gives NaN velocities when eccentricity == 1. (also in the original implementation in mercury code)

    // Input
    // gm  = grav const * (central + secondary mass)
    // q  = perihelion distance
    // e  = eccentricity
    // i  = inclination (degrees)
    // p  = longitude of perihelion !!! 
    // n  = longitude of the ascending node (degrees)
    // l  = mean anomaly 
      
    // Output
    let x: f64; // Cartesian positions  ( units the same as a )
    let y: f64; // Cartesian positions  ( units the same as a )
    let z: f64; // Cartesian positions  ( units the same as a )
    let u: f64; // Cartesian velocities ( units the same as sqrt(gm/a) )
    let v: f64; // Cartesian velocities ( units the same as sqrt(gm/a) )
    let w: f64; // Cartesian velocities ( units the same as sqrt(gm/a) )
  
    // Change from longitude of perihelion to argument of perihelion
    let g0 = p - n0;

    // Rotation factors
    let si = i0.sin();
    let ci = i0.cos();
    let sg = g0.sin();
    let cg = g0.cos();
    let sn = n0.sin();
    let cn = n0.cos();
    //let (i, si, ci) = mco_sine(i0);
    //let (g, sg, cg) = mco_sine(g0);
    //let (n, sn, cn) = mco_sine(n0);
    let mut z1 = cg * cn;
    let mut z2 = cg * sn;
    let mut z3 = sg * cn;
    let mut z4 = sg * sn;
    let d11 =  z1 - z4*ci;
    let d12 =  z2 + z3*ci;
    let d13 = sg * si;
    let d21 = -z3 - z2*ci;
    let d22 = -z4 + z1*ci;
    let d23 = cg * si;

    // Semi-major axis
    let a = q / (1. - e);

    let se: f64;
    let ce: f64;

    // Ellipse
    if e < 1. {
       let romes = (1. - e*e).sqrt();
       let temp0 = kepler_solution_for_eccentrities_smaller_than_one(e,l);
       se = temp0.sin();
       ce = temp0.cos();
       //let (_, se, ce) = mco_sine(temp0);
       z1 = a * (ce - e);
       z2 = a * romes * se;
       let temp = (gm/a).sqrt() / (1. - e*ce);
       z3 = -se * temp;
       z4 = romes * ce * temp;
    } else {
       // Parabola
       if e == 1. {
          let tmp = kepler_solution_for_a_parabola(l);
          ce = tmp.1;
          z1 = q * (1. - ce*ce);
          z2 = 2. * q * ce;
          z4 = (2.*gm/q).sqrt() / (1. + ce*ce);
          z3 = -ce * z4;
       } else {
          // Hyperbola
          let romes = (e*e - 1.).sqrt();
          let mut temp = kepler_solution_for_a_hyperbola(e,l);
          se = temp.sinh();
          ce = temp.cosh();
          z1 = a * (ce - e);
          z2 = -a * romes * se;
          temp = (gm/a.abs()).sqrt() / (e*ce - 1.);
          z3 = -se * temp;
          z4 = romes * ce * temp;
       }
    }

    x = d11 * z1  +  d21 * z2;
    y = d12 * z1  +  d22 * z2;
    z = d13 * z1  +  d23 * z2;
    u = d11 * z3  +  d21 * z4;
    v = d12 * z3  +  d22 * z4;
    w = d13 * z3  +  d23 * z4;

    return (x, y, z, u, v, w);
}




fn kepler_solution_for_eccentrities_smaller_than_one(e: f64, oldl: f64) -> f64 {
    // ! Based on the implementation of Chambers in Mercury (mco_kep)
    // Solves Kepler's equation for eccentricities less than one.
    // Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
    //
    // e = eccentricity\n
    // l = mean anomaly      (radians)\n
    // u = eccentric anomaly (   "   )\n
  
  
    //// Local
    //real(double_precision) :: l,pi,twopi,piby2,u1,u2,ome,sign
    //real(double_precision) :: x,x2,sn,dsn,z1,z2,z3,f0,f1,f2,f3
    //real(double_precision) :: p,q,p2,ss,cc
    //logical flag,big,bigg
    let ome: f64;
    let mut x: f64;
    let mut x2: f64;
    let mut dsn: f64;
    let u1: f64;
    let mut u2: f64;
    let mut f0: f64; 
    let mut f1: f64; 
    let mut f2: f64; 
    let f3: f64; 
    let mut z1: f64; 
    let mut z2: f64; 
    let mut z3: f64; 
    let p: f64;
    let p2: f64;
    let q: f64;
    let sn: f64;
    let ss: f64;
    let cc: f64;
  
    let mut l: f64;
    let piby2 = 0.5 * PI;
  
    // Reduce mean anomaly to lie in the range 0 < l < pi
    if oldl >= 0. {
       l = modulus(oldl, TWO_PI);
    } else {
       l = (modulus(oldl, TWO_PI)) + TWO_PI;
    }
    let mut sign = 1.;
    if l > PI {
       l = TWO_PI - l;
       sign = -1.;
    }
    
    ome = 1. - e;

    if (l >= 0.45) || (e < 0.55) {
        // Regions A,B or C in Nijenhuis
        // -----------------------------
        
        // Rough starting value for eccentric anomaly
        if l < ome {
           u1 = ome;
        } else {
            if l > (PI-1.-e) {
                u1 = (l+e*PI)/(1.+e);
            } else {
                u1 = l + e;
            }
        }
        
        // Improved value using Halley's method
        let flag = u1 > piby2;
        if flag {
           x = PI - u1;
        } else {
           x = u1;
        }
        x2 = x*x;
        sn = x*(1. + x2*(-0.16605 + x2*0.00761) );
        dsn = 1. + x2*(-0.49815 + x2*0.03805);
        if flag {
            dsn = -dsn;
        }
        f2 = e*sn;
        f0 = u1 - f2 - l;
        f1 = 1. - e*dsn;
        u2 = u1 - f0/(f1 - 0.5*f0*f2/f1);
    } else {
        // Region D in Nijenhuis
        // ---------------------
        
        // Rough starting value for eccentric anomaly
        z1 = 4.*e + 0.5;
        p = ome / z1;
        q = 0.5 * l / z1;
        p2 = p*p;
        z2 = (((p2*p + q*q ).sqrt() + q).log(std::f64::consts::E) / 1.5).exp();
        u1 = 2.*q / ( z2 + p + p2/z2 );
        
        // Improved value using Newton's method
        z2 = u1*u1;
        z3 = z2*z2;
        u2 = u1 - 0.075*u1*z3 / (ome + z1*z2 + 0.375*z3);
        u2 = l + e*u2*( 3. - 4.*u2*u2 );
    }
    
    // Accurate value using 3rd-order version of Newton's method
    // N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy
    
    // First get accurate values for u2 - sin(u2) and 1 - cos(u2)
    let bigg = u2 > piby2;
    if bigg {
        z3 = PI - u2;
    } else {
        z3 = u2;
    }
    
    let big = z3 > (0.5*piby2);
    if big {
       x = piby2 - z3;
    } else {
       x = z3;
    }
    
    x2 = x*x;
    
    ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. - x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))));
    cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. - x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. - x2/306.))))))));
    
    if big {
       z1 = cc + z3 - 1.;
       z2 = ss + z3 + 1. - piby2;
    } else {
       z1 = ss;
       z2 = cc;
    }
    
    if bigg {
       z1 = 2.*u2 + z1 - PI;
       z2 = 2. - z2;
    }
    
    f0 = l - u2*ome - e*z1;
    f1 = ome + e*z2;
    f2 = 0.5*e*(u2-z1);
    f3 = e/6.*(1.-z2);
    z1 = f0/f1;
    z2 = f0/(f2*z1+f1);
    return sign*( u2 + f0/((f3*z1+f2)*z2+f1) );
}


fn kepler_solution_for_a_parabola(q0: f64) -> (f64, f64) {
    // ! Based on the implementation of Duncan in Mercury (orbel_zget)
    // Solves the equivalent of Kepler's eqn. for a parabola 
    // given Q (Fitz. notation.)
    //
    // ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
    //
    // @remarks For a parabola we can solve analytically.

    // Input:
    // q0 = parabola mean anomaly
    
    // Output
    let mut eccentric_anomaly: f64; // eccentric anomaly

    // We copy the input parameter to be able to modify it.
    let mut q = q0;
    let x: f64;
    let tmp: f64;
    let iflag: i32;

    //----
    if q < 0. {
       iflag = 1;
       q = -q
    } else {
        iflag = 0;
    }

    if q < 1.0e-3 {
       eccentric_anomaly = q*(1. - (q*q/3.)*(1. -q*q));
    } else {
       x = 0.5*(3.*q + (9.* q.powf(2.) + 4.).sqrt());
       tmp = x.powf(1./3.);
       eccentric_anomaly = tmp - 1./tmp;
    }

    if iflag == 1 {
       eccentric_anomaly = -eccentric_anomaly;
       q = -q;
    }

    return (q, eccentric_anomaly);
}


fn kepler_solution_for_a_hyperbola(e: f64, n: f64) -> f64 {
    // ! Based on the implementation of Duncan in Mercury (orbel_fhybrid)
    // Solves Kepler's eqn. for hyperbola using hybrid approach. \n\n
    // ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
    //            For larger N, uses FGET

    // Inputs
    // e = eccentircity
    // n = hyperbola mean anomaly
    
    // Output
    let orbel_fhybrid: f64; // eccentric anomaly
    
    //...  Internals:
    let mut abn: f64;
  
    //----
    //...  Executable code 
  
    abn = n;
    if n <0. {
        abn = -abn;
    }
  
    //// TODO:
    if abn < 0.636*e -0.6 {
       let tmp = kepler_solution_for_a_hyperbola_hybrid_approach_for_low_n(e,n);
       orbel_fhybrid = tmp.1;
    } else {
       orbel_fhybrid = kepler_solution_for_a_hyperbola_hybrid_approach(e,n);
    }

    return orbel_fhybrid;
}


fn kepler_solution_for_a_hyperbola_hybrid_approach_for_low_n(e: f64, capn0: f64) -> (f64, f64) {
    // ! Based on the implementation of Duncan in Mercury (orbel_flon)
    // Solves Kepler's eqn. for hyperbola using hybrid approach. \n\n
    // ALGORITHM: Uses power series for N in terms of F and Newton,s method
    // REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)

    // Inputs:
    // e = eccentricty
    // capn = hyperbola mean anomaly
  
    // Output
    let mut orbel_flon;
  
    // copy of the input capn0 that is not modified
    let mut capn = capn0;
    let mut iflag : i32;
    let a: f64;
    let b: f64;
    let sq: f64;
    let biga: f64;
    let bigb: f64;
    let mut x: f64;
    let mut x2: f64;
    let mut f: f64;
    let mut fp: f64;
    let mut dx: f64;
    let diff: f64;
    let a0: f64;
    let a1: f64;
    let b1: f64;
    
    let imax: i32 = 10;
    let tiny: f64 = 4.0e-15;  // A small number
    let a3: f64 = 1037836800.;
    let a5: f64 = 51891840.;
    let a7: f64 = 1235520.;
    let a9: f64 = 17160.;
    let a11: f64 = 156.;
    
    let b3 = 3.*a3;
    let  b5 = 5.*a5;
    let  b7 = 7.*a7;
    let  b9 = 9.*a9;
    let  b11 = 11.*a11;
  
    // Function to solve "Kepler's eqn" for F (here called
    // x) for given e and CAPN. Only good for smallish CAPN 
  
    iflag = 0;
    if capn < 0. {
       iflag = 1;
       capn = -capn;
    }
  
    a1 = 6227020800. * (1. - 1./e);
    a0 = -6227020800.*capn/e;
    b1 = a1;
  
    //  Set iflag nonzero if capn < 0., in which case solve for -capn
    //  and change the sign of the final answer for F.
    //  Begin with a reasonable guess based on solving the cubic for small F  
  
  
    a = 6.*(e-1.)/e;
    b = -6.*capn/e;
    sq = (0.25*b*b +a*a*a/27.).sqrt();
    biga = (-0.5*b + sq).powf(0.3333333333333333);
    bigb = -(0.5*b + sq).powf(0.3333333333333333);
    x = biga + bigb;
    //  write(6,*) 'cubic = ',x**3 +a*x +b
    orbel_flon = x;
    // If capn is tiny (or zero) no need to go further than cubic even for
    // e =1.
    if capn >= tiny {
        let mut converge : bool = false;
        for _ in 1..imax {
           x2 = x*x;
           f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))));
           fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.*x2)))));
           dx = -f/fp;
           orbel_flon = x + dx;
           //   If we have converged here there's no point in going on
           if dx.abs() < tiny {
                converge = true;
                break;
           }
           x = orbel_flon;
        }
      
        if !converge {
            // Abnormal return here - we've gone thru the loop 
            // IMAX times without convergence
            if iflag == 1 {
               orbel_flon = -orbel_flon;
               capn = -capn;
            }
            println!("[WARNING {} UTC] FLON : RETURNING WITHOUT COMPLETE CONVERGENCE", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap());
            diff = e * orbel_flon.sinh()  - orbel_flon - capn;
            println!("N, F, ecc * F.sinh() - F - N : ");
            println!("{} {} {}", capn,orbel_flon,diff);
            return (orbel_flon, capn);
        }
    }
  
    //  Normal return here, but check if capn was originally negative
    if iflag == 1 {
       orbel_flon = -orbel_flon;
       capn = -capn;
    }
  
    return (orbel_flon, capn);
}


fn kepler_solution_for_a_hyperbola_hybrid_approach(e: f64, capn: f64) -> f64 {
    // ! Based on the implementation of Duncan in Mercury (orbel_fget)
    // Solves Kepler's eqn. for hyperbola using hybrid approach. 
    // ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
    //              Cel. Mech. ".  Quartic convergence from Danby's book.

    // Inputs:
    // e = eccentricty
    // capn = hyperbola mean anomaly
    
    // Output
    let mut orbel_fget: f64; // eccentric anomaly
  
    //...  Internals:
    let imax: i32 = 10;
    let tiny: f64 = 4.0e-15;  // A small number

    let tmp: f64;
    let mut x: f64;
    let mut shx: f64;
    let mut chx: f64;
    let mut esh: f64;
    let mut ech: f64;
    let mut f: f64;
    let mut fp: f64;
    let mut fpp: f64;
    let mut fppp: f64;
    let mut dx: f64;
  
    //----
    // Function to solve "Kepler's eqn" for F (here called
    // x) for given e and CAPN. 
  
    //  begin with a guess proposed by Danby  
    if capn < 0. {
       tmp = -2.*capn/e + 1.8;
       x = - tmp.log(std::f64::consts::E);
    } else {
       tmp = 2.*capn/e + 1.8;
       x = tmp.log(std::f64::consts::E);
    }
  
    orbel_fget = x;
  
    for _ in 1..imax {
       shx = x.sinh();
       chx = x.cosh();
       esh = e*shx;
       ech = e*chx;
       f = esh - x - capn;
       fp = ech - 1.;
       fpp = esh;
       fppp = ech;
       dx = -f/fp;
       dx = -f/(fp + dx*fpp/2.);
       dx = -f/(fp + dx*fpp/2. + dx*dx*fppp/6.);
       orbel_fget = x + dx;
       //   If we have converged here there's no point in going on
       if dx.abs() <= tiny {
        return orbel_fget;
       }
       x = orbel_fget;
    }
  
    println!("[WARNING {} UTC] FGET : RETURNING WITHOUT COMPLETE CONVERGENCE", OffsetDateTime::now_utc().format(&format_description::parse("[year].[month].[day] [hour]:[minute]:[second]").unwrap()).unwrap());
    return orbel_fget;
}


#[allow(dead_code)]
fn mco_sine(x0: f64) -> (f64, f64, f64) {
    // ! Based on the implementation of Chambers in Mercury
    // Calculates sin and cos of an angle X (in radians).
    // Chambers says: "x is modified by the routine. Without this part, outputs of mercury are different, don't ask me why."

    let mut x = x0;
    
    // TODO why results of this routine are different from simple calls of intrinsec cos() and sin()
    if x > 0. {
       x = modulus(x, TWO_PI);
    } else {
       x = modulus(x, TWO_PI + TWO_PI);
    }
    
    let cx = x.cos();
    let sx: f64;
    
    if x > PI {
       sx = -(1.0 - cx*cx).sqrt();
    } else {
       sx =  (1.0 - cx*cx).sqrt();
    }
    
    return (x, sx, cx);
}

fn modulus(a: f64, b: f64) -> f64 {
    //println!("Modulus: {}", -21.0f64 % 4.0f64);         // -1 because -21 divided by 4 gives -5 with a remainder of -1.
    //println!("Modulus: {}", 21.0f64 % 4.0f64);          //  1
    //println!("Modulus: {}", modulus(-21.0f64, 4.0f64)); //  3 because -21 + 4 x 6 is 3.
    //println!("Modulus: {}", modulus(21.0f64, 4.0f64));  //  1
    return a - (a / b).floor() * b;
}

pub fn linear_interpolation(target_x: f64, x: &[f64], y: &[f64]) -> (f64, usize) {
    let target_y : f64;
    
    // Find the nearest interval [ x(LEFT), x(RIGHT) ] to target_x.
    let (left, right) = find_indices_around_target_value(&x, target_x);

    if left == right {
        // Target value out of range, use limit values
        target_y = y[left];
    } else {
        //// Interpolate
        // Linear
        //target_y = (y[left] * (x[right] - target_x) + y[right] * (target_x - x[left])) / (x[right] - x[left])
        // Linear (alternative)
        let x_left = x[left];
        let target_percent = (target_x - x_left)/(x[right] - x_left); // Transform target to percent as in transforming x[left]..x[right] to 0..1
        target_y = y[left] * (1. - target_percent) + y[right] * target_percent;
    }

    return (target_y, left);
}

pub fn cosine_interpolation(target_x: f64, x: &[f64], y: &[f64]) -> (f64, usize) {
    let target_y : f64;
    
    // Find the nearest interval [ x(LEFT), x(RIGHT) ] to target_x.
    let (left, right) = find_indices_around_target_value(&x, target_x);

    if left == right {
        // Target value out of range, use limit values
        target_y = y[left];
    } else {
        //// Interpolate
        // Cosine interpolate (http://paulbourke.net/miscellaneous/interpolation/)
        // - Smooth around the real data points (contrary to the linear interpolation)
        let x_left = x[left];
        let mut target_percent = (target_x - x_left)/(x[right] - x_left); // Transform target to percent as in transforming x[left]..x[right] to 0..1
        target_percent = (1. - (target_percent*PI).cos()) / 2.; // Transform target percent so that it gets smoothed when close to 0 or 1 (i.e., closer to x[left] or x[right])
        target_y = y[left] * (1. - target_percent) + y[right] * target_percent;
    }

    return (target_y, left);
}

fn find_indices_around_target_value(data: &[f64], target_value: f64) -> (usize, usize) {
    // Find the nearest interval [ x(LEFT), x(RIGHT) ] to XVAL.
    let ndata = data.len();
    let last_idx = ndata.wrapping_sub(1);
    let (left_idx, right_idx) = match data.iter().position(|&r| r > target_value) {
        None => {
                    if data[last_idx] > target_value {
                        (0, 0)
                    } else {
                        (last_idx, last_idx)
                    }
                }
        Some(i) => {
                if i == 0 {
                    (0, 0)
                } else if i == ndata { 
                    (last_idx, last_idx)
                } else {
                    (i-1, i)
                }
            },
    };
    return (left_idx, right_idx)
}

pub fn calculate_pseudo_synchronization_period(semi_major_axis: f64, eccentricity: f64, star_mass: f64, planet_mass: f64) -> f64 {
    let alpha = (1.+15./2.*eccentricity.powi(2)+45./8.*eccentricity.powi(4)+5./16.
                 * eccentricity.powi(6))*1./(1.+3.*eccentricity.powi(2)+3./8.
                 * eccentricity.powi(4))*1./(1.-eccentricity.powi(2)).powf(1.5);
    let pseudo_rot = alpha * (G_SI*M_SUN*(star_mass+planet_mass)).sqrt();
    let angular_frequency = pseudo_rot * (semi_major_axis*AU).powf(-3./2.) * HOUR * 24.; // days^-1
    let pseudo_synchronization_period = TWO_PI/(angular_frequency); // days
    pseudo_synchronization_period
}

pub fn calculate_spin(angular_frequency: f64, inclination: f64, obliquity: f64, longitude_ascending_node: f64) -> Axes {
    let mut spin = Axes{x:0., y:0., z:0. };
    // Spin taking into consideration the inclination:
    spin.x = angular_frequency * (obliquity+inclination).sin() * (longitude_ascending_node).sin();
    spin.y = -angular_frequency * (obliquity+inclination).sin() * (longitude_ascending_node).cos();
    spin.z = angular_frequency * (obliquity+inclination).cos();

    spin
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calculate_keplerian_orbital_elements() {
        //---- Star (central body)
        let star_mass: f64 = 0.08; // Solar masses
        let planet_mass: f64 = 1.0 * M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)

        ////////// Specify initial position and velocity for a stable orbit
        ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
        let a: f64 = 0.018;                             // semi-major axis (in AU)
        let e: f64 = 0.1;                               // eccentricity
        let i: f64 = 5. * DEG2RAD;                      // inclination (degrees)
        let mut p: f64 = 0.;                            // argument of pericentre (degrees)
        let n: f64 = 0. * DEG2RAD;                      // longitude of the ascending node (degrees)
        let l: f64 = 0. * DEG2RAD;                      // mean anomaly (degrees)
        p = (p + n) * DEG2RAD;                          // Convert to longitude of perihelion !!
        let q = a * (1.0 - e);                          // perihelion distance
        let gm: f64 = G*(planet_mass+star_mass);
        let (x, y, z, vx, vy, vz) = calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

        assert_eq!(x, 0.0162);
        assert_eq!(y, 0.000000000000006571920583533768);
        assert_eq!(z, 0.0000000000000005749685486518799);
        assert_eq!(vx, -0.000000000000014818017716591765);
        assert_eq!(vy, 0.03987438104619194);
        assert_eq!(vz, 0.003488556306654768);
    }
}
