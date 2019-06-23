extern crate time;
use std;
use super::constants::*;
use super::particles::Axes;
//use std::cmp::Ord;

pub fn calculate_keplerian_orbital_elements(gm: f64, position: Axes, velocity: Axes) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
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
    let a: f64; // semi-major axis (in AU)
    let q: f64; // perihelion distance
    let eccentricity: f64; // eccentricity
    let mut i: f64 = 0.; // inclination
    let mut p: f64; // longitude of perihelion (NOT argument of perihelion!!)
    let mut n: f64; // longitude of ascending node
    let mut l: f64; // mean anomaly (or mean longitude if eccentricity < 1.e-8)
    let orbital_period: f64; // orbital_period (in days)
    // Local
    let hx = y * w  -  z * v;
    let hy = z * u  -  x * w;
    let hz = x * v  -  y * u;
    let h2 = hx.powf(2.) + hy.powf(2.) + hz.powf(2.);
    let v2 = u * u  +  v * v  +  w * w;
    let rv = x * u  +  y * v  +  z * w;
    let r = (x*x + y*y + z*z).sqrt();
    let h = h2.sqrt();
    let s = h2 / gm;

    // Semi-major axis
    a  = gm * r / (2.0 * gm  -  r * v2);

  
    // Inclination and node
    let ci = hz / h;
    if ci.abs() < 1. {
        i = ci.acos();
        n = hx.atan2(-hy);
        if n < 0. {
            n = n + TWO_PI;
        }
    } else {
        if ci > 0. {
            i = 0.;
        }
        if ci < 0. { 
            i = PI;
        }
        n = 0.;
    }
  
    // Eccentricity and perihelion distance
    let temp = 1.  +  s * (v2 / gm  -  2. / r);
    if temp <= 0. {
        eccentricity = 0.;
    } else {
        eccentricity = temp.sqrt();
    }
    q = s / (1. + eccentricity);
  
    // True longitude
    let mut true_longitude: f64;
    if hy != 0. {
        let to = -hx/hy;
        let temp = (1. - ci) * to;
        let tmp2 = to * to;
        true_longitude = (y*(1.+tmp2*ci)-x*temp).atan2(x*(tmp2+ci)-y*temp);
    } else {
        true_longitude = (y * ci).atan2(x);
    }
    if ci < 0. {
        true_longitude = true_longitude + PI;
    }
  
    if eccentricity < 3.0e-8 {
        p = 0.;
        l = true_longitude;
    } else {
        let mut ce = (v2*r - gm) / (eccentricity*gm);

        // Mean anomaly for ellipse
        if eccentricity < 1. {
            if ce.abs() > 1. {
                ce = ce.signum();
            }
            let mut bige = ce.acos();
            if rv < 0. {
                bige = TWO_PI - bige;
            }
            l = bige - eccentricity*bige.sin();
        } else {
            // Mean anomaly for hyperbola
            if ce < 1. {
                ce = 1.
            }
            let mut bige = (ce + (ce*ce-1.).sqrt()).log(std::f64::consts::E);
            if rv < 0. {
                bige = - bige;
            }
            l = eccentricity * bige.sinh() - bige
        }

        // Longitude of perihelion
        let mut cf = (s - r) / (eccentricity*r);
        if cf.abs() > 1. {
            cf = cf.signum();
        }
        let mut f: f64 = cf.acos();
        if rv < 0. {
            f = TWO_PI - f;
        }
        p = true_longitude - f;
        p = modulus(p + TWO_PI + TWO_PI, TWO_PI);
    }
  
    if l < 0. {
        l = l + TWO_PI;
    }
    if l > TWO_PI {
        l = modulus(l, TWO_PI);
    }

    // Given relative coordinates and velocities (of the body 1 respect to body 2),
    // and GM = G times the sum of the masses (body 1 + body 2)
    orbital_period = (TWO_PI / gm.sqrt()) * a.powf(3./2.); // in days

    return (a, q, eccentricity, i, p, n, l, orbital_period)
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
            println!("[WARNING {} UTC] FLON : RETURNING WITHOUT COMPLETE CONVERGENCE", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
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
  
    println!("[WARNING {} UTC] FGET : RETURNING WITHOUT COMPLETE CONVERGENCE", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
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

pub fn interpolate_b_spline(tdata: &[f64], ydata: &[f64], tval: f64) -> (f64, usize) {
    //// OPTIMIZATION: Return last left to reduce future searches
    // Based on: 
    //
    //    SPLINE_B_VAL evaluates a cubic B spline approximant.
    //
    //  Discussion:
    //
    //    The cubic B spline will approximate the data, but is not
    //    designed to interpolate it.
    //
    //    In effect, two "phantom" data values are appended to the data,
    //    so that the spline will interpolate the first and last data values.
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    14 August 2005
    //
    //  Author:
    //
    //    John Burkardt
    //
    //  Reference:
    //
    //    Carl de Boor,
    //    A Practical Guide to Splines,
    //    Springer Verlag, 1978.

    let mut yval : f64 = 0.;
    
    // Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    let (left, right) = find_indices_around_target_value(&tdata, tval);

    if left == right {
        // Target value out of range, use limit values
        yval = ydata[left];
    } else {
        // Evaluate the 5 nonzero B spline basis functions in the interval,
        // weighted by their corresponding data values.
        let u = (tval - tdata[left]) / (tdata[right] - tdata[left]);
        
        // B function associated with node LEFT - 1, (or "phantom node"),
        // evaluated in its 4th interval.
        let bval = ( 1.0 - 3.0 * u + 3.0 * u.powi(2) - u.powi(3) ) / 6.0;
        if left > 0 {
           yval = yval + ydata[left - 1] * bval;
        } else {
           yval = yval + ( 2.0 * ydata[0] - ydata[1] ) * bval;
        }
        
        // B function associated with node LEFT,
        // evaluated in its third interval.
        let bval = ( 4.0 - 6.0 * u.powi(2) + 3.0 * u.powi(3) ) / 6.0;
        yval = yval + ydata[left] * bval;
        
        // B function associated with node RIGHT,
        // evaluated in its second interval.
        let bval = ( 1.0 + 3.0 * u + 3.0 * u.powi(2) - 3.0 * u.powi(3) ) / 6.0;
        yval = yval + ydata[right] * bval;
        
        // B function associated with node RIGHT+1, (or "phantom node"),
        // evaluated in its first interval.
        let bval = u.powi(3) / 6.0;
        let ndata = ydata.len();
        if right + 1 < ndata {
           yval = yval + ydata[right + 1] * bval;
        } else if ndata >= 2 {
           yval = yval + ( 2.0 * ydata[ndata - 1] - ydata[ndata - 2] ) * bval;
        }

    }

    return (yval, left);
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

pub fn calculate_spin(angular_frequency: f64, inclination: f64, obliquity: f64) -> Axes {
    let mut spin = Axes{x:0., y:0., z:0. };
    // Spin taking into consideration the inclination:
    spin.x = 0.;
    spin.y = -angular_frequency * (obliquity+inclination).sin();
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

fn calculate_temperature_disk(planet_distance: f64, star_mass: f64) -> f64 {
    // planet_distance in AU, star_mass in M_SUN
    let t_ref = 280.0; //K
    let tdisk = t_ref * planet_distance.powf(-0.5) * star_mass; // K
    tdisk
}


fn calculate_disk_surface_density(time: f64, planet_distance: f64, disk_surface_density_normalization: f64, disk_inner_edge_distance: f64, disk_outer_edge_distance: f64, disk_lifetime: f64) -> f64{
    let initial_disk_surface_density = disk_surface_density_normalization * planet_distance.powi(-1)
        * (-1.0 * planet_distance/disk_outer_edge_distance).exp() * (1.0 - (disk_inner_edge_distance/planet_distance).sqrt()); // Unit of disk_surface_density_normalization

    let mut disk_surface_density = 0.0;
    if planet_distance > disk_inner_edge_distance {
        disk_surface_density = initial_disk_surface_density * (-time/disk_lifetime).exp();
    }
    disk_surface_density
}

fn calculate_keplerian_velocity(planet_distance: f64, star_mass_g: f64) -> f64 {
    let keplerian_frequency = (star_mass_g / (planet_distance).powi(3)).sqrt(); // days^-1
    keplerian_frequency
}

fn calculate_disk_viscosity(planet_distance: f64, star_mass: f64, alpha_disk: f64, disk_mean_molecular_weight: f64) -> f64{
    let star_mass_g = G * star_mass;
    let keplerian_frequency = calculate_keplerian_velocity(planet_distance, star_mass_g);
    let disk_temperature = calculate_temperature_disk(planet_distance, star_mass);
    let speed_of_sound_squared = BOLTZMANN_CONSTANT * disk_temperature / (disk_mean_molecular_weight * MASS_HYDROGEN_ATOM);
    let viscosity = alpha_disk * speed_of_sound_squared / keplerian_frequency; // AU^2.days^-1
    viscosity
}

fn calculate_gas_radial_velocity(planet_distance: f64, alpha_disk: f64, star_mass: f64, disk_mean_molecular_weight: f64, planet_semi_major_axis: f64) -> f64{
    let viscosity = calculate_disk_viscosity(planet_distance, star_mass, alpha_disk, disk_mean_molecular_weight);
    let gas_radial_velocity = -3.0 * viscosity / (2.0 * planet_semi_major_axis); // AU.days^-1
    gas_radial_velocity
}

pub fn calculate_migration_timescale(time: f64, planet_distance: f64, planet_semi_major_axis: f64, disk_surface_density_normalization: f64
                                     , disk_inner_edge_distance: f64, disk_outer_edge_distance: f64, disk_lifetime: f64
                                     , alpha_disk: f64, disk_mean_molecular_weight: f64, planet_mass: f64, star_mass: f64) -> f64 {
    let gas_radial_velocity = calculate_gas_radial_velocity(planet_distance, alpha_disk, star_mass, disk_mean_molecular_weight, planet_semi_major_axis);
    let disk_surface_density = calculate_disk_surface_density(time, planet_distance, disk_surface_density_normalization, disk_inner_edge_distance, disk_outer_edge_distance, disk_lifetime);

    let x = 1.0_f64;
    let a_dot = gas_radial_velocity * x.min(2.0 * disk_surface_density * planet_semi_major_axis.powi(2) / planet_mass);
    let timescale = -planet_semi_major_axis / a_dot;
    timescale
}
