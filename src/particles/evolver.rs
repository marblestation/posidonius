use super::super::tools::{interpolate_b_spline};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum EvolutionType {
    GalletBolmont2017(f64), // SolarLike EvolvingDissipation Evolving dissipation
    BolmontMathis2016(f64), // SolarLike EvolvingDissipation Evolving dissipation
    Baraffe2015(f64), // NEW
    Leconte2011(f64), // BrownDwarf
    Baraffe1998(f64), // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    LeconteChabrier2013, // Jupiter
    NonEvolving,

}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct Evolver {
    pub evolution_type: EvolutionType,
    pub time: Vec<f64>,
    pub radius: Vec<f64>,
    pub radius_of_gyration_2: Vec<f64>,
    pub love_number: Vec<f64>,
    pub inverse_tidal_q_factor: Vec<f64>, // Bolmont & Mathis 2016
    left_index: usize,
}

// NOTE: This is a big optimization to reduce the cost of cloning
//       Cloned particles are downgraded to a NonEvolving status
//       to avoid cloning model's data which can be big and expensive.
impl Clone for Evolver {
    fn clone(&self) -> Self {
        Evolver {
            evolution_type:EvolutionType::NonEvolving,
            time:vec![],
            radius:vec![],
            radius_of_gyration_2:vec![],
            love_number:vec![],
            inverse_tidal_q_factor:vec![],
            left_index: 0,
        }

    }
}

impl Evolver {
    // OPTIMIZATION: Skip first N elements which belong to the past
    fn idx(&self) -> usize {
        if self.left_index > 0 {
            return self.left_index - 1;
        } else {
            return 0;
        }
    }

    pub fn radius(&mut self, current_time: f64, current_radius: f64) -> f64 {
        let (new_radius, left_index) = match self.evolution_type {
            EvolutionType::NonEvolving => { (current_radius, 0) },
            _ => { interpolate_b_spline(&self.time[self.idx()..], &self.radius[self.idx()..], current_time) }
        };
        self.left_index += left_index;
        return new_radius;
    }

    pub fn radius_of_gyration_2(&mut self, current_time: f64, current_radius_of_gyration_2: f64) -> f64 {
        let (new_radius_of_gyration_2, left_index) = match self.evolution_type {
            EvolutionType::Baraffe2015(_) => { interpolate_b_spline(&self.time[self.idx()..], &self.radius_of_gyration_2[self.idx()..], current_time) },
            EvolutionType::Leconte2011(_) => { interpolate_b_spline(&self.time[self.idx()..], &self.radius_of_gyration_2[self.idx()..], current_time) },
            EvolutionType::LeconteChabrier2013 => { interpolate_b_spline(&self.time[self.idx()..], &self.radius_of_gyration_2[self.idx()..], current_time) },
            _ => { (current_radius_of_gyration_2, 0) }
        };
        self.left_index += left_index;
        return new_radius_of_gyration_2;
    }

    pub fn love_number(&mut self, current_time: f64, current_love_number: f64) -> f64 {
        let (new_love_number, left_index) = match self.evolution_type {
            EvolutionType::LeconteChabrier2013 => { interpolate_b_spline(&self.time[self.idx()..], &self.love_number[self.idx()..], current_time) },
            _ => { (current_love_number, 0) }
        };
        self.left_index += left_index;
        return new_love_number;
    }

    pub fn inverse_tidal_q_factor(&mut self, current_time: f64, current_inverse_tidal_q_factor: f64) -> f64 {
        // The planet excite/induce waves in the convective region of the star
        // and the Q factor is a quality factor that describes that loss of energy.
        //
        // The inverse_tidal_q_factor was derived from analytical formulas using 
        // parameters from stellar models. This factor is used to derive the 
        // dissipation factor of the star, which depends on the excitation 
        // frequency of the planet (the frequency at what we would see the planet
        // if we were sitting at the surface of the star).
        //
        //  - If the star does not rotate, the excitation frequency is the orbital period (i.e.,
        //  frequency)
        //  - If the planet is very far away, the excitation frequency is the spin
        //  - Planets are treated as source points and their rotation is not important here
        //
        // The tidal theory (and this code) assumes circular orbits because they are easier.
        // Excentric orbits needs more than one frequency to be described and it will be
        // included in future versions of this code.
        let (new_inverse_tidal_q_factor, left_index) = match self.evolution_type {
            EvolutionType::BolmontMathis2016(_) | EvolutionType::GalletBolmont2017(_) => {
                interpolate_b_spline(&self.time[self.idx()..], &self.inverse_tidal_q_factor[self.idx()..], current_time)
            },
            _ => (current_inverse_tidal_q_factor, 0),
        };
        self.left_index += left_index;
        return new_inverse_tidal_q_factor;
    }
}


