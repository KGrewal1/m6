use rand::{Rng, SeedableRng};
use rand_distr::{uniform::SampleUniform, Distribution, Standard, Uniform};
use rand_xoshiro::Xoshiro256StarStar;

pub struct MonadicRng(Xoshiro256StarStar);

impl MonadicRng {
    pub fn new(seed: u64) -> Self {
        MonadicRng(Xoshiro256StarStar::seed_from_u64(seed))
    }

    pub fn gen_val<T>(mut self) -> (T, MonadicRng)
    where
        Standard: Distribution<T>,
    {
        let val = self.0.gen();
        (val, self)
    }
}

pub struct UniformMonad<X: SampleUniform>(Uniform<X>);

impl<X: SampleUniform> UniformMonad<X> {
    pub fn new(dist: Uniform<X>) -> Self {
        UniformMonad(dist)
    }

    pub fn sample(&self, mut rng: MonadicRng) -> (X, MonadicRng) {
        let res = self.0.sample(&mut rng.0);
        (res, rng)
    }
}
