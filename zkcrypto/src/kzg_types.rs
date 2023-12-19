use crate::consts::{G1_GENERATOR, G1_IDENTITY, G1_NEGATIVE_GENERATOR, G2_GENERATOR, G2_NEGATIVE_GENERATOR,
                    MODULUS, R2, SCALE2_ROOT_OF_UNITY};
use crate::fft_g1::g1_linear_combination;
use crate::kzg_proofs::{
    expand_root_of_unity, pairings_verify, FFTSettings as ZFFTSettings, KZGSettings as ZKZGSettings,
};
use crate::poly::PolyData;
use crate::utils::{blst_fr_into_pc_fr, blst_p1_into_pc_g1projective, blst_p2_into_pc_g2projective,
                   montgomery_reduce, pc_fr_into_blst_fr, pc_g1projective_into_blst_p1, pc_g2projective_into_blst_p2};
// use bls12_381::{G1Affine, G1Projective, G2Affine, G2Projective, Scalar, MODULUS, R2};
use pairing_ce::bls12_381::{G1Affine, G1Prepared, G1 as G1Projective, G2Affine, G2 as G2Projective, Fr as Scalar};
use pairing_ce::bls12_381::fr::{FrRepr};
use pairing_ce::ff::{Field, PrimeField, PrimeFieldRepr, SqrtField};

use blst::{blst_fr, blst_p1};
// use ff::Field;   // imported from pairing_ce::ff
use kzg::common_utils::reverse_bit_order;
use kzg::eip_4844::{BYTES_PER_FIELD_ELEMENT, BYTES_PER_G1, BYTES_PER_G2};
use kzg::{
    FFTFr, FFTSettings, Fr as KzgFr, G1Mul, G2Mul, KZGSettings, PairingVerify, Poly, G1, G2,
};
use std::ops::{Mul, Sub};
use ff::derive::bitvec::macros::internal::funty::Fundamental;

use ff::derive::sbb;
use subtle::{ConstantTimeEq, CtOption};
use pairing_ce::{CurveProjective, GenericCurveProjective};
use rand::Rng;

fn to_scalar(zfr: &ZFr) -> Scalar {
    zfr.fr
}

fn bigint_check_mod_256(a: &[u64; 4]) -> bool {
    let (_, overflow) = a[0].overflowing_sub(MODULUS.into_repr().0[0]);
    let (_, overflow) = a[1].overflowing_sub(MODULUS.into_repr().0[1] + overflow as u64);
    let (_, overflow) = a[2].overflowing_sub(MODULUS.into_repr().0[2] + overflow as u64);
    let (_, overflow) = a[3].overflowing_sub(MODULUS.into_repr().0[3] + overflow as u64);
    overflow
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Default)]
pub struct ZFr {
    pub fr: Scalar,
}

impl ZFr {
    pub fn from_blst_fr(fr: blst_fr) -> Self {
        Self {
            fr: blst_fr_into_pc_fr(fr),
        }
    }
    pub fn to_blst_fr(&self) -> blst_fr {
        pc_fr_into_blst_fr(self.fr)
    }

    pub fn converter(points: &[ZFr]) -> Vec<Scalar> {
        let mut result = Vec::new();

        for zg1 in points {
            result.push(zg1.fr);
        }
        result
    }
}

impl KzgFr for ZFr {
    fn null() -> Self {
        Self {
            // fr: Scalar([u64::MAX, u64::MAX, u64::MAX, u64::MAX]),
            fr: Scalar::zero(),
        }
    }
    fn zero() -> Self {
        Self::from_u64(0)
    }

    fn one() -> Self {
        Self::from_u64(1)
    }

    #[cfg(feature = "rand")]
    fn rand() -> Self {
        let tmp: i64 = rand::thread_rng().gen();
        let result = FrRepr::from(tmp.as_u64());
        Self { fr: Scalar::from_repr(result).unwrap() }
    }
    // fn rand() -> Self {
    //     let rng = rand::thread_rng();
    //     let result = ff::Field::random(rng);
    //     Self { fr: result }
    // }
    #[allow(clippy::bind_instead_of_map)]
    // rewrite function BE format referencing pairing_ce function
    fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        bytes
            .try_into()
            .map_err(|_| {
                format!(
                    "Invalid byte length. Expected {}, got {}",
                    BYTES_PER_FIELD_ELEMENT,
                    bytes.len()
                )
            })
            .and_then(|bytes: &[u8; BYTES_PER_FIELD_ELEMENT]| {
                let mut tmp = FrRepr([0, 0, 0, 0]);

                tmp.0[0] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
                tmp.0[1] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
                tmp.0[2] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
                tmp.0[3] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());

                // Try to subtract the modulus
                let (_, borrow) = sbb(tmp.0[0], MODULUS.into_repr().0[0], 0);
                let (_, borrow) = sbb(tmp.0[1], MODULUS.into_repr().0[1], borrow);
                let (_, borrow) = sbb(tmp.0[2], MODULUS.into_repr().0[2], borrow);
                let (_, _borrow) = sbb(tmp.0[3], MODULUS.into_repr().0[3], borrow);

                let mut tmp2 = FrRepr::default();

                tmp2.0[0] = tmp.0[3];
                tmp2.0[1] = tmp.0[2];
                tmp2.0[2] = tmp.0[1];
                tmp2.0[3] = tmp.0[0];

                let is_zero: bool = tmp2.is_zero().into();
                if !is_zero && !bigint_check_mod_256(&tmp2.0) {
                    return Err("Invalid scalar".to_string());
                }

                let mut res = Scalar::from_repr(tmp2).unwrap();
                res.mul_assign(&R2);
                Ok(Self { fr: res })
            })
    }

    // rewrite function BE format referencing pairing_ce function
    fn from_bytes_unchecked(bytes: &[u8]) -> Result<Self, String> {
        bytes
            .try_into()
            .map_err(|_| {
                format!(
                    "Invalid byte length. Expected {}, got {}",
                    BYTES_PER_FIELD_ELEMENT,
                    bytes.len()
                )
            })
            .map(|bytes: &[u8; BYTES_PER_FIELD_ELEMENT]| {
                // let mut tmp = Scalar([0, 0, 0, 0]);
                let mut tmp = FrRepr([0, 0, 0, 0]);

                tmp.0[0] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap());
                tmp.0[1] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap());
                tmp.0[2] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap());
                tmp.0[3] = u64::from_be_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap());

                // Try to subtract the modulus
                let (_, borrow) = sbb(tmp.0[0], MODULUS.into_repr().0[0], 0);
                let (_, borrow) = sbb(tmp.0[1], MODULUS.into_repr().0[1], borrow);
                let (_, borrow) = sbb(tmp.0[2], MODULUS.into_repr().0[2], borrow);
                let (_, _borrow) = sbb(tmp.0[3], MODULUS.into_repr().0[3], borrow);

                // let mut tmp2 = Scalar::default();
                let mut tmp2 = FrRepr::default();

                tmp2.0[0] = tmp.0[3];
                tmp2.0[1] = tmp.0[2];
                tmp2.0[2] = tmp.0[1];
                tmp2.0[3] = tmp.0[0];

                let mut res = Scalar::from_repr(tmp2).unwrap();
                res.mul_assign(&R2);
                Self { fr: res }
            })
    }

    fn from_hex(hex: &str) -> Result<Self, String> {
        let bytes = hex::decode(&hex[2..]).unwrap();
        Self::from_bytes(&bytes)
    }

    // NOTE: untested
    fn from_u64_arr(u: &[u64; 4]) -> Self {
        let tmp = FrRepr(*u);
        Self {
            // fr: Scalar::from_raw(*u),
            fr: Scalar::from_repr(tmp).unwrap(),
        }
    }

    // NOTE: untested
    fn from_u64(val: u64) -> Self {
        let tmp = FrRepr::from(val);
        Self {
            // fr: Scalar::try_from(val).unwrap(),
            fr: Scalar::from_repr(tmp).unwrap(),
        }
    }

    fn to_bytes(&self) -> [u8; 32] {
        // let scalar = self.fr;
        let scalar = self.fr;
        let tmp = montgomery_reduce(
            scalar.into_repr().0[0],
            scalar.into_repr().0[1],
            scalar.into_repr().0[2],
            scalar.into_repr().0[3],
            0,
            0,
            0,
            0,
        );
        let mut res = [0; 32];
        res[0..8].copy_from_slice(&tmp.into_repr().0[3].to_be_bytes());
        res[8..16].copy_from_slice(&tmp.into_repr().0[2].to_be_bytes());
        res[16..24].copy_from_slice(&tmp.into_repr().0[1].to_be_bytes());
        res[24..32].copy_from_slice(&tmp.into_repr().0[0].to_be_bytes());
        res
    }

    //testuoti
    fn to_u64_arr(&self) -> [u64; 4] {
        let bytes = self.to_bytes();
        [
            u64::from_be_bytes(bytes[24..32].try_into().unwrap()),
            u64::from_be_bytes(bytes[16..24].try_into().unwrap()),
            u64::from_be_bytes(bytes[8..16].try_into().unwrap()),
            u64::from_be_bytes(bytes[0..8].try_into().unwrap()),
        ]
    }

    // Note: untested
    fn is_one(&self) -> bool {
        self.fr.eq(&ZFr::one().fr) == true
        // self.fr.ct_eq(&ZFr::one().fr).unwrap_u8() == 1
    }

    fn is_zero(&self) -> bool {
        self.fr.is_zero() == true
    }

    // NOTE: untested
    fn is_null(&self) -> bool {
        self.fr.eq(&ZFr::null().fr) == true
        // self.fr.ct_eq(&ZFr::null().fr).unwrap_u8() == 1
    }

    // NOTE: untested
    fn sqr(&self) -> Self {
        Self {
            fr: self.fr.sqrt().unwrap(),
        }
    }

    // NOTE: untested
    fn mul(&self, b: &Self) -> Self {
        let mut c = *self;
        c.fr.mul_assign(&b.fr);
        // Self {
        //     fr: Scalar::mul(&to_scalar(self), &to_scalar(b)),
        // }
        c
    }

    // NOTE: untested
    fn add(&self, b: &Self) -> Self {
        let mut c = *self;
        c.fr.add_assign(&b.fr);
        // Self { fr: self.fr + b.fr }
        c
    }

    // NOTE: untested
    fn sub(&self, b: &Self) -> Self {
        let mut c = *self;
        c.fr.sub_assign(&b.fr);
        // Self { fr: self.fr - b.fr }
        c
    }

    // NOTE: untested
    fn eucl_inverse(&self) -> Self {
        Self {
            fr: self.fr.inverse().unwrap(),
        }
        // Self {
        //     fr: self.fr.invert().unwrap(),
        // }
    }

    // NOTE: untested
    fn negate(&self) -> Self {
        let mut c = *self;
        c.fr.negate();
        // Self { fr: self.fr.neg() }
        c
    }

    // NOTE: untested
    fn inverse(&self) -> Self {
        Self {
            fr: self.fr.inverse().unwrap(),
        }
        // Self {
        //     fr: self.fr.invert().unwrap(),
        // }
    }

    fn pow(&self, n: usize) -> Self {
        let mut tmp = *self;
        let mut out = Self::one();
        let mut n2 = n;

        loop {
            if n2 & 1 == 1 {
                out = out.mul(&tmp);
            }
            n2 >>= 1;
            if n2 == 0 {
                break;
            }
            tmp = tmp.sqr();
        }

        out
    }

    fn div(&self, b: &Self) -> Result<Self, String> {
        if <ZFr>::is_zero(b) {
            return Err("Cannot divide by zero".to_string());
        }
        let tmp = b.eucl_inverse();
        let out = self.mul(&tmp);
        Ok(out)
    }

    fn equals(&self, b: &Self) -> bool {
        self.fr == b.fr
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct ZG1 {
    pub proj: G1Projective,
}

// NOTE: random decision to make default a zero
impl Default for ZG1 {
    fn default() -> Self {
        Self { proj: <pairing_ce::bls12_381::G1 as CurveProjective>::zero(), }
    }
}

impl ZG1 {
    pub const fn from_blst_p1(p1: blst_p1) -> Self {
        Self {
            proj: blst_p1_into_pc_g1projective(&p1),
        }
    }

    pub const fn to_blst_p1(&self) -> blst_p1 {
        pc_g1projective_into_blst_p1(self.proj)
    }
    pub const fn from_g1_projective(proj: G1Projective) -> Self {
        Self { proj }
    }

    fn affine_to_projective(p: G1Affine) -> Self {
        Self {
            proj: G1Projective::from(p),
        }
    }
    pub fn converter(points: &[ZG1]) -> Vec<G1Projective> {
        let mut result = Vec::new();

        for zg1 in points {
            result.push(zg1.proj);
        }
        result
    }
}

impl From<blst_p1> for ZG1 {
    fn from(p1: blst_p1) -> Self {
        let proj = blst_p1_into_pc_g1projective(&p1);
        Self { proj }
    }
}

impl G1 for ZG1 {
    fn identity() -> Self {
        G1_IDENTITY
    }

    fn generator() -> Self {
        G1_GENERATOR
    }

    fn negative_generator() -> Self {
        G1_NEGATIVE_GENERATOR
    }

    // NOTE: untested
    #[cfg(feature = "rand")]
    fn rand() -> Self {
        let mut rng = rand::thread_rng();
        // Self {
        //     proj: G1Projective::random(&mut rng),
        // }

        let tmp: i64 = rand::thread_rng().gen();
        let result = FrRepr::from(tmp.as_u64());
        Self { proj: G1Projective::rand(&mut rng) }
        // Self {
        //     proj: G1Projective::from(rng_affine)
        // }
    }

    #[allow(clippy::bind_instead_of_map)]
    fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        bytes
            .try_into()
            .map_err(|_| {
                format!(
                    "Invalid byte length. Expected {}, got {}",
                    BYTES_PER_G1,
                    bytes.len()
                )
            })
            .and_then(|bytes: &[u8; BYTES_PER_G1]| {
                let affine: CtOption<G1Affine> = G1Affine::from_compressed(bytes);
                match affine.into() {
                    Some(x) => Ok(ZG1::affine_to_projective(x)),
                    None => Err("Failed to deserialize G1: Affine not available".to_string()),
                }
            })
    }

    fn from_hex(hex: &str) -> Result<Self, String> {
        let bytes = hex::decode(&hex[2..]).unwrap();
        Self::from_bytes(&bytes)
    }

    fn to_bytes(&self) -> [u8; 48] {
        let g1_affine = G1Affine::from(self.proj);
        g1_affine.to_compressed()
    }

    // NOTE: untested
    fn add_or_dbl(&mut self, b: &Self) -> Self {
        self.proj.double();
        Self
    }

    // NOTE: unchanged, should not work
    fn is_inf(&self) -> bool {
        bool::from(self.proj.is_identity())
    }
    // NOTE: unchanged, should not work
    fn is_valid(&self) -> bool {
        bool::from(self.proj.is_on_curve())
    }

    // NOTE: untested
    fn dbl(&self) -> Self {
        self.proj.doube();
        Self
    }

    // NOTE: untested
    fn add(&self, b: &Self) -> Self {
        let mut c = *self;
        c.proj.add_assign(&b.proj);
        c
    }

    // NOTE: untested
    fn sub(&self, b: &Self) -> Self {
        let mut c = *self;
        c.proj.sub_assign(&b.proj);
        c
    }

    // NOTE: unchanged, untested
    fn equals(&self, b: &Self) -> bool {
        self.proj.eq(&b.proj)
    }
}

impl G1Mul<ZFr> for ZG1 {
    // NOTE: untested
    fn mul(&self, b: &ZFr) -> Self {
        let mut c = *self;
        c.proj.mul_assign(&b.fr);
        c
    }

    fn g1_lincomb(points: &[Self], scalars: &[ZFr], len: usize) -> Self {
        let mut out = Self::default();
        g1_linear_combination(&mut out, points, scalars, len);
        out
    }
}

impl PairingVerify<ZG1, ZG2> for ZG1 {
    fn verify(a1: &ZG1, a2: &ZG2, b1: &ZG1, b2: &ZG2) -> bool {
        pairings_verify(a1, a2, b1, b2)
    }
}

#[derive(Debug, Default, PartialEq, Eq, Clone)]
pub struct ZG2 {
    pub proj: G2Projective,
}

impl ZG2 {
    pub const fn from_blst_p2(p2: blst::blst_p2) -> Self {
        Self {
            proj: blst_p2_into_pc_g2projective(&p2),
        }
    }
    pub const fn from_g2_projective(proj: G2Projective) -> Self {
        Self { proj }
    }
    pub const fn to_blst_p2(&self) -> blst::blst_p2 {
        pc_g2projective_into_blst_p2(self.proj)
    }
}

impl G2 for ZG2 {
    fn generator() -> Self {
        G2_GENERATOR
    }

    fn negative_generator() -> Self {
        G2_NEGATIVE_GENERATOR
    }

    #[allow(clippy::bind_instead_of_map)]
    fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        bytes
            .try_into()
            .map_err(|_| {
                format!(
                    "Invalid byte length. Expected {}, got {}",
                    BYTES_PER_G2,
                    bytes.len()
                )
            })
            .and_then(|bytes: &[u8; BYTES_PER_G2]| {
                let affine = G2Affine::from_compressed(bytes).unwrap();
                Ok(ZG2::from_g2_projective(G2Projective::from(affine)))
            })
    }

    fn to_bytes(&self) -> [u8; 96] {
        let g2_affine = G2Affine::from(self.proj);
        g2_affine.to_compressed()
    }

    // NOTE: untested
    fn add_or_dbl(&mut self, b: &Self) -> Self {
        self.proj.add_assign(&b.proj);
        Self
    }

    // NOTE: untested
    fn dbl(&self) -> Self {
        self.proj.doube();
        Self
    }

    // Note: untested
    fn sub(&self, b: &Self) -> Self {
        self.proj.sub_assign(&b.proj);
        // Self {
        //     proj: self.proj - b.proj,
        // }
        Self
    }

    fn equals(&self, b: &Self) -> bool {
        self.proj.eq(&b.proj)
    }
}

impl G2Mul<ZFr> for ZG2 {
    // Note: untested
    fn mul(&self, b: &ZFr) -> Self {
        // FIXME: Is this right?
        self.proj.mul_assign(&b.fr);
        // Self {
        //     proj: self.proj.mul(b.fr),
        // }
        Self
    }
}

impl Default for ZFFTSettings {
    fn default() -> Self {
        Self {
            max_width: 0,
            root_of_unity: ZFr::zero(),
            expanded_roots_of_unity: Vec::new(),
            reverse_roots_of_unity: Vec::new(),
            roots_of_unity: Vec::new(),
        }
    }
}

impl FFTSettings<ZFr> for ZFFTSettings {
    fn new(scale: usize) -> Result<Self, String> {
        if scale >= SCALE2_ROOT_OF_UNITY.len() {
            return Err(String::from(
                "Scale is expected to be within root of unity matrix row size",
            ));
        }

        // max_width = 2 ^ max_scale
        let max_width: usize = 1 << scale;
        let root_of_unity = ZFr::from_u64_arr(&SCALE2_ROOT_OF_UNITY[scale]);

        // create max_width of roots & store them reversed as well
        let expanded_roots_of_unity = expand_root_of_unity(&root_of_unity, max_width).unwrap();
        let mut reverse_roots_of_unity = expanded_roots_of_unity.clone();
        reverse_roots_of_unity.reverse();

        // Permute the roots of unity
        let mut roots_of_unity = expanded_roots_of_unity.clone();
        roots_of_unity.pop();
        reverse_bit_order(&mut roots_of_unity)?;

        Ok(Self {
            max_width,
            root_of_unity,
            expanded_roots_of_unity,
            reverse_roots_of_unity,
            roots_of_unity,
        })
    }

    fn get_max_width(&self) -> usize {
        self.max_width
    }

    fn get_expanded_roots_of_unity_at(&self, i: usize) -> ZFr {
        self.expanded_roots_of_unity[i]
    }

    fn get_expanded_roots_of_unity(&self) -> &[ZFr] {
        &self.expanded_roots_of_unity
    }

    fn get_reverse_roots_of_unity_at(&self, i: usize) -> ZFr {
        self.reverse_roots_of_unity[i]
    }

    fn get_reversed_roots_of_unity(&self) -> &[ZFr] {
        &self.reverse_roots_of_unity
    }

    fn get_roots_of_unity_at(&self, i: usize) -> ZFr {
        self.roots_of_unity[i]
    }

    fn get_roots_of_unity(&self) -> &[ZFr] {
        &self.roots_of_unity
    }
}

impl KZGSettings<ZFr, ZG1, ZG2, ZFFTSettings, PolyData> for ZKZGSettings {
    fn new(
        secret_g1: &[ZG1],
        secret_g2: &[ZG2],
        _length: usize,
        fft_settings: &ZFFTSettings,
    ) -> Result<ZKZGSettings, String> {
        Ok(Self {
            secret_g1: secret_g1.to_vec(),
            secret_g2: secret_g2.to_vec(),
            fs: fft_settings.clone(),
        })
    }

    fn commit_to_poly(&self, p: &PolyData) -> Result<ZG1, String> {
        if p.coeffs.len() > self.secret_g1.len() {
            return Err(String::from("Polynomial is longer than secret g1"));
        }

        let mut out = ZG1::default();
        g1_linear_combination(&mut out, &self.secret_g1, &p.coeffs, p.coeffs.len());

        Ok(out)
    }

    fn compute_proof_single(&self, p: &PolyData, x: &ZFr) -> Result<ZG1, String> {
        if p.coeffs.is_empty() {
            return Err(String::from("Polynomial must not be empty"));
        }

        // `-(x0^n)`, where `n` is `1`
        let divisor_0 = x.negate();

        // Calculate `q = p / (x^n - x0^n)` for our reduced case (see `compute_proof_multi` for
        // generic implementation)
        let mut out_coeffs = Vec::from(&p.coeffs[1..]);
        for i in (1..out_coeffs.len()).rev() {
            let tmp = out_coeffs[i].mul(&divisor_0);
            out_coeffs[i - 1] = out_coeffs[i - 1].sub(&tmp);
        }

        let q = PolyData { coeffs: out_coeffs };
        let ret = self.commit_to_poly(&q)?;
        Ok(ret)
    }

    fn check_proof_single(&self, com: &ZG1, proof: &ZG1, x: &ZFr, y: &ZFr) -> Result<bool, String> {
        let x_g2 = G2_GENERATOR.mul(x);
        let s_minus_x: ZG2 = self.secret_g2[1].sub(&x_g2);
        let y_g1 = G1_GENERATOR.mul(y);
        let commitment_minus_y: ZG1 = com.sub(&y_g1);

        Ok(pairings_verify(
            &commitment_minus_y,
            &G2_GENERATOR,
            proof,
            &s_minus_x,
        ))
    }

    fn compute_proof_multi(&self, p: &PolyData, x: &ZFr, n: usize) -> Result<ZG1, String> {
        if p.coeffs.is_empty() {
            return Err(String::from("Polynomial must not be empty"));
        }

        if !n.is_power_of_two() {
            return Err(String::from("n must be a power of two"));
        }

        // Construct x^n - x0^n = (x - x0.w^0)(x - x0.w^1)...(x - x0.w^(n-1))
        let mut divisor = PolyData {
            coeffs: Vec::with_capacity(n + 1),
        };

        // -(x0^n)
        let x_pow_n = x.pow(n);

        divisor.coeffs.push(x_pow_n.negate());

        // Zeros
        for _ in 1..n {
            divisor.coeffs.push(ZFr { fr: Scalar::zero() });
        }

        // x^n
        divisor.coeffs.push(ZFr { fr: Scalar::one() });

        let mut new_polina = p.clone();

        // Calculate q = p / (x^n - x0^n)
        // let q = p.div(&divisor).unwrap();
        let q = new_polina.div(&divisor)?;
        let ret = self.commit_to_poly(&q)?;
        Ok(ret)
    }

    fn check_proof_multi(
        &self,
        com: &ZG1,
        proof: &ZG1,
        x: &ZFr,
        ys: &[ZFr],
        n: usize,
    ) -> Result<bool, String> {
        if !n.is_power_of_two() {
            return Err(String::from("n is not a power of two"));
        }

        // Interpolate at a coset.
        let mut interp = PolyData {
            coeffs: self.fs.fft_fr(ys, true)?,
        };

        let inv_x = x.inverse(); // Not euclidean?
        let mut inv_x_pow = inv_x;
        for i in 1..n {
            interp.coeffs[i] = interp.coeffs[i].mul(&inv_x_pow);
            inv_x_pow = inv_x_pow.mul(&inv_x);
        }

        // [x^n]_2
        let x_pow = inv_x_pow.inverse();

        let xn2 = G2_GENERATOR.mul(&x_pow);

        // [s^n - x^n]_2
        let xn_minus_yn = self.secret_g2[n].sub(&xn2);

        // [interpolation_polynomial(s)]_1
        let is1 = self.commit_to_poly(&interp).unwrap();

        // [commitment - interpolation_polynomial(s)]_1 = [commit]_1 - [interpolation_polynomial(s)]_1
        let commit_minus_interp = com.sub(&is1);
        let ret = pairings_verify(&commit_minus_interp, &G2_GENERATOR, proof, &xn_minus_yn);

        Ok(ret)
    }

    fn get_expanded_roots_of_unity_at(&self, i: usize) -> ZFr {
        self.fs.get_expanded_roots_of_unity_at(i)
    }

    fn get_roots_of_unity_at(&self, i: usize) -> ZFr {
        self.fs.get_roots_of_unity_at(i)
    }

    fn get_fft_settings(&self) -> &ZFFTSettings {
        &self.fs
    }

    fn get_g1_secret(&self) -> &[ZG1] {
        &self.secret_g1
    }

    fn get_g2_secret(&self) -> &[ZG2] {
        &self.secret_g2
    }
}
