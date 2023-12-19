use super::P1;
use crate::P2;
// use bls12_381::{Fp as ZFp, Fp2 as ZFp2, G1Projective, G2Projective, Scalar};
use pairing_ce::bls12_381::{G1 as G1Projective, G2 as G2Projective, Fr as Scalar, FqRepr};
use pairing_ce::bls12_381::{Fq as ZFp, Fq2 as ZFp2};
use pairing_ce::bls12_381::FrRepr;
use pairing_ce::ff::{Field, PrimeField};
use blst::{blst_fp, blst_fp2, blst_fr, blst_p1, blst_p2};
use ff::derive::{adc, mac};
use pairing_ce::CurveProjective;
use crate::consts::{INV, MODULUS};

#[derive(Debug, PartialEq, Eq)]
pub struct Error;

pub const fn blst_fr_into_pc_fr(fr: blst_fr) -> Scalar {
    // Scalar(fr.l)
    let tmp = FrRepr(fr.l);
    Scalar::from_repr(tmp).unwrap()
}
pub const fn pc_fr_into_blst_fr(scalar: Scalar) -> blst_fr {
    // blst_fr { l: scalar.0 }
    let tmp = scalar.into_repr();
    blst_fr { l: tmp.0 }
}
pub const fn blst_fp2_into_pc_fq2(fp: &blst_fp2) -> ZFp2 {
    // let c0 = ZFp(fp.fp[0].l);
    // let c1 = ZFp(fp.fp[1].l);
    let tmp0 = FqRepr(fp.fp[0].l);
    let tmp1 = FqRepr(fp.fp[1].l);
    let c0 = ZFp::from_repr(tmp0).unwrap();
    let c1 = ZFp::from_repr(tmp1).unwrap();
    ZFp2 { c0, c1 }
}

pub const fn blst_p1_into_pc_g1projective(p1: &P1) -> G1Projective {
    // let x = ZFp(p1.x.l);
    // let y = ZFp(p1.y.l);
    // let z = ZFp(p1.z.l);
    // G1Projective { x, y, z }
    let x = ZFp::from_repr(FqRepr(p1.x.l)).unwrap();
    let y = ZFp::from_repr(FqRepr(p1.y.l)).unwrap();
    let z = ZFp::from_repr(FqRepr(p1.z.l)).unwrap();
    G1Projective::from_xyz_checked(x, y, z).unwrap()
}

pub const fn pc_g1projective_into_blst_p1(p1: G1Projective) -> blst_p1 {
    // let x = blst_fp { l: p1.x.0 };
    // let y = blst_fp { l: p1.y.0 };
    // let z = blst_fp { l: p1.z.0 };

    let (tmp_x, tmp_y, tmp_z) = p1.as_xyz();
    let x = blst_fp { l: tmp_x.into_repr().0 };
    let y = blst_fp { l: tmp_y.into_repr().0 };
    let z = blst_fp { l: tmp_z.into_repr().0 };

    blst_p1 { x, y, z }
}

pub const fn blst_p2_into_pc_g2projective(p2: &P2) -> G2Projective {
    // G2Projective {
    //     x: blst_fp2_into_pc_fq2(&p2.x),
    //     y: blst_fp2_into_pc_fq2(&p2.y),
    //     z: blst_fp2_into_pc_fq2(&p2.z),
    // }
    let x = blst_fp2_into_pc_fq2(&p2.x);
    let y = blst_fp2_into_pc_fq2(&p2.y);
    let z = blst_fp2_into_pc_fq2(&p2.z);
    G2Projective::from_xyz_checked(x, y, z).unwrap()
}

pub const fn pc_g2projective_into_blst_p2(p2: G2Projective) -> blst_p2 {
    // let x = blst_fp2 {
    //     fp: [blst_fp { l: p2.x.c0.0 }, blst_fp { l: p2.x.c1.0 }],
    // };
    //
    // let y = blst_fp2 {
    //     fp: [blst_fp { l: p2.y.c0.0 }, blst_fp { l: p2.y.c1.0 }],
    // };
    //
    // let z = blst_fp2 {
    //     fp: [blst_fp { l: p2.z.c0.0 }, blst_fp { l: p2.z.c1.0 }],
    // };

    let (tmp_x, tmp_y, tmp_z) = p2.as_xyz();

    let x = blst_fp2 {
        fp: [blst_fp { l: tmp_x.c0.into_repr().0 }, blst_fp { l: tmp_x.c1.into_repr().0 }],
    };

    let y = blst_fp2 {
        fp: [blst_fp { l: tmp_y.c0.into_repr().0 }, blst_fp { l: tmp_y.c1.into_repr().0 }],
    };

    let z = blst_fp2 {
        fp: [blst_fp { l: tmp_z.c0.into_repr().0 }, blst_fp { l: tmp_z.c1.into_repr().0 }],
    };

    blst_p2 { x, y, z }
}

#[inline(always)]
pub const fn montgomery_reduce(
    r0: u64,
    r1: u64,
    r2: u64,
    r3: u64,
    r4: u64,
    r5: u64,
    r6: u64,
    r7: u64,
) -> Scalar {
    // The Montgomery reduction here is based on Algorithm 14.32 in
    // Handbook of Applied Cryptography
    // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

    let k = r0.wrapping_mul(INV);
    let (_, carry) = mac(r0, k, MODULUS.into_repr().0[0], 0);
    let (r1, carry) = mac(r1, k, MODULUS.into_repr().0[1], carry);
    let (r2, carry) = mac(r2, k, MODULUS.into_repr().0[2], carry);
    let (r3, carry) = mac(r3, k, MODULUS.into_repr().0[3], carry);
    let (r4, carry2) = adc(r4, 0, carry);

    let k = r1.wrapping_mul(INV);
    let (_, carry) = mac(r1, k, MODULUS.into_repr().0[0], 0);
    let (r2, carry) = mac(r2, k, MODULUS.into_repr().0[1], carry);
    let (r3, carry) = mac(r3, k, MODULUS.into_repr().0[2], carry);
    let (r4, carry) = mac(r4, k, MODULUS.into_repr().0[3], carry);
    let (r5, carry2) = adc(r5, carry2, carry);

    let k = r2.wrapping_mul(INV);
    let (_, carry) = mac(r2, k, MODULUS.into_repr().0[0], 0);
    let (r3, carry) = mac(r3, k, MODULUS.into_repr().0[1], carry);
    let (r4, carry) = mac(r4, k, MODULUS.into_repr().0[2], carry);
    let (r5, carry) = mac(r5, k, MODULUS.into_repr().0[3], carry);
    let (r6, carry2) = adc(r6, carry2, carry);

    let k = r3.wrapping_mul(INV);
    let (_, carry) = mac(r3, k, MODULUS.into_repr().0[0], 0);
    let (r4, carry) = mac(r4, k, MODULUS.into_repr().0[1], carry);
    let (r5, carry) = mac(r5, k, MODULUS.into_repr().0[2], carry);
    let (r6, carry) = mac(r6, k, MODULUS.into_repr().0[3], carry);
    let (r7, _) = adc(r7, carry2, carry);

    // Result may be within MODULUS of the correct value
    // (&Scalar([r4, r5, r6, r7])).sub(&MODULUS)
    let tmp = FrRepr([r4, r5, r6, r7]);
    let mut res = Scalar::from_repr(tmp).unwrap();
    res.sub_assign(&MODULUS);
    res
}
