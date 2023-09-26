use criterion::{criterion_group, criterion_main, Criterion};
use kzg_bench::benches::eip_4844::bench_eip_4844;
use rust_kzg_zkcrypto::eip_4844::{
    blob_to_kzk_commitment_zkc, bytes_to_blob, compute_blob_kzg_proof_zkc, compute_kzg_proof_zkc,
    load_trusted_setup_file_zkc, verify_blob_kzg_proof_zkc, verify_blob_kzg_proof_batch_zkc, verify_kzg_proof_zkc,
};
use rust_kzg_zkcrypto::fftsettings::ZkFFTSettings;
use rust_kzg_zkcrypto::kzg_proofs::KZGSettings;
use rust_kzg_zkcrypto::kzg_types::ZkG2Projective;
use rust_kzg_zkcrypto::poly::KzgPoly;
use rust_kzg_zkcrypto::utils::ZkG1Projective;
use rust_kzg_zkcrypto::zkfr::blsScalar;

fn bench_eip_4844_(c: &mut Criterion) {
    bench_eip_4844::<blsScalar, ZkG1Projective, ZkG2Projective, KzgPoly, ZkFFTSettings, KZGSettings>(
        c,
        &load_trusted_setup_file_zkc,
        &blob_to_kzk_commitment_zkc,
        &bytes_to_blob,
        &compute_kzg_proof_zkc,
        &verify_kzg_proof_zkc,
        &compute_blob_kzg_proof_zkc,
        &verify_blob_kzg_proof_zkc,
        &verify_blob_kzg_proof_batch_zkc,
    );
}

criterion_group!(benches, bench_eip_4844_,);
criterion_main!(benches);
