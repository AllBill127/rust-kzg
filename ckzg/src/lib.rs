pub mod consts;
pub mod eip_4844;
pub mod fftsettings;
pub mod fftsettings4844;
pub mod finite;
pub mod fk20settings;
pub mod kzgsettings;
pub mod kzgsettings4844;
pub mod poly;
pub mod utils;

#[cfg(feature = "parallel")]
const RUN_PARALLEL: bool = true;
#[cfg(not(feature = "parallel"))]
const RUN_PARALLEL: bool = false;
