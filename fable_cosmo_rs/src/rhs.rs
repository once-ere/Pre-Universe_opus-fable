#![forbid(unsafe_code)]
#![deny(warnings)]
#![allow(non_snake_case)]

//! The two-equation right-hand side handed to CVODE.
//!
//! The independent variable is the number of e-folds `N = ln a`, and the state
//! is `y = [ln S, s]`:
//!
//! ```text
//!   d(ln S)/dN = -nu(s) = -(3 - xi (1 - h(s)))
//!   ds/dN      = gamma s (1 - s)
//! ```

use std::any::Any;

use cvode_rs::prelude::*;

use crate::model::{bridge_flow, dilution_exponent, Params, IDX_BRIDGE, IDX_LOG_S};

/// The `CVRhsFn` given to `CVodeInit`. Returns 0 on success and -1 when the
/// parameter block is missing or a vector pointer cannot be taken.
pub fn rhs(
    _n: f64,
    y: &N_Vector,
    ydot: &N_Vector,
    user_data: &mut Option<Box<dyn Any>>,
) -> i32 {
    let p = match user_data.as_mut().and_then(|b| b.downcast_mut::<Params>()) {
        Some(p) => *p,
        None => return -1,
    };
    let state = {
        let d = match N_VGetArrayPointer(y) {
            Some(d) => d,
            None => return -1,
        };
        [d[IDX_LOG_S], d[IDX_BRIDGE]]
    };
    let mut out = match N_VGetArrayPointer(ydot) {
        Some(d) => d,
        None => return -1,
    };
    out[IDX_LOG_S] = -dilution_exponent(state[IDX_BRIDGE], p.xi);
    out[IDX_BRIDGE] = bridge_flow(state[IDX_BRIDGE], p.gamma_flow);
    0
}
