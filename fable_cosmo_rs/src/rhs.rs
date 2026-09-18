//! The two-equation right-hand side handed to CVODE.
//!
//! The independent variable is the number of e-folds `N = ln a`, and the state
//! is `y = [ln S, s]`:
//!
//! ```text
//!   d(ln S)/dN = -nu(s) = -(3 - xi (1 - h(s)))
//!   ds/dN      = gamma s (1 - s)
//! ```
//!
//! Neither equation depends on `ln S` itself, so only the bridge field is read.

use std::any::Any;

use cvode_rs::prelude::*;

use crate::model::{bridge_flow, dilution_exponent, Params, IDX_BRIDGE, IDX_LOG_S};

/// The `CVRhsFn` given to `CVodeInit`. Returns `CV_SUCCESS` on success and -1
/// when the parameter block is missing or a vector pointer cannot be taken.
pub fn rhs(
    _n: f64,
    y: &N_Vector,
    ydot: &N_Vector,
    user_data: &mut Option<Box<dyn Any>>,
) -> i32 {
    // The parameters are read, never written, so a shared downcast suffices.
    let Some(p) = user_data.as_ref().and_then(|b| b.downcast_ref::<Params>()) else {
        return -1;
    };
    // Scope the borrow of `y` so it is released before `ydot` is taken.
    let s = {
        let Some(d) = N_VGetArrayPointer(y) else {
            return -1;
        };
        d[IDX_BRIDGE]
    };
    let Some(mut out) = N_VGetArrayPointer(ydot) else {
        return -1;
    };
    out[IDX_LOG_S] = -dilution_exponent(s, p.xi);
    out[IDX_BRIDGE] = bridge_flow(s, p.gamma_flow);
    CV_SUCCESS
}
