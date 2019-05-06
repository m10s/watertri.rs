//! An implementation of the ray-triangle intersection algorithm described in:
//!
//! > Sven Woop, Carsten Benthin, and Ingo Wald. "Watertight ray/triangle intersection."
//! > Journal of Computer Graphics Techniques (JCGT) 2.1 (2013): 65-82.
//!
//! Does not perform backface culling.
#![allow(non_snake_case)]
// Variable names from the paper (appendix A) are not snake_case
extern crate cgmath;

// use cgmath::num_traits::Signed;
use cgmath::{BaseFloat, Vector3};

/// Precomputed data depending only on the ray.
#[derive(Clone, Debug)]
pub struct RayData<S> {
    kx: usize,
    ky: usize,
    kz: usize,
    Sx: S,
    Sy: S,
    Sz: S,
    org: Vector3<S>,
}

impl<S> RayData<S>
where
    S: BaseFloat,
{
    /// Pre-compute the transformation that is applied to all triangles.
    pub fn new(org: Vector3<S>, dir: Vector3<S>) -> RayData<S> {
        // The paper swaps kx and ky if dir[kz] is negative, to preserve winding order.
        // But winding order is only relevant for backface culling, which we don't perform.
        let kz = max_dim(dir);
        let kx = (kz + 1) % 3;
        let ky = (kz + 2) % 3;

        // S::div(1.0, dir[kz]);

        RayData {
            kx: kx,
            ky: ky,
            kz: kz,
            Sx: dir[kx] / dir[kz],
            Sy: dir[ky] / dir[kz],
            Sz: S::one() / dir[kz],
            org: org,
        }
    }

    /// Perform the intersection calculation.
    pub fn intersect(
        &self,
        A: Vector3<S>,
        B: Vector3<S>,
        C: Vector3<S>,
    ) -> Option<Intersection<S>> {
        let (Sx, Sy, Sz, org) = (self.Sx, self.Sy, self.Sz, self.org);
        let (kx, ky, kz) = (self.kx, self.ky, self.kz);
        let (A, B, C) = (A - org, B - org, C - org);
        let Ax = A[kx] - Sx * A[kz];
        let Ay = A[ky] - Sy * A[kz];
        let Bx = B[kx] - Sx * B[kz];
        let By = B[ky] - Sy * B[kz];
        let Cx = C[kx] - Sx * C[kz];
        let Cy = C[ky] - Sy * C[kz];

        let mut U = Cx * By - Cy * Bx;
        let mut V = Ax * Cy - Ay * Cx;
        let mut W = Bx * Ay - By * Ax;

        if U == S::zero() || V == S::zero() || W == S::zero() {
            let CxBy = Cx * By;
            let CyBx = Cy * Bx;
            U = CxBy - CyBx;
            let AxCy = Ax * Cy;
            let AyCx = Ay * Cx;
            V = AxCy - AyCx;
            let BxAy = Bx * Ay;
            let ByAx = By * Ax;
            W = BxAy - ByAx;
        }

        if (U < S::zero() || V < S::zero() || W < S::zero())
            && (U > S::zero() || V > S::zero() || W > S::zero())
        {
            return None;
        }

        let det = U + V + W;
        if det == S::zero() {
            return None;
        }

        let Az = Sz * A[kz];
        let Bz = Sz * B[kz];
        let Cz = Sz * C[kz];
        let T = U * Az + V * Bz + W * Cz;

        let rcpDet = S::one() / det;
        Some(Intersection {
            t: T * rcpDet,
            u: U * rcpDet,
            v: V * rcpDet,
            w: W * rcpDet,
        })
    }
}

/// Geometric information about a ray-triangle intersection.
pub struct Intersection<S> {
    /// Parametric distance from the ray origin to the intersection.
    pub t: S,
    /// Barycentric coordinate.
    pub u: S,
    /// Barycentric coordinate.
    pub v: S,
    /// Barycentric coordinate.
    pub w: S,
}

fn max_dim<S>(v: Vector3<S>) -> usize
where
    S: BaseFloat,
{
    let (x, y, z) = (v.x.abs(), v.y.abs(), v.z.abs());
    if x > y {
        // y isn't the maximum, so it's either x or z
        if x > z {
            0
        } else {
            2
        }
    } else if y > z {
        // x isn't the maximum, so it's either y or z
        1
    } else {
        2
    }
}
