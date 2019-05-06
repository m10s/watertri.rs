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
    sx: S,
    sy: S,
    sz: S,
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
            sx: dir[kx] / dir[kz],
            sy: dir[ky] / dir[kz],
            sz: S::one() / dir[kz],
            org: org,
        }
    }

    /// Perform the intersection calculation.
    pub fn intersect(
        &self,
        a: Vector3<S>,
        b: Vector3<S>,
        c: Vector3<S>,
    ) -> Option<Intersection<S>> {
        let (sx, sy, sz, org) = (self.sx, self.sy, self.sz, self.org);
        let (kx, ky, kz) = (self.kx, self.ky, self.kz);
        let (a, b, c) = (a - org, b - org, c - org);
        let ax = a[kx] - sx * a[kz];
        let ay = a[ky] - sy * a[kz];
        let bx = b[kx] - sx * b[kz];
        let by = b[ky] - sy * b[kz];
        let cx = c[kx] - sx * c[kz];
        let cy = c[ky] - sy * c[kz];

        let mut u = cx * by - cy * bx;
        let mut v = ax * cy - ay * cx;
        let mut w = bx * ay - by * ax;

        if u == S::zero() || v == S::zero() || w == S::zero() {
            let cxby = cx * by;
            let cybx = cy * bx;
            u = cxby - cybx;
            let axcy = ax * cy;
            let aycx = ay * cx;
            v = axcy - aycx;
            let bxay = bx * ay;
            let byax = by * ax;
            w = bxay - byax;
        }

        if (u < S::zero() || v < S::zero() || w < S::zero())
            && (u > S::zero() || v > S::zero() || w > S::zero())
        {
            return None;
        }

        let det = u + v + w;
        if det == S::zero() {
            return None;
        }

        let az = sz * a[kz];
        let bz = sz * b[kz];
        let cz = sz * c[kz];
        let t = u * az + v * bz + w * cz;

        let rcp_det = S::one() / det;
        Some(Intersection {
            t: t * rcp_det,
            u: u * rcp_det,
            v: v * rcp_det,
            w: w * rcp_det,
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
