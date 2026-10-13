//! Vector maths and the Kabsch superposition, following Bio.SVDSuperimposer (row vectors: x' = x·rot + tran).

pub type V3 = [f64; 3];
pub type M3 = [[f64; 3]; 3];

pub fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
pub fn length(a: V3) -> f64 {
    dot(a, a).sqrt()
}
pub fn dist(a: V3, b: V3) -> f64 {
    length(sub(a, b))
}

/// round to f32 precision, where the Python code stores float32 arrays
pub fn f32r(a: V3) -> V3 {
    [a[0] as f32 as f64, a[1] as f32 as f64, a[2] as f32 as f64]
}

/// row vector times matrix: v·m
pub fn vm(v: V3, m: &M3) -> V3 {
    [
        v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
        v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
        v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2],
    ]
}

pub fn mm(a: &M3, b: &M3) -> M3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    r
}

pub fn transpose(a: &M3) -> M3 {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[j][i];
        }
    }
    r
}

pub fn det(a: &M3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// angle between vectors in degrees (distances.vector_angle)
pub fn angle(a: V3, b: V3) -> f64 {
    let c = (dot(a, b) / (dot(a, a) * dot(b, b)).sqrt()).clamp(-1.0, 1.0);
    c.acos() * 180.0 / std::f64::consts::PI
}

/// Superimpose `coords` onto `reference` (Bio.SVDSuperimposer): returns (rot, tran, rms)
pub fn kabsch(reference: &[V3], coords: &[V3]) -> (M3, V3, f64) {
    let n = coords.len() as f64;
    let mut av1 = [0.0; 3];
    let mut av2 = [0.0; 3];
    for (c, r) in coords.iter().zip(reference) {
        av1 = add(av1, *c);
        av2 = add(av2, *r);
    }
    av1 = scale(av1, 1.0 / n);
    av2 = scale(av2, 1.0 / n);
    // correlation matrix a = coords^T · reference
    let mut a = nalgebra::Matrix3::<f64>::zeros();
    for (c, r) in coords.iter().zip(reference) {
        let c = sub(*c, av1);
        let r = sub(*r, av2);
        for i in 0..3 {
            for j in 0..3 {
                a[(i, j)] += c[i] * r[j];
            }
        }
    }
    let svd = a.svd(true, true);
    let u = svd.u.unwrap();
    let mut vt = svd.v_t.unwrap();
    let to_m3 = |m: &nalgebra::Matrix3<f64>| -> M3 {
        [[m[(0, 0)], m[(0, 1)], m[(0, 2)]], [m[(1, 0)], m[(1, 1)], m[(1, 2)]], [m[(2, 0)], m[(2, 1)], m[(2, 2)]]]
    };
    // rot = (vt^T u^T)^T = u · vt
    let mut rot = to_m3(&(u * vt));
    if det(&rot) < 0.0 {
        // numpy sorts singular values descending and flips vt[2], i.e. the row of the smallest one
        let k = (0..3)
            .min_by(|&i, &j| svd.singular_values[i].partial_cmp(&svd.singular_values[j]).unwrap())
            .unwrap();
        for j in 0..3 {
            vt[(k, j)] = -vt[(k, j)];
        }
        rot = to_m3(&(u * vt));
    }
    let tran = sub(av2, vm(av1, &rot));
    let mut s = 0.0;
    for (c, r) in coords.iter().zip(reference) {
        let d = sub(add(vm(*c, &rot), tran), *r);
        s += dot(d, d);
    }
    (rot, tran, (s / n).sqrt())
}
