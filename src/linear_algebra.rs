pub fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    /*
    i j k
    x y z
    a b c

    (yc - bz)i - (xc - az)j + (xb - ay)k
    x = a[0]
    y = a[1]
    z = a[2]
    a = b[0]
    b = b[1]
    c = b[2]
    */

    [
        (a[1] * b[2] - b[1] * a[2]),
        -1.0 * (a[0] * b[2] - b[0] * a[2]),
        (a[0] * b[1] - a[1] * b[0]),
    ]
}

pub fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn normalize(a: &[f64; 3]) -> [f64; 3] {
    let mag = (a[0].powi(2) + a[1].powi(2) + a[2].powi(2)).sqrt();
    [a[0] / mag, a[1] / mag, a[2] / mag]
}

pub fn multiply_4_by_3(vector: &[f64; 3], matrix: &[[f64; 4]; 4]) -> [f64; 3] {
    [
        matrix[0][0] * vector[0]
            + matrix[0][1] * vector[1]
            + matrix[0][2] * vector[2]
            + matrix[0][3] * 1.0,
        matrix[1][0] * vector[0]
            + matrix[1][1] * vector[1]
            + matrix[1][2] * vector[2]
            + matrix[1][3] * 1.0,
        matrix[2][0] * vector[0]
            + matrix[2][1] * vector[1]
            + matrix[2][2] * vector[2]
            + matrix[2][3] * 1.0,
        // matrix[3][0] * vector[0]
        //     + matrix[3][1] * vector[1]
        //     + matrix[3][2] * vector[2]
        //     + matrix[3][3] * 1.0,
    ]
}
