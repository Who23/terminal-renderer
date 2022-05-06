use crate::{color::Color, linear_algebra::multiply_4_by_3};

#[derive(Debug)]
pub struct Triangle {
    pub points: [[f64; 3]; 3],
    pub color: Color,
}

impl Triangle {
    pub fn new(point1: [f64; 3], point2: [f64; 3], point3: [f64; 3], color: Color) -> Self {
        Self {
            points: [point1, point2, point3],
            color,
        }
    }
    pub fn camera_coords(&self, w_to_c: &[[f64; 4]; 4]) -> [[f64; 3]; 3] {
        let mut camera_coords = [[0.0; 3]; 3];

        for (i, point) in (&self.points).iter().enumerate() {
            camera_coords[i] = multiply_4_by_3(&point, &w_to_c);
            camera_coords[i][2] *= -1.0;
        }

        camera_coords
    }
    pub fn project(&self, camera_coords: &[[f64; 3]; 3], k1: f64) -> [[f64; 2]; 3] {
        let mut projected = [[0.0; 2]; 3];
        for (i, point) in (&camera_coords).iter().enumerate() {
            projected[i][0] = (k1 * point[0]) / (point[2]);
            projected[i][1] = (k1 * point[1]) / (point[2]);
        }

        projected
    }
    pub fn apply_transformation(&mut self, matrix: &[[f64; 3]; 3]) {
        for point in &mut self.points {
            *point = [
                matrix[0][0] * point[0] + matrix[0][1] * point[1] + matrix[0][2] * point[2],
                matrix[1][0] * point[0] + matrix[1][1] * point[1] + matrix[1][2] * point[2],
                matrix[2][0] * point[0] + matrix[2][1] * point[1] + matrix[2][2] * point[2],
            ]
        }
    }
}
