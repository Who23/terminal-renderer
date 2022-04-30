// https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
// https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/rasterization-stage
// https://www.a1k0n.net/2011/07/20/donut-math.html
// https://www.cse.psu.edu/~rtc12/CSE486/lecture12.pdf

use std::{
    f64::consts::{PI, SQRT_2},
    io::{stdout, BufWriter, StdoutLock, Write},
};
use termion::{self, terminal_size};

const THETA: f64 = (2.0 * PI) / 480.0;
const BRIGHTNESS_CHARS: [u8; 9] = *b".:-=+*#%@";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Color {
    Red,
    Green,
    Blue,
    Magenta,
    Cyan,
    Yellow,
    None,
}

#[derive(Debug)]
struct Triangle {
    points: [[f64; 3]; 3],
    color: Color,
}

impl Triangle {
    fn new(point1: [f64; 3], point2: [f64; 3], point3: [f64; 3], color: Color) -> Self {
        Self {
            points: [point1, point2, point3],
            color,
        }
    }
    fn camera_coords(&self, w_to_c: &[[f64; 4]; 4]) -> [[f64; 3]; 3] {
        let mut camera_coords = [[0.0; 3]; 3];

        for (i, point) in (&self.points).iter().enumerate() {
            camera_coords[i] = multiply_4_by_3(&point, &w_to_c);
            camera_coords[i][2] *= -1.0;
        }

        camera_coords
    }
    fn project(&self, camera_coords: &[[f64; 3]; 3], k1: f64) -> [[f64; 2]; 3] {
        let mut projected = [[0.0; 2]; 3];
        for (i, point) in (&camera_coords).iter().enumerate() {
            projected[i][0] = (k1 * point[0]) / (point[2]);
            projected[i][1] = (k1 * point[1]) / (point[2]);
        }

        projected
    }
    fn apply_transformation(&mut self, matrix: &[[f64; 3]; 3]) {
        for point in &mut self.points {
            *point = [
                matrix[0][0] * point[0] + matrix[0][1] * point[1] + matrix[0][2] * point[2],
                matrix[1][0] * point[0] + matrix[1][1] * point[1] + matrix[1][2] * point[2],
                matrix[2][0] * point[0] + matrix[2][1] * point[1] + matrix[2][2] * point[2],
            ]
        }
    }
}

fn main() {
    let r_x: [[f64; 3]; 3] = [
        [1.0, 0.0, 0.0],
        [0.0, THETA.cos(), -1.0 * THETA.sin()],
        [0.0, THETA.sin(), THETA.cos()],
    ];

    let r_y: [[f64; 3]; 3] = [
        [THETA.cos(), 0.0, THETA.sin()],
        [0.0, 1.0, 0.0],
        [-1.0 * THETA.sin(), 0.0, THETA.cos()],
    ];

    let r_z: [[f64; 3]; 3] = [
        [THETA.cos(), -1.0 * THETA.sin(), 0.0],
        [THETA.sin(), THETA.cos(), 0.0],
        [0.0, 0.0, 1.0],
    ];

    let k1: f64 = 1.5;
    let k2: f64 = 10.0;
    let (w, h) = {
        let size = terminal_size().unwrap();
        (size.0.into(), size.1.into())
    };
    let scaling = std::cmp::max(w, h);

    let stdout = stdout();
    let mut stdout = BufWriter::new(stdout.lock());
    let mut buffer: Vec<Vec<u8>>; // vec![vec![32u8; w]; h];
    let mut zbuffer: Vec<Vec<f64>>; // vec![vec![32u8; w]; h];
    let mut color_buffer: Vec<Vec<Color>>;

    // cube
    let mut triangles = vec![
        // top face (+z)
        // Triangle::new(
        //     [-1.0, 1.0, 1.0],
        //     [1.0, 1.0, 1.0],
        //     [-1.0, -1.0, 1.0],
        //     Color::Cyan,
        // ),
        Triangle::new(
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            Color::Cyan,
        ),
        // bottom face (-z)
        Triangle::new(
            [-1.0, 1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, -1.0, -1.0],
            Color::Red,
        ),
        Triangle::new(
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, -1.0, -1.0],
            Color::Red,
        ),
        // far face (+y)
        Triangle::new(
            [-1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, -1.0],
            Color::Green,
        ),
        Triangle::new(
            [1.0, 1.0, -1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, -1.0],
            Color::Green,
        ),
        // front face (-y)
        Triangle::new(
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [-1.0, -1.0, -1.0],
            Color::Blue,
        ),
        Triangle::new(
            [1.0, -1.0, -1.0],
            [1.0, -1.0, 1.0],
            [-1.0, -1.0, -1.0],
            Color::Blue,
        ),
        // right face (+x)
        Triangle::new(
            [1.0, 1.0, -1.0],
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            Color::Magenta,
        ),
        Triangle::new(
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            Color::Magenta,
        ),
        // left face (-x)
        Triangle::new(
            [-1.0, 1.0, -1.0],
            [-1.0, 1.0, 1.0],
            [-1.0, -1.0, -1.0],
            Color::Yellow,
        ),
        Triangle::new(
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, 1.0],
            [-1.0, -1.0, -1.0],
            Color::Yellow,
        ),
    ];

    // camera aligned with xyz axes (not rotated) at (0, 0, 10)
    let world_to_camera = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, -10.0],
        [0.0, 0.0, 0.0, 1.0],
    ];

    // light pointing in +z
    let light_direction = [0.0, 0.0, 1.0];
    let light_direction = normalize(&multiply_4_by_3(&light_direction, &world_to_camera));

    let w_f = w as f64;
    let h_f = h as f64;
    let scaling_f = scaling as f64;
    loop {
        buffer = vec![vec![b' '; w]; h];
        zbuffer = vec![vec![f64::MAX; w]; h];
        color_buffer = vec![vec![Color::None; w]; h];
        for triangle in &mut triangles {
            let camera = triangle.camera_coords(&world_to_camera);
            let mut projected = triangle.project(&camera, k1);

            // scale projected points to screen size
            for point in &mut projected {
                point[0] = point[0] * scaling_f + w_f / 2.0;
                // characters are about half as wide as they are tall, so we
                // scale the image down vertically by half to compensate
                point[1] = point[1] * scaling_f / 2.0 + h_f / 2.0;
            }

            // double the area of the triangle, and the denominator for the
            // barycentric coordinates formula
            let triangle_area_2 = edge(&projected[0], &projected[1], &projected[2]);

            // loop over pixels and check if they are contained in the triangle
            for buf_x in 0..w {
                for buf_y in 0..h {
                    let p = [buf_x as f64, buf_y as f64];
                    let E_01 = edge(&projected[0], &projected[1], &p);
                    let E_12 = edge(&projected[1], &projected[2], &p);
                    let E_20 = edge(&projected[2], &projected[0], &p);
                    if contains_point(&[E_01, E_12, E_20]) {
                        // barymetric coordinate coefficients
                        let l0 = E_12 / triangle_area_2;
                        let l1 = E_20 / triangle_area_2;
                        let l2 = E_01 / triangle_area_2;
                        // 1/P.z = 1/v0.z * l0 + 1/v1.z * l1 + 1/v2.z * l2
                        let inverse_z: f64 =
                            l0 / camera[0][2] + l1 / camera[1][2] + l2 / camera[2][2];

                        let z = 1.0 / inverse_z;

                        // if this point is closer to the camera than the previous points
                        if z < zbuffer[buf_y][buf_x] {
                            // calculate brightness / luminance of pixel
                            // dot product is equal to ||a|| * ||b|| * cos T
                            // making magnitude of both vectors -> result is cosT
                            // range is -1 to 1.
                            // tranform this range from [0, 8] to choose which character
                            // to draw to the screen.
                            let lumininance = dot(
                                &normalize(&cross(
                                    // vector from v0 to v1
                                    &[
                                        camera[1][0] - camera[0][0],
                                        camera[1][1] - camera[0][1],
                                        camera[1][2] - camera[0][2],
                                    ],
                                    // vector from v0 to v2
                                    &[
                                        camera[2][0] - camera[0][0],
                                        camera[2][1] - camera[0][1],
                                        camera[2][2] - camera[0][2],
                                    ],
                                )),
                                &light_direction,
                            );

                            let luminance = lumininance.abs() * (BRIGHTNESS_CHARS.len() - 1) as f64;
                            zbuffer[buf_y][buf_x] = z;
                            buffer[buf_y][buf_x] = BRIGHTNESS_CHARS[luminance as usize];
                            color_buffer[buf_y][buf_x] = triangle.color;
                        }
                    }
                }
            }

            // if 0 < y && y < h && 0 < x && x < w {
            //     buffer[y][x] = b'@';
            // }

            triangle.apply_transformation(&r_x);
            triangle.apply_transformation(&r_x);
            triangle.apply_transformation(&r_x);

            triangle.apply_transformation(&r_y);
            triangle.apply_transformation(&r_y);
            triangle.apply_transformation(&r_y);
            triangle.apply_transformation(&r_y);
            triangle.apply_transformation(&r_y);

            triangle.apply_transformation(&r_z);
            triangle.apply_transformation(&r_z);
        }

        write_to_screen(&buffer, &color_buffer, &mut stdout);
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
}

fn edge(a: &[f64; 2], b: &[f64; 2], c: &[f64; 2]) -> f64 {
    (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0])
}

fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
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

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalize(a: &[f64; 3]) -> [f64; 3] {
    let mag = (a[0].powi(2) + a[1].powi(2) + a[2].powi(2)).sqrt();
    [a[0] / mag, a[1] / mag, a[2] / mag]
}

fn multiply_4_by_3(vector: &[f64; 3], matrix: &[[f64; 4]; 4]) -> [f64; 3] {
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

// 0 1 p
// (p[0] - t[0][0]) * (t[1][1] - t[0][1]) - (p[1] - t[0][1]) * (t[1][0] - t[0][0])
// (p[0] - t[0][0]) * (t[1][1] - t[0][1]) - (p[1] - t[0][1]) * (t[1][0] - t[0][0])
// fn original_edge_all(p: &[f64; 2], t: &[[f64; 2]; 3]) -> bool {
//     let E_01 = (p[0] - t[0][0]) * (t[1][1] - t[0][1]) - (p[1] - t[0][1]) * (t[1][0] - t[0][0]);
//     let E_12 = (p[0] - t[1][0]) * (t[2][1] - t[1][1]) - (p[1] - t[1][1]) * (t[2][0] - t[1][0]);
//     let E_20 = (p[0] - t[2][0]) * (t[0][1] - t[2][1]) - (p[1] - t[2][1]) * (t[0][0] - t[2][0]);

//     if (E_01 >= 0.0 && E_12 >= 0.0 && E_20 >= 0.0) || (E_01 <= 0.0 && E_12 <= 0.0 && E_20 <= 0.0) {
//         return true;
//     } else {
//         return false;
//     }
// }

// edge: the output of the edge function for the point p
fn contains_point(edge: &[f64; 3]) -> bool {
    // we have to check if they are all the same sign - this accounts for the
    // normal facing either towards or away from the camera depending on the
    // winding. We don't need that info anyway unless we're doing backface
    // culling, which we are not. Outside the shape, there cannot be a point
    // where all three edge calculations result in the same sign. i don't think.
    if (edge[0] >= 0.0 && edge[1] >= 0.0 && edge[2] >= 0.0)
        || (edge[0] <= 0.0 && edge[1] <= 0.0 && edge[2] <= 0.0)
    {
        return true;
    } else {
        return false;
    }
}

fn write_to_screen(
    buffer: &Vec<Vec<u8>>,
    color_buffer: &Vec<Vec<Color>>,
    stdout: &mut BufWriter<StdoutLock>,
) {
    let colors_map = [
        format!("{}", termion::color::Fg(termion::color::Red)),
        format!("{}", termion::color::Fg(termion::color::Green)),
        format!("{}", termion::color::Fg(termion::color::Blue)),
        format!("{}", termion::color::Fg(termion::color::Magenta)),
        format!("{}", termion::color::Fg(termion::color::Cyan)),
        format!("{}", termion::color::Fg(termion::color::Yellow)),
        format!("{}", termion::color::Fg(termion::color::Reset)),
    ]
    .map(|c| c.into_bytes());

    write!(stdout, "{}", termion::clear::All).unwrap();
    for (characters, colors) in buffer.iter().zip(color_buffer.iter()) {
        for (character, color) in characters.iter().zip(colors.iter()) {
            let idx = match color {
                &Color::Red => 0,
                &Color::Green => 1,
                &Color::Blue => 2,
                &Color::Magenta => 3,
                &Color::Cyan => 4,
                &Color::Yellow => 5,
                &Color::None => 6,
            };

            stdout.write_all(&colors_map[idx]).unwrap();
            stdout.write_all(&[*character]).unwrap();
        }
    }
    stdout.flush().unwrap();
}
