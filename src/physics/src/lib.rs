use wasm_bindgen::prelude::*;
use serde::{ Serialize, Deserialize };

// Field type constants
pub const U_FIELD: u8 = 0;
pub const V_FIELD: u8 = 1;
pub const S_FIELD: u8 = 2;

// Vec2 for positions and velocities
#[wasm_bindgen]
#[derive(Serialize, Deserialize, Clone, Copy, Default)]
pub struct Vec2 {
    pub x: f32,
    pub y: f32,
}

#[wasm_bindgen]
impl Vec2 {
    #[wasm_bindgen(constructor)]
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }
}

// Fluid particle for rendering
#[wasm_bindgen]
#[derive(Serialize, Deserialize, Clone)]
pub struct FluidParticle {
    pub pos: Vec2,
    pub vel: Vec2,
    pub color: u32,
}

#[wasm_bindgen]
impl FluidParticle {
    #[wasm_bindgen(constructor)]
    pub fn new(x: f32, y: f32, vx: f32, vy: f32, color: u32) -> Self {
        Self {
            pos: Vec2::new(x, y),
            vel: Vec2::new(vx, vy),
            color,
        }
    }
}

// Implementation enum for solver type
#[wasm_bindgen]
#[derive(Serialize, Deserialize, Clone, Copy, PartialEq)]
pub enum Implementation {
    Euler = 0,
    Flip = 1,
}

// Cell type constants for FLIP
const FLUID_CELL: i32 = 0;
const AIR_CELL: i32 = 1;
const SOLID_CELL: i32 = 2;

// FLIP Fluid implementation - particle-based fluid simulation
#[wasm_bindgen]
#[derive(Serialize, Deserialize, Clone)]
pub struct FlipFluid {
    // Grid properties
    density: f32,
    num_x: usize,
    num_y: usize,
    h: f32,
    inv_spacing: f32,
    num_cells: usize,

    // Grid velocities (staggered MAC grid)
    u: Vec<f32>,
    v: Vec<f32>,
    du: Vec<f32>,
    dv: Vec<f32>,
    prev_u: Vec<f32>,
    prev_v: Vec<f32>,
    p: Vec<f32>,
    s: Vec<f32>, // 0 = solid, 1 = fluid
    cell_type: Vec<i32>,

    // Particles
    max_particles: usize,
    num_particles: usize,
    particle_pos: Vec<f32>, // x, y pairs
    particle_vel: Vec<f32>, // vx, vy pairs
    particle_color: Vec<f32>, // r, g, b triples
    particle_density: Vec<f32>,
    particle_rest_density: f32,

    // Particle grid for neighbor search
    particle_radius: f32,
    p_inv_spacing: f32,
    p_num_x: usize,
    p_num_y: usize,
    p_num_cells: usize,
    num_cell_particles: Vec<i32>,
    first_cell_particle: Vec<i32>,
    cell_particle_ids: Vec<i32>,
}

impl FlipFluid {
    pub fn new(
        density: f32,
        width: f32,
        height: f32,
        spacing: f32,
        particle_radius: f32,
        max_particles: usize
    ) -> FlipFluid {
        let num_x = ((width / spacing).floor() as usize) + 1;
        let num_y = ((height / spacing).floor() as usize) + 1;
        let h = (width / (num_x as f32)).max(height / (num_y as f32));
        let inv_spacing = 1.0 / h;
        let num_cells = num_x * num_y;

        // Particle grid
        let p_inv_spacing = 1.0 / (2.2 * particle_radius);
        let p_num_x = ((width * p_inv_spacing).floor() as usize) + 1;
        let p_num_y = ((height * p_inv_spacing).floor() as usize) + 1;
        let p_num_cells = p_num_x * p_num_y;

        // Initialize particle colors to blue
        let mut particle_color = vec![0.0; 3 * max_particles];
        for i in 0..max_particles {
            particle_color[3 * i + 2] = 1.0; // Blue
        }

        FlipFluid {
            density,
            num_x,
            num_y,
            h,
            inv_spacing,
            num_cells,
            u: vec![0.0; num_cells],
            v: vec![0.0; num_cells],
            du: vec![0.0; num_cells],
            dv: vec![0.0; num_cells],
            prev_u: vec![0.0; num_cells],
            prev_v: vec![0.0; num_cells],
            p: vec![0.0; num_cells],
            s: vec![0.0; num_cells],
            cell_type: vec![AIR_CELL; num_cells],
            max_particles,
            num_particles: 0,
            particle_pos: vec![0.0; 2 * max_particles],
            particle_vel: vec![0.0; 2 * max_particles],
            particle_color,
            particle_density: vec![0.0; num_cells],
            particle_rest_density: 0.0,
            particle_radius,
            p_inv_spacing,
            p_num_x,
            p_num_y,
            p_num_cells,
            num_cell_particles: vec![0; p_num_cells],
            first_cell_particle: vec![0; p_num_cells + 1],
            cell_particle_ids: vec![0; max_particles],
        }
    }

    fn clamp(x: f32, min: f32, max: f32) -> f32 {
        x.max(min).min(max)
    }

    fn clamp_usize(x: i32, min: usize, max: usize) -> usize {
        (x.max(min as i32) as usize).min(max)
    }

    pub fn integrate_particles(&mut self, dt: f32, gravity: f32) {
        for i in 0..self.num_particles {
            self.particle_vel[2 * i + 1] += dt * gravity;
            self.particle_pos[2 * i] += self.particle_vel[2 * i] * dt;
            self.particle_pos[2 * i + 1] += self.particle_vel[2 * i + 1] * dt;
        }
    }

    pub fn push_particles_apart(&mut self, num_iters: usize) {
        let color_diffusion_coeff = 0.001;

        // Count particles per cell
        self.num_cell_particles.fill(0);

        for i in 0..self.num_particles {
            let x = self.particle_pos[2 * i];
            let y = self.particle_pos[2 * i + 1];

            let xi = Self::clamp_usize(
                (x * self.p_inv_spacing).floor() as i32,
                0,
                self.p_num_x - 1
            );
            let yi = Self::clamp_usize(
                (y * self.p_inv_spacing).floor() as i32,
                0,
                self.p_num_y - 1
            );
            let cell_nr = xi * self.p_num_y + yi;
            self.num_cell_particles[cell_nr] += 1;
        }

        // Partial sums
        let mut first = 0i32;
        for i in 0..self.p_num_cells {
            first += self.num_cell_particles[i];
            self.first_cell_particle[i] = first;
        }
        self.first_cell_particle[self.p_num_cells] = first;

        // Fill particles into cells
        for i in 0..self.num_particles {
            let x = self.particle_pos[2 * i];
            let y = self.particle_pos[2 * i + 1];

            let xi = Self::clamp_usize(
                (x * self.p_inv_spacing).floor() as i32,
                0,
                self.p_num_x - 1
            );
            let yi = Self::clamp_usize(
                (y * self.p_inv_spacing).floor() as i32,
                0,
                self.p_num_y - 1
            );
            let cell_nr = xi * self.p_num_y + yi;
            self.first_cell_particle[cell_nr] -= 1;
            self.cell_particle_ids[self.first_cell_particle[cell_nr] as usize] = i as i32;
        }

        // Push particles apart
        let min_dist = 2.0 * self.particle_radius;
        let min_dist2 = min_dist * min_dist;

        for _ in 0..num_iters {
            for i in 0..self.num_particles {
                let px = self.particle_pos[2 * i];
                let py = self.particle_pos[2 * i + 1];

                let pxi = (px * self.p_inv_spacing).floor() as i32;
                let pyi = (py * self.p_inv_spacing).floor() as i32;
                let x0 = (pxi - 1).max(0) as usize;
                let y0 = (pyi - 1).max(0) as usize;
                let x1 = ((pxi + 1) as usize).min(self.p_num_x - 1);
                let y1 = ((pyi + 1) as usize).min(self.p_num_y - 1);

                for xi in x0..=x1 {
                    for yi in y0..=y1 {
                        let cell_nr = xi * self.p_num_y + yi;
                        let first = self.first_cell_particle[cell_nr] as usize;
                        let last = self.first_cell_particle[cell_nr + 1] as usize;

                        for j in first..last {
                            let id = self.cell_particle_ids[j] as usize;
                            if id == i {
                                continue;
                            }

                            let qx = self.particle_pos[2 * id];
                            let qy = self.particle_pos[2 * id + 1];

                            let dx = qx - px;
                            let dy = qy - py;
                            let d2 = dx * dx + dy * dy;

                            if d2 > min_dist2 || d2 == 0.0 {
                                continue;
                            }

                            let d = d2.sqrt();
                            let s = (0.5 * (min_dist - d)) / d;
                            let dx = dx * s;
                            let dy = dy * s;

                            self.particle_pos[2 * i] -= dx;
                            self.particle_pos[2 * i + 1] -= dy;
                            self.particle_pos[2 * id] += dx;
                            self.particle_pos[2 * id + 1] += dy;

                            // Diffuse colors
                            for k in 0..3 {
                                let color0 = self.particle_color[3 * i + k];
                                let color1 = self.particle_color[3 * id + k];
                                let color = (color0 + color1) * 0.5;
                                self.particle_color[3 * i + k] =
                                    color0 + (color - color0) * color_diffusion_coeff;
                                self.particle_color[3 * id + k] =
                                    color1 + (color - color1) * color_diffusion_coeff;
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn handle_particle_collisions(
        &mut self,
        obstacle_x: f32,
        obstacle_y: f32,
        obstacle_radius: f32,
        obstacle_vel_x: f32,
        obstacle_vel_y: f32
    ) {
        let h = self.h;
        let r = self.particle_radius;
        let or = obstacle_radius;
        let min_dist = or + r;
        let n = self.num_y;

        // Boundary limits (keep particles away from edges based on wall thickness)
        let wall_thickness = 2.0; // Match the basin walls
        let min_x = wall_thickness * h + r;
        let max_x = ((self.num_x as f32) - wall_thickness) * h - r;
        let min_y = wall_thickness * h + r;
        let max_y = ((self.num_y as f32) - wall_thickness) * h - r;

        for i in 0..self.num_particles {
            let mut x = self.particle_pos[2 * i];
            let mut y = self.particle_pos[2 * i + 1];
            let mut vx = self.particle_vel[2 * i];
            let mut vy = self.particle_vel[2 * i + 1];

            // Obstacle collision - push particles out of the ball
            if or > 0.0 {
                let dx = x - obstacle_x;
                let dy = y - obstacle_y;
                let d = (dx * dx + dy * dy).sqrt();

                if d < min_dist && d > 0.001 {
                    // Normalize direction from obstacle center to particle
                    let nx = dx / d;
                    let ny = dy / d;

                    // Push particle out to the surface of the obstacle
                    x = obstacle_x + nx * min_dist;
                    y = obstacle_y + ny * min_dist;

                    // Reflect velocity off the obstacle surface
                    let vn = vx * nx + vy * ny;
                    if vn < 0.0 {
                        // Remove velocity component going into obstacle
                        vx -= vn * nx;
                        vy -= vn * ny;
                        // Add obstacle velocity
                        vx += obstacle_vel_x;
                        vy += obstacle_vel_y;
                    }
                }
            }

            // Wall collisions - clamp position and zero velocity
            if x < min_x {
                x = min_x;
                vx = 0.0;
            }
            if x > max_x {
                x = max_x;
                vx = 0.0;
            }
            if y < min_y {
                y = min_y;
                vy = 0.0;
            }
            if y > max_y {
                y = max_y;
                vy = 0.0;
            }

            // Solid cell collisions (user-placed walls)
            // Check the cell the particle is in and neighboring cells
            let xi = (x * self.inv_spacing).floor() as i32;
            let yi = (y * self.inv_spacing).floor() as i32;

            for dxi in -1..=1_i32 {
                for dyi in -1..=1_i32 {
                    let cx = xi + dxi;
                    let cy = yi + dyi;

                    if cx >= 0 && cx < (self.num_x as i32) && cy >= 0 && cy < (self.num_y as i32) {
                        let cell_idx = (cx as usize) * n + (cy as usize);

                        // If this cell is solid (s == 0.0)
                        if self.s[cell_idx] == 0.0 {
                            // Cell boundaries
                            let cell_min_x = (cx as f32) * h;
                            let cell_max_x = ((cx as f32) + 1.0) * h;
                            let cell_min_y = (cy as f32) * h;
                            let cell_max_y = ((cy as f32) + 1.0) * h;

                            // Check if particle overlaps with this solid cell
                            // Find closest point on cell to particle
                            let closest_x = x.max(cell_min_x).min(cell_max_x);
                            let closest_y = y.max(cell_min_y).min(cell_max_y);

                            let pdx = x - closest_x;
                            let pdy = y - closest_y;
                            let dist_sq = pdx * pdx + pdy * pdy;

                            if dist_sq < r * r {
                                if dist_sq > 0.001 {
                                    // Particle is outside but overlapping
                                    let dist = dist_sq.sqrt();
                                    let overlap = r - dist;
                                    let nx = pdx / dist;
                                    let ny = pdy / dist;

                                    // Push particle out
                                    x += nx * overlap * 1.01;
                                    y += ny * overlap * 1.01;

                                    // Zero velocity component into wall
                                    let vn = vx * nx + vy * ny;
                                    if vn < 0.0 {
                                        vx -= vn * nx;
                                        vy -= vn * ny;
                                    }
                                } else {
                                    // Particle center is inside the cell - push to nearest edge
                                    let cell_cx = ((cx as f32) + 0.5) * h;
                                    let cell_cy = ((cy as f32) + 0.5) * h;
                                    let to_center_x = x - cell_cx;
                                    let to_center_y = y - cell_cy;

                                    // Push out along the axis with greater distance to center
                                    if to_center_x.abs() > to_center_y.abs() {
                                        if to_center_x > 0.0 {
                                            x = cell_max_x + r;
                                        } else {
                                            x = cell_min_x - r;
                                        }
                                        vx = 0.0;
                                    } else {
                                        if to_center_y > 0.0 {
                                            y = cell_max_y + r;
                                        } else {
                                            y = cell_min_y - r;
                                        }
                                        vy = 0.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            self.particle_pos[2 * i] = x;
            self.particle_pos[2 * i + 1] = y;
            self.particle_vel[2 * i] = vx;
            self.particle_vel[2 * i + 1] = vy;
        }
    }

    pub fn update_particle_density(&mut self) {
        let n = self.num_y;
        let h = self.h;
        let h1 = self.inv_spacing;
        let h2 = 0.5 * h;

        self.particle_density.fill(0.0);

        for i in 0..self.num_particles {
            let mut x = self.particle_pos[2 * i];
            let mut y = self.particle_pos[2 * i + 1];

            x = Self::clamp(x, h, ((self.num_x - 1) as f32) * h);
            y = Self::clamp(y, h, ((self.num_y - 1) as f32) * h);

            let x0 = ((x - h2) * h1).floor() as usize;
            let tx = (x - h2 - (x0 as f32) * h) * h1;
            let x1 = (x0 + 1).min(self.num_x - 2);

            let y0 = ((y - h2) * h1).floor() as usize;
            let ty = (y - h2 - (y0 as f32) * h) * h1;
            let y1 = (y0 + 1).min(self.num_y - 2);

            let sx = 1.0 - tx;
            let sy = 1.0 - ty;

            if x0 < self.num_x && y0 < self.num_y {
                self.particle_density[x0 * n + y0] += sx * sy;
            }
            if x1 < self.num_x && y0 < self.num_y {
                self.particle_density[x1 * n + y0] += tx * sy;
            }
            if x1 < self.num_x && y1 < self.num_y {
                self.particle_density[x1 * n + y1] += tx * ty;
            }
            if x0 < self.num_x && y1 < self.num_y {
                self.particle_density[x0 * n + y1] += sx * ty;
            }
        }

        // Calculate rest density
        if self.particle_rest_density == 0.0 {
            let mut sum = 0.0;
            let mut num_fluid_cells = 0;

            for i in 0..self.num_cells {
                if self.cell_type[i] == FLUID_CELL {
                    sum += self.particle_density[i];
                    num_fluid_cells += 1;
                }
            }

            if num_fluid_cells > 0 {
                self.particle_rest_density = sum / (num_fluid_cells as f32);
            }
        }
    }

    pub fn transfer_velocities(&mut self, to_grid: bool, flip_ratio: f32) {
        let n = self.num_y;
        let h = self.h;
        let h1 = self.inv_spacing;
        let h2 = 0.5 * h;

        if to_grid {
            self.prev_u.copy_from_slice(&self.u);
            self.prev_v.copy_from_slice(&self.v);

            self.du.fill(0.0);
            self.dv.fill(0.0);
            self.u.fill(0.0);
            self.v.fill(0.0);

            // Set cell types
            for i in 0..self.num_cells {
                self.cell_type[i] = if self.s[i] == 0.0 { SOLID_CELL } else { AIR_CELL };
            }

            for i in 0..self.num_particles {
                let x = self.particle_pos[2 * i];
                let y = self.particle_pos[2 * i + 1];
                let xi = Self::clamp_usize((x * h1).floor() as i32, 0, self.num_x - 1);
                let yi = Self::clamp_usize((y * h1).floor() as i32, 0, self.num_y - 1);
                let cell_nr = xi * n + yi;
                if self.cell_type[cell_nr] == AIR_CELL {
                    self.cell_type[cell_nr] = FLUID_CELL;
                }
            }
        }

        for component in 0..2 {
            let dx = if component == 0 { 0.0 } else { h2 };
            let dy = if component == 0 { h2 } else { 0.0 };

            let (f, prev_f, d) = if component == 0 {
                (&mut self.u, &self.prev_u, &mut self.du)
            } else {
                (&mut self.v, &self.prev_v, &mut self.dv)
            };

            for i in 0..self.num_particles {
                let mut x = self.particle_pos[2 * i];
                let mut y = self.particle_pos[2 * i + 1];

                x = Self::clamp(x, h, ((self.num_x - 1) as f32) * h);
                y = Self::clamp(y, h, ((self.num_y - 1) as f32) * h);

                let x0 = ((x - dx) * h1).floor().min((self.num_x - 2) as f32) as usize;
                let tx = (x - dx - (x0 as f32) * h) * h1;
                let x1 = (x0 + 1).min(self.num_x - 2);

                let y0 = ((y - dy) * h1).floor().min((self.num_y - 2) as f32) as usize;
                let ty = (y - dy - (y0 as f32) * h) * h1;
                let y1 = (y0 + 1).min(self.num_y - 2);

                let sx = 1.0 - tx;
                let sy = 1.0 - ty;

                let d0 = sx * sy;
                let d1 = tx * sy;
                let d2 = tx * ty;
                let d3 = sx * ty;

                let nr0 = x0 * n + y0;
                let nr1 = x1 * n + y0;
                let nr2 = x1 * n + y1;
                let nr3 = x0 * n + y1;

                if to_grid {
                    let pv = self.particle_vel[2 * i + component];
                    f[nr0] += pv * d0;
                    d[nr0] += d0;
                    f[nr1] += pv * d1;
                    d[nr1] += d1;
                    f[nr2] += pv * d2;
                    d[nr2] += d2;
                    f[nr3] += pv * d3;
                    d[nr3] += d3;
                } else {
                    let offset = if component == 0 { n } else { 1 };
                    let valid0 = if
                        self.cell_type[nr0] != AIR_CELL ||
                        (nr0 >= offset && self.cell_type[nr0 - offset] != AIR_CELL)
                    {
                        1.0
                    } else {
                        0.0
                    };
                    let valid1 = if
                        self.cell_type[nr1] != AIR_CELL ||
                        (nr1 >= offset && self.cell_type[nr1 - offset] != AIR_CELL)
                    {
                        1.0
                    } else {
                        0.0
                    };
                    let valid2 = if
                        self.cell_type[nr2] != AIR_CELL ||
                        (nr2 >= offset && self.cell_type[nr2 - offset] != AIR_CELL)
                    {
                        1.0
                    } else {
                        0.0
                    };
                    let valid3 = if
                        self.cell_type[nr3] != AIR_CELL ||
                        (nr3 >= offset && self.cell_type[nr3 - offset] != AIR_CELL)
                    {
                        1.0
                    } else {
                        0.0
                    };

                    let d_sum = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                    if d_sum > 0.0 {
                        let pic_v =
                            (valid0 * d0 * f[nr0] +
                                valid1 * d1 * f[nr1] +
                                valid2 * d2 * f[nr2] +
                                valid3 * d3 * f[nr3]) /
                            d_sum;
                        let corr =
                            (valid0 * d0 * (f[nr0] - prev_f[nr0]) +
                                valid1 * d1 * (f[nr1] - prev_f[nr1]) +
                                valid2 * d2 * (f[nr2] - prev_f[nr2]) +
                                valid3 * d3 * (f[nr3] - prev_f[nr3])) /
                            d_sum;
                        let flip_v = self.particle_vel[2 * i + component] + corr;

                        self.particle_vel[2 * i + component] =
                            (1.0 - flip_ratio) * pic_v + flip_ratio * flip_v;
                    }
                }
            }

            if to_grid {
                // Normalize
                for idx in 0..f.len() {
                    if d[idx] > 0.0 {
                        f[idx] /= d[idx];
                    }
                }

                // Restore solid cells
                for i in 0..self.num_x {
                    for j in 0..self.num_y {
                        let solid = self.cell_type[i * n + j] == SOLID_CELL;
                        if component == 0 {
                            if solid || (i > 0 && self.cell_type[(i - 1) * n + j] == SOLID_CELL) {
                                self.u[i * n + j] = self.prev_u[i * n + j];
                            }
                        } else {
                            if solid || (j > 0 && self.cell_type[i * n + j - 1] == SOLID_CELL) {
                                self.v[i * n + j] = self.prev_v[i * n + j];
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn solve_incompressibility(
        &mut self,
        num_iters: usize,
        dt: f32,
        over_relaxation: f32,
        compensate_drift: bool
    ) {
        self.p.fill(0.0);
        self.prev_u.copy_from_slice(&self.u);
        self.prev_v.copy_from_slice(&self.v);

        let n = self.num_y;
        let cp = (self.density * self.h) / dt;

        for _ in 0..num_iters {
            for i in 1..self.num_x - 1 {
                for j in 1..self.num_y - 1 {
                    if self.cell_type[i * n + j] != FLUID_CELL {
                        continue;
                    }

                    let center = i * n + j;
                    let left = (i - 1) * n + j;
                    let right = (i + 1) * n + j;
                    let bottom = i * n + j - 1;
                    let top = i * n + j + 1;

                    let sx0 = self.s[left];
                    let sx1 = self.s[right];
                    let sy0 = self.s[bottom];
                    let sy1 = self.s[top];
                    let s = sx0 + sx1 + sy0 + sy1;

                    if s == 0.0 {
                        continue;
                    }

                    // Divergence: d = u[i+1,j] - u[i,j] + v[i,j+1] - v[i,j]
                    let mut div = self.u[right] - self.u[center] + self.v[top] - self.v[center];

                    // Apply overrelaxation
                    div *= over_relaxation;

                    if self.particle_rest_density > 0.0 && compensate_drift {
                        // Reduce divergence in dense regions to push particles apart
                        let k = 1.0;
                        let compression =
                            self.particle_density[center] - self.particle_rest_density;
                        if compression > 0.0 {
                            div -= k * compression;
                        }
                    }

                    // Pressure correction (for tracking only)
                    let p_val = div / s;
                    self.p[center] += cp * p_val;

                    // Velocity corrections using s values from neighbors
                    // u[i,j] += d * s[i-1,j] / s
                    // u[i+1,j] -= d * s[i+1,j] / s
                    // v[i,j] += d * s[i,j-1] / s
                    // v[i,j+1] -= d * s[i,j+1] / s
                    self.u[center] += sx0 * p_val;
                    self.u[right] -= sx1 * p_val;
                    self.v[center] += sy0 * p_val;
                    self.v[top] -= sy1 * p_val;
                }
            }
        }
    }

    pub fn update_particle_colors(&mut self) {
        let h1 = self.inv_spacing;

        for i in 0..self.num_particles {
            // Fade towards default blue color
            let s = 0.01;

            self.particle_color[3 * i] = Self::clamp(self.particle_color[3 * i] - s, 0.0, 1.0);
            self.particle_color[3 * i + 1] = Self::clamp(
                self.particle_color[3 * i + 1] - s,
                0.0,
                1.0
            );
            self.particle_color[3 * i + 2] = Self::clamp(
                self.particle_color[3 * i + 2] + s,
                0.0,
                1.0
            );

            let x = self.particle_pos[2 * i];
            let y = self.particle_pos[2 * i + 1];
            let xi = Self::clamp_usize((x * h1).floor() as i32, 1, self.num_x - 1);
            let yi = Self::clamp_usize((y * h1).floor() as i32, 1, self.num_y - 1);
            let cell_nr = xi * self.num_y + yi;

            let d0 = self.particle_rest_density;

            if d0 > 0.0 && cell_nr < self.particle_density.len() {
                let rel_density = self.particle_density[cell_nr] / d0;
                
                // Color based on density:
                // Low density (sparse/negative space) = green
                // High density (dense pockets) = red
                // Normal density = stays blue
                if rel_density < 0.7 {
                    // Low density - green
                    self.particle_color[3 * i] = 0.2;     // Low red
                    self.particle_color[3 * i + 1] = 0.9; // High green
                    self.particle_color[3 * i + 2] = 0.3; // Low blue
                } else if rel_density > 1.3 {
                    // High density - red
                    self.particle_color[3 * i] = 0.9;     // High red
                    self.particle_color[3 * i + 1] = 0.2; // Low green
                    self.particle_color[3 * i + 2] = 0.2; // Low blue
                }
                // Normal density (0.7-1.3) keeps fading to blue
            }
        }
    }

    pub fn simulate(
        &mut self,
        dt: f32,
        gravity: f32,
        flip_ratio: f32,
        num_pressure_iters: usize,
        num_particle_iters: usize,
        over_relaxation: f32,
        compensate_drift: bool,
        separate_particles: bool,
        obstacle_x: f32,
        obstacle_y: f32,
        obstacle_radius: f32,
        obstacle_vel_x: f32,
        obstacle_vel_y: f32
    ) {
        let num_sub_steps = 1;
        let sdt = dt / (num_sub_steps as f32);

        for _ in 0..num_sub_steps {
            self.integrate_particles(sdt, gravity);
            if separate_particles {
                self.push_particles_apart(num_particle_iters);
            }
            self.handle_particle_collisions(
                obstacle_x,
                obstacle_y,
                obstacle_radius,
                obstacle_vel_x,
                obstacle_vel_y
            );
            self.transfer_velocities(true, flip_ratio);
            self.update_particle_density();
            self.solve_incompressibility(
                num_pressure_iters,
                sdt,
                over_relaxation,
                compensate_drift
            );
            self.transfer_velocities(false, flip_ratio);
        }

        self.update_particle_colors();
    }
}

#[wasm_bindgen]
#[derive(Serialize, Deserialize, Clone)]
pub struct Fluid {
    density: f32,
    num_x: usize,
    num_y: usize,
    num_cells: usize,
    h: f32,
    u: Vec<f32>,
    v: Vec<f32>,
    new_u: Vec<f32>,
    new_v: Vec<f32>,
    p: Vec<f32>,
    s: Vec<f32>,
    m: Vec<f32>,
    new_m: Vec<f32>,
    // Age tracking for particles (frame number when smoke was added, 0 = no smoke)
    particle_age: Vec<u32>,
}

impl Fluid {
    pub fn new(density: f32, num_x: usize, num_y: usize, h: f32) -> Fluid {
        let actual_num_x = num_x + 2;
        let actual_num_y = num_y + 2;
        let num_cells = actual_num_x * actual_num_y;

        let mut m = vec![0.0_f32; num_cells];
        m.fill(1.0);

        Fluid {
            density,
            num_x: actual_num_x,
            num_y: actual_num_y,
            num_cells,
            h,
            u: vec![0.0; num_cells],
            v: vec![0.0; num_cells],
            new_u: vec![0.0; num_cells],
            new_v: vec![0.0; num_cells],
            p: vec![0.0; num_cells],
            s: vec![0.0; num_cells],
            m,
            new_m: vec![0.0; num_cells],
            particle_age: vec![0; num_cells],
        }
    }

    pub fn integrate(&mut self, dt: f32, gravity: f32) {
        let n = self.num_y;
        for i in 1..self.num_x {
            for j in 1..self.num_y - 1 {
                if self.s[i * n + j] != 0.0 && self.s[i * n + j - 1] != 0.0 {
                    self.v[i * n + j] += gravity * dt;
                }
            }
        }
    }

    pub fn solve_incompressibility(&mut self, num_iters: usize, dt: f32, over_relaxation: f32) {
        let n = self.num_y;
        let cp = (self.density * self.h) / dt;

        for _ in 0..num_iters {
            for i in 1..self.num_x - 1 {
                for j in 1..self.num_y - 1 {
                    if self.s[i * n + j] == 0.0 {
                        continue;
                    }

                    let sx0 = self.s[(i - 1) * n + j];
                    let sx1 = self.s[(i + 1) * n + j];
                    let sy0 = self.s[i * n + j - 1];
                    let sy1 = self.s[i * n + j + 1];
                    let s = sx0 + sx1 + sy0 + sy1;
                    if s == 0.0 {
                        continue;
                    }

                    let div =
                        self.u[(i + 1) * n + j] -
                        self.u[i * n + j] +
                        self.v[i * n + j + 1] -
                        self.v[i * n + j];

                    let mut p = -div / s;
                    p *= over_relaxation;
                    self.p[i * n + j] += cp * p;

                    self.u[i * n + j] -= sx0 * p;
                    self.u[(i + 1) * n + j] += sx1 * p;
                    self.v[i * n + j] -= sy0 * p;
                    self.v[i * n + j + 1] += sy1 * p;
                }
            }
        }
    }

    pub fn extrapolate(&mut self) {
        let n = self.num_y;
        for i in 0..self.num_x {
            self.u[i * n + 0] = self.u[i * n + 1];
            self.u[i * n + self.num_y - 1] = self.u[i * n + self.num_y - 2];
        }
        for j in 0..self.num_y {
            self.v[0 * n + j] = self.v[1 * n + j];
            self.v[(self.num_x - 1) * n + j] = self.v[(self.num_x - 2) * n + j];
        }
    }

    pub fn sample_field(&self, x: f32, y: f32, field: u8) -> f32 {
        let n = self.num_y;
        let h = self.h;
        let h1 = 1.0 / h;
        let h2 = 0.5 * h;

        let x = x.max(h).min((self.num_x as f32) * h);
        let y = y.max(h).min((self.num_y as f32) * h);

        let (dx, dy, f): (f32, f32, &Vec<f32>) = match field {
            U_FIELD => (0.0, h2, &self.u),
            V_FIELD => (h2, 0.0, &self.v),
            S_FIELD => (h2, h2, &self.m),
            _ => (h2, h2, &self.m),
        };

        let x0 = (((x - dx) * h1).floor() as usize).min(self.num_x - 1);
        let tx = (x - dx - (x0 as f32) * h) * h1;
        let x1 = (x0 + 1).min(self.num_x - 1);

        let y0 = (((y - dy) * h1).floor() as usize).min(self.num_y - 1);
        let ty = (y - dy - (y0 as f32) * h) * h1;
        let y1 = (y0 + 1).min(self.num_y - 1);

        let sx = 1.0 - tx;
        let sy = 1.0 - ty;

        sx * sy * f[x0 * n + y0] +
            tx * sy * f[x1 * n + y0] +
            tx * ty * f[x1 * n + y1] +
            sx * ty * f[x0 * n + y1]
    }

    pub fn avg_u(&self, i: usize, j: usize) -> f32 {
        let n = self.num_y;
        (self.u[i * n + j - 1] +
            self.u[i * n + j] +
            self.u[(i + 1) * n + j - 1] +
            self.u[(i + 1) * n + j]) *
            0.25
    }

    pub fn avg_v(&self, i: usize, j: usize) -> f32 {
        let n = self.num_y;
        (self.v[(i - 1) * n + j] +
            self.v[i * n + j] +
            self.v[(i - 1) * n + j + 1] +
            self.v[i * n + j + 1]) *
            0.25
    }

    pub fn advect_vel(&mut self, dt: f32) {
        self.new_u.copy_from_slice(&self.u);
        self.new_v.copy_from_slice(&self.v);

        let n = self.num_y;
        let h = self.h;
        let h2 = 0.5 * h;

        for i in 1..self.num_x {
            for j in 1..self.num_y {
                // u component
                if self.s[i * n + j] != 0.0 && self.s[(i - 1) * n + j] != 0.0 && j < self.num_y - 1 {
                    let mut x = (i as f32) * h;
                    let mut y = (j as f32) * h + h2;
                    let u = self.u[i * n + j];
                    let v = self.avg_v(i, j);
                    x = x - dt * u;
                    y = y - dt * v;
                    let u = self.sample_field(x, y, U_FIELD);
                    self.new_u[i * n + j] = u;
                }
                // v component
                if self.s[i * n + j] != 0.0 && self.s[i * n + j - 1] != 0.0 && i < self.num_x - 1 {
                    let mut x = (i as f32) * h + h2;
                    let mut y = (j as f32) * h;
                    let u = self.avg_u(i, j);
                    let v = self.v[i * n + j];
                    x = x - dt * u;
                    y = y - dt * v;
                    let v = self.sample_field(x, y, V_FIELD);
                    self.new_v[i * n + j] = v;
                }
            }
        }

        self.u.copy_from_slice(&self.new_u);
        self.v.copy_from_slice(&self.new_v);
    }

    pub fn advect_smoke(&mut self, dt: f32) {
        self.new_m.copy_from_slice(&self.m);
        let mut new_particle_age = self.particle_age.clone();

        let n = self.num_y;
        let h = self.h;
        let h2 = 0.5 * h;

        for i in 1..self.num_x - 1 {
            for j in 1..self.num_y - 1 {
                if self.s[i * n + j] != 0.0 {
                    let u = (self.u[i * n + j] + self.u[(i + 1) * n + j]) * 0.5;
                    let v = (self.v[i * n + j] + self.v[i * n + j + 1]) * 0.5;
                    let x = (i as f32) * h + h2 - dt * u;
                    let y = (j as f32) * h + h2 - dt * v;

                    self.new_m[i * n + j] = self.sample_field(x, y, S_FIELD);

                    // Also advect particle age - get the age from the source cell
                    // Find the source cell indices
                    let src_i = ((x / h).floor() as usize).clamp(1, self.num_x - 2);
                    let src_j = ((y / h).floor() as usize).clamp(1, self.num_y - 2);
                    let src_idx = src_i * n + src_j;

                    // If the source cell has smoke with a valid age, transfer it
                    // Otherwise keep the current age
                    if self.m[src_idx] < 0.99 && self.particle_age[src_idx] > 0 {
                        new_particle_age[i * n + j] = self.particle_age[src_idx];
                    }
                }
            }
        }
        self.m.copy_from_slice(&self.new_m);
        self.particle_age = new_particle_age;
    }

    // Enforce boundary conditions - zero velocities at solid boundaries
    pub fn enforce_boundaries(&mut self) {
        let n = self.num_y;

        for i in 0..self.num_x {
            for j in 0..self.num_y {
                if self.s[i * n + j] == 0.0 {
                    // This cell is solid - zero all velocities touching it
                    self.u[i * n + j] = 0.0;
                    self.v[i * n + j] = 0.0;
                    if i + 1 < self.num_x {
                        self.u[(i + 1) * n + j] = 0.0;
                    }
                    if j + 1 < self.num_y {
                        self.v[i * n + j + 1] = 0.0;
                    }
                }
            }
        }
    }

    pub fn simulate(&mut self, dt: f32, gravity: f32, num_iters: usize, over_relaxation: f32) {
        // Enforce boundaries before simulation
        self.enforce_boundaries();

        self.integrate(dt, gravity);

        self.p.fill(0.0);
        self.solve_incompressibility(num_iters, dt, over_relaxation);

        self.extrapolate();
        self.advect_vel(dt);
        self.advect_smoke(dt);

        // Enforce boundaries after simulation
        self.enforce_boundaries();
    }

    pub fn get_u(&self, i: usize, j: usize) -> f32 {
        let n = self.num_y;
        self.u[i * n + j]
    }

    pub fn get_v(&self, i: usize, j: usize) -> f32 {
        let n = self.num_y;
        self.v[i * n + j]
    }

    pub fn get_p(&self, i: usize, j: usize) -> f32 {
        let n = self.num_y;
        self.p[i * n + j]
    }

    pub fn get_s(&self, i: usize, j: usize) -> f32 {
        let n = self.num_y;
        self.s[i * n + j]
    }

    pub fn get_m(&self, i: usize, j: usize) -> f32 {
        let n = self.num_y;
        self.m[i * n + j]
    }

    /// Remove the oldest particles (smoke cells) from the simulation.
    /// This is called when the particle count exceeds the maximum.
    pub fn remove_oldest_particles(&mut self, count: usize) {
        let n = self.num_y;

        // Collect all cells that have smoke (m < 0.99) with their ages
        // Only include cells that have a valid age (> 0)
        let mut smoke_cells: Vec<(usize, usize, u32)> = Vec::new();

        for i in 1..self.num_x - 1 {
            for j in 1..self.num_y - 1 {
                if self.m[i * n + j] < 0.99 && self.s[i * n + j] != 0.0 {
                    let age = self.particle_age[i * n + j];
                    // Only include particles with valid age (age > 0)
                    // age = 0 means uninitialized, treat as newest (don't remove)
                    if age > 0 {
                        smoke_cells.push((i, j, age));
                    }
                }
            }
        }

        if smoke_cells.is_empty() {
            return;
        }

        // Sort by age (oldest first - lowest frame number means oldest)
        // Particles created earlier have lower frame numbers
        smoke_cells.sort_by_key(|&(_, _, age)| age);

        // Remove the oldest particles (those with lowest frame numbers)
        let to_remove = count.min(smoke_cells.len());
        for (i, j, _) in smoke_cells.iter().take(to_remove) {
            let idx = i * n + j;
            self.m[idx] = 1.0; // Clear smoke
            self.particle_age[idx] = 0; // Reset age
        }
    }

    /// Remove particles that have reached the border of the grid
    /// border_size is the thickness of the border in cells
    pub fn remove_border_particles(&mut self, border_size: usize) {
        let n = self.num_y;
        let bs = border_size.max(1); // At least 1 cell border

        // Clear smoke at border cells (within border_size cells of edges)
        for i in 0..self.num_x {
            for b in 0..bs {
                if b < self.num_y {
                    // Bottom rows
                    self.m[i * n + b] = 1.0;
                    self.particle_age[i * n + b] = 0;
                    // Top rows
                    if self.num_y - 1 - b < self.num_y {
                        self.m[i * n + (self.num_y - 1 - b)] = 1.0;
                        self.particle_age[i * n + (self.num_y - 1 - b)] = 0;
                    }
                }
            }
        }

        for j in 0..self.num_y {
            for b in 0..bs {
                if b < self.num_x {
                    // Left columns
                    self.m[b * n + j] = 1.0;
                    self.particle_age[b * n + j] = 0;
                    // Right columns
                    if self.num_x - 1 - b < self.num_x {
                        self.m[(self.num_x - 1 - b) * n + j] = 1.0;
                        self.particle_age[(self.num_x - 1 - b) * n + j] = 0;
                    }
                }
            }
        }
    }
}

// Helper function to get scientific color mapping
#[wasm_bindgen]
pub fn get_sci_color(val: f32, min_val: f32, max_val: f32) -> Vec<u8> {
    let val = val.max(min_val).min(max_val - 0.0001);
    let d = max_val - min_val;
    let val = if d == 0.0 { 0.5 } else { (val - min_val) / d };
    let m = 0.25_f32;
    let num = (val / m).floor() as i32;
    let s = (val - (num as f32) * m) / m;

    let (r, g, b) = match num {
        0 => (0.0, s, 1.0),
        1 => (0.0, 1.0, 1.0 - s),
        2 => (s, 1.0, 0.0),
        3 => (1.0, 1.0 - s, 0.0),
        _ => (1.0, 0.0, 0.0),
    };

    vec![(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8, 255]
}

// Main Universe struct - all simulation logic integrated here
#[wasm_bindgen]
#[derive(Serialize, Deserialize, Clone)]
pub struct Universe {
    // Fluid grid (Euler implementation)
    fluid: Option<Fluid>,

    // FLIP fluid (particle-based implementation)
    flip_fluid: Option<FlipFluid>,

    // Simulation parameters
    gravity: f32,
    dt: f32,
    num_iters: usize,
    frame_nr: u32,
    over_relaxation: f32,

    // Wind
    wind_x: f32,
    wind_y: f32,

    // Emitter
    emitter_enabled: bool,
    emitter_x: f32,
    emitter_y: f32,
    emitter_radius: f32,
    max_particles: usize,

    // Obstacle
    obstacle_x: f32,
    obstacle_y: f32,
    obstacle_radius: f32,

    // State
    paused: bool,
    scene_nr: u8,

    // Visualization flags
    show_obstacle: bool,
    show_streamlines: bool,
    show_velocities: bool,
    show_pressure: bool,
    show_smoke: bool,
    show_density: bool,

    // Rendering parameters
    particle_radius: f32,
    cell_size: f32,
    border_size: usize, // Border thickness in cells
    grid_size: usize, // Grid resolution (number of cells per side)

    // Speed control
    speed: f32,

    // Implementation type
    implementation: Implementation,

    // FLIP-specific parameters
    flip_ratio: f32, // 0 = PIC, 1 = FLIP (default 0.9)
    compensate_drift: bool, // Drift compensation in pressure solver
    separate_particles: bool, // Push particles apart to avoid clumping
}

#[wasm_bindgen]
impl Universe {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Universe {
        let mut universe = Universe {
            fluid: None,
            flip_fluid: None,
            gravity: 0.0, // Wind tunnel default (no gravity)
            dt: 1.0 / 60.0, // 60 FPS
            num_iters: 40, // Solver iterations
            frame_nr: 0,
            over_relaxation: 1.9, // SOR acceleration
            wind_x: 2.0, // Default wind from left
            wind_y: 0.0,
            emitter_enabled: true,
            emitter_x: 0.5, // Center of grid
            emitter_y: 0.5, // Center of grid
            emitter_radius: 5.0, // Fixed radius in cells
            max_particles: 10000,
            obstacle_x: 0.0,
            obstacle_y: 0.0,
            obstacle_radius: 0.15,
            paused: false,
            scene_nr: 0,
            show_obstacle: false,
            show_streamlines: false,
            show_velocities: false,
            show_pressure: false,
            show_smoke: true,
            show_density: false,
            particle_radius: 4.0,
            cell_size: 10.0,
            border_size: 1, // 1 cell thick border
            grid_size: 100, // Default grid size (resolution)
            speed: 1.0,
            implementation: Implementation::Euler,
            flip_ratio: 0.9, // 90% FLIP, 10% PIC (like flip.html default)
            compensate_drift: false, // Disable drift compensation by default
            separate_particles: true, // Enable particle separation (like flip.html default)
        };
        // Initialize with blank canvas
        universe.setup_blank(16.0, 9.0);
        universe
    }

    pub fn setup_scene(&mut self, scene_nr: u8, sim_width: f32, sim_height: f32) {
        self.scene_nr = scene_nr;
        self.obstacle_radius = 0.15;
        self.over_relaxation = 1.9;
        self.dt = 1.0 / 60.0;
        self.num_iters = 40;

        let res: usize = match scene_nr {
            0 => 50, // Tank
            3 => 200, // Hires Tunnel
            _ => 100, // Wind Tunnel, Paint
        };

        let domain_height = 1.0_f32;
        let domain_width = (domain_height / sim_height) * sim_width;
        let h = domain_height / (res as f32);

        let num_x = (domain_width / h).floor() as usize;
        let num_y = (domain_height / h).floor() as usize;

        let density = 1000.0_f32;

        let mut f = Fluid::new(density, num_x, num_y, h);
        let n = f.num_y;

        match scene_nr {
            0 => {
                // Tank
                for i in 0..f.num_x {
                    for j in 0..f.num_y {
                        let mut s = 1.0_f32;
                        if i == 0 || i == f.num_x - 1 || j == 0 {
                            s = 0.0;
                        }
                        f.s[i * n + j] = s;
                    }
                }
                self.gravity = -9.81;
                self.show_pressure = true;
                self.show_smoke = false;
                self.show_streamlines = false;
                self.show_velocities = false;
            }
            1 | 3 => {
                // Wind Tunnel / Hires Tunnel
                let in_vel = 2.0_f32;
                for i in 0..f.num_x {
                    for j in 0..f.num_y {
                        let mut s = 1.0_f32;
                        if i == 0 || j == 0 || j == f.num_y - 1 {
                            s = 0.0;
                        }
                        f.s[i * n + j] = s;

                        if i == 1 {
                            f.u[i * n + j] = in_vel;
                        }
                    }
                }

                let pipe_h = (0.1 * (f.num_y as f32)) as usize;
                let min_j = (0.5 * (f.num_y as f32) - 0.5 * (pipe_h as f32)) as usize;
                let max_j = (0.5 * (f.num_y as f32) + 0.5 * (pipe_h as f32)) as usize;

                for j in min_j..max_j {
                    f.m[j] = 0.0;
                }

                self.fluid = Some(f);
                self.set_obstacle(0.4, 0.5, true);

                self.gravity = 0.0;
                self.show_pressure = false;
                self.show_smoke = true;
                self.show_streamlines = false;
                self.show_velocities = false;

                if scene_nr == 3 {
                    self.dt = 1.0 / 120.0;
                    self.num_iters = 100;
                    self.show_pressure = true;
                }
                return;
            }
            2 => {
                // Paint
                self.gravity = 0.0;
                self.over_relaxation = 1.0;
                self.show_pressure = false;
                self.show_smoke = true;
                self.show_streamlines = false;
                self.show_velocities = false;
                self.obstacle_radius = 0.1;
            }
            _ => {}
        }

        self.fluid = Some(f);
    }

    // Setup a blank canvas with fixed square size and colored border
    pub fn setup_blank(&mut self, _sim_width: f32, _sim_height: f32) {
        self.scene_nr = 255; // Blank canvas marker
        self.obstacle_radius = 0.15;
        self.over_relaxation = 1.9;
        self.dt = 1.0 / 60.0;
        self.num_iters = 40;

        // Reset FLIP settings to defaults
        self.flip_ratio = 0.9;
        self.compensate_drift = false;
        self.separate_particles = true;

        // Reset emitter to upper area for FLIP (will pour down), center for Euler
        match self.implementation {
            Implementation::Flip => {
                // FLIP: Gravity enabled, no wind, emitter at top
                self.gravity = -9.81;
                self.wind_x = 0.0;
                self.wind_y = 0.0;
                self.emitter_x = 0.5;
                self.emitter_y = 0.8; // Near top so particles fall down
                self.show_obstacle = true;
                self.obstacle_x = 0.5;
                self.obstacle_y = 0.3; // Ball in lower half
            }
            Implementation::Euler => {
                // Euler: No gravity, wind enabled, emitter at left-center
                self.gravity = 0.0;
                self.wind_x = 2.0;
                self.wind_y = 0.0;
                self.emitter_x = 0.15; // Left side
                self.emitter_y = 0.5; // Center height
                self.show_obstacle = true;
                self.obstacle_x = 0.5;
                self.obstacle_y = 0.5; // Ball in center
            }
        }

        // Fixed square grid using grid_size
        let res: usize = self.grid_size;
        let domain_size = 1.0_f32;
        let h = domain_size / (res as f32);

        let num_x = res;
        let num_y = res;

        let density = 1000.0_f32;

        let mut f = Fluid::new(density, num_x, num_y, h);
        let n = f.num_y;

        // All cells are fluid - no walls (border is just visual)
        for i in 0..f.num_x {
            for j in 0..f.num_y {
                f.s[i * n + j] = 1.0; // All fluid
            }
        }

        self.show_pressure = false;
        self.show_smoke = true;
        self.show_streamlines = false;
        self.show_velocities = false;

        self.fluid = Some(f);

        // Also initialize FLIP fluid
        let flip_width = (num_x as f32) * h;
        let flip_height = (num_y as f32) * h;
        let flip_spacing = h;
        let flip_particle_radius = 0.3 * flip_spacing;
        let flip_max_particles = 100000;

        let mut flip = FlipFluid::new(
            density,
            flip_width,
            flip_height,
            flip_spacing,
            flip_particle_radius,
            flip_max_particles
        );

        // Set up basin walls for FLIP - solid walls around edges to contain particles
        let wall_thickness = 2; // 2 cells thick walls
        for i in 0..flip.num_x {
            for j in 0..flip.num_y {
                let idx = i * flip.num_y + j;
                // Create thick walls around all edges (basin)
                if
                    i < wall_thickness ||
                    i >= flip.num_x - wall_thickness ||
                    j < wall_thickness ||
                    j >= flip.num_y - wall_thickness
                {
                    flip.s[idx] = 0.0; // Solid
                } else {
                    flip.s[idx] = 1.0; // Fluid
                }
            }
        }

        self.flip_fluid = Some(flip);
        
        // Initialize obstacle walls in the grid (so they appear immediately)
        if self.show_obstacle {
            let ox = self.obstacle_x;
            let oy = self.obstacle_y;
            self.set_obstacle(ox, oy, true);
        }
    }

    // Get the bounding box of all particles (smoke cells)
    pub fn get_particle_bounds(&self) -> Vec<f32> {
        if let Some(ref f) = self.fluid {
            let n = f.num_y;
            let mut min_x = f32::MAX;
            let mut min_y = f32::MAX;
            let mut max_x = f32::MIN;
            let mut max_y = f32::MIN;
            let mut has_particles = false;

            for i in 1..f.num_x - 1 {
                for j in 1..f.num_y - 1 {
                    if f.m[i * n + j] < 0.99 && f.s[i * n + j] != 0.0 {
                        let x = ((i as f32) + 0.5) * self.cell_size;
                        let y = ((j as f32) + 0.5) * self.cell_size;
                        min_x = min_x.min(x);
                        min_y = min_y.min(y);
                        max_x = max_x.max(x);
                        max_y = max_y.max(y);
                        has_particles = true;
                    }
                }
            }

            if has_particles {
                vec![min_x, min_y, max_x, max_y]
            } else {
                // Return grid center if no particles
                let cx = ((f.num_x as f32) * self.cell_size) / 2.0;
                let cy = ((f.num_y as f32) * self.cell_size) / 2.0;
                vec![cx - 100.0, cy - 100.0, cx + 100.0, cy + 100.0]
            }
        } else {
            vec![0.0, 0.0, 100.0, 100.0]
        }
    }

    // Get total grid dimensions in pixels
    pub fn get_grid_dimensions(&self) -> Vec<f32> {
        if let Some(ref f) = self.fluid {
            vec![(f.num_x as f32) * self.cell_size, (f.num_y as f32) * self.cell_size]
        } else {
            vec![0.0, 0.0]
        }
    }

    pub fn set_obstacle(&mut self, x: f32, y: f32, reset: bool) {
        let mut vx = 0.0_f32;
        let mut vy = 0.0_f32;

        if !reset {
            vx = (x - self.obstacle_x) / self.dt;
            vy = (y - self.obstacle_y) / self.dt;
        }

        self.obstacle_x = x;
        self.obstacle_y = y;
        let r = self.obstacle_radius;

        if let Some(ref mut f) = self.fluid {
            let n = f.num_y;

            for i in 1..f.num_x - 2 {
                for j in 1..f.num_y - 2 {
                    f.s[i * n + j] = 1.0;

                    let dx = ((i as f32) + 0.5) * f.h - x;
                    let dy = ((j as f32) + 0.5) * f.h - y;

                    if dx * dx + dy * dy < r * r {
                        f.s[i * n + j] = 0.0;
                        if self.scene_nr == 2 {
                            f.m[i * n + j] = 0.5 + 0.5 * (0.1 * (self.frame_nr as f32)).sin();
                        } else {
                            f.m[i * n + j] = 1.0;
                        }
                        f.u[i * n + j] = vx;
                        f.u[(i + 1) * n + j] = vx;
                        f.v[i * n + j] = vy;
                        f.v[i * n + j + 1] = vy;
                    }
                }
            }
        }

        self.show_obstacle = true;
    }

    pub fn simulate(&mut self) {
        if !self.paused {
            match self.implementation {
                Implementation::Euler => self.simulate_euler(),
                Implementation::Flip => self.simulate_flip(),
            }
        }
        self.frame_nr += 1;
    }

    fn simulate_euler(&mut self) {
        // Get current particle count before mutable borrow
        let current_particles = self.get_particle_count();
        let emitter_enabled = self.emitter_enabled;
        let max_particles = self.max_particles;
        let emitter_x = self.emitter_x;
        let emitter_y = self.emitter_y;
        let emitter_radius = self.emitter_radius;
        let wind_x = self.wind_x;
        let wind_y = self.wind_y;
        let current_frame = self.frame_nr;
        let border_size = self.border_size;

        if let Some(ref mut f) = self.fluid {
            let n = f.num_y;

            // Note: The obstacle is now rendered separately and doesn't modify the grid
            // This prevents the ball from erasing walls when moved

            // Apply wind as a global force - gently push all fluid cells
            // This makes wind changes more responsive
            let wind_strength = 0.1; // How quickly wind affects velocity
            for i in 1..f.num_x - 1 {
                for j in 1..f.num_y - 1 {
                    if f.s[i * n + j] != 0.0 {
                        // Blend current velocity towards wind velocity
                        f.u[i * n + j] += (wind_x - f.u[i * n + j]) * wind_strength;
                        f.v[i * n + j] += (wind_y - f.v[i * n + j]) * wind_strength;
                    }
                }
            }

            // Also set wind at inlet boundaries for inflow
            for j in 0..f.num_y {
                if f.s[1 * n + j] != 0.0 {
                    f.u[1 * n + j] = wind_x;
                    f.v[1 * n + j] = wind_y;
                }
            }

            // Remove oldest particles if we're at or over the limit
            if current_particles >= max_particles {
                let particles_to_remove =
                    current_particles - max_particles + (max_particles / 20).max(1);
                f.remove_oldest_particles(particles_to_remove);
            }

            // Emit particles from emitter if enabled
            if emitter_enabled {
                let emitter_i = ((emitter_x * (f.num_x as f32)) as usize).min(f.num_x - 2).max(1);
                let emitter_j = ((emitter_y * (f.num_y as f32)) as usize).min(f.num_y - 2).max(1);
                let emit_radius = (emitter_radius as i32).max(1);

                for di in -emit_radius..=emit_radius {
                    for dj in -emit_radius..=emit_radius {
                        let i = ((emitter_i as i32) + di) as usize;
                        let j = ((emitter_j as i32) + dj) as usize;

                        if i >= 1 && i < f.num_x - 1 && j >= 1 && j < f.num_y - 1 {
                            let dist_sq = (di * di + dj * dj) as f32;
                            let r_sq = (emit_radius * emit_radius) as f32;

                            if dist_sq <= r_sq && f.s[i * n + j] != 0.0 {
                                // Only set age if this is a new smoke cell
                                if f.m[i * n + j] > 0.5 {
                                    f.particle_age[i * n + j] = current_frame;
                                }
                                f.m[i * n + j] = 0.0; // Emit smoke/fluid
                            }
                        }
                    }
                }
            }

            // Apply speed multiplier to dt
            let effective_dt = self.dt * self.speed.abs();
            if self.speed != 0.0 {
                f.simulate(effective_dt, self.gravity, self.num_iters, self.over_relaxation);
            }

            // Remove particles that reach the border
            f.remove_border_particles(border_size);
        }
    }

    fn simulate_flip(&mut self) {
        let gravity = self.gravity;
        let dt = self.dt * self.speed.abs();
        let emitter_enabled = self.emitter_enabled;
        let emitter_x = self.emitter_x;
        let emitter_y = self.emitter_y;
        let emitter_radius = self.emitter_radius;
        let obstacle_x = self.obstacle_x;
        let obstacle_y = self.obstacle_y;
        let obstacle_radius = self.obstacle_radius;

        if dt == 0.0 {
            return;
        }

        if let Some(ref mut flip) = self.flip_fluid {
            // Emit particles from emitter if enabled
            if emitter_enabled && flip.num_particles < flip.max_particles {
                // Convert emitter position (0-1) to world coordinates
                let world_x = emitter_x * (flip.num_x as f32) * flip.h;
                let world_y = emitter_y * (flip.num_y as f32) * flip.h;
                let emit_radius = emitter_radius * flip.h;

                // Minimum distance between particles (2x particle radius)
                let min_dist = flip.particle_radius * 2.0;

                // Build a quick occupancy grid for the emitter region
                // Use cells of size min_dist for fast neighbor lookup
                let cell_size = min_dist;
                let grid_min_x = world_x - emit_radius - cell_size;
                let grid_min_y = world_y - emit_radius - cell_size;
                let grid_size =
                    (((emit_radius * 2.0 + cell_size * 2.0) / cell_size).ceil() as usize) + 1;

                // Mark cells that have particles
                let mut occupied: Vec<bool> = vec![false; grid_size * grid_size];

                for i in 0..flip.num_particles {
                    let px = flip.particle_pos[2 * i];
                    let py = flip.particle_pos[2 * i + 1];

                    // Check if particle is near the emitter region
                    let dx = px - world_x;
                    let dy = py - world_y;
                    if dx * dx + dy * dy <= (emit_radius + cell_size * 2.0).powi(2) {
                        let cx = ((px - grid_min_x) / cell_size).floor() as i32;
                        let cy = ((py - grid_min_y) / cell_size).floor() as i32;
                        if
                            cx >= 0 &&
                            cy >= 0 &&
                            (cx as usize) < grid_size &&
                            (cy as usize) < grid_size
                        {
                            occupied[(cx as usize) * grid_size + (cy as usize)] = true;
                        }
                    }
                }

                // Add particles in a circle around the emitter
                let spacing = min_dist;
                let num_per_axis = (emit_radius / spacing).ceil() as i32;

                // Limit how many particles we add per frame
                let max_new_particles_per_frame = 20;
                let mut particles_added = 0;

                'outer: for di in -num_per_axis..=num_per_axis {
                    for dj in -num_per_axis..=num_per_axis {
                        if particles_added >= max_new_particles_per_frame {
                            break 'outer;
                        }
                        if flip.num_particles >= flip.max_particles {
                            break 'outer;
                        }

                        let px = world_x + (di as f32) * spacing;
                        let py = world_y + (dj as f32) * spacing;

                        let dx = px - world_x;
                        let dy = py - world_y;

                        // Check if within emitter radius
                        if dx * dx + dy * dy > emit_radius * emit_radius {
                            continue;
                        }

                        // Check occupancy grid - look at this cell and neighbors
                        let cx = ((px - grid_min_x) / cell_size).floor() as i32;
                        let cy = ((py - grid_min_y) / cell_size).floor() as i32;

                        let mut cell_occupied = false;
                        for dcx in -1..=1 {
                            for dcy in -1..=1 {
                                let ncx = cx + dcx;
                                let ncy = cy + dcy;
                                if
                                    ncx >= 0 &&
                                    ncy >= 0 &&
                                    (ncx as usize) < grid_size &&
                                    (ncy as usize) < grid_size
                                {
                                    if occupied[(ncx as usize) * grid_size + (ncy as usize)] {
                                        cell_occupied = true;
                                        break;
                                    }
                                }
                            }
                            if cell_occupied {
                                break;
                            }
                        }

                        if cell_occupied {
                            continue;
                        }

                        // Add particle
                        let idx = flip.num_particles;
                        flip.particle_pos[2 * idx] = px;
                        flip.particle_pos[2 * idx + 1] = py;
                        flip.particle_vel[2 * idx] = 0.0; // Start with zero velocity, gravity will take over
                        flip.particle_vel[2 * idx + 1] = 0.0;
                        flip.particle_color[3 * idx] = 0.0;
                        flip.particle_color[3 * idx + 1] = 0.4;
                        flip.particle_color[3 * idx + 2] = 1.0;
                        flip.num_particles += 1;
                        particles_added += 1;

                        // Mark this cell as occupied for subsequent checks
                        if
                            cx >= 0 &&
                            cy >= 0 &&
                            (cx as usize) < grid_size &&
                            (cy as usize) < grid_size
                        {
                            occupied[(cx as usize) * grid_size + (cy as usize)] = true;
                        }
                    }
                }
            }

            // Scale down gravity for FLIP (it's more sensitive to forces)
            // Using a smaller scale to prevent particles from going through walls
            let scaled_gravity = gravity * 0.3;

            // No wind for FLIP - particles just fall with gravity

            // Use sub-stepping for stability with fast particles
            let num_sub_steps = 4; // More sub-steps for better collision handling
            let sub_dt = dt / (num_sub_steps as f32);

            // Run FLIP simulation with sub-stepping using configurable parameters
            let flip_ratio = self.flip_ratio;
            let num_pressure_iters = self.num_iters;
            let num_particle_iters = 2;
            let over_relaxation = self.over_relaxation;
            let compensate_drift = self.compensate_drift;
            let separate_particles = self.separate_particles;

            // Convert obstacle position from normalized (0-1) to world coordinates
            let world_obstacle_x = obstacle_x * (flip.num_x as f32) * flip.h;
            let world_obstacle_y = obstacle_y * (flip.num_y as f32) * flip.h;
            let world_obstacle_radius = obstacle_radius * (flip.num_x as f32) * flip.h;

            for _ in 0..num_sub_steps {
                flip.simulate(
                    sub_dt,
                    scaled_gravity,
                    flip_ratio,
                    num_pressure_iters,
                    num_particle_iters,
                    over_relaxation,
                    compensate_drift,
                    separate_particles,
                    world_obstacle_x,
                    world_obstacle_y,
                    world_obstacle_radius,
                    0.0, // obstacle_vel_x (TODO: track obstacle velocity)
                    0.0 // obstacle_vel_y
                );
            }

            // Remove particles at boundaries
            let h = flip.h;
            let border = ((self.border_size as f32) + 1.0) * h;
            let max_x = (flip.num_x as f32) * h - border;
            let max_y = (flip.num_y as f32) * h - border;

            // Remove particles outside bounds (swap with last particle)
            let mut i = 0;
            while i < flip.num_particles {
                let x = flip.particle_pos[2 * i];
                let y = flip.particle_pos[2 * i + 1];

                if x < border || x > max_x || y < border || y > max_y {
                    // Swap with last particle
                    let last = flip.num_particles - 1;
                    if i != last {
                        flip.particle_pos[2 * i] = flip.particle_pos[2 * last];
                        flip.particle_pos[2 * i + 1] = flip.particle_pos[2 * last + 1];
                        flip.particle_vel[2 * i] = flip.particle_vel[2 * last];
                        flip.particle_vel[2 * i + 1] = flip.particle_vel[2 * last + 1];
                        flip.particle_color[3 * i] = flip.particle_color[3 * last];
                        flip.particle_color[3 * i + 1] = flip.particle_color[3 * last + 1];
                        flip.particle_color[3 * i + 2] = flip.particle_color[3 * last + 2];
                    }
                    flip.num_particles -= 1;
                } else {
                    i += 1;
                }
            }
        }
    }

    pub fn time_step(&mut self, _dt: f32) -> u8 {
        if self.paused {
            return 0;
        }
        self.simulate();
        0
    }

    pub fn reset(&mut self) {
        self.setup_blank(16.0, 9.0);
    }

    // Drawing methods
    pub fn paint_fluid(&mut self, x: f32, y: f32, brush_size: f32) {
        match self.implementation {
            Implementation::Euler => {
                let current_frame = self.frame_nr;
                if let Some(ref mut f) = self.fluid {
                    // Convert pixel coordinates to grid cell indices
                    // brush_size is in pixels, cell_size is pixels per cell
                    let grid_radius = (brush_size / self.cell_size).max(1.0).ceil() as i32;

                    // Use floor for proper pixel-to-cell mapping
                    let gi = (x / self.cell_size).floor() as i32;
                    let gj = (y / self.cell_size).floor() as i32;

                    for di in -grid_radius..=grid_radius {
                        for dj in -grid_radius..=grid_radius {
                            let i = (gi + di) as usize;
                            let j = (gj + dj) as usize;

                            if i >= 1 && i < f.num_x - 1 && j >= 1 && j < f.num_y - 1 {
                                let n = f.num_y;
                                let dx = di as f32;
                                let dy = dj as f32;
                                let dist_sq = dx * dx + dy * dy;
                                let r_sq = (grid_radius as f32) * (grid_radius as f32);

                                if dist_sq <= r_sq {
                                    if f.s[i * n + j] != 0.0 {
                                        // Only set age if this is a new smoke cell
                                        if f.m[i * n + j] > 0.5 {
                                            f.particle_age[i * n + j] = current_frame;
                                        }
                                        f.m[i * n + j] = 0.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            Implementation::Flip => {
                if let Some(ref mut flip) = self.flip_fluid {
                    // Convert pixel coordinates to world coordinates
                    let scale = flip.h / self.cell_size;
                    let world_x = x * scale;
                    let world_y = y * scale;
                    let radius = brush_size * scale;

                    // Minimum distance between particles (2x particle radius)
                    let min_dist = flip.particle_radius * 2.0;
                    let min_dist_sq = min_dist * min_dist;

                    // Add particles in a circle around brush position
                    let spacing = min_dist;
                    let num_per_axis = (radius / spacing).ceil() as i32;

                    for di in -num_per_axis..=num_per_axis {
                        for dj in -num_per_axis..=num_per_axis {
                            let px = world_x + (di as f32) * spacing;
                            let py = world_y + (dj as f32) * spacing;

                            let dx = px - world_x;
                            let dy = py - world_y;

                            if dx * dx + dy * dy <= radius * radius {
                                // Check bounds
                                if
                                    px > flip.h &&
                                    px < ((flip.num_x - 1) as f32) * flip.h &&
                                    py > flip.h &&
                                    py < ((flip.num_y - 1) as f32) * flip.h
                                {
                                    // Check if there's already a particle too close
                                    let mut too_close = false;
                                    for k in 0..flip.num_particles {
                                        let ex = flip.particle_pos[2 * k];
                                        let ey = flip.particle_pos[2 * k + 1];
                                        let edx = ex - px;
                                        let edy = ey - py;
                                        if edx * edx + edy * edy < min_dist_sq {
                                            too_close = true;
                                            break;
                                        }
                                    }

                                    if too_close {
                                        continue;
                                    }

                                    // Check if we have space for more particles
                                    if flip.num_particles < flip.max_particles {
                                        let idx = flip.num_particles;
                                        flip.particle_pos[2 * idx] = px;
                                        flip.particle_pos[2 * idx + 1] = py;
                                        flip.particle_vel[2 * idx] = 0.0; // Start with zero velocity
                                        flip.particle_vel[2 * idx + 1] = 0.0;
                                        flip.particle_color[3 * idx] = 0.0;
                                        flip.particle_color[3 * idx + 1] = 0.4;
                                        flip.particle_color[3 * idx + 2] = 1.0;
                                        flip.num_particles += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn paint_wall(&mut self, x: f32, y: f32, brush_size: f32) {
        // Works the same for both Euler and FLIP (modifies the grid)
        if let Some(ref mut f) = self.fluid {
            // Convert pixel coordinates to grid cell indices
            let grid_radius = (brush_size / self.cell_size).max(1.0).ceil() as i32;

            let gi = (x / self.cell_size).floor() as i32;
            let gj = (y / self.cell_size).floor() as i32;

            for di in -grid_radius..=grid_radius {
                for dj in -grid_radius..=grid_radius {
                    let i = (gi + di) as usize;
                    let j = (gj + dj) as usize;

                    if i >= 1 && i < f.num_x - 1 && j >= 1 && j < f.num_y - 1 {
                        let n = f.num_y;
                        let dx = di as f32;
                        let dy = dj as f32;
                        let dist_sq = dx * dx + dy * dy;
                        let r_sq = (grid_radius as f32) * (grid_radius as f32);

                        if dist_sq <= r_sq {
                            // Set cell as solid
                            f.s[i * n + j] = 0.0;
                            // Set smoke to 1.0 (no smoke in solid)
                            f.m[i * n + j] = 1.0;
                            // Zero velocities at all cell boundaries (staggered grid)
                            f.u[i * n + j] = 0.0;
                            if i + 1 < f.num_x {
                                f.u[(i + 1) * n + j] = 0.0;
                            }
                            f.v[i * n + j] = 0.0;
                            if j + 1 < f.num_y {
                                f.v[i * n + j + 1] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        // Also set solid cells in FLIP fluid
        if let Some(ref mut flip) = self.flip_fluid {
            let grid_radius = (brush_size / self.cell_size).max(1.0).ceil() as i32;

            let gi = (x / self.cell_size).floor() as i32;
            let gj = (y / self.cell_size).floor() as i32;

            for di in -grid_radius..=grid_radius {
                for dj in -grid_radius..=grid_radius {
                    let i = (gi + di) as usize;
                    let j = (gj + dj) as usize;

                    if i >= 1 && i < flip.num_x - 1 && j >= 1 && j < flip.num_y - 1 {
                        let n = flip.num_y;
                        let dx = di as f32;
                        let dy = dj as f32;
                        let dist_sq = dx * dx + dy * dy;
                        let r_sq = (grid_radius as f32) * (grid_radius as f32);

                        if dist_sq <= r_sq {
                            flip.s[i * n + j] = 0.0;
                        }
                    }
                }
            }
        }
    }

    pub fn erase(&mut self, x: f32, y: f32, brush_size: f32) {
        match self.implementation {
            Implementation::Euler => {
                if let Some(ref mut f) = self.fluid {
                    // Convert pixel coordinates to grid cell indices
                    let grid_radius = (brush_size / self.cell_size).max(1.0).ceil() as i32;

                    let gi = (x / self.cell_size).floor() as i32;
                    let gj = (y / self.cell_size).floor() as i32;

                    for di in -grid_radius..=grid_radius {
                        for dj in -grid_radius..=grid_radius {
                            let i = (gi + di) as usize;
                            let j = (gj + dj) as usize;

                            if i >= 1 && i < f.num_x - 1 && j >= 1 && j < f.num_y - 1 {
                                let n = f.num_y;
                                let dx = di as f32;
                                let dy = dj as f32;
                                let dist_sq = dx * dx + dy * dy;
                                let r_sq = (grid_radius as f32) * (grid_radius as f32);

                                if dist_sq <= r_sq {
                                    f.s[i * n + j] = 1.0;
                                    f.m[i * n + j] = 1.0;
                                    f.particle_age[i * n + j] = 0;
                                    f.u[i * n + j] = 0.0;
                                    f.v[i * n + j] = 0.0;
                                }
                            }
                        }
                    }
                }
            }
            Implementation::Flip => {
                if let Some(ref mut flip) = self.flip_fluid {
                    // Convert pixel coordinates to world coordinates
                    let scale = flip.h / self.cell_size;
                    let world_x = x * scale;
                    let world_y = y * scale;
                    let radius = brush_size * scale;

                    // Remove particles within brush radius
                    let radius_sq = radius * radius;
                    let mut i = 0;
                    while i < flip.num_particles {
                        let px = flip.particle_pos[2 * i];
                        let py = flip.particle_pos[2 * i + 1];
                        let dx = px - world_x;
                        let dy = py - world_y;

                        if dx * dx + dy * dy <= radius_sq {
                            // Swap with last particle
                            let last = flip.num_particles - 1;
                            if i != last {
                                flip.particle_pos[2 * i] = flip.particle_pos[2 * last];
                                flip.particle_pos[2 * i + 1] = flip.particle_pos[2 * last + 1];
                                flip.particle_vel[2 * i] = flip.particle_vel[2 * last];
                                flip.particle_vel[2 * i + 1] = flip.particle_vel[2 * last + 1];
                                flip.particle_color[3 * i] = flip.particle_color[3 * last];
                                flip.particle_color[3 * i + 1] = flip.particle_color[3 * last + 1];
                                flip.particle_color[3 * i + 2] = flip.particle_color[3 * last + 2];
                            }
                            flip.num_particles -= 1;
                        } else {
                            i += 1;
                        }
                    }

                    // Also clear solid cells in the area
                    let grid_radius = (brush_size / self.cell_size).max(1.0).ceil() as i32;
                    let gi = (x / self.cell_size).floor() as i32;
                    let gj = (y / self.cell_size).floor() as i32;

                    for di in -grid_radius..=grid_radius {
                        for dj in -grid_radius..=grid_radius {
                            let ii = (gi + di) as usize;
                            let jj = (gj + dj) as usize;

                            if ii >= 1 && ii < flip.num_x - 1 && jj >= 1 && jj < flip.num_y - 1 {
                                let n = flip.num_y;
                                let dx = di as f32;
                                let dy = dj as f32;
                                let dist_sq = dx * dx + dy * dy;
                                let r_sq = (grid_radius as f32) * (grid_radius as f32);

                                if dist_sq <= r_sq {
                                    flip.s[ii * n + jj] = 1.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn clear_particles(&mut self) {
        match self.implementation {
            Implementation::Euler => {
                if let Some(ref mut f) = self.fluid {
                    f.m.fill(1.0);
                }
            }
            Implementation::Flip => {
                if let Some(ref mut flip) = self.flip_fluid {
                    flip.num_particles = 0;
                }
            }
        }
    }

    pub fn clear_walls(&mut self) {
        if let Some(ref mut f) = self.fluid {
            let n = f.num_y;
            for i in 1..f.num_x - 1 {
                for j in 1..f.num_y - 1 {
                    f.s[i * n + j] = 1.0;
                }
            }
        }

        if let Some(ref mut flip) = self.flip_fluid {
            let n = flip.num_y;
            for i in 1..flip.num_x - 1 {
                for j in 1..flip.num_y - 1 {
                    flip.s[i * n + j] = 1.0;
                }
            }
        }
    }

    // Get particles for rendering - now with BLUE color
    // Handles both Euler and FLIP implementations
    pub fn get_particles(&self) -> JsValue {
        let mut particles: Vec<FluidParticle> = Vec::new();

        match self.implementation {
            Implementation::Euler => {
                if let Some(ref f) = self.fluid {
                    let n = f.num_y;

                    for i in 1..f.num_x - 1 {
                        for j in 1..f.num_y - 1 {
                            let smoke = f.m[i * n + j];
                            if smoke < 0.99 && f.s[i * n + j] != 0.0 {
                                let x = ((i as f32) + 0.5) * self.cell_size;
                                let y = ((j as f32) + 0.5) * self.cell_size;

                                let u = (f.u[i * n + j] + f.u[(i + 1) * n + j]) * 0.5;
                                let v = (f.v[i * n + j] + f.v[i * n + j + 1]) * 0.5;

                                // Blue color based on smoke density
                                // More smoke (lower m) = darker blue
                                let intensity = smoke;
                                let r = (intensity * 100.0) as u32; // Low red
                                let g = (intensity * 150.0) as u32; // Medium green
                                let b = (200.0 + intensity * 55.0) as u32; // High blue (200-255)
                                let color = (r << 16) | (g << 8) | b | 0xff000000;

                                particles.push(
                                    FluidParticle::new(x, y, u * 100.0, v * 100.0, color)
                                );
                            }
                        }
                    }
                }
            }
            Implementation::Flip => {
                if let Some(ref flip) = self.flip_fluid {
                    // Scale factor from world coordinates to pixel coordinates
                    let scale = self.cell_size / flip.h;

                    for i in 0..flip.num_particles {
                        let x = flip.particle_pos[2 * i] * scale;
                        let y = flip.particle_pos[2 * i + 1] * scale;
                        let vx = flip.particle_vel[2 * i];
                        let vy = flip.particle_vel[2 * i + 1];

                        // Get color from particle (RGB floats 0-1)
                        let r = (flip.particle_color[3 * i] * 255.0) as u32;
                        let g = (flip.particle_color[3 * i + 1] * 255.0) as u32;
                        let b = (flip.particle_color[3 * i + 2] * 255.0) as u32;
                        let color = (r << 16) | (g << 8) | b | 0xff000000;

                        particles.push(FluidParticle::new(x, y, vx * 100.0, vy * 100.0, color));
                    }
                }
            }
        }

        serde_wasm_bindgen::to_value(&particles).unwrap()
    }

    // Get walls for rendering - returns (x, y, is_border) tuples
    // Border is visual only (not a physics wall), user walls are actual obstacles
    pub fn get_walls(&self) -> JsValue {
        let mut walls: Vec<(i32, i32, bool)> = Vec::new();

        // Get the correct fluid's s array based on implementation
        let (num_x, num_y, s_array): (usize, usize, Option<&Vec<f32>>) = match self.implementation {
            Implementation::Euler => {
                if let Some(ref f) = self.fluid {
                    (f.num_x, f.num_y, Some(&f.s))
                } else {
                    return serde_wasm_bindgen::to_value(&walls).unwrap();
                }
            }
            Implementation::Flip => {
                if let Some(ref flip) = self.flip_fluid {
                    (flip.num_x, flip.num_y, Some(&flip.s))
                } else {
                    return serde_wasm_bindgen::to_value(&walls).unwrap();
                }
            }
        };

        if let Some(s) = s_array {
            let n = num_y;
            let bs = self.border_size.max(1);

            // Add visual border (not physics walls) with configurable thickness
            for b in 0..bs {
                for i in 0..num_x {
                    // Bottom rows
                    if b < num_y {
                        walls.push((i as i32, b as i32, true));
                    }
                    // Top rows
                    if num_y - 1 - b < num_y {
                        walls.push((i as i32, (num_y - 1 - b) as i32, true));
                    }
                }
                for j in bs..num_y - bs {
                    // Left columns
                    if b < num_x {
                        walls.push((b as i32, j as i32, true));
                    }
                    // Right columns
                    if num_x - 1 - b < num_x {
                        walls.push(((num_x - 1 - b) as i32, j as i32, true));
                    }
                }
            }

            // Add user-placed walls (actual physics obstacles)
            for i in bs..num_x - bs {
                for j in bs..num_y - bs {
                    if s[i * n + j] == 0.0 {
                        walls.push((i as i32, j as i32, false));
                    }
                }
            }
        }

        serde_wasm_bindgen::to_value(&walls).unwrap()
    }

    pub fn get_particle_count(&self) -> usize {
        match self.implementation {
            Implementation::Euler => {
                if let Some(ref f) = self.fluid {
                    let n = f.num_y;
                    let mut count = 0;

                    for i in 1..f.num_x - 1 {
                        for j in 1..f.num_y - 1 {
                            if f.m[i * n + j] < 0.99 && f.s[i * n + j] != 0.0 {
                                count += 1;
                            }
                        }
                    }
                    count
                } else {
                    0
                }
            }
            Implementation::Flip => {
                if let Some(ref flip) = self.flip_fluid { flip.num_particles } else { 0 }
            }
        }
    }

    // Pressure range for color mapping
    pub fn get_pressure_range(&self) -> Vec<f32> {
        match &self.fluid {
            Some(f) => {
                let mut min_p = f.p[0];
                let mut max_p = f.p[0];
                for &p in &f.p {
                    if p < min_p {
                        min_p = p;
                    }
                    if p > max_p {
                        max_p = p;
                    }
                }
                vec![min_p, max_p]
            }
            None => vec![0.0, 0.0],
        }
    }

    // Toggle methods
    pub fn toggle_pause(&mut self) {
        self.paused = !self.paused;
    }
    pub fn toggle_streamlines(&mut self) {
        self.show_streamlines = !self.show_streamlines;
    }
    pub fn toggle_velocities(&mut self) {
        self.show_velocities = !self.show_velocities;
    }
    pub fn toggle_pressure(&mut self) {
        self.show_pressure = !self.show_pressure;
    }
    pub fn toggle_smoke(&mut self) {
        self.show_smoke = !self.show_smoke;
    }
    pub fn toggle_density(&mut self) {
        self.show_density = !self.show_density;
    }
    pub fn toggle_over_relaxation(&mut self) {
        if self.over_relaxation == 1.0 {
            self.over_relaxation = 1.9;
        } else {
            self.over_relaxation = 1.0;
        }
    }

    // Getters
    pub fn get_gravity(&self) -> f32 {
        self.gravity
    }
    pub fn get_dt(&self) -> f32 {
        self.dt
    }
    pub fn get_num_iters(&self) -> usize {
        self.num_iters
    }
    pub fn get_over_relaxation(&self) -> f32 {
        self.over_relaxation
    }
    pub fn get_is_paused(&self) -> bool {
        self.paused
    }
    pub fn get_scene_nr(&self) -> u8 {
        self.scene_nr
    }
    pub fn get_cell_size(&self) -> f32 {
        self.cell_size
    }
    pub fn get_particle_radius(&self) -> f32 {
        self.particle_radius
    }
    pub fn get_speed(&self) -> f32 {
        self.speed
    }
    pub fn get_implementation(&self) -> Implementation {
        self.implementation
    }
    pub fn get_show_pressure(&self) -> bool {
        self.show_pressure
    }
    pub fn get_show_smoke(&self) -> bool {
        self.show_smoke
    }
    pub fn get_show_density(&self) -> bool {
        self.show_density
    }
    pub fn get_show_streamlines(&self) -> bool {
        self.show_streamlines
    }
    pub fn get_show_velocities(&self) -> bool {
        self.show_velocities
    }

    // Obstacle getters
    pub fn get_obstacle_x(&self) -> f32 {
        self.obstacle_x
    }
    pub fn get_obstacle_y(&self) -> f32 {
        self.obstacle_y
    }
    pub fn get_obstacle_radius(&self) -> f32 {
        self.obstacle_radius
    }
    pub fn get_show_obstacle(&self) -> bool {
        self.show_obstacle
    }

    // Setters
    pub fn set_gravity(&mut self, gravity: f32) {
        self.gravity = gravity;
    }
    pub fn set_dt(&mut self, dt: f32) {
        self.dt = dt;
    }
    pub fn set_num_iters(&mut self, num_iters: usize) {
        self.num_iters = num_iters;
    }
    pub fn set_over_relaxation(&mut self, over_relaxation: f32) {
        self.over_relaxation = over_relaxation;
    }
    pub fn set_is_paused(&mut self, paused: bool) {
        self.paused = paused;
    }
    pub fn set_cell_size(&mut self, size: f32) {
        self.cell_size = size;
    }
    pub fn set_particle_radius(&mut self, radius: f32) {
        self.particle_radius = radius;
    }
    pub fn set_speed(&mut self, speed: f32) {
        self.speed = speed;
    }
    pub fn set_implementation(&mut self, implementation: Implementation) {
        if self.implementation != implementation {
            self.implementation = implementation;
            // Reset simulation when switching implementations
            self.setup_blank(16.0, 9.0);
        }
    }
    pub fn set_show_density(&mut self, show: bool) {
        self.show_density = show;
    }

    // Obstacle setters
    pub fn set_obstacle_radius(&mut self, radius: f32) {
        self.obstacle_radius = radius;
    }
    pub fn set_show_obstacle(&mut self, show: bool) {
        self.show_obstacle = show;
    }

    // Wind getters/setters
    pub fn get_wind_x(&self) -> f32 {
        self.wind_x
    }
    pub fn get_wind_y(&self) -> f32 {
        self.wind_y
    }
    pub fn set_wind_x(&mut self, wind_x: f32) {
        self.wind_x = wind_x;
    }
    pub fn set_wind_y(&mut self, wind_y: f32) {
        self.wind_y = wind_y;
    }

    // Emitter getters/setters
    pub fn get_emitter_enabled(&self) -> bool {
        self.emitter_enabled
    }
    pub fn get_emitter_x(&self) -> f32 {
        self.emitter_x
    }
    pub fn get_emitter_y(&self) -> f32 {
        self.emitter_y
    }
    pub fn get_emitter_radius(&self) -> f32 {
        self.emitter_radius
    }
    pub fn get_max_particles(&self) -> usize {
        self.max_particles
    }
    pub fn set_emitter_enabled(&mut self, enabled: bool) {
        self.emitter_enabled = enabled;
    }
    pub fn set_emitter_x(&mut self, x: f32) {
        self.emitter_x = x;
    }
    pub fn set_emitter_y(&mut self, y: f32) {
        self.emitter_y = y;
    }
    pub fn set_emitter_radius(&mut self, radius: f32) {
        self.emitter_radius = radius;
    }
    pub fn set_max_particles(&mut self, max: usize) {
        self.max_particles = max;
    }
    pub fn toggle_emitter(&mut self) {
        self.emitter_enabled = !self.emitter_enabled;
    }

    // Border getters/setters
    pub fn get_border_size(&self) -> usize {
        self.border_size
    }
    pub fn set_border_size(&mut self, size: usize) {
        self.border_size = size.max(1); // Minimum 1 cell border
    }

    // Grid size getters/setters
    pub fn get_grid_size(&self) -> usize {
        self.grid_size
    }
    pub fn set_grid_size(&mut self, size: usize) {
        let new_size = size.clamp(20, 500); // Clamp between 20 and 500 cells
        if new_size != self.grid_size {
            self.grid_size = new_size;
            // Reinitialize the grid with the new size
            self.setup_blank(16.0, 9.0);
        }
    }

    // FLIP-specific getters/setters
    pub fn get_flip_ratio(&self) -> f32 {
        self.flip_ratio
    }
    pub fn set_flip_ratio(&mut self, ratio: f32) {
        self.flip_ratio = ratio.clamp(0.0, 1.0); // 0 = PIC, 1 = FLIP
    }
    pub fn get_compensate_drift(&self) -> bool {
        self.compensate_drift
    }
    pub fn set_compensate_drift(&mut self, compensate: bool) {
        self.compensate_drift = compensate;
    }
    pub fn get_separate_particles(&self) -> bool {
        self.separate_particles
    }
    pub fn set_separate_particles(&mut self, separate: bool) {
        self.separate_particles = separate;
    }

    // Fluid data access
    pub fn get_fluid_num_x(&self) -> usize {
        match &self.fluid {
            Some(f) => f.num_x,
            None => 0,
        }
    }

    pub fn get_fluid_num_y(&self) -> usize {
        match &self.fluid {
            Some(f) => f.num_y,
            None => 0,
        }
    }

    pub fn get_fluid_h(&self) -> f32 {
        match &self.fluid {
            Some(f) => f.h,
            None => 0.0,
        }
    }

    pub fn sample_field(&self, x: f32, y: f32, field: u8) -> f32 {
        match &self.fluid {
            Some(f) => f.sample_field(x, y, field),
            None => 0.0,
        }
    }

    pub fn get_fluid_u(&self, i: usize, j: usize) -> f32 {
        match &self.fluid {
            Some(f) => f.get_u(i, j),
            None => 0.0,
        }
    }

    pub fn get_fluid_v(&self, i: usize, j: usize) -> f32 {
        match &self.fluid {
            Some(f) => f.get_v(i, j),
            None => 0.0,
        }
    }

    pub fn get_fluid_p(&self, i: usize, j: usize) -> f32 {
        match &self.fluid {
            Some(f) => f.get_p(i, j),
            None => 0.0,
        }
    }

    pub fn get_fluid_s(&self, i: usize, j: usize) -> f32 {
        match &self.fluid {
            Some(f) => f.get_s(i, j),
            None => 0.0,
        }
    }

    pub fn get_fluid_m(&self, i: usize, j: usize) -> f32 {
        match &self.fluid {
            Some(f) => f.get_m(i, j),
            None => 0.0,
        }
    }

    pub fn get_data(&self) -> JsValue {
        serde_wasm_bindgen::to_value(&self).unwrap()
    }
}
