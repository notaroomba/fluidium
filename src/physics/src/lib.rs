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
    // Fluid grid
    fluid: Option<Fluid>,

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
}

#[wasm_bindgen]
impl Universe {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Universe {
        let mut universe = Universe {
            fluid: None,
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
        self.gravity = 0.0;
        self.wind_x = 2.0; // Default wind from left
        self.wind_y = 0.0;

        // Reset emitter to center
        self.emitter_x = 0.5;
        self.emitter_y = 0.5;

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
                    let emitter_i = ((emitter_x * (f.num_x as f32)) as usize)
                        .min(f.num_x - 2)
                        .max(1);
                    let emitter_j = ((emitter_y * (f.num_y as f32)) as usize)
                        .min(f.num_y - 2)
                        .max(1);
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
        self.frame_nr += 1;
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

    pub fn paint_wall(&mut self, x: f32, y: f32, brush_size: f32) {
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
                            // u is stored at left face of cell, so u[i] and u[i+1] bound this cell
                            // v is stored at bottom face of cell, so v[j] and v[j+1] bound this cell
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
    }

    pub fn erase(&mut self, x: f32, y: f32, brush_size: f32) {
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

    pub fn clear_particles(&mut self) {
        if let Some(ref mut f) = self.fluid {
            f.m.fill(1.0);
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
    }

    // Get particles for rendering - now with BLUE color
    pub fn get_particles(&self) -> JsValue {
        let mut particles: Vec<FluidParticle> = Vec::new();

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

                        particles.push(FluidParticle::new(x, y, u * 100.0, v * 100.0, color));
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

        if let Some(ref f) = self.fluid {
            let n = f.num_y;
            let bs = self.border_size.max(1);

            // Add visual border (not physics walls) with configurable thickness
            for b in 0..bs {
                for i in 0..f.num_x {
                    // Bottom rows
                    if b < f.num_y {
                        walls.push((i as i32, b as i32, true));
                    }
                    // Top rows
                    if f.num_y - 1 - b < f.num_y {
                        walls.push((i as i32, (f.num_y - 1 - b) as i32, true));
                    }
                }
                for j in bs..f.num_y - bs {
                    // Left columns
                    if b < f.num_x {
                        walls.push((b as i32, j as i32, true));
                    }
                    // Right columns
                    if f.num_x - 1 - b < f.num_x {
                        walls.push(((f.num_x - 1 - b) as i32, j as i32, true));
                    }
                }
            }

            // Add user-placed walls (actual physics obstacles)
            for i in bs..f.num_x - bs {
                for j in bs..f.num_y - bs {
                    if f.s[i * n + j] == 0.0 {
                        walls.push((i as i32, j as i32, false));
                    }
                }
            }
        }

        serde_wasm_bindgen::to_value(&walls).unwrap()
    }

    pub fn get_particle_count(&self) -> usize {
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
        self.implementation = implementation;
    }
    pub fn set_show_density(&mut self, show: bool) {
        self.show_density = show;
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
