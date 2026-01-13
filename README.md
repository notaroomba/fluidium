<h1 align="center">
    <br>
    <a href="https://fluidium.notaroomba.dev"><img src="public/logo.png" alt="Fluidium" width="200"></a>
    <br>
    Fluidium
    <br>
</h1>

<p align="center">
  Interactive 2D fluid dynamics simulation with Euler and FLIP solvers, powered by a Rust physics engine compiled to WebAssembly.
</p>

<div align="center">

![Rust](https://img.shields.io/badge/rust-%23000000.svg?style=for-the-badge&logo=rust&logoColor=white)
![WASM](https://img.shields.io/badge/WebAssembly-%2300969C.svg?style=for-the-badge&logo=webassembly&logoColor=white)
![React](https://img.shields.io/badge/react-%2320232a.svg?style=for-the-badge&logo=react&logoColor=%2361DAFB)
![PixiJS](https://img.shields.io/badge/Pixi.js-%23FF66CC.svg?style=for-the-badge)

</div>

An interactive fluid simulation sandbox that combines a high-performance physics core (written in Rust and compiled to WebAssembly) with a responsive web UI (React + TypeScript) and PixiJS for high-performance 2D rendering.

<p align="center">
  <a href="https://fluidium.notaroomba.dev">Live demo</a>
</p>

## Key features

- Physics engine implemented in Rust and compiled to WebAssembly for accurate, real-time simulation
- Two fluid solvers: Euler (grid-based) and FLIP (particle-based hybrid)
- Interactive drawing tools: paint fluid, walls, and erase
- Configurable wind vector and particle emitter position
- Draggable obstacle (ball) that interacts with fluid
- Adjustable simulation parameters: gravity, time step, solver iterations, over-relaxation
- FLIP-specific settings: PIC/FLIP ratio, drift compensation, particle separation
- Visual overlays: grid, velocity vectors, density visualization, FPS counter
- Responsive UI with draggable property editor panels

## Preview

Open the app locally with `bun run dev` and try:
- Draw fluid and walls using the pencil tool from the sidebar
- Drag the obstacle ball to interact with the fluid
- Toggle between Euler and FLIP implementations
- Adjust wind direction using the vector canvas in settings
- Configure emitter position and simulation parameters

## Credits

- Rust + wasm-bindgen for the physics engine
- React + TypeScript for the UI, TailwindCSS for styling, PixiJS for rendering

## License

MIT

---

> [notaroomba.dev](https://notaroomba.dev) &nbsp;&middot;&nbsp;
> GitHub [@NotARoomba](https://github.com/NotARoomba)