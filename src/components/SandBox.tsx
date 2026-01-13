import { extend, useTick } from "@pixi/react";
import type { Universe, FluidParticle } from "physics-engine";
import { Container, Graphics } from "pixi.js";
import { useCallback, useRef, useState } from "react";
import { useSimulation } from "../contexts/SimulationContext";

extend({
  Container,
  Graphics,
});

interface SandBoxProps {
  universe: Universe;
}

export default function SandBox({ universe }: SandBoxProps) {
  const {
    render,
    setRender,
    showGrid,
    showVelocity,
    toolMode,
    brushType,
    brushSize,
    setFps,
  } = useSimulation();

  const [isDrawing, setIsDrawing] = useState(false);
  const fpsCounterRef = useRef({ frames: 0, lastTime: performance.now() });
  const pixiContainerRef = useRef<any>(null);
  const lastDrawPosRef = useRef<{ x: number; y: number } | null>(null);

  // Apply Y-flip to the container
  if (pixiContainerRef.current && pixiContainerRef.current.scale.y !== -1) {
    pixiContainerRef.current.scale.y = -1;
  }

  const handleDraw = (x: number, y: number) => {
    if (toolMode === "none") return;

    if (toolMode === "pencil") {
      switch (brushType) {
        case "water":
          universe.paint_fluid(x, y, brushSize);
          break;
        case "wall":
          universe.paint_wall(x, y, brushSize);
          break;
        case "air":
          universe.erase(x, y, brushSize);
          break;
      }
    } else if (toolMode === "eraser") {
      universe.erase(x, y, brushSize);
    }

    setRender((prev) => prev + 1);
  };

  const handlePointerDown = (event: any) => {
    if (toolMode === "none") return;

    setIsDrawing(true);
    const localPos = pixiContainerRef.current?.toLocal(event.data.global);
    if (localPos) {
      handleDraw(localPos.x, localPos.y);
      lastDrawPosRef.current = { x: localPos.x, y: localPos.y };
    }
  };

  const handlePointerMove = (event: any) => {
    if (!isDrawing || toolMode === "none") return;

    const localPos = pixiContainerRef.current?.toLocal(event.data.global);
    if (localPos) {
      // Interpolate between last position and current for smooth drawing
      if (lastDrawPosRef.current) {
        const dx = localPos.x - lastDrawPosRef.current.x;
        const dy = localPos.y - lastDrawPosRef.current.y;
        const dist = Math.sqrt(dx * dx + dy * dy);
        const step = brushSize / 4;

        if (dist > step) {
          const steps = Math.ceil(dist / step);
          for (let i = 1; i <= steps; i++) {
            const t = i / steps;
            const ix = lastDrawPosRef.current.x + dx * t;
            const iy = lastDrawPosRef.current.y + dy * t;
            handleDraw(ix, iy);
          }
        } else {
          handleDraw(localPos.x, localPos.y);
        }
      }
      lastDrawPosRef.current = { x: localPos.x, y: localPos.y };
    }
  };

  const handlePointerUp = () => {
    setIsDrawing(false);
    lastDrawPosRef.current = null;
  };

  const drawCallback = useCallback(
    (graphics: any) => {
      graphics.clear();

      const particles = universe.get_particles() as FluidParticle[];
      const walls = universe.get_walls() as [number, number][];
      const cellSize = universe.get_cell_size();
      const particleRadius = universe.get_particle_radius();

      // Draw grid if enabled
      if (showGrid) {
        const gridRange = 2000; // visible grid range
        const gridStep = cellSize;

        graphics.setStrokeStyle({
          width: 1,
          color: 0xcccccc,
          alpha: 0.3,
        });

        for (let x = -gridRange; x <= gridRange; x += gridStep) {
          graphics.moveTo(x, -gridRange);
          graphics.lineTo(x, gridRange);
          graphics.stroke();
        }

        for (let y = -gridRange; y <= gridRange; y += gridStep) {
          graphics.moveTo(-gridRange, y);
          graphics.lineTo(gridRange, y);
          graphics.stroke();
        }
      }

      // Draw walls
      for (const [gx, gy] of walls) {
        const wx = (gx + 0.5) * cellSize;
        const wy = (gy + 0.5) * cellSize;
        const halfSize = cellSize / 2;

        graphics.rect(wx - halfSize, wy - halfSize, cellSize, cellSize);
        graphics.fill({ color: 0x4a5568, alpha: 1 });
      }

      // Draw fluid particles
      for (const particle of particles) {
        graphics.circle(particle.pos.x, particle.pos.y, particleRadius);
        // Convert ARGB (u32) to RGB by masking off the alpha channel
        const rgbColor = particle.color & 0x00ffffff;
        graphics.fill({ color: rgbColor, alpha: 0.9 });
      }

      // Draw velocity vectors if enabled
      if (showVelocity) {
        for (const particle of particles) {
          const velMag = Math.sqrt(
            particle.vel.x * particle.vel.x + particle.vel.y * particle.vel.y
          );

          if (velMag > 1) {
            const scale = 0.05;
            const vx = particle.vel.x * scale;
            const vy = particle.vel.y * scale;

            graphics.setStrokeStyle({
              width: 1,
              color: 0xff0000,
              alpha: 0.7,
            });
            graphics.moveTo(particle.pos.x, particle.pos.y);
            graphics.lineTo(particle.pos.x + vx, particle.pos.y + vy);
            graphics.stroke();
          }
        }
      }

      // Draw brush preview
      if (toolMode !== "none" && pixiContainerRef.current) {
        // We'll draw this in the render loop based on mouse position
      }
    },
    [universe, render, showGrid, showVelocity, toolMode, brushSize]
  );

  useTick((delta) => {
    // Calculate FPS
    const counter = fpsCounterRef.current;
    counter.frames++;
    const now = performance.now();
    if (now - counter.lastTime >= 1000) {
      setFps(Math.round((counter.frames * 1000) / (now - counter.lastTime)));
      counter.frames = 0;
      counter.lastTime = now;
    }

    if (universe.get_is_paused()) return;

    setRender((prev) => prev + 1);
    universe.time_step(delta.deltaTime / 60);
  });

  return (
    <>
      <pixiContainer
        x={0}
        y={0}
        ref={pixiContainerRef}
        interactive
        onPointerDown={handlePointerDown}
        onGlobalPointerMove={handlePointerMove}
        onPointerUp={handlePointerUp}
        onPointerUpOutside={handlePointerUp}
      >
        <pixiGraphics draw={drawCallback} />
      </pixiContainer>
    </>
  );
}
