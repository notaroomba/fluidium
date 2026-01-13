import { extend, useTick } from "@pixi/react";
import type { Universe, FluidParticle } from "physics-engine";
import { Container, Graphics, Rectangle } from "pixi.js";
import { useCallback, useRef, useState, useEffect } from "react";
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
  const [isDraggingObstacle, setIsDraggingObstacle] = useState(false);
  const fpsCounterRef = useRef({ frames: 0, lastTime: performance.now() });
  const pixiContainerRef = useRef<any>(null);
  const lastDrawPosRef = useRef<{ x: number; y: number } | null>(null);

  // Apply Y-flip to the container and set up hit area for pointer events
  useEffect(() => {
    if (pixiContainerRef.current) {
      if (pixiContainerRef.current.scale.y !== -1) {
        pixiContainerRef.current.scale.y = -1;
      }
      // Set a large hit area so the container receives pointer events on empty space
      // The hit area needs to be large enough to cover the visible canvas
      pixiContainerRef.current.hitArea = new Rectangle(
        -10000,
        -10000,
        20000,
        20000
      );
    }
  }, []);

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

  // Check if point is inside the obstacle (in pixel coordinates)
  const isInsideObstacle = (x: number, y: number): boolean => {
    const u = universe as any;
    if (!u.get_obstacle_x || !u.get_show_obstacle || !u.get_show_obstacle())
      return false;

    const gridDims = universe.get_grid_dimensions();
    const gridWidth = gridDims[0];
    const gridHeight = gridDims[1];
    const cellSize = universe.get_cell_size();

    // Obstacle position is in normalized coordinates (0-1), convert to pixels
    // Offset by 1 cell to match visual rendering
    const obstacleX = u.get_obstacle_x() * gridWidth - cellSize;
    const obstacleY = u.get_obstacle_y() * gridHeight - cellSize;
    const obstacleRadius = u.get_obstacle_radius() * gridWidth;

    const dx = x - obstacleX;
    const dy = y - obstacleY;
    return dx * dx + dy * dy <= obstacleRadius * obstacleRadius;
  };

  // Update obstacle position (in pixel coordinates, converts to normalized)
  const updateObstacle = (x: number, y: number) => {
    const u = universe as any;
    if (!u.set_obstacle) return;

    const gridDims = universe.get_grid_dimensions();
    const gridWidth = gridDims[0];
    const gridHeight = gridDims[1];
    const cellSize = universe.get_cell_size();

    // Convert pixel position to normalized coordinates (0-1)
    // Add 1 cell offset to account for visual offset
    const normalizedX = Math.max(
      0.1,
      Math.min(0.9, (x + cellSize) / gridWidth)
    );
    const normalizedY = Math.max(
      0.1,
      Math.min(0.9, (y + cellSize) / gridHeight)
    );

    u.set_obstacle(normalizedX, normalizedY, false);
    setRender((prev) => prev + 1);
  };

  const handlePointerDown = (event: any) => {
    const localPos = pixiContainerRef.current?.toLocal(event.data.global);
    if (!localPos) return;

    // Check if clicking on obstacle first
    if (isInsideObstacle(localPos.x, localPos.y)) {
      setIsDraggingObstacle(true);
      return;
    }

    if (toolMode === "none") return;

    setIsDrawing(true);
    handleDraw(localPos.x, localPos.y);
    lastDrawPosRef.current = { x: localPos.x, y: localPos.y };
  };

  const handlePointerMove = (event: any) => {
    const localPos = pixiContainerRef.current?.toLocal(event.data.global);
    if (!localPos) return;

    // If dragging obstacle, update its position
    if (isDraggingObstacle) {
      updateObstacle(localPos.x, localPos.y);
      return;
    }

    if (!isDrawing || toolMode === "none") return;

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
  };

  const handlePointerUp = () => {
    setIsDrawing(false);
    setIsDraggingObstacle(false);
    lastDrawPosRef.current = null;
  };

  const drawCallback = useCallback(
    (graphics: any) => {
      graphics.clear();

      const particles = universe.get_particles() as FluidParticle[];
      const walls = universe.get_walls() as [number, number, boolean][];
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

      // Draw walls - border walls are cyan, user walls are gray
      for (const [gx, gy, isBorder] of walls) {
        const wx = (gx + 0.5) * cellSize;
        const wy = (gy + 0.5) * cellSize;
        const halfSize = cellSize / 2;

        graphics.rect(wx - halfSize, wy - halfSize, cellSize, cellSize);
        // Border walls are cyan (0x00bcd4), user-placed walls are gray (0x4a5568)
        const wallColor = isBorder ? 0x00bcd4 : 0x4a5568;
        graphics.fill({ color: wallColor, alpha: 1 });
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

      // Draw obstacle (interactive ball) - solid, not transparent
      const u = universe as any;
      if (u.get_show_obstacle && u.get_show_obstacle()) {
        const gridDims = universe.get_grid_dimensions();
        const gridWidth = gridDims[0];
        const gridHeight = gridDims[1];

        // Obstacle position is in normalized coordinates (0-1), convert to pixels
        // Offset by 1 cell to align visual circle with physics grid cells
        const obstacleX = u.get_obstacle_x() * gridWidth - cellSize;
        const obstacleY = u.get_obstacle_y() * gridHeight - cellSize;
        const obstacleRadius = u.get_obstacle_radius() * gridWidth;

        // Draw solid obstacle circle
        graphics.circle(obstacleX, obstacleY, obstacleRadius);
        graphics.fill({ color: 0x555555, alpha: 1.0 }); // Solid, not transparent

        // Add a highlight/outline for visibility
        graphics.setStrokeStyle({
          width: 3,
          color: 0x888888,
          alpha: 1,
        });
        graphics.circle(obstacleX, obstacleY, obstacleRadius);
        graphics.stroke();

        // Inner highlight for 3D effect
        graphics.circle(
          obstacleX - obstacleRadius * 0.25,
          obstacleY + obstacleRadius * 0.25,
          obstacleRadius * 0.25
        );
        graphics.fill({ color: 0x999999, alpha: 0.8 });
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
