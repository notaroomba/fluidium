import { Link } from "react-router-dom";
import { Application, extend } from "@pixi/react";
import SandBox from "../components/SandBox";
import { useWindowDimension } from "../utils/useWindowDimension";
import Transitions from "../utils/Transitions";
import { Universe } from "physics-engine";
import { useEffect, useRef, memo } from "react";
import { Viewport } from "../utils/Viewport";
import {
  SimulationProvider,
  useSimulation,
} from "../contexts/SimulationContext";
import { Container } from "pixi.js";

extend({ Container });

// Memoize the canvas to prevent re-renders from context changes
const PixiCanvas = memo(function PixiCanvas({
  universe,
}: {
  universe: Universe;
}) {
  const [width, height] = useWindowDimension();
  
  // Get grid dimensions for centering
  const gridDims = universe.get_grid_dimensions();
  const gridWidth = gridDims[0];
  const gridHeight = gridDims[1];

  return (
    <Application
      background={"#ffffff"}
      width={width}
      height={height}
      className="overflow-hidden fixed top-0 left-0"
      antialias
      onInit={(app) => {
        // Store app reference globally for viewport access
        (window as any).pixiApp = app;
      }}
    >
      <Viewport centerOnGrid={{ width: gridWidth, height: gridHeight }}>
        <SandBox universe={universe} />
      </Viewport>
    </Application>
  );
});

function HomeContent({ universe }: { universe: Universe }) {
  const { setShowGrid, setShowVelocity } = useSimulation();
  useEffect(() => {
    setShowGrid(false);
    setShowVelocity(false);
  }, [setShowGrid, setShowVelocity]);
  return (
    <Transitions>
      <PixiCanvas universe={universe} />
      <div className="fixed -translate-x-1/2 -translate-y-1/3 left-1/2 top-1/3 z-10">
        <h1 className="text-4xl md:text-6xl font-semibold mb-8 text-gray-900">
          Fluidium
        </h1>
        <div className="flex flex-col gap-3 mx-auto justify-center">
          <Link
            to="/simulation"
            className="no-underline w-full mx-auto bg-black text-white text-xl lg:text-3xl font-bold px-8 py-3 rounded-lg shadow-lg text-center transition-transform hover:scale-105"
          >
            Simulation
          </Link>
        </div>
      </div>
    </Transitions>
  );
}

export default function Home() {
  const universeRef = useRef<Universe>(new Universe());

  useEffect(() => {
    // Initialize universe with Euler wind tunnel
    const universe = universeRef.current;
    
    // Setup blank canvas first to ensure Euler mode
    universe.setup_blank(16.0, 9.0);
    
    // Hide the obstacle/ball on the home screen
    (universe as any).set_show_obstacle(false);
    
    // Get grid info for placing random walls
    const cellSize = universe.get_cell_size();
    const gridDims = universe.get_grid_dimensions();
    const gridWidth = gridDims[0];
    const gridHeight = gridDims[1];
    
    // Add random walls around the simulation area
    const numWalls = 8 + Math.floor(Math.random() * 5); // 8-12 random walls
    for (let i = 0; i < numWalls; i++) {
      // Place walls in random positions, avoiding the center and edges
      const margin = gridWidth * 0.15;
      const x = margin + Math.random() * (gridWidth - margin * 2);
      const y = margin + Math.random() * (gridHeight - margin * 2);
      
      // Random brush size for variety
      const brushSize = cellSize * (2 + Math.random() * 4);
      
      // Paint the wall
      universe.paint_wall(x, y, brushSize);
    }
  }, []);

  return (
    <SimulationProvider universe={universeRef.current}>
      <HomeContent universe={universeRef.current} />
    </SimulationProvider>
  );
}
