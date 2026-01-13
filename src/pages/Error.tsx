import { Link } from "react-router-dom";
import { Application } from "@pixi/react";
import SandBox from "../components/SandBox";
import { useWindowDimension } from "../utils/useWindowDimension";
import Transitions from "../utils/Transitions";
import { Universe } from "physics-engine";
import { useEffect, useRef } from "react";
import { Viewport } from "../utils/Viewport";
import { SimulationProvider } from "../contexts/SimulationContext";

function ErrorContent({ universe }: { universe: Universe }) {
  const [width, height] = useWindowDimension();
  
  // Get grid dimensions for centering
  const gridDims = universe.get_grid_dimensions();
  const gridWidth = gridDims[0];
  const gridHeight = gridDims[1];

  return (
    <Transitions>
      <Application
        background={"#ffffff"}
        width={width}
        height={height}
        className="overflow-hidden fixed top-0 left-0"
        antialias
      >
        <Viewport centerOnGrid={{ width: gridWidth, height: gridHeight }}>
          <SandBox universe={universe} />
        </Viewport>
      </Application>
      <div className="fixed -translate-x-1/2 -translate-y-1/3 left-1/2 top-1/3 z-10 text-center">
        <h1 className="text-6xl md:text-9xl font-bold mb-4 text-red-600">
          404
        </h1>
        <p className="text-2xl md:text-4xl font-semibold mb-12 text-gray-900">
          Page Not Found
        </p>
        <div className="flex flex-col gap-3 mx-auto justify-center">
          <Link
            to="/"
            className="no-underline w-1/2 mx-auto bg-black text-white text-xl lg:text-3xl font-bold px-8 py-3 rounded-lg shadow-lg text-center transition-transform hover:scale-105"
          >
            Home
          </Link>
        </div>
      </div>
    </Transitions>
  );
}

export default function Error() {
  const universeRef = useRef<Universe>(new Universe());

  useEffect(() => {
    // Initialize universe with Euler wind tunnel and random walls
    const universe = universeRef.current;
    
    // Setup blank canvas to ensure Euler mode
    universe.setup_blank(16.0, 9.0);
    
    // Hide the obstacle/ball on the error screen
    (universe as any).set_show_obstacle(false);
    
    // Get grid info for placing random walls
    const cellSize = universe.get_cell_size();
    const gridDims = universe.get_grid_dimensions();
    const gridWidth = gridDims[0];
    const gridHeight = gridDims[1];
    
    // Add random walls around the simulation area
    const numWalls = 10 + Math.floor(Math.random() * 6); // 10-15 random walls for error page
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
      <ErrorContent universe={universeRef.current} />
    </SimulationProvider>
  );
}
