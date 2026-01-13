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
      <Viewport>
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
            className="no-underline w-1/2 mx-auto bg-black text-white text-xl lg:text-3xl font-bold px-8 py-3 rounded-lg shadow-lg text-center transition-transform hover:scale-105"
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
    // Initialize universe with a fluid scene
    const universe = universeRef.current;
    // Setup wind tunnel scene (scene 1) for the home page background
    universe.setup_scene(1, 16.0, 9.0);
  }, []);

  return (
    <SimulationProvider universe={universeRef.current}>
      <HomeContent universe={universeRef.current} />
    </SimulationProvider>
  );
}
