import {
  Rewind,
  Pause,
  Play,
  FastForward,
  RotateCcw,
  Calculator,
  Waves,
  Maximize2,
  Settings,
  Grid3x3,
  Navigation,
} from "lucide-react";
import { Implementation } from "physics-engine";
import { useEffect, useState } from "react";
import { useSimulation } from "../contexts/SimulationContext";

export default function SettingsBar() {
  const {
    universe,
    setRender,
    isUniverseEditorOpen,
    setIsUniverseEditorOpen,
    showGrid,
    setShowGrid,
    showVelocity,
    setShowVelocity,
    isPaused,
    setIsPaused,
  } = useSimulation();
  
  const [implementation, setImplementation] = useState<Implementation>(
    universe.get_implementation()
  );

  // Speed multipliers
  const multipliers = [0.25, 0.5, 1, 2, 4];
  const [multiplier, setMultiplier] = useState(() => {
    const current = universe.get_speed();
    return current;
  });

  useEffect(() => {
    universe.set_implementation(implementation);
    console.log("Implementation set to", implementation);
  }, [implementation]);

  useEffect(() => {
    universe.set_is_paused(isPaused);
    setRender((prev) => prev + 1);
  }, [isPaused]);

  const rewind = () => {
    const currentIndex = multipliers.indexOf(multiplier);
    const newIndex = Math.max(0, currentIndex - 1);
    const newMultiplier = multipliers[newIndex];
    setMultiplier(newMultiplier);
    universe.set_speed(newMultiplier);
    setRender((prev) => prev + 1);
  };

  const fastForward = () => {
    const currentIndex = multipliers.indexOf(multiplier);
    const newIndex = Math.min(multipliers.length - 1, currentIndex + 1);
    const newMultiplier = multipliers[newIndex];
    setMultiplier(newMultiplier);
    universe.set_speed(newMultiplier);
    setRender((prev) => prev + 1);
  };

  const reset = () => {
    universe.reset();
    setMultiplier(1);
    universe.set_speed(1.0);
    setRender((prev) => prev + 1);
  };

  const resetView = () => {
    // Access the viewport through the stage children
    const app = (window as any).pixiApp;
    if (app && app.stage && app.stage.children) {
      const viewport = app.stage.children.find(
        (child: any) => child.constructor.name === "ViewportWrapper"
      );
      if (viewport) {
        viewport.position.set(app.canvas.width / 2, app.canvas.height / 2);
        viewport.scale.set(1, 1);
        viewport.rotation = 0;
      }
    }
  };

  return (
    <div className="flex flex-col items-center gap-2">
      {/* Settings Bar */}
      <div className="w-fit max-w-full py-2 px-2 sm:px-4 bg-white border pointer-events-auto border-gray-200 rounded-lg shadow-xl overflow-x-auto">
        <div className="flex items-center justify-center divide-x divide-gray-300 flex-wrap sm:flex-nowrap gap-y-2">
          {/* Integration Method */}
          <div className="flex items-center gap-1 sm:gap-2 px-2 sm:px-4">
            <button
              onClick={() => setImplementation(Implementation.Euler)}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                implementation === Implementation.Euler
                  ? "bg-blue-100"
                  : "hover:bg-gray-100"
              }`}
              title="Euler Method (Simple)"
            >
              <Calculator className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
            <button
              onClick={() => setImplementation(Implementation.FLIP)}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                implementation === Implementation.FLIP
                  ? "bg-blue-100"
                  : "hover:bg-gray-100"
              }`}
              title="FLIP Method (More Realistic)"
            >
              <Waves className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
          </div>

          {/* Visualization Options */}
          <div className="flex items-center gap-1 sm:gap-2 px-2 sm:px-4">
            <button
              onClick={() => {
                const newValue = !showGrid;
                setShowGrid(newValue);
                universe.set_show_grid(newValue);
                setRender((prev) => prev + 1);
              }}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                showGrid ? "bg-blue-100" : "hover:bg-gray-100"
              }`}
              title={showGrid ? "Hide Grid" : "Show Grid"}
            >
              <Grid3x3 className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
            <button
              onClick={() => {
                const newValue = !showVelocity;
                setShowVelocity(newValue);
                universe.set_show_velocity(newValue);
                setRender((prev) => prev + 1);
              }}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                showVelocity ? "bg-blue-100" : "hover:bg-gray-100"
              }`}
              title={
                showVelocity
                  ? "Hide Velocity Vectors"
                  : "Show Velocity Vectors"
              }
            >
              <Navigation className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
          </div>

          {/* Playback Controls */}
          <div className="flex items-center gap-1 sm:gap-2 px-2 sm:px-4">
            <div className="w-8 sm:w-12 text-center text-xs sm:text-base">
              <p>{multiplier}x</p>
            </div>
            <button
              onClick={rewind}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 active:bg-blue-200 hover:bg-gray-100`}
              title="Slower"
            >
              <Rewind className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
            <button
              onClick={() => setIsPaused(!isPaused)}
              className="p-1.5 sm:p-2 hover:bg-gray-100 rounded cursor-pointer transition-all duration-200"
              title={isPaused ? "Play" : "Pause"}
            >
              {isPaused ? (
                <Play className="w-4 h-4 sm:w-5 sm:h-5" />
              ) : (
                <Pause className="w-4 h-4 sm:w-5 sm:h-5" />
              )}
            </button>
            <button
              onClick={fastForward}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 active:bg-blue-200 hover:bg-gray-100`}
              title="Faster"
            >
              <FastForward className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
            <button
              onClick={reset}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 active:bg-blue-200 hover:bg-gray-100`}
              title="Reset Simulation"
            >
              <RotateCcw className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
          </div>

          {/* Settings */}
          <div className="flex items-center gap-1 sm:gap-2 px-2 sm:px-4">
            <button
              onClick={() => {
                setIsUniverseEditorOpen(!isUniverseEditorOpen);
              }}
              className={`p-1.5 sm:p-2 hover:bg-gray-100 rounded cursor-pointer transition-all duration-200 ${
                isUniverseEditorOpen ? "bg-blue-100" : ""
              }`}
              title={
                isUniverseEditorOpen
                  ? "Close Settings"
                  : "Open Settings"
              }
            >
              <Settings className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
            <button
              onClick={resetView}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 active:bg-blue-200 hover:bg-gray-100`}
              title="Reset View"
            >
              <Maximize2 className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
          </div>
        </div>
      </div>
    </div>
  );
}
