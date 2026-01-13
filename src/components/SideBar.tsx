import { useSimulation } from "../contexts/SimulationContext";
import type { ToolMode, BrushType } from "../contexts/SimulationContext";
import { useState, useEffect } from "react";
import {
  Pencil,
  Eraser,
  MousePointer2,
  Droplets,
  Square,
  Wind,
  Minus,
  Plus,
  Circle,
} from "lucide-react";

export default function SideBar() {
  const {
    universe,
    setRender,
    toolMode,
    setToolMode,
    brushType,
    setBrushType,
    brushSize,
    setBrushSize,
  } = useSimulation();

  // Ball (obstacle) state
  const [showBall, setShowBall] = useState(() => {
    const u = universe as any;
    return u.get_show_obstacle ? u.get_show_obstacle() : false;
  });
  const [ballSize, setBallSize] = useState(() => {
    const u = universe as any;
    return u.get_obstacle_radius
      ? Math.round(u.get_obstacle_radius() * 100)
      : 15;
  });

  // Sync ball state with universe
  useEffect(() => {
    const u = universe as any;
    if (u.set_show_obstacle) {
      u.set_show_obstacle(showBall);
      setRender((prev) => prev + 1);
    }
  }, [showBall]);

  useEffect(() => {
    const u = universe as any;
    if (u.set_obstacle_radius) {
      u.set_obstacle_radius(ballSize / 100);
      setRender((prev) => prev + 1);
    }
  }, [ballSize]);

  const handleToolChange = (mode: ToolMode) => {
    setToolMode(mode);
  };

  const handleBrushTypeChange = (type: BrushType) => {
    setBrushType(type);
  };

  const adjustBrushSize = (delta: number) => {
    const newSize = Math.max(10, Math.min(200, brushSize + delta));
    setBrushSize(newSize);
  };

  return (
    <div className="w-fit max-w-full py-2 px-2 sm:px-4 bg-white border pointer-events-auto border-gray-200 rounded-lg shadow-xl overflow-x-auto">
      <div className="flex flex-col items-center justify-center divide-y divide-gray-300 gap-y-2">
        {/* Tool Mode Selection */}
        <div className="flex flex-col items-center gap-2 py-2">
          <span className="text-xs text-gray-500 font-medium">Tool</span>
          <div className="flex items-center gap-1">
            <button
              onClick={() => handleToolChange("none")}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                toolMode === "none" ? "bg-blue-100" : "hover:bg-gray-100"
              }`}
              title="Move/Pan (No drawing)"
            >
              <MousePointer2 className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
            <button
              onClick={() => handleToolChange("pencil")}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                toolMode === "pencil" ? "bg-blue-100" : "hover:bg-gray-100"
              }`}
              title="Pencil (Draw)"
            >
              <Pencil className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
            <button
              onClick={() => handleToolChange("eraser")}
              className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                toolMode === "eraser" ? "bg-blue-100" : "hover:bg-gray-100"
              }`}
              title="Eraser"
            >
              <Eraser className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
          </div>
        </div>

        {/* Brush Type Selection (only when pencil is selected) */}
        {toolMode === "pencil" && (
          <div className="flex flex-col items-center gap-2 py-2">
            <span className="text-xs text-gray-500 font-medium">Brush</span>
            <div className="flex items-center gap-1">
              <button
                onClick={() => handleBrushTypeChange("water")}
                className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                  brushType === "water" ? "bg-blue-100" : "hover:bg-gray-100"
                }`}
                title="Water"
              >
                <Droplets className="w-4 h-4 sm:w-5 sm:h-5 text-blue-500" />
              </button>
              <button
                onClick={() => handleBrushTypeChange("wall")}
                className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                  brushType === "wall" ? "bg-blue-100" : "hover:bg-gray-100"
                }`}
                title="Wall"
              >
                <Square className="w-4 h-4 sm:w-5 sm:h-5 text-gray-700" />
              </button>
              <button
                onClick={() => handleBrushTypeChange("air")}
                className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
                  brushType === "air" ? "bg-blue-100" : "hover:bg-gray-100"
                }`}
                title="Air (Erase fluid only)"
              >
                <Wind className="w-4 h-4 sm:w-5 sm:h-5 text-gray-400" />
              </button>
            </div>
          </div>
        )}

        {/* Brush Size */}
        {(toolMode === "pencil" || toolMode === "eraser") && (
          <div className="flex flex-col items-center gap-2 py-2">
            <span className="text-xs text-gray-500 font-medium">Size</span>
            <div className="flex items-center gap-1">
              <button
                onClick={() => adjustBrushSize(-10)}
                className="p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 hover:bg-gray-100"
                title="Decrease brush size"
              >
                <Minus className="w-4 h-4 sm:w-5 sm:h-5" />
              </button>
              <span className="w-8 text-center text-sm font-medium">
                {brushSize}
              </span>
              <button
                onClick={() => adjustBrushSize(10)}
                className="p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 hover:bg-gray-100"
                title="Increase brush size"
              >
                <Plus className="w-4 h-4 sm:w-5 sm:h-5" />
              </button>
            </div>
            {/* Brush size slider */}
            <input
              type="range"
              min="10"
              max="200"
              value={brushSize}
              onChange={(e) => setBrushSize(parseInt(e.target.value))}
              className="w-20 h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer"
              title="Brush size"
            />
          </div>
        )}

        {/* Quick Actions */}
        <div className="flex flex-col items-center gap-2 py-2">
          <span className="text-xs text-gray-500 font-medium">Clear</span>
          <div className="flex items-center gap-1">
            <button
              onClick={() => {
                universe.clear_particles();
                setRender((prev) => prev + 1);
              }}
              className="p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 hover:bg-red-100 text-red-500"
              title="Clear all water"
            >
              <Droplets className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
            <button
              onClick={() => {
                universe.clear_walls();
                setRender((prev) => prev + 1);
              }}
              className="p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 hover:bg-red-100 text-red-500"
              title="Clear all walls"
            >
              <Square className="w-4 h-4 sm:w-5 sm:h-5" />
            </button>
          </div>
        </div>

        {/* Ball/Obstacle Controls */}
        <div className="flex flex-col items-center gap-2 py-2">
          <span className="text-xs text-gray-500 font-medium">Ball</span>
          <button
            onClick={() => setShowBall(!showBall)}
            className={`p-1.5 sm:p-2 rounded cursor-pointer transition-all duration-200 ${
              showBall ? "bg-blue-100" : "hover:bg-gray-100"
            }`}
            title={showBall ? "Hide Ball" : "Show Ball"}
          >
            <Circle
              className={`w-4 h-4 sm:w-5 sm:h-5 ${
                showBall ? "text-blue-600" : "text-gray-500"
              }`}
            />
          </button>
          {showBall && (
            <div className="flex flex-col items-center gap-1">
              <span className="text-xs text-gray-400">Size</span>
              <div className="flex items-center gap-1">
                <button
                  onClick={() => setBallSize(Math.max(5, ballSize - 2))}
                  className="p-1 rounded cursor-pointer transition-all duration-200 hover:bg-gray-100"
                  title="Decrease ball size"
                >
                  <Minus className="w-3 h-3" />
                </button>
                <span className="w-6 text-center text-xs font-medium">
                  {ballSize}
                </span>
                <button
                  onClick={() => setBallSize(Math.min(40, ballSize + 2))}
                  className="p-1 rounded cursor-pointer transition-all duration-200 hover:bg-gray-100"
                  title="Increase ball size"
                >
                  <Plus className="w-3 h-3" />
                </button>
              </div>
              <input
                type="range"
                min="5"
                max="40"
                value={ballSize}
                onChange={(e) => setBallSize(parseInt(e.target.value))}
                className="w-16 h-1.5 bg-gray-200 rounded-lg appearance-none cursor-pointer"
                title="Ball size"
              />
            </div>
          )}
        </div>

        {/* Particle Count */}
        <div className="flex flex-col items-center gap-1 py-2">
          <span className="text-xs text-gray-500 font-medium">Particles</span>
          <span className="text-sm font-bold text-blue-600">
            {universe.get_particle_count()}
          </span>
        </div>
      </div>
    </div>
  );
}
