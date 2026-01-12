import { X, GripVertical } from "lucide-react";
import { useState, useRef, useEffect } from "react";
import { useSimulation } from "../contexts/SimulationContext";

export default function UniversePropertyEditor() {
  const {
    universe,
    isUniverseEditorOpen,
    setIsUniverseEditorOpen,
    setRender,
  } = useSimulation();

  const [position, setPosition] = useState({
    x: window.innerWidth - 420,
    y: 10,
  });
  const [isDragging, setIsDragging] = useState(false);
  const [dragOffset, setDragOffset] = useState({ x: 0, y: 0 });
  const [isEditing, setIsEditing] = useState<string | null>(null);
  const editorRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);

  // Local state for input values
  const [gravity, setGravity] = useState(universe.get_gravity().toFixed(2));
  const [windX, setWindX] = useState(universe.get_wind_x());
  const [windY, setWindY] = useState(universe.get_wind_y());
  const [windXValue, setWindXValue] = useState(universe.get_wind_x().toFixed(2));
  const [windYValue, setWindYValue] = useState(universe.get_wind_y().toFixed(2));
  const [speed, setSpeed] = useState(universe.get_speed().toFixed(2));
  const [flipRatio, setFlipRatio] = useState(universe.get_flip_ratio().toFixed(2));
  const [particleRadius, setParticleRadius] = useState(universe.get_particle_radius().toFixed(1));
  const [cellSize, setCellSize] = useState(universe.get_cell_size().toFixed(0));

  // Canvas interaction state
  const [isDraggingWind, setIsDraggingWind] = useState(false);

  // Update local values when universe changes
  useEffect(() => {
    if (!isEditing) {
      setGravity(universe.get_gravity().toFixed(2));
      setWindXValue(universe.get_wind_x().toFixed(2));
      setWindYValue(universe.get_wind_y().toFixed(2));
      setSpeed(universe.get_speed().toFixed(2));
      setFlipRatio(universe.get_flip_ratio().toFixed(2));
      setParticleRadius(universe.get_particle_radius().toFixed(1));
      setCellSize(universe.get_cell_size().toFixed(0));
    }
    setWindX(universe.get_wind_x());
    setWindY(universe.get_wind_y());
  }, [universe, isEditing]);

  // Draw wind vector canvas
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    const centerX = canvas.width / 2;
    const centerY = canvas.height / 2;
    const scale = 0.5; // pixels per unit

    // Clear canvas
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Draw background
    ctx.fillStyle = "#f9fafb";
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // Draw grid
    ctx.strokeStyle = "#e5e7eb";
    ctx.lineWidth = 1;
    for (let i = 0; i <= canvas.width; i += 20) {
      ctx.beginPath();
      ctx.moveTo(i, 0);
      ctx.lineTo(i, canvas.height);
      ctx.stroke();
    }
    for (let i = 0; i <= canvas.height; i += 20) {
      ctx.beginPath();
      ctx.moveTo(0, i);
      ctx.lineTo(canvas.width, i);
      ctx.stroke();
    }

    // Draw center point
    ctx.fillStyle = "#3b82f6";
    ctx.beginPath();
    ctx.arc(centerX, centerY, 4, 0, Math.PI * 2);
    ctx.fill();

    // Draw wind vector
    const windVectorX = windX * scale;
    const windVectorY = -windY * scale; // Flip Y for canvas coordinates

    if (Math.abs(windVectorX) > 0.1 || Math.abs(windVectorY) > 0.1) {
      // Draw line
      ctx.strokeStyle = "#3b82f6";
      ctx.lineWidth = 3;
      ctx.beginPath();
      ctx.moveTo(centerX, centerY);
      ctx.lineTo(centerX + windVectorX, centerY + windVectorY);
      ctx.stroke();

      // Draw arrowhead
      const angle = Math.atan2(windVectorY, windVectorX);
      const arrowLength = 10;
      const arrowAngle = Math.PI / 6;

      ctx.fillStyle = "#3b82f6";
      ctx.beginPath();
      ctx.moveTo(centerX + windVectorX, centerY + windVectorY);
      ctx.lineTo(
        centerX + windVectorX - arrowLength * Math.cos(angle - arrowAngle),
        centerY + windVectorY - arrowLength * Math.sin(angle - arrowAngle)
      );
      ctx.lineTo(
        centerX + windVectorX - arrowLength * Math.cos(angle + arrowAngle),
        centerY + windVectorY - arrowLength * Math.sin(angle + arrowAngle)
      );
      ctx.closePath();
      ctx.fill();
    }
  }, [windX, windY]);

  const handleCanvasMouseDown = (e: React.MouseEvent<HTMLCanvasElement>) => {
    setIsDraggingWind(true);
    updateWindFromMouse(e);
  };

  const handleCanvasMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (isDraggingWind) {
      updateWindFromMouse(e);
    }
  };

  const handleCanvasMouseUp = () => {
    setIsDraggingWind(false);
  };

  const updateWindFromMouse = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;

    const centerX = canvas.width / 2;
    const centerY = canvas.height / 2;
    const scale = 0.5;

    const newWindX = (x - centerX) / scale;
    const newWindY = -(y - centerY) / scale; // Flip Y back

    setWindX(newWindX);
    setWindY(newWindY);
    setWindXValue(newWindX.toFixed(2));
    setWindYValue(newWindY.toFixed(2));
    universe.set_wind_x(newWindX);
    universe.set_wind_y(newWindY);
    setRender((prev) => prev + 1);
  };

  useEffect(() => {
    const handleMouseMove = (e: MouseEvent) => {
      if (isDragging) {
        setPosition({
          x: e.clientX - dragOffset.x,
          y: e.clientY - dragOffset.y,
        });
      }
    };

    const handleMouseUp = () => {
      setIsDragging(false);
    };

    if (isDragging) {
      document.addEventListener("mousemove", handleMouseMove);
      document.addEventListener("mouseup", handleMouseUp);
    }

    return () => {
      document.removeEventListener("mousemove", handleMouseMove);
      document.removeEventListener("mouseup", handleMouseUp);
    };
  }, [isDragging, dragOffset]);

  const handleDragStart = (e: React.MouseEvent) => {
    if (editorRef.current) {
      const rect = editorRef.current.getBoundingClientRect();
      setDragOffset({
        x: e.clientX - rect.left,
        y: e.clientY - rect.top,
      });
      setIsDragging(true);
    }
  };

  const handlePropertyUpdate = (property: string, value: number) => {
    switch (property) {
      case "gravity":
        universe.set_gravity(value);
        setGravity(value.toFixed(2));
        break;
      case "windX":
        universe.set_wind_x(value);
        setWindX(value);
        setWindXValue(value.toFixed(2));
        break;
      case "windY":
        universe.set_wind_y(value);
        setWindY(value);
        setWindYValue(value.toFixed(2));
        break;
      case "speed":
        universe.set_speed(value);
        setSpeed(value.toFixed(2));
        break;
      case "flipRatio":
        universe.set_flip_ratio(value);
        setFlipRatio(value.toFixed(2));
        break;
      case "particleRadius":
        universe.set_particle_radius(value);
        setParticleRadius(value.toFixed(1));
        break;
      case "cellSize":
        universe.set_cell_size(value);
        setCellSize(value.toFixed(0));
        break;
    }
    setRender((prev) => prev + 1);
  };

  if (!isUniverseEditorOpen) return null;

  return (
    <div
      ref={editorRef}
      className="fixed bg-white border border-gray-200 rounded-lg shadow-2xl pointer-events-auto max-h-[90vh] overflow-y-auto"
      style={{
        left:
          typeof window !== "undefined" && window.innerWidth < 640
            ? "50%"
            : `${position.x}px`,
        top:
          typeof window !== "undefined" && window.innerWidth < 640
            ? "50%"
            : `${position.y}px`,
        transform:
          typeof window !== "undefined" && window.innerWidth < 640
            ? "translate(-50%, -50%)"
            : "none",
        width:
          typeof window !== "undefined" && window.innerWidth < 640
            ? "90vw"
            : "400px",
        maxWidth: "400px",
        zIndex: 1000,
      }}
    >
      {/* Header */}
      <div className="flex sticky top-0 items-center justify-between p-2 sm:p-3 border-b border-gray-200 bg-gray-50 rounded-t-lg">
        <div
          className="flex items-center gap-2 cursor-move flex-1"
          onMouseDown={handleDragStart}
        >
          <GripVertical className="w-4 h-4 sm:w-5 sm:h-5 text-gray-400 hidden sm:block" />
          <GripVertical className="w-4 h-4 sm:w-5 sm:h-5 text-gray-400 -ml-4 hidden sm:block" />
          <h2 className="text-base sm:text-lg font-semibold text-gray-800 sm:ml-2">
            Fluid Settings
          </h2>
        </div>
        <button
          onClick={() => setIsUniverseEditorOpen(false)}
          className="p-1 hover:bg-gray-200 cursor-pointer rounded transition-all duration-200"
        >
          <X className="w-5 h-5 text-gray-600" />
        </button>
      </div>

      {/* Property fields */}
      <div className="p-3 sm:p-4 space-y-3 sm:space-y-4">
        {/* Wind Vector Canvas */}
        <div className="space-y-1">
          <label className="text-sm font-medium text-gray-700">
            Wind Vector:
          </label>
          <canvas
            ref={canvasRef}
            width={352}
            height={200}
            className="border border-gray-300 rounded cursor-crosshair"
            onMouseDown={handleCanvasMouseDown}
            onMouseMove={handleCanvasMouseMove}
            onMouseUp={handleCanvasMouseUp}
            onMouseLeave={handleCanvasMouseUp}
          />
          <p className="text-xs text-gray-500">
            Click and drag to set wind direction and speed
          </p>
        </div>

        {/* Wind X */}
        <div className="space-y-1">
          <label className="text-sm font-medium text-gray-700">
            Wind X:
          </label>
          <input
            type="text"
            value={windXValue}
            onChange={(e) => {
              setIsEditing("windX");
              setWindXValue(e.target.value);
              const val = parseFloat(e.target.value);
              if (!isNaN(val)) {
                handlePropertyUpdate("windX", val);
              }
            }}
            onBlur={() => setIsEditing(null)}
            className="w-full px-3 py-2 text-base border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </div>

        {/* Wind Y */}
        <div className="space-y-1">
          <label className="text-sm font-medium text-gray-700">
            Wind Y:
          </label>
          <input
            type="text"
            value={windYValue}
            onChange={(e) => {
              setIsEditing("windY");
              setWindYValue(e.target.value);
              const val = parseFloat(e.target.value);
              if (!isNaN(val)) {
                handlePropertyUpdate("windY", val);
              }
            }}
            onBlur={() => setIsEditing(null)}
            className="w-full px-3 py-2 text-base border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </div>

        {/* Gravity */}
        <div className="space-y-1">
          <label className="text-sm font-medium text-gray-700">
            Gravity:
          </label>
          <input
            type="text"
            value={gravity}
            onChange={(e) => {
              setIsEditing("gravity");
              setGravity(e.target.value);
              const val = parseFloat(e.target.value);
              if (!isNaN(val)) {
                handlePropertyUpdate("gravity", val);
              }
            }}
            onBlur={() => setIsEditing(null)}
            className="w-full px-3 py-2 text-base border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
          <span className="text-xs text-gray-500">Default: -200 (downward)</span>
        </div>

        {/* Simulation Speed */}
        <div className="space-y-1">
          <label className="text-sm font-medium text-gray-700">
            Simulation Speed:
          </label>
          <input
            type="text"
            value={speed}
            onChange={(e) => {
              setIsEditing("speed");
              setSpeed(e.target.value);
              const val = parseFloat(e.target.value);
              if (!isNaN(val) && val > 0) {
                handlePropertyUpdate("speed", val);
              }
            }}
            onBlur={() => setIsEditing(null)}
            className="w-full px-3 py-2 text-base border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </div>

        {/* FLIP Ratio */}
        <div className="space-y-1">
          <label className="text-sm font-medium text-gray-700">
            FLIP Ratio (0 = PIC, 1 = FLIP):
          </label>
          <input
            type="range"
            min="0"
            max="1"
            step="0.05"
            value={flipRatio}
            onChange={(e) => {
              const val = parseFloat(e.target.value);
              handlePropertyUpdate("flipRatio", val);
            }}
            className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer"
          />
          <div className="flex justify-between text-xs text-gray-500">
            <span>PIC (Smooth)</span>
            <span>{flipRatio}</span>
            <span>FLIP (Splashy)</span>
          </div>
        </div>

        {/* Particle Radius */}
        <div className="space-y-1">
          <label className="text-sm font-medium text-gray-700">
            Particle Radius:
          </label>
          <input
            type="text"
            value={particleRadius}
            onChange={(e) => {
              setIsEditing("particleRadius");
              setParticleRadius(e.target.value);
              const val = parseFloat(e.target.value);
              if (!isNaN(val) && val > 0) {
                handlePropertyUpdate("particleRadius", val);
              }
            }}
            onBlur={() => setIsEditing(null)}
            className="w-full px-3 py-2 text-base border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
        </div>

        {/* Cell Size */}
        <div className="space-y-1">
          <label className="text-sm font-medium text-gray-700">
            Grid Cell Size:
          </label>
          <input
            type="text"
            value={cellSize}
            onChange={(e) => {
              setIsEditing("cellSize");
              setCellSize(e.target.value);
              const val = parseFloat(e.target.value);
              if (!isNaN(val) && val > 0) {
                handlePropertyUpdate("cellSize", val);
              }
            }}
            onBlur={() => setIsEditing(null)}
            className="w-full px-3 py-2 text-base border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
          />
          <span className="text-xs text-gray-500">Affects simulation resolution</span>
        </div>
      </div>
    </div>
  );
}
