import { createContext, useContext, useEffect, useState } from "react";
import type { ReactNode } from "react";
import { Universe } from "physics-engine";

// Tool modes for drawing
export type ToolMode = "none" | "pencil" | "eraser";
export type BrushType = "water" | "wall" | "air";

interface SimulationContextType {
  universe: Universe;
  render: number;
  setRender: React.Dispatch<React.SetStateAction<number>>;
  isUniverseEditorOpen: boolean;
  setIsUniverseEditorOpen: (open: boolean) => void;
  isPaused: boolean;
  setIsPaused: (paused: boolean) => void;
  fps: number;
  setFps: React.Dispatch<React.SetStateAction<number>>;
  // Drawing tools
  toolMode: ToolMode;
  setToolMode: (mode: ToolMode) => void;
  brushType: BrushType;
  setBrushType: (type: BrushType) => void;
  brushSize: number;
  setBrushSize: (size: number) => void;
  // Visualization
  showGrid: boolean;
  setShowGrid: (show: boolean) => void;
  showVelocity: boolean;
  setShowVelocity: (show: boolean) => void;
}

const SimulationContext = createContext<SimulationContextType | undefined>(
  undefined
);

export function useSimulation() {
  const context = useContext(SimulationContext);
  if (!context) {
    throw new Error("useSimulation must be used within SimulationProvider");
  }
  return context;
}

interface SimulationProviderProps {
  children: ReactNode;
  universe: Universe;
}

export function SimulationProvider({
  children,
  universe,
}: SimulationProviderProps) {
  const [render, setRender] = useState(0);
  const [isUniverseEditorOpen, setIsUniverseEditorOpen] = useState(false);
  const [isPaused, setIsPaused] = useState(universe.get_is_paused());
  const [fps, setFps] = useState(0);

  // Drawing tools
  const [toolMode, setToolMode] = useState<ToolMode>("pencil");
  const [brushType, setBrushType] = useState<BrushType>("water");
  const [brushSize, setBrushSize] = useState(30);

  // Visualization
  const [showGrid, setShowGrid] = useState(false);
  const [showVelocity, setShowVelocity] = useState(false);

  useEffect(() => {
    universe.set_is_paused(isPaused);

    setRender((prev) => prev + 1);
  }, [isPaused]);

  return (
    <SimulationContext.Provider
      value={{
        universe,
        render,
        setRender,
        isUniverseEditorOpen,
        setIsUniverseEditorOpen,
        isPaused,
        setIsPaused,
        fps,
        setFps,
        toolMode,
        setToolMode,
        brushType,
        setBrushType,
        brushSize,
        setBrushSize,
        showGrid,
        setShowGrid,
        showVelocity,
        setShowVelocity,
      }}
    >
      {children}
    </SimulationContext.Provider>
  );
}
