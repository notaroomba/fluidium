import { Viewport as BaseViewport, type IViewportOptions } from "pixi-viewport";
import { extend, useApplication } from "@pixi/react";
import { type PropsWithChildren, useEffect, useRef } from "react";
import { Application } from "pixi.js";

type ViewportProps = Omit<IViewportOptions, "events"> & {
  centerOnGrid?: { width: number; height: number };
};

class ViewportWrapper extends BaseViewport {
  constructor(
    options: Omit<ViewportProps, "centerOnGrid"> & { app: Application }
  ) {
    const { app, ...rest } = options;
    super({
      ...rest,
      // events is the only required argument to the constructor.
      events: app.renderer.events,
    });
    // Only enable drag/pan on simulation page, use middle mouse button
    const url = new URL(window.location.href);
    if (url.pathname !== "/") {
      this.drag({ mouseButtons: "middle" }).pinch().wheel().decelerate();
    }
  }
}

extend({ ViewportWrapper });

function Viewport(props: PropsWithChildren<ViewportProps>) {
  const { children, centerOnGrid, ...rest } = props;
  const { app } = useApplication();
  const viewportRef = useRef<BaseViewport | null>(null);

  // Center viewport on grid when centerOnGrid changes or on mount
  useEffect(() => {
    if (viewportRef.current && centerOnGrid && app) {
      const vp = viewportRef.current;
      // Move viewport to center the grid in the screen
      // The grid's center is at (width/2, height/2)
      // We need to position the viewport so this point appears at the screen center
      vp.moveCenter(centerOnGrid.width / 2, -centerOnGrid.height / 2);
    }
  }, [centerOnGrid, app]);

  return (
    app?.renderer && (
      <pixiViewportWrapper
        ref={(ref: BaseViewport | null) => {
          viewportRef.current = ref;
          // Store viewport reference globally
          (window as any).pixiViewport = ref;
          // Center immediately when ref is set
          if (ref && centerOnGrid) {
            ref.moveCenter(centerOnGrid.width / 2, -centerOnGrid.height / 2);
          }
        }}
        position={{ x: app.canvas.width / 2, y: -app.canvas.height / 2 }}
        app={app}
        {...rest}
      >
        {children}
      </pixiViewportWrapper>
    )
  );
}

export { Viewport };
