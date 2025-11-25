# Integrating Adaptive Composite Projections into Observable Versor Zooming Notebook

This guide shows how to add adaptive composite projection logic to the existing [Observable versor-zooming notebook](https://observablehq.com/@d3/versor-zooming).

## Overview

The adaptive composite projection automatically switches between different map projections based on:
- **Zoom level** (scale)
- **Center latitude** (location on the globe)

This provides optimal projection for each zoom level and location.

## Step-by-Step Integration

### Step 1: Add the Adaptive Projection Function

Create a new cell and add the `getAdaptiveProjection` function from `observable-adaptive-projection.js`. This function determines which projection to use based on zoom and center.

### Step 2: Add the Wrapper Projection

Add the `adaptiveProjection` function. This creates a wrapper projection that:
- Works with versor zooming
- Automatically switches projection types
- Maintains rotation/center state during switches

### Step 3: Replace the Zoom Function

Replace the existing `zoom` function with the modified version that supports adaptive projections. The key change is checking for `projection.updateAdaptive()` and handling projection switches.

### Step 4: Update the Projection Cell

Replace:
```javascript
projection = d3[projectionName]().precision(0.1)
```

With:
```javascript
projection = adaptiveProjection(width, 'naturalEarth')
projection.precision(0.1)
```

### Step 5: Update Height Calculation

The height calculation should work as-is, but you may need to call `projection.updateHeight(height)` after calculating height.

### Step 6: (Optional) Add World Projection Selector

Add a selector to choose between different world map projections (Natural Earth, Robinson, Mollweide, etc.) for the low-zoom levels.

## Projection Switching Rules

The adaptive system switches projections based on these rules:

1. **Mercator** (zoom ≥ 14): All latitudes, high zoom
2. **Albers Conic** (zoom 6-14, lat 20°-75°): Mid-latitudes, medium-high zoom
3. **Lambert Cylindrical** (zoom 6-14, lat 0°-20°): Low latitudes, medium-high zoom
4. **Lambert Azimuthal Polar** (zoom 2-13, lat ≥ 75°): Polar regions
5. **Lambert Azimuthal Oblique** (zoom 2-6, lat < 75°): Mid-zoom, non-polar
6. **Natural Earth/Robinson/etc.** (zoom 1-2.5, lat ≤ 68°): World view, low zoom

## Key Features

- **Seamless switching**: Projections switch automatically as you zoom/pan
- **Versor zooming compatible**: Works with the existing versor zooming implementation
- **State preservation**: Rotation and center are preserved during projection switches
- **Multiple world projections**: Choose from Natural Earth, Robinson, Mollweide, etc. for world view

## Testing

1. Start at world view - should show Natural Earth (or selected world projection)
2. Zoom in gradually - should switch to Lambert Azimuthal around zoom 2
3. Continue zooming - should switch to Albers/Mercator at higher zooms
4. Pan to polar regions - should use polar projections
5. Pan to equator - should use cylindrical projections

## Notes

- The projection switching happens automatically during zoom/pan
- The versor zooming continues to work smoothly
- Projection name is available via `projection.getProjectionName()`
- You can change the world projection via `projection.updateWorldProjection(name)`

