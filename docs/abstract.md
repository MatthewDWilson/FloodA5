# FloodA5: A Pentagonal Discrete Global Grid System for 2D Surface Water Modelling

## Abstract

Two-dimensional flood inundation models discretise the landscape into computational cells and route water between them. The choice of grid geometry has direct consequences for numerical accuracy, computational efficiency, and ease of use. Most operational models use structured rectangular grids, which introduce directional bias in flow routing and produce artefacts at mesh boundaries when the domain is irregular. Unstructured triangular or quadrilateral grids avoid these biases but require complex mesh generation and neighbourhood traversal. We present FloodA5, a 2D surface water model built on the A5 pentagonal Discrete Global Grid System (DGGS), and describe its key architectural advantages for urban and catchment-scale flood simulation.

### The A5 Grid

The A5 DGGS tessellates the sphere with pentagons arranged in a hierarchical, multi-resolution structure. At any given resolution level each cell has exactly five edge-sharing neighbours, a property that eliminates the directional bias inherent in rectilinear grids while maintaining the topological simplicity of a regular lattice. Cells at adjacent resolution levels share a strict parent–child relationship, allowing adaptive refinement without re-meshing. For flood modelling at urban scales, resolution level 14 yields cells of approximately 12.6 ha (0.126 km²) — comparable to the spatial scale at which hydrological connectivity is resolved by 1 m LiDAR — while level 17 provides ~1,976 m² cells suitable for detailed drainage analysis.

Cell geometry is intrinsically geographic: boundaries and centres are expressed as geodetic coordinates, and all geometric operations (area, edge length, adjacency) are computed on the sphere using haversine distances. This avoids the distortion introduced by projecting large domains onto a plane, and means a single model configuration can be applied at any location on Earth without coordinate system selection.

Mesh generation over an arbitrary area of interest is fully automated: the user provides a GeoJSON polygon, and the A5 DGGS library (pya5) returns all cells whose centres fall within the domain, together with their boundaries, centre coordinates, and topological adjacency via the `grid_disk` function. For the 175 km² Christchurch test domain at resolution 14, this produces 2953 cells in under 10 seconds.

### Sub-Grid Sampling

A key limitation of storage-cell 2D models is the assumption that each cell has a single representative bed elevation. In practice, a 12.6 ha pentagon (resolution 14) contains substantial internal topographic variability — ridges, channels, roads, and built structures — that controls whether and how water accumulates and routes. FloodA5 addresses this through Sub-Grid Sampling (SGS), which replaces the single cell elevation with a hypsometric curve: a pre-computed lookup table relating water surface elevation (WSE) to stored volume and wetted plan area.

For each cell, 512 quasi-random sample points are distributed within the cell polygon using a Halton low-discrepancy sequence and filtered through a point-in-polygon test. The 1 m LiDAR DEM is sampled at each point using bilinear interpolation; the resulting elevation distribution is binned into 100 quantile-spaced levels from which cumulative volume and area curves are derived. These tables are computed once as a pre-processing step and stored alongside the mesh in the GeoParquet file as Arrow list columns, so subsequent model runs load them directly without re-processing the DEM.

At each timestep, the current stored volume is converted to WSE via inverse interpolation of the volume curve. Flow between adjacent cells is governed by Manning's equation applied across the shared edge, with the sill elevation set to the minimum terrain elevation along the shared boundary — itself sampled from the DEM during pre-processing. This means that water in a cell fills the topographic lows first, routing through sub-grid channels even when those channels are narrower than the cell diameter. Cells transition smoothly from dry to wet rather than switching instantaneously, eliminating the numerical oscillation at the wetting front that affects flat-bottomed storage models.

### Flow Routing

The diffusive wave approximation is used to compute inter-cell fluxes. The volumetric flow rate across a shared edge of width $W$ is:

$$Q = \frac{A_\text{wet}}{n} R^{2/3} S^{1/2}$$

where $A_\text{wet}$ is the cross-sectional wetted area at the sill, $n$ is Manning's roughness, $R$ is the hydraulic radius (approximated as the depth above the sill for a wide channel), and $S$ is the water surface slope. Manning's $n$ can be specified as a global constant or sampled from a per-cell friction raster. The timestep is determined adaptively by the diffusive-wave CFL stability criterion, with a user-specified upper bound.

### Implementation

FloodA5 is implemented in Julia, leveraging multithreaded CPU parallelism and optional GPU acceleration (CUDA.jl) for the polygon-in-polygon sampling step. DEM ingestion uses ArchGDAL.jl, supporting any GDAL-readable raster format and handling CRS reprojection automatically. Simulation output is written to HDF5, with per-frame datasets for depth, stored volume, scalar velocity, and the saturation fraction (wetted area / total cell area) — the last being a physically meaningful output variable that has no equivalent in conventional flat-bottomed models.

The complete pipeline — mesh generation, DEM ingestion, SGS pre-processing, and a one-hour simulation on the 2953-cell Christchurch domain — runs in under 15 seconds on a desktop workstation with a consumer GPU.

### Conclusions

The A5 pentagonal DGGS offers a principled foundation for 2D surface water modelling: uniform five-connectivity eliminates directional bias, the hierarchical resolution structure supports multi-scale analysis, and the geographic coordinate system avoids projection distortion. Combined with the SGS pre-processing approach, which captures sub-cell topographic variability from high-resolution LiDAR without requiring fine-grid resolution, FloodA5 achieves physically realistic wetting-front behaviour at modest computational cost. Future work will focus on validation against observed flood events, implementation of upstream inflow boundary conditions, and extension to spatially variable rainfall forcing.
