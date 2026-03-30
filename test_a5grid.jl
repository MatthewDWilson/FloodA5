"""
test_a5grid.jl
--------------
Basic tests for A5Grid.jl. Run with:
    julia test_a5grid.jl
"""

#push!(LOAD_PATH, joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "A5Grid.jl"))
using .A5Grid
using Test

@testset "A5Grid Tests" begin

    @testset "lon_lat_to_cell / cell_to_lon_lat round-trip" begin
        lon, lat = 151.2093, -33.8688   # Sydney
        res = 10
        cell_id = lon_lat_to_cell(lon, lat, res)
        @test !isempty(cell_id)
        @test length(cell_id) > 0

        center = cell_to_lon_lat(cell_id)
        @test length(center) == 2
        # Centre should be within ~1° of the input at res 10
        @test abs(center[1] - lon) < 1.0
        @test abs(center[2] - lat) < 1.0
    end

    @testset "cell_to_boundary" begin
        cell_id = lon_lat_to_cell(0.0, 0.0, 8)
        boundary = cell_to_boundary(cell_id)
        @test length(boundary) >= 5    # pentagon = 5 vertices (+ closing vertex)
        @test all(length(v) == 2 for v in boundary)
    end

    @testset "get_resolution" begin
        for res in [5, 8, 12]
            cell_id = lon_lat_to_cell(10.0, 50.0, res)
            @test get_resolution(cell_id) == res
        end
    end

    @testset "cell_to_children / cell_to_parent" begin
        cell_id = lon_lat_to_cell(10.0, 50.0, 8)
        children = cell_to_children(cell_id, 9)
        @test length(children) > 0
        @test all(!isempty(c) for c in children)

        parent = cell_to_parent(cell_id)
        @test !isempty(parent)
        @test get_resolution(parent) == 7
    end

    @testset "mesh_for_aoi — small bounding box" begin
        # Small 0.2° × 0.2° box — should produce a handful of cells at res 8
        mini_aoi = """
        {
          "type": "Feature",
          "geometry": {
            "type": "Polygon",
            "coordinates": [[[10.0, 50.0], [10.2, 50.0], [10.2, 50.2], [10.0, 50.2], [10.0, 50.0]]]
          },
          "properties": {}
        }
        """
        mesh = mesh_for_aoi(mini_aoi, 8)
        @test length(mesh) > 0
        @test mesh.resolution == 8
        @test all(!isempty(c.id) for c in mesh.cells)
    end

    @testset "save_mesh_geojson" begin
        mini_aoi = """
        {
          "type": "Feature",
          "geometry": {
            "type": "Polygon",
            "coordinates": [[[0.0, 0.0], [0.2, 0.0], [0.2, 0.2], [0.0, 0.2], [0.0, 0.0]]]
          },
          "properties": {}
        }
        """
        mesh = mesh_for_aoi(mini_aoi, 7)
        out = tempname() * ".geojson"
        save_mesh_geojson(mesh, out)
        @test isfile(out)
        content = read(out, String)
        @test occursin("FeatureCollection", content)
        rm(out)
    end

end

println("\nAll tests passed ✓")
