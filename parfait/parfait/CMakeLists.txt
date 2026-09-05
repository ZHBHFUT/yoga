cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

PROJECT(parfait-core-unit-tests)
set(CMAKE_CXX_STANDARD 14)

find_package(tracer REQUIRED)

set(StringToolsFiles
        Dictionary.cpp
        DictionaryRules.cpp
        StringTools.cpp
        ToString.cpp
        )

set(StringToolsHeaders
        Throw.h
        Dictionary.h
        DictionaryRules.h
        FortranSyntax.h
        FortranUnformatted.h
        DictionaryEntry.h
        DictionarySchemaValidator.h
        JsonParser.h
        JsonPrinters.h
        Lexer.h
        OmlParser.h
        NamelistParser.h
        Namelist.h
        PomlParser.h
        StringTools.h
        ToString.h
        RegexMatchers.h
        )

set(CoreFiles
        FileTools.cpp
        MetricDecomposition.cpp
        Sphere.cpp
        DataFrame.cpp
        )

set(CoreHeaders
        Adt.h
        Adt.hpp
        Adt3dExtent.h
        Adt3dExtent.hpp
        Adt3dPoint.h
        Adt3dPoint.hpp
        Barycentric.h
        ByteSwap.h
        CellTesselator.h
        CGNSElements.h
        CGNSFaceExtraction.h
        CommandLine.h
        CommandLineMenu.h
        CRS.h
        CartBlock.h
        CartBlock.hpp
        CartBlockVisualize.h
        CellWindingConverters.h
        CellContainmentChecker.h
        Checkpoint.h
        ConsoleRedirect.h
        ContiguousMatrixMultiplication.h
        DataFrame.h
        DataFrameTools.h
        DenseMatrix.h
        Ellipsoid.h
        FloatingPointExceptions.h
        MetricDecomposition.h
        HermiteSpline.h
        Decompositions.h
        DistanceTree.h
        DistanceTree.hpp
        DualMetrics.h
        EdgeBuilder.h
        EdgeSegmentProjection.h
        Intersections.h
        Edges.h
        Extent.h
        Extent.hpp
        ExtentBuilder.h
        ExtentMeshWrapper.h
        ExtentWriter.h
        ExtentWriter.hpp
        Exception.h
        Facet.h
        Facet.hpp
        FacetIntersection.h
        Facetize.h
        FileTools.h
        Flatten.h
        GraphPlotter.h
        GraphDistance.h
        GreedyColoring.h
        ImportedUgrid.h
        ImportedUgrid.hpp
        ImportedUgridFactory.h
        ImpliedMetric.h
        Interpolation.h
        JsonCommon.h
        KMeans.h
        LagrangeTriangle.h
        LagrangeTriangle.hpp
        LaplacianSmoothing.h
        LeastSquaresReconstruction.h
        LinearAlgebra.h
        LinearPartitioner.h
        LinearPartitioner.hpp
        LinePlotters.h
        LinePlotters.hpp
        Linspace.h
        MetricTensor.h
        MapbcParser.h
        MapbcReader.h
        MapbcReader.hpp
        MapbcWriter.h
        MapTools.h
        MarchingCubes.h
        MatrixMarketReader.h
        MinimumStorageSet.h
        MatrixBase.h
        MatrixStorage.h
        MatrixOperations.h
        Morton.h
        MotionMatrix.h
        MotionMatrix.hpp
        MovingAverages.h
        NodeToCell.h
        NodeToNode.h
        NodeToNode.hpp
        Normals.h
        ObjReader.h
        Octree.h
        OctreeRasterizer.h
        OctreeStorage.h
        ParfaitMacros.h
        Point.h
        Point.hpp
        PointwiseSources.h
        PointWriter.h
        PolygonClipping.h
        PointGenerator.h
        QuadTree.h
        Region.h
        RegionFactory.h
        RewindLeftHandedElements.h
        RingBuffer.h
        SCurve.h
        SetTools.h
        Shell.h
        Sphere.h
        StackTrace.h
        STL.h
        STL.hpp
        STLFactory.h
        STLReader.h
        STLReader.hpp
        SparseMatrix.h
        SparseMatrix.hpp
        StitchMesh.h
        TaubinSmoothing.h
        TecplotBinaryHelpers.h
        TecplotBrickCollapse.h
        TecplotBinaryReader.h
        TecplotBinaryWriter.h
        TecplotBinaryWriter.hpp
        TecplotWriter.h
        TecplotWriter.hpp
        TetGenWriter.h
        Timing.h
        Timing.hpp
        Topology.h
        TriReader.h
        UgridReader.h
        UgridReader.hpp
        UgridWriter.h
        UgridWriter.hpp
        UniquePoints.h
        UnitTransformer.h
        VTKElements.h
        VTKWriter.h
        VTKWriter.hpp
        VectorTools.h
        VolumeMeshAdapter.h
        VolumeMeshInterface.h
        VizualizationWriter.h
        Plane.h
        Plot3DReader.h
        Plot3DWriter.h
        VTKWriter.h
        Wiggle.h
        )

set(CoreMPIHeaders
        GraphColoring.h
        GraphOrderings.h
        Inspector.h
        Metrics.h
        ParallelExtent.h
        PointWriter.h
        RecursiveBisection.h
        RecursiveBisection.hpp
        SyncField.h
        SyncPattern.h
        )

add_library(StringTools_obj OBJECT ${StringToolsFiles})
set_target_properties(StringTools_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(StringTools_obj PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
        )

add_library(StringTools SHARED $<TARGET_OBJECTS:StringTools_obj>)
add_library(parfait::StringTools ALIAS StringTools)
#target_link_libraries(StringTools PUBLIC MessagePasser::MessagePasser)
set_target_properties(StringTools PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(StringTools PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
        )

add_library(StringTools_static STATIC $<TARGET_OBJECTS:StringTools_obj>)
if(WIN32)
  set_target_properties(StringTools_static PROPERTIES OUTPUT_NAME StringTools_static)
else()
  set_target_properties(StringTools_static PROPERTIES OUTPUT_NAME StringTools)
endif()
add_library(parfait::StringTools_static ALIAS StringTools_static)
set_target_properties(StringTools_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
#target_link_libraries(StringTools_static PUBLIC MessagePasser::MessagePasser)
target_include_directories(StringTools_static PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
        )

set_standard_ring_rpath(StringTools)

install(EXPORT parfait
        NAMESPACE parfait::
        DESTINATION cmake)

install(FILES ${StringToolsHeaders} DESTINATION include/parfait)
install(TARGETS StringTools StringTools_static
        EXPORT parfait
        DESTINATION lib
        ARCHIVE DESTINATION lib
        COMPONENT StringTools)

add_library(parfait_obj ${CoreFiles})
add_library(parfait_static STATIC $<TARGET_OBJECTS:parfait_obj>)
add_library(parfait SHARED $<TARGET_OBJECTS:parfait_obj>)
target_include_directories(parfait_static PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
        )
target_include_directories(parfait PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:include>
        )

add_library(parfait::parfait ALIAS parfait)
add_library(parfait::parfait_static ALIAS parfait_static)
if(WIN32)
  set_target_properties(parfait_static PROPERTIES OUTPUT_NAME "parfait_static")
else()
  set_target_properties(parfait_static PROPERTIES OUTPUT_NAME "parfait")
endif()
target_link_libraries(parfait_obj PUBLIC StringTools_obj)
target_link_libraries(parfait PUBLIC StringTools)
target_link_libraries(parfait_static PUBLIC StringTools_static)
ring_precompiled_headers(parfait PUBLIC StringTools.h)

set_standard_ring_rpath(parfait)

install(TARGETS parfait_static parfait
        EXPORT parfait
        DESTINATION lib
        ARCHIVE DESTINATION lib
        )

if (MessagePasser_FOUND)
    target_link_libraries(parfait_static PUBLIC MessagePasser::MessagePasser)
    target_link_libraries(parfait PUBLIC MessagePasser::MessagePasser)
    install(FILES ${CoreMPIHeaders} ${CoreHeaders} DESTINATION include/parfait)
else()
    install(FILES ${CoreHeaders} DESTINATION include/parfait)
endif ()

