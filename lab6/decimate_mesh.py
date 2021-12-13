## \file decimate_mesh.py
#  Some simple mesh decimation routines.
#  Use data structure HALF_EDGE_MESH_DCMT_BASE (DCMT = decimate).

import math
from math import sqrt
from math import acos
import sys
import numpy as np

import half_edge_mesh
import half_edge_mesh_DCMT
import half_edge_mesh_IO
from half_edge_mesh_DCMT import HALF_EDGE_MESH_DCMT_BASE
from half_edge_mesh_DCMT\
    import VERTEX_DCMT_BASE, HALF_EDGE_DCMT_BASE, CELL_DCMT_BASE


def main(argv):
    global input_filename, output_filename
    global flag_allow_non_manifold, flag_fail_on_non_manifold

    # Initialize
    input_filename = None
    output_filename = None
    InitFlags()

    parse_command_line(sys.argv)
    mesh = HALF_EDGE_MESH_DCMT_BASE(VERTEX_DCMT_BASE, HALF_EDGE_DCMT_BASE, CELL_DCMT_BASE)
    half_edge_mesh_IO.open_and_read_off_file(input_filename, mesh)

    try:
        split_all_cells(mesh)
        old_triangles = mesh.NumCells()

        iteration = 0
        while True:
            join_triangle_with_large_angle(mesh)
            split_all_cells(mesh)
            collapse_shortest_edge_for_small_angle_cell(mesh)
            if old_triangles == mesh.NumCells():
                break
            old_triangles = mesh.NumCells()
            iteration += 1
            print("finish stage 1 iter {:d}".format(iteration))
            if iteration >= 5:
                break

        iteration = 0
        while True:
            split_longest_edge_for_large_angle_cell(mesh)
            split_all_cells(mesh)
            iteration += 1
            print("finish stage 1 iter {:d}".format(iteration))
            if iteration >= 5:
                break

        passed_check = check_mesh(mesh)
        print("Mesh data structure passed check.")
        print_mesh_info(mesh)

    except Exception as e:
        print(e)
        sys.stderr.write("Exiting.")
        exit(-1)

    if (output_filename is None):
        output_filename = "out.off"
    if (output_filename == input_filename):
        output_filename = "out2.off"

    print()
    print("Writing file: " + output_filename + ".")

    half_edge_mesh_IO.open_and_write_file(output_filename, mesh)


# *** Collapse edge routines ***
## Collapse edge
def collapse_edge(mesh, half_edge):
    global flag_allow_non_manifold

    flag = check_edge_collapse(mesh, half_edge)

    if mesh.IsIllegalEdgeCollapseH(half_edge):
        return

    if flag or flag_allow_non_manifold:
        print("Collapsing edge (" + half_edge.EndpointsStr(",") + ").")

        vnew = mesh.CollapseEdge(half_edge.Index())
        if (vnew is None):
            print("Skipped illegal collapse of edge (" + half_edge.EndpointsStr(",") + ").")
    else:
        print("Skipped collapse of edge (" + half_edge.EndpointsStr(",") + ").")

## Collapse shortest cell edge.
def collapse_shortest_cell_edge(mesh, icell):
    cell = mesh.Cell(icell)
    if (cell is None):
        return

    minL, maxL, ihalf_edge_min, ihalf_edge_max = cell.ComputeMinMaxEdgeLengthSquared()

    half_edge_min = mesh.HalfEdge(ihalf_edge_min)
    collapse_edge(mesh, half_edge_min)


## Collapse shortest edge in each cell.
def collapse_shortest_edge_for_small_angle_cell(mesh):
    n = mesh.NumCells()

    # Create a list of the cell indices.
    cell_indices_list = list(mesh.CellIndices())
    cell_indices_list.sort()

    # DO NOT iterate over CellIndices() directly, since collapse/split/join
    #   may delete cells.
    for icell in cell_indices_list:
        # Check if cell index is valid.
        #   cell may be none if it has been deleted from the cell dictionary
        cell = mesh.Cell(icell)
        if (cell is None):
            continue

        cos_minA, cos_maxA, ihalf_edge_min, ihalf_edge_max = cell.ComputeCosMinMaxAngle()
        if cos_minA > np.cos(30.*math.pi/180.) and cos_maxA > np.cos(150.*math.pi/180.):
            # print(icell, math.degrees(acos(cos_minA)), math.degrees(acos(cos_maxA)))
            collapse_shortest_cell_edge(mesh, icell)

# *** Split edge routines. ***
## Split edge.
def split_edge(mesh, half_edge):
    print("Splitting edge (" + half_edge.EndpointsStr(",") + ").")

    vnew = mesh.SplitEdge(half_edge.Index())
    if (vnew is None):
        print("Split of edge (" + half_edge.EndpointsStr(",") + ") failed.")


## Split longest cell edge.
def split_longest_cell_edge(mesh, icell):
    cell = mesh.Cell(icell)
    if (cell is None):
        return

    minL, maxL, ihalf_edge_min, ihalf_edge_max = cell.ComputeMinMaxEdgeLengthSquared()
    half_edge_max = mesh.HalfEdge(ihalf_edge_max)
    split_edge(mesh, half_edge_max)


## Split longest edge in each cell.
def split_longest_edge_for_large_angle_cell(mesh):
    # Create a list of the cell indices.
    cell_indices_list = list(mesh.CellIndices())
    cell_indices_list.sort()

    # DO NOT iterate over CellIndices() directly, since collapse/split/join
    #   may delete cells.
    for icell in cell_indices_list:
        # Check if cell index is valid.
        #   cell may be none if it has been deleted from the cell dictionary
        cell = mesh.Cell(icell)
        if (cell is None):
            continue

        cos_minA, cos_maxA, ihalf_edge_min, ihalf_edge_max = cell.ComputeCosMinMaxAngle()
        if cos_maxA < np.cos(120.*math.pi/180.):
            # min_angle = math.degrees(acos(cos_minA))
            # max_angle = math.degrees(acos(cos_maxA))
            # print(min_angle, 180 - min_angle - max_angle, max_angle)
            split_longest_cell_edge(mesh, icell,)


# *** Split cell routines. ***

## Split cell with diagonal (half_edgeA.FromVertex(), half_edgeB.FromVertex())
#  - Returns split edge.
#  - Returns None if split fails.
def split_cell(mesh, half_edgeA, half_edgeB):
    ihalf_edgeA = half_edgeA.Index()
    ihalf_edgeB = half_edgeB.Index()
    vA = half_edgeA.FromVertex()
    vB = half_edgeB.FromVertex()
    ivA = vA.Index()
    ivB = vB.Index()
    icell = half_edgeA.CellIndex()

    flag = check_split_cell(mesh, half_edgeA, half_edgeB)

    if (mesh.IsIllegalSplitCell(half_edgeA, half_edgeB)):
        return None

    if (flag or flag_allow_non_manifold):
        print("Splitting cell {:d} at diagonal ({:d},{:d}).".format(icell, ivA, ivB))
        split_edge = mesh.SplitCell(ihalf_edgeA, ihalf_edgeB)
        if (split_edge is None):
            print("Split of cell {:d} at diagonal (""{:d},{:d}) failed.".format(icell, ivA, ivB))
        return split_edge
    else:
        print("Skipping split of cell {:d} at diagonal ({:d},{:d}).".format(icell, ivA, ivB))
        return None


## Split cell at largest angle.
#  - Split cell at vertex forming the largest angle.
#  - Split as many times as necessary to triangulate.
def split_cell_at_largest_angle(mesh, cell):
    cos_minA, cos_maxA, ihalf_edge_min, ihalf_edge_max = cell.ComputeCosMinMaxAngle()

    half_edgeA = mesh.HalfEdge(ihalf_edge_max)

    while not half_edgeA.Cell().IsTriangle():
        half_edgeB = half_edgeA.PrevHalfEdgeInCell().PrevHalfEdgeInCell()
        vA = half_edgeA.FromVertex()
        vB = half_edgeB.FromVertex()
        ivA = vA.Index()
        ivB = vB.Index()

        split_edge = split_cell(mesh, half_edgeA, half_edgeB)
        if (split_edge is None):
            # Cannot split cell at largest angle.
            return

        # Get largest angle in remaining cell.
        cellA = split_edge.Cell()
        cos_minA, cos_maxA, ihalf_edge_min, ihalf_edge_max =\
            cell.ComputeCosMinMaxAngle()
        half_edgeA = mesh.HalfEdge(ihalf_edge_max)


## Split all cells.
def split_all_cells(mesh):
    n = mesh.MaxCellIndex()
    flag_check = False

    # Create a list of the cell indices.
    cell_indices_list = list(mesh.CellIndices())
    cell_indices_list.sort()

    # DO NOT iterate over CellIndices() directly, since collapse/split/join
    #   may delete cells.
    kount = 0
    for icell in cell_indices_list:
        # Check if cell index is valid.
        #   cell may be none if it has been deleted from the cell dictionary
        cell = mesh.Cell(icell)
        if (cell is None):
            continue

        split_cell_at_largest_angle(mesh, cell)
        kount = kount + 1


# *** Join cell routines ***
def join_two_cells(mesh, half_edge):
    half_edgeX = half_edge.NextHalfEdgeAroundEdge()
    icell = half_edge.CellIndex()
    icellX = half_edgeX.CellIndex()

    flag = check_join_cell(mesh, half_edge)

    if (mesh.IsIllegalJoinCells(half_edge)):
        return

    if flag:
        print("Joining cell {:d} to cell {:d} by deleting edge (".format(icell, icellX) +
                half_edge.EndpointsStr(",") + ").")

        half_edgeB = mesh.JoinTwoCells(half_edge.Index())
        if half_edgeB is None:
            print("Join of cell {:d} to cell {:d} failed.".format(icell, icellX))

        return
    else:
        print("Skipping join of cell {:d} with cell {:d}.".format(icell, icellX))


## Attempt to join each cell by deleting longest edge.
def join_triangle_with_large_angle(mesh):
    # Create a list of the cell indices.
    cell_indices_list = list(mesh.CellIndices())
    cell_indices_list.sort()

    # DO NOT iterate over CellIndices() directly, since collapse/split/join
    #   may delete cells.
    for icell in cell_indices_list:
        # Check if cell index is valid.
        #   cell may be none if it has been deleted from the cell dictionary
        cell = mesh.Cell(icell)
        if (cell is None):
            continue

        if (cell.NumVertices() > 3):
            continue

        cos_minA, cos_maxA, ihalf_edge_min, ihalf_edge_max = cell.ComputeCosMinMaxAngle()
        if cos_maxA < np.cos(120.*math.pi/180):
            # min_angle = math.degrees(acos(cos_minA))
            # max_angle = math.degrees(acos(cos_maxA))
            # print(min_angle, 180 - min_angle - max_angle, max_angle)
            half_edge = mesh.HalfEdge(ihalf_edge_max).NextHalfEdgeInCell()
            join_two_cells(mesh, half_edge)

# *** Check routines ***

def check_oriented_manifold(mesh):
    flag_non_manifold_vertex, flag_non_manifold_edge, iv, ihalf_edgeA = mesh.CheckManifold()

    if flag_non_manifold_edge:
        half_edgeA = mesh.HalfEdge(ihalf_edgeA)
        sys.stderr.write("Warning: Non-manifold edge (" + half_edgeA.EndpointsStr(",") + ").\n")

        # Non-manifold edge automatically implies inconsistent orientations.
        return False

    flag_orientation, ihalf_edgeB = mesh.CheckOrientation()

    if flag_orientation:
        if flag_non_manifold_vertex:
            sys.stderr.write("Warning: Non-manifold vertex {iv}.")

            return False
    else:
        if flag_non_manifold_vertex:
            sys.stderr.write\
                ("Warning: Non-manifold vertex or inconsistent orientations in cells incident on vertex {:d}.\n".format(iv))
        else:
            sys.stderr.write\
                ("Warning: Inconsistent orientation of cells incident on edge (" + half_edgeB.EndpointsStr(",") + ").")
        return False

    return True


def check_mesh(mesh):
    global flag_fail_on_non_manifold

    flag, error_msg = mesh.CheckAll()
    if not(flag):
        sys.stderr.write("Error detected in mesh data structure.\n")
        if not(error_msg is None):
            sys.stderr.write(error_msg + "\n")
        exit(-1)

    if (flag_fail_on_non_manifold):
        flag_oriented_manifold = check_oriented_manifold(mesh)
        if (flag_fail_on_non_manifold and not(flag_oriented_manifold)):
            sys.stderr.write("Detected non-manifold or inconsistent orientations")
            sys.stderr.write("Exiting.")
            exit(-1)
        return flag_oriented_manifold
    else:
        return True


## Print a warning message if collapsing half_edge is illegal or
#    will change mesh topology.
#  - Return True if collapse is not illegal and does not change
#    mesh topology.
def check_edge_collapse(mesh, half_edge):
    icell = half_edge.CellIndex()
    return_flag = True

    if (mesh.IsIllegalEdgeCollapseH(half_edge)):
        print("Collapse of edge  (" + half_edge.EndpointsStr(",") +") is illegal.")
        print("  Some cell contains vertices " +\
                half_edge.EndpointsStr(" and ") + " but not edge (" +\
                half_edge.EndpointsStr(",") + ").")

        return_flag = False

    flag, ivC = mesh.FindTriangleHole(half_edge)
    if (flag):
        print("Collapsing edge (" + half_edge.EndpointsStr(",") +") will change the mesh topology.")
        print("  Vertices (" + half_edge.EndpointsStr(", ") +", " + str(ivC) + ") form a triangle hole.")

        return_flag = False

    if not(half_edge.IsBoundary()):
        if half_edge.FromVertex().IsBoundary() and half_edge.ToVertex().IsBoundary():
            print("Collapsing edge (" + half_edge.EndpointsStr(",") + ") merges two non-adjacent boundary vertices.")
            return_flag = False

    if mesh.IsIsolatedTriangle(icell):
        print("Collapsing edge(" + half_edge.EndpointsStr(",") +") will delete isolated cell {:d}.".format(icell))
        return_flag = False

    if mesh.IsInTetrahedron(icell):
        print("Collapsing edge (" + half_edge.EndpointsStr(",") +") will collapse a tetrahedron.")
        return_flag = False

    return return_flag

## Print a warning message if splitting cell at diagonal
#    (half_edgeA.FromVertex(), half_edgeB.FromVertex())
#    will change the mesh topology.
#  - Return True if split does not change manifold topology.
def check_split_cell(mesh, half_edgeA, half_edgeB):
    vA = half_edgeA.FromVertex()
    vB = half_edgeB.FromVertex()
    ivA = vA.Index()
    ivB = vB.Index()
    icell = half_edgeA.CellIndex()
    half_edgeC = mesh.FindEdge(vA, vB)

    flag_cell_edge = False
    return_flag = True

    if (mesh.IsIllegalSplitCell(half_edgeA, half_edgeB)):
        if (vA is half_edgeB.ToVertex()) or (vB is half_edgeA.FromVertex()):
            flag_cell_edge = True

        if flag_cell_edge:
            print("({:d},{:d}) is a cell edge, not a cell diagonal.".format(ivA, ivB))
        else:
            print("Illegal split of cell {:d} with diagonal ({:d}{:d}).".format(icell, ivA, ivB))
        return_flag = False

    if not(half_edgeC is None) and not(flag_cell_edge):
        sys.stdout.write("Splitting cell {:d} with diagonal ({:d},{:d})".format(icell, ivA, ivB))
        sys.stdout.write(" creates an edge incident on three or nmore cells.\n")
        return_flag = False

    return return_flag


## Print a warning if joining cells separated by half_edge is illegal.
#  - Return true if join is legal.
def check_join_cell(mesh, half_edge):
    TWO = 2
    return_flag = True

    if (mesh.IsIllegalJoinCells(half_edge)):
        half_edgeX = half_edge.NextHalfEdgeAroundEdge()

        if (half_edge.IsBoundary()):
            print("Only one cell contains edge (" + \
                  half_edge.EndpointsStr(",") + ").")
        elif not (half_edge.FromVertex().IsIncidentOnMoreThanTwoEdges()):
            print("Half edge endpoint {half_edge.FromVertexIndex()} is incident on only two edges.")
        elif not (half_edge.ToVertex().IsIncidentOnMoreThanTwoEdges()):
            print("Half edge endpoint {half_edge.ToVertexIndex()} is incident on only two edges.")
        elif not (half_edge is half_edgeX.NextHalfEdgeAroundEdge()):
            print("More than two cells are incident on edge (" + \
                  half_edge.EndpointsStr(",") + ").")
        else:
            cell = half_edge.Cell()
            cellX = half_edgeX.Cell()
            num_shared_vertices = \
                mesh.CountNumVerticesSharedByTwoCells(cell, cellX)
            if (num_shared_vertices > TWO):
                print("Cells {cell.Index()} and {cellX.Index()} share {num_shared_vertices} vertices.")
            else:
                print("Join of two cells incident on edge (" + \
                      half_edge.EndpointsStr(",") + ") is illegal.")

        return_flag = False

    return return_flag


# *** Init/parse/print/prompt functions. ***

## Initialize global flags.
def InitFlags():
        global input_filename, output_filename
        global flag_allow_non_manifold, flag_fail_on_non_manifold

        # Initialize
        input_filename = None
        output_filename = None
        flag_allow_non_manifold = False
        flag_fail_on_non_manifold = False


def parse_command_line(argv):
    global input_filename, output_filename
    global flag_allow_non_manifold, flag_fail_on_non_manifold

    iarg = 1
    while iarg < len(argv) and argv[iarg][0] == '-':
        s = argv[iarg]
        if (s == "-allow_non_manifold"):
            flag_allow_non_manifold = True
        elif (s == "-fail_on_non_manifold"):
            flag_fail_on_non_manifold = True
        else:
            sys.stderr.write("Usage error. Option " + s + " is undefined.\n")
            usage_error()

        iarg = iarg + 1

    if (iarg >= len(argv) or iarg+2 < len(argv)):
        usage_error()

    input_filename = argv[iarg]

    if (iarg+1 < len(argv)):
        output_filename = argv[iarg+1]

## Print mesh information, such as number of vertices, edges, cells
#    min and max edge lengths and cell angles.
def print_mesh_info(mesh):
    FIVE = 5
    num_vertices = mesh.NumVertices()
    num_edges = mesh.CountNumEdges()
    num_boundary_edges = mesh.CountNumBoundaryEdges()
    num_cells = mesh.NumCells()
    num_triangles = mesh.CountNumTriangles()
    num_quads = mesh.CountNumQuads()
    num_large_cells = mesh.CountNumCellsOfSizeGE(FIVE)

    minL_squared, maxL_squared, ihmin, ihmax =\
        mesh.ComputeMinMaxEdgeLengthSquared()
    min_ratio_squared, icmin, Lmin, Lmax, ihmin, ihmax =\
        mesh.ComputeMinCellEdgeLengthRatioSquared()
    cos_minA, cos_maxA, ihmin, ihmax = mesh.ComputeCosMinMaxAngle()

    flag_non_manifoldV, flag_non_manifoldE, iv, ie = mesh.CheckManifold()
    is_oriented, iv = mesh.CheckOrientation()

    print("Number of vertices: ", num_vertices)
    print("Number of mesh edges: ", num_edges)
    print("Number of boundary edges: ", num_boundary_edges)
    print("Number of mesh cells: ", num_cells)
    print("  Number of mesh triangles: ", num_triangles)
    print("  Number of mesh quadrilaterals: ", num_quads)
    print("  Number of mesh cells with > 4 vertices: ", num_large_cells)
    print("Min edge length: {:.4f}".format(sqrt(minL_squared)))
    print("Max edge length: {:.4f}".format(sqrt(maxL_squared)))
    print("Min cell edge length ratio: {:.4f}".format(sqrt(min_ratio_squared)))
    print("Minimum cell angle: {:.4f}".format(math.degrees(acos(cos_minA))))
    print("Maximum cell angle: {:.4f}".format(math.degrees(acos(cos_maxA))))

    if not flag_non_manifoldV and not flag_non_manifoldE and is_oriented:
        print("Mesh is an oriented manifold.")
    else:
        print("Mesh is non-manifold or has inconsistent cell orientations.")

def usage_error():
    sys.stderr.flush()
    exit(-1)

if __name__ == '__main__':
    main(sys.argv)
