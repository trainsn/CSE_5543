##   file test_half_edge_mesh.py
#    Test program for half_edge_mesh.py.
#
#  - Reads a .off file into HALF_EDGE_MESH_BASE.
#  - Checks mesh, manifold and orientations properties of the mesh.
#  - Writes the mesh to a .off file.

import argparse
import sys
import half_edge_mesh
import half_edge_mesh_IO
from half_edge_mesh import VERTEX_BASE, HALF_EDGE_BASE, CELL_BASE
from half_edge_mesh import HALF_EDGE_MESH_BASE

# parse arguments
def parse_args():
  parser = argparse.ArgumentParser(description="Deep Learning Model")

  parser.add_argument("--num_iter", type=int, default=1,
                      help="the number of the Catmull-Clark subdivision.")
  parser.add_argument("--input", type=str, required=True,
                      help="the input off file.")

  return parser.parse_args()

def main(args):
    # log hyperparameters
    print(args)
    input_filename = args.input

    mesh = HALF_EDGE_MESH_BASE(VERTEX_BASE, HALF_EDGE_BASE, CELL_BASE)
    half_edge_mesh_IO.open_and_read_off_file(input_filename, mesh)
    flag, error_msg = mesh.CheckAll()
    if not flag:
        sys.stderr.write("Error detected in mesh data structure.\n")
        if not(error_msg is None):
            sys.stderr.write(error_msg + "\n")
        exit(-1)

    for iter in range(args.num_iter):
        meshNew = HALF_EDGE_MESH_BASE(VERTEX_BASE, HALF_EDGE_BASE, CELL_BASE)
        idx = 0

        # generate face points
        for icell in mesh.CellIndices():
            cell = mesh.Cell(icell)
            e_start = cell.HalfEdge()
            e = e_start
            num = 0
            X, Y, Z = 0., 0., 0.
            while True:
                num += 1
                vertex = e.fromVertex()
                X += vertex.KthCoord(0)
                Y += vertex.KthCoord(1)
                Z += vertex.KthCoord(2)
                e = e.NextHalfEdgeInCell()
                if e.Index() == e_start.Index():
                    break
            coord = [X / num, Y / num, Z / num]
            meshNew.AddVertex(idx)
            meshNew.SetCoord(idx, coord)
            cell.face_point = meshNew.Cell(idx)
            idx += 1

        # generate edge points
        for ihalf_edge in mesh.HalfEdgeIndices():
            edge = mesh.HalfEdge(ihalf_edge)
            if edge.edge_point is not None:
                continue
            if edge.IsBoundary():
                st, en = edge.FromVertex(), edge.ToVertex()
                coord = [(st.KthCoord(0) + en.KthCoord(0)) / 2.0,
                         (st.KthCoord(1) + en.KthCoord(1)) / 2.0,
                         (st.KthCoord(2) + en.KthCoord(2)) / 2.0]
                meshNew.AddVertex(idx)
                meshNew.SetCoord(idx, coord)
                edge.edge_point = meshNew.Cell(idx)
            else:
                st, en = edge.FromVertex(), edge.ToVertex()
                face1, face2 = edge.Cell().face_point, edge.NextHalfEdgeAroundEdge().Cell().face_point
                coord = [(st.KthCoord(0) + en.KthCoord(0) + face1.KthCoord(0) + face2.KthCoord(0)) / 4.0,
                         (st.KthCoord(1) + en.KthCoord(1) + face1.KthCoord(1) + face2.KthCoord(1)) / 4.0,
                         (st.KthCoord(2) + en.KthCoord(2) + face1.KthCoord(2) + face2.KthCoord(2)) / 4.0]
                meshNew.AddVertex(idx)
                meshNew.SetCoord(idx, coord)
                edge.edge_point = meshNew.Cell(idx)
                edge.NextHalfEdgeAroundEdge().edge_point = meshNew.Cell(idx)

            idx += 1


    output_filename = "out.off"
    half_edge_mesh_IO.open_and_write_file(output_filename, mesh)

if __name__ == '__main__':
    main(parse_args())
