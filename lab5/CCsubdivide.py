#  - Reads a .off file into HALF_EDGE_MESH_BASE.
#  - Checks mesh, manifold and orientations properties of the mesh.
#  - Apply Catmull-Clark subdivision algorithm
#  - Writes the output mesh to a .off file.

import argparse
import sys
import half_edge_mesh_IO
from half_edge_mesh import VERTEX_BASE, HALF_EDGE_BASE, CELL_BASE
from half_edge_mesh import HALF_EDGE_MESH_BASE

# parse arguments
def parse_args():
  parser = argparse.ArgumentParser(description="Catmull-Clark subdivision")

  parser.add_argument("--num_iter", type=int, default=1,
                      help="the number of iterations of the Catmull-Clark subdivision.")
  parser.add_argument("--input", type=str, required=True,
                      help="the input off file.")

  return parser.parse_args()

def main(args):
    # log hyperparameters
    print(args)
    input_filename = args.input

    initial_mesh = HALF_EDGE_MESH_BASE(VERTEX_BASE, HALF_EDGE_BASE, CELL_BASE)
    half_edge_mesh_IO.open_and_read_off_file(input_filename, initial_mesh)
    flag, error_msg = initial_mesh.CheckAll()
    if not flag:
        sys.stderr.write("Error detected in mesh data structure.\n")
        if not(error_msg is None):
            sys.stderr.write(error_msg + "\n")
        exit(-1)

    meshes = [initial_mesh]
    for iter in range(args.num_iter):
        mesh = meshes[iter]
        mesh_new = HALF_EDGE_MESH_BASE(VERTEX_BASE, HALF_EDGE_BASE, CELL_BASE)
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
                vertex = e.FromVertex()
                X += vertex.KthCoord(0)
                Y += vertex.KthCoord(1)
                Z += vertex.KthCoord(2)
                e = e.NextHalfEdgeInCell()
                if e.Index() == e_start.Index():
                    break
            coord = [X / num, Y / num, Z / num]
            mesh_new.AddVertex(idx)
            mesh_new.SetCoord(idx, coord)
            cell.face_point_idx = idx
            idx += 1

        # generate edge points
        for ihalf_edge in mesh.HalfEdgeIndices():
            edge = mesh.HalfEdge(ihalf_edge)
            if edge.edge_point_idx != -1:
                continue
            if edge.IsBoundary():
                st, en = edge.FromVertex(), edge.ToVertex()
                coord = [(st.KthCoord(0) + en.KthCoord(0)) / 2.0,
                         (st.KthCoord(1) + en.KthCoord(1)) / 2.0,
                         (st.KthCoord(2) + en.KthCoord(2)) / 2.0]
                mesh_new.AddVertex(idx)
                mesh_new.SetCoord(idx, coord)
                edge.edge_point_idx = idx
            else:
                st, en = edge.FromVertex(), edge.ToVertex()
                face1, face2 = \
                    mesh_new.Vertex(edge.Cell().face_point_idx), \
                    mesh_new.Vertex(edge.NextHalfEdgeAroundEdge().Cell().face_point_idx)
                coord = [(st.KthCoord(0) + en.KthCoord(0) + face1.KthCoord(0) + face2.KthCoord(0)) / 4.0,
                         (st.KthCoord(1) + en.KthCoord(1) + face1.KthCoord(1) + face2.KthCoord(1)) / 4.0,
                         (st.KthCoord(2) + en.KthCoord(2) + face1.KthCoord(2) + face2.KthCoord(2)) / 4.0]
                mesh_new.AddVertex(idx)
                mesh_new.SetCoord(idx, coord)
                edge.edge_point_idx = idx
                edge.NextHalfEdgeAroundEdge().edge_point_idx = idx

            idx += 1

        # generate new vertices
        for iv in mesh.VertexIndices():
            v = mesh.Vertex(iv)
            valence = v.NumHalfEdgesFrom()

            is_boundary = False
            for k in range(valence):
                half_edge = v.KthHalfEdgeFrom(k)
                if half_edge.IsBoundary():
                    is_boundary = True
                    break

            coord = [0., 0., 0.]
            if is_boundary:     # old vertex point
                coord[0] = v.KthCoord(0)
                coord[1] = v.KthCoord(1)
                coord[2] = v.KthCoord(2)
            else:
                # find the face average
                face_avg = [0., 0., 0.]
                for k in range(valence):
                    half_edge = v.KthHalfEdgeFrom(k)
                    face_point = mesh_new.Vertex(half_edge.Cell().face_point_idx)
                    face_avg[0] += face_point.KthCoord(0) / valence
                    face_avg[1] += face_point.KthCoord(1) / valence
                    face_avg[2] += face_point.KthCoord(2) / valence

                # find the edge_point average
                edge_avg = [0., 0., 0.]
                for k in range(valence):
                    half_edge_point = mesh_new.Vertex(v.KthHalfEdgeFrom(k).edge_point_idx)
                    edge_avg[0] += half_edge_point.KthCoord(0) / valence
                    edge_avg[1] += half_edge_point.KthCoord(1) / valence
                    edge_avg[2] += half_edge_point.KthCoord(2) / valence

                coord[0] = (face_avg[0] + edge_avg[0] * 2 + v.KthCoord(0) * (valence - 3)) / valence
                coord[1] = (face_avg[1] + edge_avg[1] * 2 + v.KthCoord(1) * (valence - 3)) / valence
                coord[2] = (face_avg[2] + edge_avg[2] * 2 + v.KthCoord(2) * (valence - 3)) / valence
            mesh_new.AddVertex(idx)
            mesh_new.SetCoord(idx, coord)
            v.new_point_idx = idx
            idx += 1

        # connect new mesh
        num_faces = 0
        for icell in mesh.CellIndices():
            cell = mesh.Cell(icell)
            e_start = cell.HalfEdge()
            e = e_start
            while True:
                v1 = cell.face_point_idx
                v2 = e.edge_point_idx
                v3 = e.ToVertex().new_point_idx
                v4 = e.NextHalfEdgeInCell().edge_point_idx
                mesh_new.AddCell(num_faces, [v1, v2, v3, v4])
                num_faces += 1
                e = e.NextHalfEdgeInCell()
                if e.Index() == e_start.Index():
                    break

        meshes.append(mesh_new)

    output_filename = "out.off"
    half_edge_mesh_IO.open_and_write_file(output_filename, meshes[args.num_iter])

if __name__ == '__main__':
    main(parse_args())
