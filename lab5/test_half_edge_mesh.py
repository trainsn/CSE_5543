## \file test_half_edge_mesh.py
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

    output_filename = "out.off"
    half_edge_mesh_IO.open_and_write_file(output_filename, mesh)

if __name__ == '__main__':
    main(parse_args())
