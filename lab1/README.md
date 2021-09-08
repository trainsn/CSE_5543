# CSE_5543 Lab1

run the command:
``` 
python ibezier.py --inputA {param}
``` 
If the input parameter is ab interger n, the user can interactively specify the location of the n control points by drawing on
a window region with a mouse, and the Bezier curve has degree n − 1.

Using the mouse, the user is able to interactively select any control point and move it. When any control point is moved, the program can automatically recompute and redraw
the Bezier curve.

When clicking the “Subdivide” button, the program would split the Bezier curve of degree
n − 1 into two Bezier curves at u = 1/2. Further calls to “Subdivide” would split each current Bezier curve into two, creating 4, 8 Bezier curves, and so on.

When clicking the “Save” button, the program allows the user to save the control
points to a file.

When input to the program is the name of a previously saved file, the program would
read the curve degree and control points from the file and use those as the starting point.

Sample file text:
``` 
BEZIER
# Sample file containing Bezier control points
# {dimension} {number of points} {degree}
2 4 3
0.2 0.2
0.3 0.8
0.6 0.7
0.8 0.3
``` 