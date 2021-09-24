# CSE_5543 Lab2

run the command:
``` 
python ibezier.py --mode {`file` or `integer`} --inputA {param}
``` 
If the mode parameter is integer, then "inputA" parameter should be two intergers n,d (without blank). 
The user can interactively specify the location of the n control points by drawing on a window region with a mouse, and the B-spline curve has degree d.

Using the mouse, the user is able to interactively select any control point and move it. When any control point is moved, the program can automatically recompute and redraw
the B-spline curve.

When clicking the “Add CP” button, the program would prompt the user to select a control point q<sub>h</sub>, where h ≥ d. 
The program would create a new knot vector with a new knot t = (t<sub>h</sub> +t<sub>h+1</sub>)/2 between knots t<sub>h</sub> and t<sub>h+1</sub>.
The program replaces the d−1 control points q<sub>h−d+1</sub>, q<sub>h−d+2</sub>, ..., q<sub>h−1</sub> with d new control
points using the Boehm algorithm.

When clicking the “Save” button, the program allows the user to save the control
points and the knot vector to a file.

When input to the program is the name of a previously saved file, the program would
read the curve degree, knot vector control points from the file and use those as the starting point.

Sample file text:
``` 
BSPLINE
# Sample file containing Bspline control points
# {dimension} {number of points} {degree}
2 7 2
0 0 0 1 2 3 4 5 5 5
0.15 0.2
0.20 0.2
0.30 0.8
0.34 0.6
0.42 0.7
0.63 0.5
0.7 0.3
``` 