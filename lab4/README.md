# CSE_5543 Lab3

run the command:
``` 
python bsplinesurf.py --infile1 <infile1> --infile2 <infile2> --outfile <outfile>
``` 

I use the format for the B-splines in files infile1 and infile2 in the instructions for Lab 2.
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

The output is a file in Geomview (ASCII) OFF format. 
You can use the program meshlab to read and view the output GEOMVIEW OFF file.

The operating system I used to test my program is Windows 10. 