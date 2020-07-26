using Gnuplot

"""
View
The set view command sets the viewing angle for splots.
It controls how the 3-d coordinates of the plot are mapped into the 2-d screen space.
It provides controls for both rotation and scaling of the plotted data, but supports
orthographic projections only. It supports both 3D projection or orthogonal 2D projection
into a 2D plot-like map.

Syntax:

     set view { <rot_x>{,{<rot_z>}{,{<scale>}{,<scale_z>}}} | map }
     show view

where rot_x and rot_z control the rotation angles (in degrees) in a virtual 3-d
coordinate system aligned with the screen such that initially (that is, before
the rotations are performed) the screen horizontal axis is x, screen vertical
axis is y, and the axis perpendicular to the screen is z. The first rotation
applied is rot_x around the x axis. The second rotation applied is rot_z
around the new z axis.

Command set view map is used to represent the drawing as a map. It can be used
for contour plots, or for color pm3d maps. In the latter, take care that you properly
use zrange and cbrange for input data point filtering and color range scaling, respectively.

rot_x is bounded to the [0:180] range with a default of 60 degrees,
while rot_z is bounded to the [0:360] range with a default of 30 degrees.
scale controls the scaling of the entire splot, while scale_z scales the z axis only. Both scales default to 1.0.

Examples:

     set view 60, 30, 1, 1
     set view ,,0.5

The first sets all the four default values. The second changes only scale, to 0.5. 
"""
function example_01()
    x = y = -15:0.33:15
    fz(x,y) = sin.(sqrt.(x.^2 + y.^2))./sqrt.(x.^2+y.^2)
    fxy = [fz(x,y) for x in x, y in y]

    @gp "set term pdfcairo size 13cm,12cm" :-
    @gp :- "set view 60, 30, 1, 1" :-
    @gp :- "set output 'IMG_test_gnuplot.pdf'" :-
    @gsp :- x y fxy "w pm3d" "set ylabel 'y' textcolor rgb 'red'" :-
    @gsp :- "set xlabel 'x' textcolor rgb 'blue'"
end

example_01()