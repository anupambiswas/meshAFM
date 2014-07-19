This project case is about triangulating a polygonal region centered about origin (0, 0). The distance of the vertices are in the range of 0.65 to 1. These values may be changed in the test() function in proj_2.cpp. They are stored in variables rmin and rmax. The total number of points (n) and the number of boundary points are input by the user at the command line. The default values are n = 100 and nb = 40. This project also makes use of the class conv2D the details of which can be obtained from the conv2D section in my guthub account. It is basically used to generate internal points inside the arbitrary region. The final mesh is of poor quality near the boundary because very high aspect ratio triangles are formed. This will be rectified in the next project.

Instructions to build and run:

$ make

$ ./run
# This will run the case for default total number of points = 100 and boundary points = 40.

$ run 200 60
# This will run the case for total number of points = 200 and boundary points = 60.

In order to delete all files and folders created after a run, use the following.
$ make resetf
# deletes the results folder and other files

In order to vizualize the mesh generation process, start gnuplot.

$ gnuplot
> load 'results/dyn'
# This will start the animation and show how the faces are generated and also the advancing front.
