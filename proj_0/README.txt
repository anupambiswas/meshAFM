This project case is about triangulating a circular region. The radius is 1. 

Instructions to build and run:

$ make

$ ./run
# This will run the case for default total number of points = 100 and boundary points = 40.

$ run 200 65
# This will run the case for total number of points = 200 and boundary points = 65.

In order to vizualize the mesh generation process, go to the results folder and start gnuplot.

$ gnuplot
> load 'dyn'
# This will start the animation and show how the faces are generated and also the advancing front.
