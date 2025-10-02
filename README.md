This is a basic implementation of the semi-implicit method for pressure linked equations (SIMPLE) algorithm. The code also includes a solver for species flow.
The algorithm works well for low Reynolds numbers (<100), some tuning of the under-relaxiation factors may be needed to get the simulation to converge.
The plotting and simulation driver where originaly done with OpenGL. Since the code for rendering the data in OpenGL was developed by a third party, it has not been included here.
If you intend to use this algorithm, you will have to adapt the code to use your own method of displaying the results.
