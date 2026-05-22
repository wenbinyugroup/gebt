
# Outputs

GEBT outputs 1D beam displacements, rotations, sectional forces and sectional moments for each point defined in the input file and each element in a file named input file name.out. The element is identified using the coordinates of its mid-point.

First, we output the coordinates in the global coordinate system arranged as \(x _ { 1 } , x _ { 2 } , x _ { 3 }\) . Immediately following, we output the results for displacement/rotation variables. For each line, there are 6 values arranged as u1 u2 u3 \(\theta _ { 1 } ~ \theta _ { 2 } ~ \theta _ { 3 }\)

Then we output the results for force/monent variables. For each line, there are 6 values arranged as \(F _ { 1 } ~ F _ { 2 } ~ F _ { 3 } ~ M _ { 1 } ~ M _ { 2 } ~ M _ { 3 }\)

For dynamic analysis, we also output the linear and angular momenta for each element and there are 6 values values arranged as \(P _ { 1 } P _ { 2 } P _ { 3 } H _ { 1 } H _ { 2 } H _ { 3 }\)

For eigenvalue analysis, GEBT outputs the steady state solution first. Then it outputs the eigenvalues (in Hz) and eigenvectors corresponding to this state. The eigenvalues are listed from smallest to the largest in magnitude.

It is noted that the displacements and rotations are expressed in the global coordinate system (ai) while the forces/moments and momenta are expressed in the deformed beam coordinate system (Bi) for each element. The forces/moments are expressed in the global coordinate system (ai) for each boundary point. As the forces/moments of each connection point could be different for different elements it is associated with, which can be obtained from the elements connected by this point, such results are not reported but the places are replaced with zeroes. Note the calculated forces and moments at a boundary point are the internal forces and moments evaluated at this point. If the point is the ending point of the member, the internal forces/moments will be equal to applied forces/moments. If the point is the starting point of the member, the magnitude of internal forces/moments will be equal to applied forces/moments but with a different sign. For boundary points with applied forces/moments, the prescribed values are directly copied in the output file. The output file is in pure text format and can be opened by any text editor.

