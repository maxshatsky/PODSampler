Hello! PODSampler collects data at the selected surface inside the mesh.
First of all you need to make a surface to collect POD-data on.

0) go to system and open ControlDict
1) set endTime 0.1 sec
2) uncomment ControlDict line:  libs ("libsampledEnsightSurfaceMesh.so");
3) uncomment ControlDict line:  #include "samplePlane"
4) if you dare, feel free to change surface parameters in "system/samplePlane" file
5) run laplacianFoam

You made a surface "postProcessing/makePlane/s1/s1.case".
You can open and inspect it via paraView.
Now you can collect POD-data at the surface.

6) set endTime 65 secs
7) comment ControlDict line  #include "samplePlane"
8) uncomment ControlDict line  #include "functionPOD_T"
9) if you dare, feel free to change PODSampler settings in "system/funtionPOD_Global" file
10) run laplacianFoam

You have collected POD-data at the surface. Congrats.
