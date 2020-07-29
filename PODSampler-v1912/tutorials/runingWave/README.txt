You need to make a surface to collect POD-data on.

0) go to system and open ControlDict
1) set endTime 0.1 sec
2) uncomment ControlDict line  libs ("libsampledEnsightSurfaceMesh.so");
3) uncomment ControlDict line  #include "samplePlane"
4) run laplacianFoam

You made a surface. Surface is located at "postProcessing/makePlane/s1/s1.case". You can open and inspect it via paraView.

Now you can collect POD-data at the surface.

5) set endTime 65 secs
6) comment ControlDict line  #include "samplePlane"
7) uncomment ControlDict line  #include "functionPOD_T"
8) run laplacianFoam

You have collected POD-data at the surface. Congrats!
