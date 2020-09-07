This is a POD function object directory.
The function object works with OpenFOAM-v1812.

Current version of POD is able to read any field on the given surface every <timeIndexPOD>-th iteration
and create a covariance matrix, acoeffs, timelist, basis and map of coordinates of this data after <nSnapshots> interpolations.

To let the function object work User should:
- run "wmake" in current directory.
- run "wmake" in thirdparty/sampleEnsightSurfaceMesh/ .
- include libs ("libsampledEnsightSurfaceMesh.so") to controlDict file.
- include functionPOD_Global instructions in the end of each functionPOD_* file.
- include functionPOD_* instructions int the end of controlDict file for each scalar field like:
	functions
	{
	        #include "functionUcmpts"
		#include "functionPOD_T"
		#include "functionPOD_p"
		#include "functionPOD_Ux"
		#include "functionPOD_Uy"
		#include "functionPOD_Uz"
	}
- Add interpolation surface (*.stl) in directory <casefoler>/constant/triSurface/ .
- Check surface directory to be set in functionPOD_Global file.
- Check settings in controlDict file.
- Double-check settings in conrolDict file.
- Also read comments in functionPOD_Global file and make preset to evade mistakes.

Good luck!
