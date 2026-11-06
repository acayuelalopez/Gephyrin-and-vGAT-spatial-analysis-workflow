import qupath.lib.objects.PathObjects
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import qupath.lib.io.PathIO
import qupath.lib.gui.scripting.QPEx
import qupath.ext.biop.cellpose.Cellpose2D


// Get nucleus diameter in microns
def diameterNuc = 1


// Get the current image data
def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
// Get the image name
def imageName = imageData.getServer().getMetadata().getName()
// Get the project path
def projectPath =  QPEx.getProject().getPath()
// Get the pixel size in microns
def pixelSize = server.getPixelCalibration().getPixelWidthMicrons()

// Get all annotations
def annotations = getAnnotationObjects()

// Filter by name (case-insensitive)
def target = annotations.findAll { it.getPathClass()?.toString() == "Cortex_ipsilateral" }

// Select them in the viewer
selectObjects(target)

print "Selected ${target.size()} annotation(s)."

def pathModel = 'model_telomer.901241'
def cellpose = Cellpose2D.builder(pathModel)
        .pixelSize(pixelSize)                  // Resolution for detection in um
         .channels( 'Channel 2' )          // Overwrites the logic of this plugin with these two values. These will be sent directly to --chan and --chan2
        .diameter(diameterNuc/pixelSize)                    // Median object diameter. Set to 0.0 for the `bact_omni` model or for automatic computation
        .measureIntensity() 
        .flowThreshold( 0.6)
        //.cellExpansion(0) // Add cell measurements (in all compartments)
        .cellprobThreshold(-3.0)
        .build()

// Run detection for the selected objects
def pathObjects = QPEx.getSelectedObjects() // To process only selected annotations, useful while testing
if (pathObjects.isEmpty()) {
    Dialogs.showErrorMessage("Cellpose", "Please select a parent object!")
    return
}



cellpose.detectObjects( imageData, pathObjects )

saveAnnotationMeasurements(projectPath.toString().replaceAll("project.qpproj","output")+File.separator+imageName.toString().replaceAll("/","_")+"_Gephyrin_ipsilateral_summary.csv")

detectionToAnnotationDistances(true)

saveDetectionMeasurements(projectPath.toString().replaceAll("project.qpproj","output")+File.separator+imageName.toString().replaceAll("/","_")+"_Gephyrin_ipsilateral_single_cell.csv")

// You could do some post-processing here, e.g. to remove objects that are too small, but it is usually better to
// do this in a separate script so you can see the results before deleting anything.

println 'Cellpose detection script done'

import qupath.ext.biop.cellpose.Cellpose2D


