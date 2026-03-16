'''
Groovy script to read skeleton ImageJ ROIs derived from QuPath cell anotations back into QuPath for quality control. 
Also adding Sholl metrics from ImageJ analysis to cell data.

James Crichton, RILD, University of Exeter. 
Written for Alex Mellor: Kate Ellacott lab (RILD)
'''


// Import necessary packages
import qupath.imagej.tools.IJTools
import qupath.lib.objects.PathObjects
import qupath.lib.measurements.MeasurementList

def plane = ImagePlane.getDefaultPlane()

// Find cell data 
def pathExportForSholl = buildFilePath(PROJECT_BASE_DIR, 'exported_for_sholl')
def imageName = getCurrentImageData().getServer().getMetadata().getName()

getCellObjects().each { cell -> 
    def cell_name = cell.getName().toString()
    def cellExportPath = buildFilePath(pathExportForSholl, imageName, cell_name)
    def roiFilePath = buildFilePath(cellExportPath, 'rois.zip')
    def measurementsFilePath = buildFilePath(cellExportPath, 'cell_skeleton_Sholl-Profiles.csv')
    
    
    // Import Skeleton as an annotation
    def file = new File(roiFilePath)
    def ijRois = IJTools.readImageJRois(file)   
    def skeleton_ROI = ijRois.toArray().findAll{r -> r.getName().toString() == "skeleton"}[0]  // Get the first roi caled "skeleton" in the collection for this cell
    
    double downsample = 1.0
    double xOrigin = -cell.getROI().getBoundsX()
    double yOrigin = -cell.getROI().getBoundsY()
       
    def skeletonObject = IJTools.convertToAnnotation(skeleton_ROI, xOrigin, yOrigin, downsample, plane) 
    
    def skeleton_name = cell_name.replace("cell", "skeleton")
    skeletonObject.setName(skeleton_name)
    skeletonObject.setPathClass(getPathClass("Skeleton"))
    
    addObject(skeletonObject)
    
/////////////////    
// From the sholl data add the log intersection/area for each radius, as a measurement to the cell 

    // Read in the sholl measurements table
    def lines = new File(measurementsFilePath).readLines()
    def header = lines[0].split(",").toList()
    
    // Find the col positions
    radius_col = header.indexOf("Radius")
    value_col = header.indexOf("log(Inters./Area)")
    
    
    for (int i = 1; i < lines.size(); i++) {
        radius =  lines[i].split(",").toList()[radius_col]
        measurement_name = "Sholl radius "+ radius.toString() +": log(inters/area)"
        value = lines[i].split(",").toList()[value_col].toDouble()
        cell.measurements.put(measurement_name, value)  // add measurements to cell        
        }    
        
    }
    




