'''
Script to segment cells with complex morphology from tissue sections, using QuPath. 
This involves three main steps:
1. Nuclear segmentation (instance segmentation)
2. Segmentation of cell body (semantic segmentation) as a single mask using a trained px classifier
3. Splitting of the cell body mask using associated nuclei as seeds, to generate cell objects 

Requires:
StarDist QuPath Extension
StarDist model (https://github.com/qupath/models/tree/main/stardist)
        
James Crichton, RILD  Bioimaging, University of Exeter
Written for Alex Mellor, Kate Ellacott lab, University of Exeter
'''

////////// USER INPUT //////////////////////////////////////////

def modelPath = "C:/Users/spoonbill/Documents/Image_Analysis_Tools/Stardist_models/dsb2018_heavy_augment.pb"  // Stardist model path. Can download from  https://github.com/qupath/models/tree/main/stardist)
def PX_CLASSIFIER = "Iba1 px classifier"  // Name of trained px_classifier saved within the QuPath Project
def MIN_AREA = 40.0  // minimum area for annotation detected
def MIN_HOLE = 3.0  // minimum hole allowed in annotation
def CELL_CLASS = "Iba1 - Microglia"  // String name of the annotations generated from px classifier 
def overlap_threshold = 80  // Minimum % area for nucleus to overlap with cell mask and be included 
def circularity_threshold = 0.6 // Minumum nuclear circularity

////////////////////////////////////////////////////////////////


////////// IMPORTS /////////////////////////////////////////////

import qupath.lib.analysis.DelaunayTools
import qupath.lib.objects.PathObjects
import qupath.ext.stardist.StarDist2D
import qupath.lib.scripting.QP

////////////////////////////////////////////////////////////////



////////// STEP 1 - NUCLEAR SEG/////////////////////////////////

removeObjects(getAllObjects()) // clear objects
createFullImageAnnotation(true) //make whole image annotaiton
frame_obj = getAnnotationObjects()[0]
    
// Segment nuclei with Stardist 
selectObjects(frame_obj)

def stardist = StarDist2D
    .builder(modelPath)
    .channels('Channel 1')            // Extract channel called 'DAPI'
    .normalizePercentiles(1, 99) // Percentile normalization
    .threshold(0.5)              // Probability (detection) threshold
    .pixelSize(0.5)              // Resolution for detection
    .cellExpansion(0)            // Expand nuclei to approximate cell boundaries
    .measureShape()              // Add shape measurements
    .measureIntensity()          // Add cell measurements (in all compartments)
    .build()
	
def pathObjects = QP.getSelectedObjects()

// Run detection for the selected objects
def imageData = QP.getCurrentImageData()
stardist.detectObjects(imageData, pathObjects)
stardist.close() // This can help clean up & regain memory

////////////////////////////////////////////////////////////////


////////// STEP 2 - CELL BODY SEG //////////////////////////////

// Segment cell body with px classifier
selectObjects(frame_obj)
createAnnotationsFromPixelClassifier(PX_CLASSIFIER, MIN_AREA, MIN_HOLE)


// Associate segmented mask with individual nuclei
resolveHierarchy()


// Remove cells not touching the cell mask
mask_collection = getAnnotationObjects().findAll{it.getPathClass() == getPathClass(CELL_CLASS)}
mask_obj = mask_collection[0]
mask_children = mask_obj.getChildObjects()

non_mask_cells = getDetectionObjects()-mask_children
removeObjects(non_mask_cells)


// Add nuclei shape features
selectObjects(getDetectionObjects())
addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY")


// Remove nuclei which are <threshold% covered by mask
mask_geom = mask_obj.getROI().getGeometry()

mask_cells = getDetectionObjects().collect{cell->
    cell_geom = cell.getROI().getGeometry()
    total_area = cell_geom.getArea()
    cell_mask_overlap = cell_geom.intersection(mask_geom).getArea()

    percentage_overlap = (cell_mask_overlap/total_area)*100
    print("Percentage overlap = " + percentage_overlap +"%")
          
    if(percentage_overlap>= overlap_threshold) {
        cell_circularity = cell.measurements.get("Circularity")
        print("circularity = "+ cell_circularity.toString())
        
        if(cell_circularity >=circularity_threshold) {
            return cell 
        }}}     
    
    
//remove non-mask cells
non_mask_cells_B = getDetectionObjects()-mask_cells
removeObjects(non_mask_cells_B)
resetSelection()


// Remove Iba1 not touching a nucleus
selectAllObjects()
runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin', '{}')

resolveHierarchy()

getAnnotationObjects().each{anno->
    if (anno.hasChildObjects()==false) {
        removeObject(anno)}
    }
 
 
////////////////////////////////////////////////////////////////


////////// STEP 3 - SPLIT CELL BODY BY ASSOCIATED NUCLEI ///////

selectObjectsByClassification(CELL_CLASS)
mergeSelectedAnnotations()
cell_mask = getAnnotationObjects().findAll {it.getPathClass() == getPathClass(CELL_CLASS)}


// Segment cell body using Voronoi
/*Script to create Voronoi Diagrams from annotations and limit the expansion to a reference Annotation.
 * This script also counts the number of touching neighbors and puts this value into the Voronoi annotations measurements. 
 * @author Isaac Vieco-Martí
 */


//select detections (nuclei) to do Voronoi Diagrams
def annoInterest = getDetectionObjects()

//annotation to limit the expansion
def reference = getAnnotationObjects().findAll {it.getPathClass() == getPathClass(CELL_CLASS)}
def referenceGeom = reference[0].getROI().getGeometry()
def plane = reference[0].getROI().getImagePlane()


// Get the Voronoi Faces, it is a map with the original annotation and the geometry of the Voronoi Polygon
def voronoiDiagram = DelaunayTools.createFromGeometryCoordinates(annoInterest,false,4.0).getVoronoiFaces()

//Create Voronoi annotations, limit the expansion and set the class of its origin annotation
voronoiDiagram = voronoiDiagram.collectEntries { annotationObj, voronoiFace ->
    objClass = annotationObj.getPathClass()
    intersect = voronoiFace.intersection(referenceGeom)
    roi = GeometryTools.geometryToROI(intersect, plane)
    [(annotationObj): PathObjects.createAnnotationObject(roi,objClass)]
} 


// Get the Voronoi annotations
def finalAnnotations = voronoiDiagram.values()

addObjects(finalAnnotations)

// Tidy up annotations:
// Remove the original annotations
removeObjects(cell_mask)

//Split all objects. Some are composed of multiple pieces currently. Splitting helps to tidy them up
selectObjects(getAllObjects())
runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin', '{}')

//Remove objects below size threshold
sizethreshold = 300 //currently in px as images are not scaled
small_anno=getAnnotationObjects().findAll {it.getROI().getGeometry().getArea()< sizethreshold}
removeObjects(small_anno)


//Fill holes in cell annotations. This simplifies their shape, so circularity can be quantified
getAnnotationObjects().each{anno -> 
    selectObjects(anno)
    runPlugin('qupath.lib.plugins.objects.FillAnnotationHolesPlugin', '{}')
    }
    
//Resolve hierarchy (parent/child relationships)
resolveHierarchy()


//Remove annotations without 1 detection i.e. their nucleus
empty_anno = getAnnotationObjects().findAll{it.nDescendants()!=1}
removeObjects(empty_anno)


// Classify detections to cells
// Make cells
getAnnotationObjects().each{anno -> anno.setName(anno.getID().toString())} //make IDs annotation names

def detections = getDetectionObjects()

for (nucleus in detections) {
    parent_ID = nucleus.getParent().toString()
    parent_ID = parent_ID.split(" ")[0] //remove extra characters added to string
    cytoplasm = getAnnotationObjects().findAll{it.getName()==parent_ID}  
    
    //Has cytoplasm been found? If so, make a cell
    if (cytoplasm.size()>0){
        cell = PathObjects.createCellObject(cytoplasm[0].getROI(), nucleus.getROI(), null) //Set cell nucleus, cytoplasm and parent        
        addObjects(cell)
        removeObjects(nucleus)
        removeObjects(cytoplasm)
    }

}

////////////////////////////////////////////////////////////////


////////// STEP 4 - ADD CALCULATIONS ///////////////////////////

// Calculate ramification 
selectObjects(getCellObjects())
addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY", "SOLIDITY", "MAX_DIAMETER", "MIN_DIAMETER", "NUCLEUS_CELL_RATIO")


getCellObjects().each{cell ->
    m_list=cell.getMeasurementList()
    cell_area = m_list.get("Cell: Area µm^2")
    cell_perimeter = m_list.get("Cell: Length µm")
    ratio = cell_area/cell_perimeter
    cell.measurements.put("Area:Perimeter Ratio", ratio)
    }

