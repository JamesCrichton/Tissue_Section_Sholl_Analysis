'''
Script to segment cells with complex morphology from tissue sections, using QuPath. 
This involves three main steps:
1. Nuclear segmentation (instance segmentation)
2. Segmentation of cell body (semantic segmentation) as a single mask using a trained px classifier
3. Splitting of the cell body mask using associated nuclei as seeds, to generate cell objects 

Requires:
Cellpose extension
        
Author: James Crichton, RILD  Bioimaging, University of Exeter
Written for Alex Mellor, Kate Ellacott lab, University of Exeter
'''

////////// USER INPUT //////////////////////////////////////////

def pathModel = "//phineas/argus/Kate Ellacott Lab/Alex Mellor/JC Analysis Support/QP Analysis/models/Custom_model_2026-03-26_14_57.cpm"  // Bespoke cellpose model path. 

def PX_CLASSIFIER = "Iba1 GFAP pixel classifier"  // Name of trained px_classifier saved within the QuPath Project
def MIN_AREA = 40.0  // minimum area for annotation detected
def MIN_HOLE = 3.0  // minimum hole allowed in annotation
//def CELL_CLASS = "Iba1 - Microglia"  // String name of the annotations generated from px classifier 
def IBA1_MASK = "Iba1 - Microglia"  // String name of the Iba1 annotations generated from px classifier 
def GFAP_MASK = "GFAP- Astrocyte"  // String name of the GFAP annotations generated from px classifier 
def OVERLAP_THRESHOLD = 20  // Minimum % area for nucleus to overlap with cell mask and be included 
def CIRCULARITY_THRESHOLD = 0.5 // Minumum nuclear circularity

////////////////////////////////////////////////////////////////


////////// IMPORTS /////////////////////////////////////////////

import qupath.lib.analysis.DelaunayTools
import qupath.lib.objects.PathObjects
import qupath.ext.biop.cellpose.Cellpose2D
import qupath.lib.scripting.QP
import qupath.lib.roi.GeometryTools

////////////////////////////////////////////////////////////////



////////// STEP 1 - NUCLEAR SEG/////////////////////////////////

removeObjects(getAllObjects()) // clear objects
createFullImageAnnotation(true) //make whole image annotaiton
frame_obj = getAnnotationObjects()[0]
    
// Segment nuclei with Cellpose
selectObjects(frame_obj)

def cellpose = Cellpose2D.builder( pathModel )
        .pixelSize( 0.25 )                      // Resolution for detection in um
        .channels( 'Channel 1' )	               // Select detection channel(s)
        .flowThreshold( 0.4 )              // Threshold for the flows, defaults to 0.4
        .diameter( 20 )                    // Median object diameter. Set to 0.0 for the `bact_omni` model or for automatic computation
        .build()

// Run detection for the selected objects
def imageData = getCurrentImageData()
def pathObjects = getSelectedObjects() // To process only selected annotations, useful while testing

cellpose.detectObjects( imageData, pathObjects )

// You could do some post-processing here, e.g. to remove objects that are too small, but it is usually better to
// do this in a separate script so you can see the results before deleting anything.

println 'Cellpose detection script done'


////////////////////////////////////////////////////////////////


////////// STEP 2 - CELL BODY SEG //////////////////////////////

// Segment cell body with px classifier
selectObjects(frame_obj)
createAnnotationsFromPixelClassifier(PX_CLASSIFIER, MIN_AREA, MIN_HOLE)


// Add nuclei shape features
selectObjects(getDetectionObjects())
addShapeMeasurements("AREA", "LENGTH", "CIRCULARITY")


// Remove low circularity nuclei
low_circularity_detections = getDetectionObjects().findAll{it.measurements.get("Circularity") <= CIRCULARITY_THRESHOLD}
removeObjects(low_circularity_detections)


// Assign parent objects to nuclei (detection objects default to those with an overlapping centroid, we want anything which overlaps above a threshold)
// Merge annotations for each mask

// Percentage and area overlap calculation for each nucleus
GFAP_collection = getAnnotationObjects().findAll{it.getPathClass() == getPathClass(GFAP_MASK)}
GFAP_geom = GFAP_collection[0].getROI().getGeometry()

IBA1_collection = getAnnotationObjects().findAll{it.getPathClass() == getPathClass(IBA1_MASK)}
IBA1_geom = IBA1_collection[0].getROI().getGeometry()

getDetectionObjects().each{ nuc -> 
    nuc_geom = nuc.getROI().getGeometry()
    percentage_GFAP = (nuc_geom.intersection(GFAP_geom).getArea() / nuc_geom.getArea()) * 100
    percentage_IBA1 = (nuc_geom.intersection(IBA1_geom).getArea() / nuc_geom.getArea()) * 100
    nuc.measurements.put("Percentage_GFAP", percentage_GFAP)
    nuc.measurements.put("Percentage_IBA1", percentage_IBA1)
    }   
 
 
// Remove nuclei with overlap below threshold
no_overlap_nuclei = getDetectionObjects().findAll{(it.measurements.get("Percentage_GFAP") < OVERLAP_THRESHOLD) && (it.measurements.get("Percentage_IBA1") < OVERLAP_THRESHOLD)}

removeObjects(no_overlap_nuclei)
resetSelection()


//// Add the nuclei to their most overlapping cell mask
// 
getDetectionObjects().each{ nuc ->
    if (nuc.measurements.get("Percentage_IBA1") > nuc.measurements.get("Percentage_GFAP")) {
        IBA1_geom = nuc.getROI().getGeometry().union(IBA1_geom)
    } else {
        GFAP_geom = nuc.getROI().getGeometry().union(GFAP_geom)
    }
}

removeObjects(IBA1_collection)
removeObjects(GFAP_collection)

IBA1_roi = GeometryTools.geometryToROI(IBA1_geom)
GFAP_roi = GeometryTools.geometryToROI(GFAP_geom)

IBA1_anno = PathObjects.createAnnotationObject(IBA1_roi, getPathClass(IBA1_MASK))
GFAP_anno = PathObjects.createAnnotationObject(GFAP_roi, getPathClass(GFAP_MASK))

addObjects(IBA1_anno)
addObjects(GFAP_anno)


//// Remove isolated regions of mask
selectObjectsByClassification(GFAP_MASK)
runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin', '{}')

selectObjectsByClassification(IBA1_MASK)
runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin', '{}')

resolveHierarchy()

non_cell_annotations = getAnnotationObjects().findAll{!it.hasChildObjects()}
removeObjects(non_cell_annotations)

    
     
////////////////////////////////////////////////////////////////

////////// STEP 3 - SPLIT CELL BODY BY ASSOCIATED NUCLEI ///////

// Make a function to split the cell body of any names cell class
def cell_split(CELL_CLASS){
'''
Function will split objects with the named classification, using a Voronoi diagram with segmented nuclei as the seed points. 
Remaining objects will become defined as cell bodies. 
The purpose of this is to separate merged cells. 

CELL_CLASS = str. Name of the annotation class being split. 

'''

    selectObjectsByClassification(CELL_CLASS)
    mergeSelectedAnnotations()    
    resolveHierarchy() 
    cell_mask = getAnnotationObjects().findAll {it.getPathClass() == getPathClass(CELL_CLASS)}
    
    
    // Segment cell body using Voronoi
    /*Script to create Voronoi Diagrams from annotations and limit the expansion to a reference Annotation.
     * This script also counts the number of touching neighbors and puts this value into the Voronoi annotations measurements. 
     * @author Isaac Vieco-Martí
     */
    
    // select detections (nuclei) to make Voronoi Diagram
    def nuc_detections = getDetectionObjects().findAll{it.getParent().getPathClass() == getPathClass(CELL_CLASS)}
    
    // annotation to limit the expansion
    def reference = getAnnotationObjects().findAll {it.getPathClass() == getPathClass(CELL_CLASS)}
    def referenceGeom = reference[0].getROI().getGeometry()
    def plane = reference[0].getROI().getImagePlane()
    
    
    // Get the Voronoi Faces, it is a map with the original annotation and the geometry of the Voronoi Polygon
    def voronoiDiagram = DelaunayTools.createFromGeometryCoordinates(nuc_detections,false,4.0).getVoronoiFaces()
    
    // Create Voronoi annotations, limit the expansion and set the class of its origin annotation
    voronoiDiagram = voronoiDiagram.collectEntries { annotationObj, voronoiFace ->
        intersect = voronoiFace.intersection(referenceGeom)
        roi = GeometryTools.geometryToROI(intersect, plane)
        [(annotationObj): PathObjects.createAnnotationObject(roi, getPathClass(CELL_CLASS))]
    } 
    
    
    // Get the Voronoi annotations
    def finalAnnotations = voronoiDiagram.values()
    
    addObjects(finalAnnotations)
    //////////
    
    
    // Tidy up annotations:
    // Remove the original annotations
    removeObjects(cell_mask)
    
    
    //Split all objects. Some are composed of multiple pieces currently. Splitting helps to tidy them up
    selectObjects(getAllObjects())
    runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin', '{}')
     
    
    //Fill holes in cell annotations. This simplifies their shape, so circularity can be quantified
    getAnnotationObjects().each{anno -> 
        selectObjects(anno)
        runPlugin('qupath.lib.plugins.objects.FillAnnotationHolesPlugin', '{}')
        }
    
    
    // Remove annotations with no children
    resolveHierarchy()
    
    def anno_collection = getAnnotationObjects().findAll{it.getPathClass() == getPathClass(CELL_CLASS)}
    no_child_collection = anno_collection.findAll{!it.hasChildObjects()}
    removeObjects(no_child_collection)
        
      
    // Make the detections into new cells. Remove any without parents (can arise from Voronoi artefacts occasionally)
    for (nucleus in nuc_detections) {        
                
        parent_anno = nucleus.getParent()
        
        if (parent_anno.getPathClass() == getPathClass(CELL_CLASS)){
            nucleus_geom = nucleus.getROI().getGeometry()
            cytoplasm_geom = parent_anno.getROI().getGeometry().union(nucleus_geom)  // Modify the cytoplasm geom to ensure it fully contains its associated nucleus
            cytoplasm_ROI = GeometryTools.geometryToROI(cytoplasm_geom, plane)
                
            cell = PathObjects.createCellObject(cytoplasm_ROI, nucleus.getROI(), getPathClass(CELL_CLASS)) //Set cell nucleus, cytoplasm and parent        
            addObjects(cell)
            removeObjects(nucleus)
            removeObjects(parent_anno)
            } else {
                removeObjects(nucleus)    
            }
        }
    }



// Split merged cells 
cell_split(IBA1_MASK)
cell_split(GFAP_MASK)


// Remove any final lone annotations
lone_annotations = getAnnotationObjects().findAll{!it.hasChildObjects()}
removeObjects(lone_annotations)


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

resetSelection()
