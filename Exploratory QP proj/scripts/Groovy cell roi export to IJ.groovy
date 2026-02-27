'''
Groovy script to export pre-segmented cell detection shapes (cell and nucleus) as ImageJ ROIS for downstream analysis. 
Initially intended for Sholl analysis using IJ Neuoranatomy plugin as a searate script. 
NB. Neuoroanatomy plugin cannot currently be run directly from QuPath. 
James Crichton, RILD, University of Exeter. 
Written for Alex Mellor: Kate Ellacott lab (RILD)
'''

// Import necessary packages
import ij.ImageJ
import qupath.lib.scripting.QP
import qupath.lib.gui.scripting.QPEx
import qupath.imagej.tools.IJTools
import ij.plugin.frame.RoiManager
import ij.gui.Roi
import qupath.imagej.gui.IJExtension
import qupath.lib.regions.RegionRequest
import ij.io.FileSaver


// Set image constants
def viewer = QPEx.getCurrentViewer()
def server = viewer.getImageData().getServer()


// Loop through cell detections
getCellObjects().each{anno ->
   
    // Get cell shape data
    cellRoi = anno.getROI()
    annoX = cellRoi.boundsX
    annoY = cellRoi.boundsY
    nucRoi = anno.getNucleusROI()
    

    // Copy the cell image into IJ 
    region = RegionRequest.createInstance(server.getPath(), 1, cellRoi);
    
    def imp = IJExtension.extractROIWithOverlay(server, anno, null, region, false, null).getImage();
    
    ImageJ ij = IJExtension.getImageJInstance();
    if (ij != null) {
        ij.setVisible(true);
        imp.show()
    }
    
    RoiManager roiMan = RoiManager.getRoiManager() 
    roiMan.setVisible(true)
    roiMan.reset()
    
    // Add cell ROI
    cellX = cellRoi.boundsX - annoX
    cellY = cellRoi.boundsY - annoY
    Roi cellRoiIJ = IJTools.convertToIJRoi(cellRoi, 0, 0, 1)
    cellRoiIJ.setLocation(cellX, cellY)
    roiMan.addRoi(cellRoiIJ)
    
   
    // Add nuclear ROI
    nucX = nucRoi.boundsX - annoX
    nucY = nucRoi.boundsY - annoY
    Roi nucRoiIJ = IJTools.convertToIJRoi(nucRoi, 0, 0, 1)
    nucRoiIJ.setLocation(nucX, nucY)
    roiMan.addRoi(nucRoiIJ)
          
    // Save and close
    def cell_ID = anno.getID().toString()
    def pathExportForSholl = buildFilePath(PROJECT_BASE_DIR, 'exported_for_sholl')
    def imageName = getCurrentImageData().getServer().getMetadata().getName()
    def cellExportPath = buildFilePath(pathExportForSholl, imageName, "cell_"+cell_ID)
    mkdirs(cellExportPath) // make a dir for the cell data
    
    def roiFilePath = buildFilePath(cellExportPath, 'rois.zip')
    def imgFilePath = buildFilePath(cellExportPath, 'cell_img.tiff')
    
    roiMan.runCommand("Save", roiFilePath)
    roiMan.reset()
    
    FileSaver fs = new FileSaver(imp)
    fs.saveAsTiff(imgFilePath)
    imp.close()
    }

println('Done!')