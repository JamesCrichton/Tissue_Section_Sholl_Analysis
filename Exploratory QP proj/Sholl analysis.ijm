/*
 * ImageJ Macro for Sholl Analysis
 * Script takes existing ROIs (nuclear and cytoplasmic) and corresponding cell images. 
 * Expected path structure:
 * dataset_dir > img_dir > cell_dir > cell_img.tiff, rois.zip
 * 
 * Requires:
 * ImageJ Neuoroanatomy plugin
 * 
 * Written by James Crichton, RILD Bioimaging, University of Exeter
 * For Alex Mellor, Kate Ellacott lab, University of Exeter
 */
 

/////////////////////////////////////////////////////////////////
//////////// STEP 1 - FIND FILE PATHS ///////////////////////////

//dataset_dir = getDirectory("Choose a Directory ");  
#@ File (label="Choose dataset directory", style = "directory") dataset_dir
#@ Boolean (label="Run in background?") background_selection

setBatchMode(background_selection); 
files = listFiles(dataset_dir);  //get full file paths

// Get dirs containing both a tiff and roi zip. 
// Make list of dirs containing rois or images, then compare. Keep the overlapwith both files

roi_dir_list = newArray();
cell_img_dir_list = newArray();

for (i = 0; i < files.length; i++) {
	file = files[i];
	if (endsWith(file, ".zip"))
		roi_dir_list = Array.concat(roi_dir_list, File.getDirectory(file));
	if (endsWith(file, ".tiff")) 
		cell_img_dir_list = Array.concat(cell_img_dir_list, File.getDirectory(file));
}

cell_dirs = getOverlap(roi_dir_list, cell_img_dir_list);

print(dataset_dir);
Array.show(getFileList(dataset_dir));
print(files.length);
print(cell_dirs.length + " cells found"); 


/////////////////////////////////////////////////////////////////
//////////// STEP 2 - OPEN FILES AND RUN SHOLL //////////////////

// Set colour variables
run("Colors...", "foreground=white background=black selection=orange");
setBackgroundColor(0, 0, 0);
setForegroundColor(255, 255, 255);


// Loop through files
for (i = 0; i < cell_dirs.length; i++) {
	
	run("Fresh Start");
	roiManager("reset");
	
	// Open image and ROIs
	img_path = cell_dirs[i] + "cell_img.tiff";
	roi_path = cell_dirs[i] + "rois.zip";
	open(img_path);
	roiManager("open", roi_path);
	
	
	// Make masks for Sholl
	selectWindow("cell_img.tiff");	
	width = getWidth();
	height = getHeight();
	
	// Remove ROIs > index 2
	if (roiManager("count")>2){
		for (j = 2; j < roiManager("count"); j++) {
			roiManager("select", j);
    		roiManager("delete");
		}
	}
	
	newImage("cell_mask", "8-bit black", width, height, 1);
	roiManager("Select", 0);
	run("Fill");
	run("Select None");
	
//	newImage("nucleus_mask", "8-bit", width, height, 1);
//	roiManager("Select", 1);
//	run("Fill");
//	run("Select None");
	
	
	// Make cell skeleton for Sholl
	selectWindow("cell_mask");
	run("Skeletonize");
	rename("cell_skeleton");
	
	selectImage("cell_skeleton"); 
	run("Create Selection");
	roiManager("Add");  // Add to ROI Manager
	run("Select None");

	selectImage("cell_skeleton");
	roiManager("select", 0); // This will set the cell ROI centroid as the centre-point for Sholl rings
	

	// Run Sholl:	
	// NB: have a look at the example scripts in Templates>Neuroanatomy> for more
	// robust ways to automate Sholl. E.g., Sholl_Extract_Profile_From_Image_Demo.py
	// exemplifies how to parse an image programmatically using API calls
	run("Sholl Analysis (From Image)...", "datamodechoice=Intersections startradius=0.0 stepsize=10.0 endradius=52.79549876399019 hemishellchoice=[None. Use full shells] previewshells=false nspans=1 nspansintchoice=N/A primarybrancheschoice=[Infer from starting radius] primarybranches=0 polynomialchoice=['Best fitting' degree] polynomialdegree=0 normalizationmethoddescription=[Automatically choose] normalizerdescription=Default plotoutputdescription=[Linear plot] tableoutputdescription=[Detailed table] annotationsdescription=[ROIs (points and 2D shells)] lutchoice=mpl-viridis.lut luttable=net.imglib2.display.ColorTable8@61827ed3 save=true savedir=["+cell_dirs[i]+"] analysisaction=[Analyze image]");
	

	// Save and close
	saveAs("tiff", cell_dirs[i]+"skeleton.tiff");  // save skeleton img
	roiManager("save", roi_path);  // Save ROI Manager with skeleton ROI
	close("*");
		
	print("Processed cell " + (i+1) + "of " + cell_dirs.length);
}


/////////////////////////////////////////////////////////////////
//////////// FUNCTIONS //////////////////////////////////////////

// Get an array of files recrsively identified from a dir
function listFiles(dir) {
    dir = dir + File.separator;
    output_list = newArray();
    list = getFileList(dir);

    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) {
            sub = listFiles(dir + list[i]);
            output_list = Array.concat(output_list, sub);
        }
        else {
            path = dir + list[i];
            output_list = Array.concat(output_list, newArray(path));
        }
    }
    return output_list;
}

// Compare two arrays and retun the overlap
function getOverlap(a, b) {
    overlap = newArray();

    for (i = 0; i < a.length; i++) {
        for (j = 0; j < b.length; j++) {
            if (a[i] == b[j]) {
                overlap = Array.concat(overlap, newArray(a[i]));
                break;  // stop checking once matched
            }
        }
    }
    return overlap;
}