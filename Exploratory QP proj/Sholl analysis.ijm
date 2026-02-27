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
 
run("Fresh Start");

/////////////////////////////////////////////////////////////////
//////////// STEP 1 - FIND FILE PATHS ///////////////////////////

dataset_dir = getDirectory("Choose a Directory ");  
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


/////////////////////////////////////////////////////////////////
//////////// STEP 2 - OPEN FILES AND ROIS ///////////////////////

// Open image and ROIs
img_path = cell_dirs[1] + "cell_img.tiff"
roi_path = cell_dirs[1] + "rois.zip"
open(img_path);
roiManager("open", roi_path);


// Make masks for Sholl
width = getWidth();
height = getHeight();

newImage("cell_mask", "8-bit", width, height, 1);
roiManager("Select", 0);
run("Clear Outside");
run("Select None");

newImage("nucleus_mask", "8-bit", width, height, 1);
roiManager("Select", 1);
run("Clear Outside");
run("Select None");

// Make cell skeleton for Sholl
selectWindow("cell_mask");
run("Skeletonize");

/////////////////////////////////////////////////////////////////
//////////// FUNCTIONS //////////////////////////////////////////

// Get an array of files recrsively identified from a dir
function listFiles(dir) {
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