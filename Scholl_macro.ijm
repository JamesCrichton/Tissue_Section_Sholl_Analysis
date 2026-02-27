
//Add cell outline to ROI Manager
roiManager("reset");
roiManager("Add");

//Make new image as cell mask
width = getWidth()  
height= getHeight()  

newImage("cell_mask", "8-bit black", width , height, 1);

roiManager("Select", 0);

setForegroundColor(0, 204, 255);
run("Fill", "slice");

run("None");

// Scholl

run("Convert to Mask");
run("Skeletonize");

