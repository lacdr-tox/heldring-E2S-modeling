input_ROIs = getDirectory("Please select directory where the ROI sets are saved");
output = getDirectory("Please select directory where the output shall be saved");

File.makeDirectory(output + "/blackmap");
output_black_map = output + "/blackmap/";
File.makeDirectory(output + "/roi_subsets");
output_roisubs = output + "/roi_subsets/";
File.makeDirectory(output + "/greyscale");
output_greyscale = output + "/greyscale/";
File.makeDirectory(output + "/cropped_greyscale");
output_cropped = output + "/cropped_greyscale/";

setBackgroundColor(0, 0, 0);
setForegroundColor(255, 255, 255);
setBatchMode(true);

// Generate binary maps from ROIs
generate_black_map ();

list_ROI_files = getFileList(input_ROIs);
Array.sort(list_ROI_files);

//print("Start making quadrants...");
for (i = 0; i < list_ROI_files.length; i++) {
	count = IJ.pad(i+1,4);
	// quad 1
	ftmp = select_corner (list_ROI_files[i],1,"0(25[0-5]|2[0-4][0-9]|[0-1][0-9][0-9])-0(25[0-5]|2[0-4][0-9]|[0-1][0-9][0-9])",count);
	crop (ftmp, 2, 1);
	
	// quad 2
	ftmp = select_corner (list_ROI_files[i],2,"0(25[0-5]|2[0-4][0-9]|[0-1][0-9][0-9])-0(25[5-9]|2[6-9][0-9]|[3-9][0-9][0-9])",count);
	crop (ftmp, 2, 2);

	// quad 3
	ftmp = select_corner (list_ROI_files[i],3,"0(25[5-9]|2[6-9][0-9]|[3-9][0-9][0-9])-0(25[0-5]|2[0-4][0-9]|[0-1][0-9][0-9])",count);
	crop (ftmp, 2, 3);

	// quad 4
	ftmp = select_corner (list_ROI_files[i],4,"0(25[5-9]|2[6-9][0-9]|[3-9][0-9][0-9])-0(25[5-9]|2[6-9][0-9]|[3-9][0-9][0-9])",count);
	crop (ftmp, 2, 4);
}
	

/* in this function quadrants are from left to right, top to bottom: 
 *  1: upper left
 *  2: upper right
 *  3: lower left
 *  4: lower right
 *  
 */
function select_corner(filename, quad, regex, count) {

	open(output_black_map + "black_map.tif");
	//open("D:\\muriel\\Projects\\NNS\\trainingsets\\trainingset_MCF7_20x_new\\segmentation\\blackmap\\black_map.tif");
	roiManager("Open", input_ROIs + filename);

	//Find ROIs in corner specified by regex
	rois = findRoisWithName(regex); 

	roiManager("Deselect");
	roiManager("Select",rois);
	roiManager("Delete");
	roiManager("Save", output_roisubs + filename + "_quad" + quad + ".zip");

	// Make the colored segmentation
	for (i=0;i<roiManager("count");i++){
		roiManager("select",i);
		setForegroundColor(i,i,i); // This will change the color value (which has to be between 0 and 255 for an 8-bit image) to the size of the region divided by 3
		roiManager("fill");
		roiManager("Show None");
	}
	saveAs("Tiff", output_greyscale + count + "_greyscale_quad" + quad);
	run("Close All");

	return count + "_greyscale_quad" + quad + ".tif"
	}


function generate_black_map () {
	waitForUser("Please select a representatory \n microscopy image to extract the pixel dimensions. \n Please ensure that the image is in the tiff file format");
	open();
	height = getHeight();
	width = getWidth();
	//newImage("blackmap", "8-bit black", width, height, 1);
	makeRectangle(0, 0, width, height);
	run("Clear", "slice");
	saveAs("Tiff", output_black_map + "black_map");		
}

/* 
 * Returns an array of indexes of ROIs that match  
 * the given regular expression 
 */ 
function findRoisWithName(roiName) { 
	nR = roiManager("Count"); 
	roiIdx = newArray(nR); 
	k=0; 
	clippedIdx = newArray(0); 
	 
	for (i=0; i<nR; i++) { 
		roiManager("Select", i); 
		rName = Roi.getName(); 
		
		if (!matches(rName, roiName) ) { 
			//print("found one!");
			//print(rName);
			roiIdx[k] = i; 
			k++; 
		} 
	} 
	if (k>0) { 
		clippedIdx = Array.trim(roiIdx,k); 
	} 
	 
	return clippedIdx; 
} 

function crop (filename, n, quad) {
	// This macro chops an image into NxN tiles, where N is the number
	// of divisions chosen by the user.
	open(output_greyscale + filename);

	if (quad==1) {
		xselect = 0;
		yselect = 0;
	} else if (quad==2) {
		xselect = 1;
		yselect = 0;
	} else if (quad==3) {
		xselect = 0;
		yselect = 1;
	} else if (quad==4) {
		xselect = 1;
		yselect = 1;
	}
	
	id = getImageID(); 
	title = getTitle(); 
	getLocationAndSize(locX, locY, sizeW, sizeH); 
	width = getWidth(); 
	height = getHeight(); 
	tileWidth = width / n; 
	tileHeight = height / n; 
	for (y = 0; y < n; y++) { 
	offsetY = y * height / n; 
	 for (x = 0; x < n; x++) { 
	offsetX = x * width / n; 
	selectImage(id); 
	 call("ij.gui.ImageWindow.setNextLocation", locX + offsetX, locY + offsetY); 
	tileTitle = title + " [" + x + "," + y + "]"; 
	 run("Duplicate...", "title=" + tileTitle); 
	makeRectangle(offsetX, offsetY, tileWidth, tileHeight); 
	 run("Crop"); 
	 rename(title+ "_x"+x+"_y"+y);
	 if (x==xselect && y==yselect){
	 	saveAs("Tiff", output_cropped + title + "_x"+x+"_y"+y);
	 } 
	} 	
	} 
	selectImage(id); 
	run("Close All");
}
