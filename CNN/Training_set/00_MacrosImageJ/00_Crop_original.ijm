input_img = getDirectory("Please select the directory that contains the original images");
output = getDirectory("Please select directory where the cropped images shall be saved");

// Generate cropped files
File.makeDirectory(output + "/cropped_original");
output_crop = output + "/cropped_original/";

list_orig = getFileList(input_img);
Array.sort(list_orig);

for (i = 0; i < list_orig.length; i++) {
	crop (list_orig[i], 2);
	//crop_quads (list_orig[i], 2, 1);
	//crop_quads (list_orig[i], 2, 2);
	//crop_quads (list_orig[i], 2, 3);
	//crop_quads (list_orig[i], 2, 4);
}

function crop_quads (filename, n, quad) {
	// This macro chops an image into NxN tiles, where N is the number
	// of divisions chosen by the user.
	open(input_img + filename);

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
	 	saveAs("Tiff", output_crop + title + "_quad" + quad + "_x"+x+"_y"+y);
	 } 
	} 	
	} 
	selectImage(id); 
	run("Close All");
}

function crop (filename, n) {
	// This macro chops an image into NxN tiles, where N is the number
	// of divisions chosen by the user.
	open(input_img + filename);
	
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
	 if (x==0 && y==0){
	 	quad = 1;
	 } else if (x==1 && y==0){
	 	quad = 2;
	 } else if (x==0 && y==1){
	 	quad = 3;
	 } else if (x==1 && y==1){
	 	quad = 4;
	 }
	 saveAs("Tiff", output_crop + title + "_quad" + quad + "_x"+x+"_y"+y);
	} 	
	} 
	selectImage(id); 
	run("Close All");
}
