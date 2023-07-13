To generate predictions, the python script nnsegmentation.py should be run with the following system arguments: 
1. the folder that contains the trained models (i.e., a model for cell distance predictions cdp_model.pt and a model for neighbor distance predictions ndp_model.pt)
2. The folder that contains the 8-bit images that need to be segmented
3. The location where the segmented images should be saved

Example:

python nnsegmentation.py "path/to/models/" "path/to/input/dir/" "path/to/output/dir/"