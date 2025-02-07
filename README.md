# Imaging_cleaning
GUI to clean calcium imaging data of 1 and 2 photon data.

It loads CNMF (Caiman) .mat or .hdf5 files and allows to user to manually curate data. After deleting, merging, or both, it saves data as a .mat file for future analyses.

Settings: allows to apply thresholds for including ROIs (spatial minimun pixels, spatial maximum pixels, SNR minimum value) and thresholds to identify potential merge pairs (distance threshold, correlation threshold)
ROIs that do not fullfil spatial and SNR will be moved to the delete list
ROI pairs that fullfil the merge pair thresholds will be accessible in the Merge panel
Plotting: allows to change visulization while running the GUI (contour threshold, scaling of plot, correlation image vs maximum projection, DF/F vs denoised)

While visualizing single ROIs in the Delete panel, the Single ROI info will update
While visualizing ROI pairs in the Merge panel, Merge pair info will update

Make sure there are no repeated ROIs in the merge list before merging. Find them by Find repeats and manually delete them from the merge list, or press Run multi-merge a few times until there are not repeats anymore


