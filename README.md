# Imaging_cleaning
GUI to clean CNMF (Caiman) calcium imaging data of 1 and 2 photon data.

It loads CNMF (Caiman) .mat or .hdf5 files and allows the user to manually curate data. After deleting, merging, or both, it saves data as a .mat file for future analyses.

Clone/download the code, add to your Matlab path, run GUI_cnmf_cleaning.

Settings: allows to apply thresholds for including ROIs (spatial minimun pixels, spatial maximum pixels, SNR minimum value) and thresholds to identify potential merge pairs (distance threshold, correlation threshold)
ROIs that do not fullfil spatial and SNR will be moved to the delete list
ROI pairs that fullfil the merge pair thresholds will be accessible in the Merge panel
Plotting: allows to change visulization while running the GUI (contour threshold, scaling of plot, correlation image vs maximum projection, DF/F vs denoised)

While visualizing single ROIs in the Delete panel, the Single ROI info will update
While visualizing ROI pairs in the Merge panel, Merge pair info will update

Delete panel: allows for single ROI visualization, adding and removing from delete list. Hit Delete (red button) to actually delete ROIs that are in the delete list.
Merge panel: go through ROI pairs that are identified by the distance and correlation threshold settings. Use Add pair to add to the merge list. Use the dropdown buttons to manually select 2 ROIs and press Manual add to add to the merge list. Use Find repeats to identify ROIs that occur more than once in the merge list. Use Run multi-merge to automatically find ideal merge pairs. Rerun multi-merge if ROIs are still repeated. Hit Merge (red button) to actually merge ROIs that are in the merge list. 

Workflow: change settings, load data, start gui, work on delete or merge panel, delete or merge the data, save data, restart gui for a next round of cleaning

Make sure there are no repeated ROIs in the merge list before merging. Find them by Find repeats and manually delete them from the merge list, or press Run multi-merge a few times until there are not repeats anymore

![screenshot_GUI](https://github.com/user-attachments/assets/cfc76e9c-7b3a-46e5-832e-26a32d7470ca)

