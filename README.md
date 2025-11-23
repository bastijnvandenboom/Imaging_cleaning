# Imaging_cleaning
GUI to clean CNMF (Caiman) calcium imaging data of 1 and 2 photon data.

It loads CNMF (Caiman) .mat or .hdf5 files and allows the user to manually curate data. After deleting, merging, or both, it saves data as a .mat file for future analyses.

Clone/download the code, add to your Matlab path, run GUI_cnmf_cleaning.

Start by changing settings:

Selection: allows to filter ROIs
- Spatial min: minimum number of pixels to include an ROI (0=include all)
- Spatial max: maximum number of pixels to include an ROI (0=include all)
- SNR min: minimum signal-to-noise (SNR) value to include ROI
- Bad ROIs: move all ROIs identified as unusable by Caiman to delete list (unchecked, don't plot at all)

Merge: filters to identify merge pairs (merge pairs will be in Merge panel)
- Dist thr: find ROIs within this distance (number of pixels)
- Correl thr: minimum correlation between ROIs

- 
- Thresholds for including ROIs (spatial minimun pixels, spatial maximum pixels, signal-to-noise (SNR) minimum value)
- Thresholds to identify potential merge pairs (distance threshold, correlation threshold)
- Checking Bad ROIs will move ROIs identified as unusable by Caiman to the delete list
- Checking Flip bg will flip the background (useful for old Matlab-based CNMF-e)
- ROIs that do not fullfil spatial and SNR will be moved to the delete list
- ROI pairs that fullfil the merge pair thresholds will be accessible in the Merge panel

Plotting: allows to change visulization while running the GUI
- Contour thr: threshold for the fraction of pixels to include to plot contours or ROIs
- Scaling plot: change the scaling of the colorbar
- Cn/Mn: plot correlation image (Cn) or maximum projection (Mn)
- C_raw/C: plot DF/F (C_raw) or denoised signal (C)
- Concat mark: plots repetitive black dotted lines on the temporal traces (useful to indicate recording duration or file size)
- Flip bg: checked will flip the background (useful for old Matlab-based CNMF-e)
- Plot all contours: during GUI, checking this will plot contours of all ROIs

Start GUI
- Load File: manual select file (.hdf5 or .mat)
- START GUI: start GUI
- restart GUI: restart the GUI (either do all deletes or all merges, then save data, then restart GUI)

Save data
- Delete: delete ROIs in delete list
- Merge: merge ROIs in merge list (make sure no ROI values are repeated)
- Save data: store cleaned data

While visualizing single ROIs in the Delete panel, the Single ROI info will update

While visualizing ROI pairs in the Merge panel, Merge pair info will update

Delete panel: allows for single ROI visualization, and adding and removing from delete list. Hit Delete (red button) to actually delete ROIs that are in the delete list.

Merge panel: go through ROI pairs that are identified by the distance and correlation threshold settings. Use Add pair to add to the merge list. Use the dropdown buttons to manually select 2 ROIs and press Manual add to add to the merge list. Use Find repeats to identify ROIs that occur more than once in the merge list. Use Run multi-merge to automatically find ideal merge pairs. Rerun multi-merge if ROIs are still repeated. Hit Merge (red button) to actually merge ROIs that are in the merge list. 

Workflow: change settings, load file, start gui, work on either delete or merge panel, delete or merge the data, save data, restart gui for a next round of cleaning (to remove backend data)

Make sure there are no repeated ROIs in the merge list before merging. Find them by Find repeats and manually delete them from the merge list, or press Run multi-merge a few times until there are not repeats anymore

![screenshot_GUI](https://github.com/bastijnvandenboom/Imaging_cleaning/blob/226f63c8961a1c7056a61777bf309bc3979b2415/GUI_example.png)

