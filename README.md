# Imaging_cleaning
Matlab GUI to clean calcium imaging dataset of 1 and 2 photon imaging experiments.

It loads CNMF-E (.mat), CNMF (Caiman) (.mat or .hdf5), and Suite2P (.mat) files. GUI allows the user to manually curate data. After deleting, merging, or both, it saves data as a .mat file for future analyses.

You load a raw session and use the GUI to either delete ROIs or merge ROIs (figure 1). You can include the "bad" ROIs (as defined by the ROI extraction software you used) and you can use multi-merge (merge more than 2 ROIs). Instead of going through all the ROIs, you can click on the spatial overview screen on an ROI which will be visualized (spatially and temporally) and can be deleted (figure 2).

Clone/download the code, add to your Matlab path, run GUI_imaging_cleaning.

Start by changing settings:

Selection: allows to filter ROIs
- Spat min: minimum number of pixels (spatial) to include an ROI (0=include all)
- Spat max: maximum number of pixels (spatial) to include an ROI (0=include all)
- SNR min: minimum signal-to-noise (SNR) value to include ROI
- Bad ROIs: move all ROIs identified as unusable by Caiman to delete list (unchecked, don't plot at all)

Merge: filters to identify merge pairs (merge pairs will be in Merge panel)
- Dist thr: find ROIs within this distance (number of pixels)
- Correl thr: minimum correlation between ROIs

Plotting: allows to change visulization while running the GUI
- Cont thr: threshold for the fraction of pixels to include to plot contours or ROIs
- Scale plot: change the scaling of the colorbar
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

While visualizing single ROIs in the Delete panel, the Single ROI info will update.

While visualizing ROI pairs in the Merge panel, Merge pair info will update.

Delete panel: allows for single ROI visualization, and adding and removing from delete list. Hit Delete (red button) to actually delete ROIs that are in the delete list.

Merge panel: go through ROI pairs that are identified by the distance and correlation threshold settings. Use Add pair to add to the merge list. Use the dropdown buttons to manually select 2 ROIs and press Manual add to add to the merge list. Use Find repeats to identify ROIs that occur more than once in the merge list. Use Run multi-merge to automatically find ideal merge pairs. Rerun multi-merge if ROIs are still repeated. Hit Merge (red button) to actually merge ROIs that are in the merge list. 

Workflow: change settings, load file, start gui, work on either delete or merge panel, delete or merge the data, save data, restart gui for a next round of cleaning (to remove backend data).

Make sure there are no repeated ROIs in the merge list before merging. Find them by Find repeats and manually delete them from the merge list, or press Run multi-merge a few times until there are not repeats anymore.

GUI made with Matlab GUIDE but can be run with the latest Matlab version (2026a).

FIGURE 1 - overview of the GUI
![screenshot_GUI](https://github.com/bastijnvandenboom/Imaging_cleaning/blob/5b715442cc151e10e633421d4c515f25ddf1d454/GUI_example.png)

FIGURE 2 - Use mouse to select an ROI

![screenshot_GUI](https://github.com/bastijnvandenboom/Imaging_cleaning/blob/5b715442cc151e10e633421d4c515f25ddf1d454/GUI_example.png)
