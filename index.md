### What is it?
It's a matlab 'toolbox' [gui + functions] for plotting group level SPM12 M/EEG datasets. Load up a bunch of SPM datasets and plot specific conditions (averaged over group), optionally with reference MNI coordinates projected onto the mesh. Visualisation options include thresholding, inflation and alpha. Export high res images. Plot specific conditions, or type in part of the condition name to find all matching and average them all (e.g. if I have conditions 'Face_1','Face_2' ... 'Face_n', then type 'Face' to partial string match and plot the average overlay to all faces). 

Also runs t-contrasts between conditions - either select condition from dropdown or type in pairs of conditions to compute contrast for. Results can visualised in a sep window with linked-rotation and high res image output option.

![Screenshot](https://s17.postimg.org/jccp6721b/Screen_Shot_2016_12_05_at_10_31_40.png)

### What does it do?
It loads up an averaged mesh and plots the group averaged overlay for a given trial(s) on top. From this, you can click an area of interest to find the coordinates of this 'peak' for each subject [within a given radius]. You can then either reconstruct the data [time by trials] for the given identified vertices and store in a matrix, or call an spm built in function which extracts those vertices into a new LFP-SPM object for each individual dataset.

### Why'd you make it?
Because I recently came to SPM from other softwares and I wanted to keep the same analysis approach that I'd used previously. 

### How to use?
For the gui: VSExtractor_gui

For a script: Click_to_extract_virtual_sensors

