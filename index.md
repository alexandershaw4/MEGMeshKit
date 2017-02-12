### What is it?
It's a set of scripts and functions for visualising group mesh's from source localised SPM MEEG objects.

See the scripts example.m to get started and mesh_group_to_video.m

Functions included

- meshmesh: plot a quick patch mesh from the data stored in the object (e.g. D.inv.forward.mesh)
- plotmesh: plot a mesh as above with MNI coords projcted on it and optional inflation
- plotmeshfo: plot mesh patch with inflation, transparency & function overlay [vector] 
- plotmesh_fo: as above but it find the overlay itself based on condition name or index of D.cond
- plotmesh_fo_grp: as above but averages overlays of a set of objects [w/ sudo-inflation of common blobs]and returns them
- plotmesh_fo_tmap: plot mesh with a t-map overlay, thresholded at a critical-t
- plotmeshfov: generate a mesh video from overlay matrix that is (nverts * ntime_pnts)

- mesh_pca1: mesh based pca for a given trial, time and freq window 
- mesh_ica: do a spatio-temporal ica for a single subjects source localisation [uses FastICA alogrithm]
- vsmooth: vertex smoothing that is essentially 'inflation'
- writeims: write out mesh and function overlay [vector] as proper nifti / gifti images
- get_xyz2: capture crosshair point [click] from 3d mesh
- ind2mni: convert figure indices from above to MNI values [using affine transfm data in tess_mni]
- gifti_inflator: inflate gifti object / mesh
- contrastmesh: perform ttests between subjects mesh points [not very good, use spm]

- loadarrayspm: from an array of spm meeg filenames, load the objects into an array
- linksubplots: link the subplots of a figure so that they rotate together
- reduce_eig_mesh: iteratively smooth mesh overlay until 90% variance explained by n-components
- plotmesh_fo_grp_pca: applied above to group


[![Mesh Vid](https://img.youtube.com/vi/VID/0.jpg)](https://www.youtube.com/watch?v=l3z34g9Zogo&feature=youtu.be)


* may also require some functions in my 'misc' repo.
