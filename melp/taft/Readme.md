# Time alignment for tiles (TAFT)

```
import melp.taft as taft

taft.select(mu3e_detector)

```

## Generate Histogram file

- Python Version:

```
options_hist = {
    "histo_options": (10000, -32, 32),  # nbins, min, max (10000, -64, 64)
    "hist_file": "hist_test2.root",  # histogram save file
    "ttree_loc": "alignment/mu3e"   # alignment/
}

melp.taft.generate_hist("input_file", **options_hist)
```

- ROOT Macro:

```
CONDOR_INPUT="inputfile1.root inputfile2.root" root

.x build_histograms.C
```

## Calibrating Detector

```
options_cal = {
    "debug_station": 2,
    "tof": "advanced_new",
    "dt_mode": "median",
    "overwrite": False,
    "hist_file": "histo_2mio_frames_1.root"  # histogram file
}

resid_z, resid_phi, cal_z, cal_phi = melp.taft.calibrate(**options_cal)
```
