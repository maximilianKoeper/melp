# Time alignment for tiles (TAFT)

___

### Example Jupyter Notebook:
[Jupyter Notebook](../../jupyter-notbooks/taft_v2-global_correction_z-Copy2.ipynb)

___

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
    "debug_station": 2,     # 1 / 2
    "tof": "simple",  # advanced_new / advanced_graf / simple / None
    "dt_mode": "median",    # MEDIAN / mean / gaus
    "overwrite": True,      # True / False
    "hist_file": 'merged.root', 
    
    "cosmic_correction": False,
    "cosmic_mc_primary": True,
    "cosmic_n_modes" : 5,  # (x2 for cos and sin)
    "ttree_loc": "alignment/mu3e",
    "cosmic_threshold": 0.05,  #m
    "cosmic_file": 'mu3e_sorted_000002_cosmic.root'
}

resid_z, resid_phi, cal_z, cal_phi = melp.taft.calibrate(**options_cal)
```
