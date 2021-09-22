# Mu3eHelper for tile analysis (for python) aka melp  

# melp.TileHitRate (class)

### Usage example:
```
test = melp.TileHitRate("sorted.root" ,"Test", "Test")

test.tileHitRateHID()

z_total, z_primary, z_secondary, z_tertiary, edep_total, edep_primary, edep_secondary, edep_tertiary = test.getResultHID()
```

### Functions:
1. **Tile hit rate and energy deposition in tiles in z-direction**
```
tileHitRate(n)
```
Where n is the number of frames. When left blank it uses all frames.
Returns two arrays:
- z_arr: contains the z-direction of the hits.
- edep_arr: contains the energy deposition in z-direction.

2. **Tile hit rate and energy deposition in tiles in z-direction split in different HIDs**
```
tileHitRateHID(n)
```
Where n is the number of frames. When left blank it uses all frames.
Returns eight arrays:
- z_total_arr: contains z-direction of the hits.
- z_primary_arr: contains z-direction of the primary hits.
- z_secondary_arr: contains z-direction of the secondary hits.
- z_tertiary_arr: contains z-direction of the tertiary hits.

- edep_total_arr: contains the energy deposition in z-direction.
- edep_primary_arr: contains the energy deposition of primary hits in z-direction.
- edep_secondary_arr: contains the energy deposition of secondary hits in z-direction.
- edep_tertiary_arr: contains the energy deposition of tertiary hits in z-direction.

3. **Get results**
```
getResult()
```

4.**Save results**
```
saveNpz()
```

# melp.TileHitAngle (class)

### Usage example:
```
test = melp.TileHitAngle("sorted.root", "outtest")

test.hitAngleRec(angle="phi")

z, angle = test.getResult()

test.saveTxt()
test.saveCompressed()
```

### Functions:
1. **Tile hit rate and energy deposition in tiles in z-direction**
Using a linear path to match tile hits and pixel hits.
```
hitAngleRec(n, angle = ["norm", "theta", "phi"])
```

Using helix reconstruction:
```
hitAngleHelix(n = 0, angle = ["norm", "theta", "phi"])
```

Using MC truth information:

Note: The _digi.json_ file for the mu3eSim need to be modified. Set both ```truth_options``` to 2.
```
hitAngleTruth(n=0, angle = ["norm", "theta", "phi"], hit_type = ["primary", "secondary", "all"], particle_type = ["electron", "positron", "all"])
```

Where n is the number of frames. When left blank it uses all frames.

For ```angle``` one of the following can be chosen:

- ```norm``` returns the angle between the normal vector of the tile and the direction of the hit. (Default)
- ```theta``` returns the polar angle.
- ```phi``` returns the azimuth angle.

For ```hit_type``` one of the following can be chosen:

- ```primary``` uses only primary hits. (Default)
- ```secondary``` uses only secondary hits.
- ```all``` uses all hits.

For ```particle_type``` one of the following can be chosen:

- ```electron``` uses only electrons.
- ```positron``` uses only positrons.
- ```all``` uses all particles. (Default)

2. **Get results**
```
getResult()  

getBinned()
```

3. **Save results**

Save .txt file.
```
saveTxt()
```
Save .npz file.
```
saveNpz()
```
Save compressed .npz file
```
saveCompressed()
```
Save binned .npz file
```
saveBinned()
```


# melp.multithreading (module)
Runs multiple instances of the same function,  therefore multiple root files need to be provided, since the functions itself can't run on multiple cores. The corresponding output is saved as _.npz_ or _.txt_. For some functions also a binned output can be saved as _.npz_. These output files can later be merged. For an example, of how this can be done, take a look at the Jupyter notebook _TileHitAngle_Analyser.ipynb_.

The default is that all but two cores are used. This can be changed in _multithreading.py_.

### Usage example:
```
from melp import multithreading as mt

args = (False, True, True, "norm")
mt.run_mt("mt_tileHitRateHID", "./testdata/sorted/sorted*.root", args)
```

### Functions:
1. **Run selected function with multithreading**
```
run_mt(function_str, src, args)
```
Arguments:

- ```function_str```: name of function that should be run, as str.

- ```src```: directory of root files to analyse. The file names should include a digit in order for the glob module to work properly (see usage example).

- ```args```: arguments for parallelized function which is called with ```function_str```.


2. **Available funtions that can be parallelized.**
```
mt_tileHitRateHID(npz)

mt_hitAngleRec(txt, npz, binned, angle)

mt_hitAngleHelix(txt, npz, binned, angle)

mt_hitAngleTruth(txt, npz, binned, angle, hit_type, particle_type)
```
Arguments:
- txt, npz, binned: bool, if true it saves in that format. 

