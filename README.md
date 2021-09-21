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

test.hitAngleTID(angle="phi")

z, angle = test.getResult()

test.saveTxt()
test.saveCompressed()
```

### Functions:
1. **Tile hit rate and energy deposition in tiles in z-direction**
```
hitAngleTID(n, angle = [norm, theta, phi])
```

Using helix reconstruction:
```
hitAngleHelix(n = 0, angle_req = [norm, theta, phi])
```

Where n is the number of frames. When left blank it uses all frames. For ```angle``` and ```angle_req``` one of the following can be chosen:


    - ```norm``` returns the angle between the normal vector of the tile and the direction of the hit.

    - ```theta``` returns the polar angle. 

    - ```phi``` returns the azimuth angle.

Returns two arrays:
- z_arr: z-component of tile position
- angle_sensor_tile: angle (depends on setting)

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

mt_hitAngleTID(txt, npz, binned, angle)

mt_hitAngleHelix(txt, npz, binned)
```
Arguments: 
- txt, npz, binned: bool, if true it saves in that format
- angle:
    - ```norm``` returns the angle between the normal vector of the tile and the direction of the hit.

    - ```theta``` returns the polar angle. 

    - ```phi``` returns the azimuth angle.
