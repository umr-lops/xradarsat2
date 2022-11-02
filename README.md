# xsar

radarSat2 Level 1 python reader for efficient xarray/dask based processor

 

# Install


```
conda install -c conda-forge radarSat2_xarray_reader
```

```pycon
>>> import radarSat2_xarray_reader
>>> folder_path = "/home/datawork-cersat-public/cache/project/sarwing/data/RS2/L1/VV/2010/288/RS2_OK72200_PK649463_DK111111_SCWA_20101015_210132_VV_SGF"
>>> radarSat2_xarray_reader.rs2_reader(folder_path)

``