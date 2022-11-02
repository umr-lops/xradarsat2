# xsar

radarSat2 Level 1 python reader for efficient xarray/dask based processor

 

# Install


```
conda install -c conda-forge radarSat2_xarray_reader
```

```python
>>> import radarSat2_xarray_reader
>>> folder_path = "/home/datawork-cersat-public/cache/project/sarwing/data/RS2/L1/VV/2010/288/RS2_OK72200_PK649463_DK111111_SCWA_20101015_210132_VV_SGF"
>>> radarSat2_xarray_reader.rs2_reader(folder_path)

datatree.DataTree

    Groups:
        orbitAndAttitude
            Groups: (0)
            Dimensions:
                timeStamp: 11
            Coordinates:
                timeStamp
                (timeStamp)
                datetime64[ns]
                2022-03-23T22:18:49.924403 ... 2...
            Data variables:
                yaw
                (timeStamp)
                float64
                3.814 3.808 3.792 ... 3.74 3.731
                roll
                (timeStamp)
                float64
                -29.8 -29.8 -29.8 ... -29.79 -29.79
                pitch
                (timeStamp)
                float64
                -0.003779 -0.004797 ... 0.002047
                xPosition
                (timeStamp)
                float64
                -2.721e+06 ... -2.545e+06
                yPosition
                (timeStamp)
                float64
                6.426e+06 6.419e+06 ... 6.334e+06
                zPosition
                (timeStamp)
                float64
                -1.666e+06 ... -2.209e+06
                xVelocity
                (timeStamp)
                float64
                2.185e+03 2.206e+03 ... 2.388e+03
                yVelocity
                (timeStamp)
                float64
                -924.9 -980.2 ... -1.476e+03
                zVelocity
                (timeStamp)
                float64
                -7.166e+03 ... -7.005e+03
            Attributes:

            attitudeDataSource :
                Downlink
            attitudeOffsetsApplied :
                true
            Description :
                Attitude Information Data Store. Orbit Information Data Store. State Vector Data Store. Earth Centered Rotating (ECR) coordinates.

        geolocationGrid
            Groups: (0)
            Dimensions:
                line: 11pixel: 11
            Coordinates:
                line
                (line)
                int64
                0 1026 2053 ... 8215 9242 10269
                pixel
                (pixel)
                int64
                0 1057 2114 ... 8457 9514 10572
            Data variables:
                latitude
                (line, pixel)
                float64
                -11.88 -12.0 ... -17.43 -17.54
                longitude
                (line, pixel)
                float64
                106.0 106.4 106.9 ... 109.1 109.6
                height
                (line, pixel)
                float64
                0.0 0.0 0.0 0.0 ... 0.0 0.0 0.0 0.0
            Attributes: (13)
        imageGenerationParameters
            Groups: (2)
            Dimensions:
            Coordinates: (0)
            Data variables: (0)
            Attributes: (0)
        radarParameters
            Groups: (0)
            Dimensions:
                beam: 4pole: 2NbOfNoiseLevelValues: 99
            Coordinates:
                beam
                (beam)
                <U3
                'W1' 'W2' 'W30' 'S7'
                pole
                (pole)
                <U2
                'VH' 'VV'
            Data variables:
                pulsesReceivedPerDwell
                (beam)
                int64
                58 58 58 58
                numberOfPulseIntervalsPerDwell
                (beam)
                int64
                65 66 65 67
                rank
                (beam)
                int64
                7 8 7 9
                settableGain
                (beam, pole)
                float64
                -1.0 -1.0 -1.0 ... -1.0 -1.0 -1.0
                pulseRepetitionFrequency
                (beam)
                float64
                1.275e+03 1.335e+03 ... 1.287e+03
                samplesPerEchoLine
                (beam)
                int64
                6768 7920 7944 7320
                noiseLevelValues_BetaNought
                (NbOfNoiseLevelValues)
                float64
                -26.41 -27.32 ... -23.75 -21.88
                noiseLevelValues_SigmaNought
                (NbOfNoiseLevelValues)
                float64
                -27.6 -28.52 ... -28.4 -26.61
                noiseLevelValues_Gamma
                (NbOfNoiseLevelValues)
                float64
                -25.73 -26.67 ... -28.13 -26.35
            Attributes: (20)
        lut
            Groups: (0)
            Dimensions:
                pixels: 10573
            Coordinates:
                pixels
                (pixels)
                int64
                0 1 2 3 ... 10569 10570 10571 10572
            Data variables:
                lutBeta
                (pixels)
                float64
                1.358e+07 1.358e+07 ... 1.358e+07
                lutGamma
                (pixels)
                float64
                1.159e+07 1.159e+07 ... 3.838e+07
                lutSigma
                (pixels)
                float64
                1.786e+07 1.786e+07 ... 4.071e+07
            Attributes:

            Description :
                RADARSAT Product LUT. (c) COPYRIGHT MacDonald Dettwiler and Associates Ltd., 2003 All Rights Reserved.Three output scaling Look-up Tables (LUTs) are included with every product. These LUTs allow one to convert the digital numbers foundin the output product to sigma-nought, beta-nought, or gamma-noughtvalues (depending on which LUT is used) by applying a constantoffset and range dependent gain to the SAR imagery.There is one entry in the gains list for each range sample in theimagery. In order to convert the digital number of a given rangesample to a calibrated value, the digital value is first squared,then the offset (B) is added and the result is divided by thegain value (A) corresponding to the range sample.i.e., calibrated value = ( digital value^2 + B ) / A

    Dimensions:
    Coordinates: (0)
    Data variables: (0)
    Attributes: (0)

``