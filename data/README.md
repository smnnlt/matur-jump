# Codebook for the data sheet

## Variables

For more details on the methodology for measuring/calculating the variables, please refer to the manuscript.

#### Anthropometrics

-   `ID`: Individual participant identifier 
-   `sex`: Biological sex, coded as `f` (female) or `m` (male) 
-   `SQ`: highest attained squad level within the German Athletics Association. Athletes are labeled as `0` if they were too young to be categorized. `1`: Regional Squad; `2`: National Squad 2; `3`: National Squad 1; `4` higher than National Squad 1.
-   `group`: discipline group (here: `J` jump for all).
-   `BM`: Body mass (in kg) 
-   `BH`: Body height (in cm) 
-   `SH`: Sitting height (in cm)

#### Age

All age variables are given in years as a continuous value. 
-   `CA`: Chronological age 
-   `dAPHV`: Difference to age at peak velocity, i.e., `CA` minus `APHV`. Used as a proxy of the biological age. 
-   `APHV`: Age at peak height velocity. Calculated from `CA`, `BM`, `BH`, and `SH` using the Mirwald method.

#### Performance measures

Sprint values do not need to add up, because the best times can be from different trials. 
For the sprint variables, lower values indicate better performance (faster sprinting times). 
-   `0-10m`: 60m sprint, first 10 meters (in s)
-   `10-30m`: 60m sprint, meter 10-30 (in s) 
-   `0-30m`: 60m sprint, first half (in s)
-   `30-60m`: 60m sprint, second half (in s)
-   `0-60m`: 60m sprint, full distance (in s) 
-   `FOST`: Forward shot throw distance (in m) 
-   `TH`: Triple hop distance (in m) 
-   `FJT`: Five-jump test distance (in m) 
-   `CMJ`: Countermovement jump height (in cm) 
-   `DJ`: Drop jump efficiency 
-   `12MR`: 12-Minute run test distance (in m)
-   `FJT`: Five-jump test for distance (in m)
-   `SUb`: swing-up exercise big (count)
-   `SUs`: swing-up exercise small (count)

## License

This dataset is released under a CC-BY license (<https://creativecommons.org/licenses/by/4.0/>).
