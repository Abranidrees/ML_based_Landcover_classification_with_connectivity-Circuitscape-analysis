# Project structure

The following is the structure of the project.

-   project/

    -   README.md

    -   project.Rproj

    -   .Rhistory

    -   .RData

    -   renv.lock

    -   .Rprofile

    -   Config

        -   paths.yml (all file paths)

        -   colors_base7.csv (color code of land cover classes)

        -   circuitscape_gen.ini

        -   circuitscape_pfl.ini

    -   logs (other useful information e.g., sessionInfo.txt, alignment_check.csv etc.)

    -   outputs

        -   figures (contain the output figuresof all the windows)

        -   rasters (contain the predicted rasters of all the windows)

        -   tables (contain the csv tables of all the windows)

        -   vectors (contain the training points of all the windows)

    -   rasters

        -   Composities (contains the input rasters e.g., all_bands_04-06.tif)

        -   reference (ValidationLandcover2023.tif)

        -   resistance (GEN and PEL resistance rasters)

    -   renv

    -   reports (methods_results.html and methods_results.Rmd)

    -   scripts (contain all the scripts e.g., 00_setup.R; 10_train_points.R; 20_classification.R etc.)

    -   tables (train_counts_seed 2023.csv)

    -   temp_tiles (temporaory tiles)

    -   vectors (input training points)

## How to Run

Install the locked packages versions from renv.lock file.

```{r}
install.packages("renv")
renv::restore()
```

Run the scripts in the following order:

1.  00_training_seed_alignment.R

2.  00_setup.R

3.  10_train_points.R (Note: change the window and the respective points_in path)

4.  20_classification.R (Note: change the window )

5.  30_accuracy.R (Note: change the window )

6.  40_change.R

7.  50_fragmentation.R (Note: change the window)

8.  55_frag_summary.R

9.  60_gdistance.R (Note: change the specie type and Period)

10. 61_gdistance_summary.R (Note: change the specie type)

11. 65_circuitscape.R (Note: change the period)

## Download the full reproducible project file
You can downlaod the whole project with this link: https://doi.org/10.5281/zenodo.18876331
## Output Maps
<img width="724" height="737" alt="Screenshot 2026-03-03 124219" src="https://github.com/user-attachments/assets/0c59122d-beab-455d-8cc1-fa418bab209c" />
<img width="719" height="740" alt="Screenshot 2026-03-03 124226" src="https://github.com/user-attachments/assets/39f25e42-292c-456b-9014-d48dab1b4840" />
<img width="728" height="746" alt="Screenshot 2026-03-03 124233" src="https://github.com/user-attachments/assets/c5604c46-839d-4758-a991-c26db00f349e" />
<img width="737" height="744" alt="Screenshot 2026-03-03 124211" src="https://github.com/user-attachments/assets/9a121c29-dd48-4b48-a440-301f2e59f053" />
 

