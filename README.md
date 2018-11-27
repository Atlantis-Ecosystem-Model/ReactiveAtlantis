[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/jporobicg/ReactiveAtlantis.svg?branch=master)](https://travis-ci.org/jporobicg/ReactiveAtlantis)
# ReactiveAtlantis
*ReactiveAtlantis* is a R-package builded using the Shiny
package as the main platform for the reactive programming approach.

*ReactiveAtlantis* has several tools that were created to help in the tuning,
parameterization and analysis of the processes and parameters most often modified
during the calibration of Atlantis (e.g. growth rate, predation, recruitment,
Audzijonyte *et al.* 2017. Among the processes performed by this
package are:
* Visualization and analysis of the input, output and initial conditions of an Atlantis model.
*  Interactive modification of Atlantis configuration files.
*  Simulation of new parameters to help in the calibration on an Atlantis model.
*  Execution of a model skill assessment, to evaluate the performance of the model
  to reflect the observed data.

## Getting Started
These instructions will give you access to use the R-package *ReactiveAtlantis*. If
you have some problem for your installation, please let me know and I will try to
solve it as soon as possible.

### Prerequisites and installation

What things you need to install To run *ReactiveAtlantis* on R.

```R
# install packages
install.packages('devtools')   ## you need to do this step just once
# running
library("devtools")
install_github('Atlantis-Ecosystem-Model/ReactiveAtlantis','Atlantis-Ecosystem-Model', force=TRUE, dependencies=TRUE)
library("ReactiveAtlantis")
```

## Running *ReactiveAtlantis*
### Compare outputs and Biomass visualization
```R
nc.current  <- 'your_current_output.nc'
nc.old      <- 'your_previous_output.nc'
grp.csv     <- 'your_groups_definition_file.csv'
bgm.file    <- 'your_spatial_configuration_file.bgm'
cum.depths  <- c(0, 20, 50, 150, 250, 400, 650, 1000, 4300) ## This should be the cummulative depth of your model
## individual file
compare(nc.current, nc.out.old = NULL, grp.csv, bgm.file, cum.depths)
## compare to previuos run
compare(nc.current, nc.old, grp.csv, bgm.file, cum.depths)
```

### Predation analysis from the Atlantis output
```R
biom        <- 'your_BiomIndx.txt'
diet.file   <- 'your_DietCheck.txt'
bio.age     <- 'your_AgeBiomIndx.txt' ## optional file. just if you want to check the predation by age
grp.csv     <- 'your_groups_definition_file.csv'
## Predation by Age
predation(biom, grp.csv, diet.file, bio.age)
## No predation by Age
predation(biom, grp.csv, diet.file, bio.age = NULL)

```

### Exploring predator-prey interactions from the initial conditions
```R
prm.file    <- 'your_prm_file.prm'
nc.file     <- 'your_current_output.nc'
grp.csv     <- 'your_groups_definition_file.csv'
bgm.file    <- 'your_spatial_configuration_file.bgm'
cum.depths  <- c(0, 20, 50, 150, 250, 400, 650, 1000, 4300) ## This should be the cummulative depth of your model
feeding.mat(prm.file, grp.file, nc.file, bgm.file, cum.depths)
```

### Atlantis food web and trophic level composition
```R
grp.csv     <- 'your_groups_definition_file.csv'
prm.file    <- 'your_prm_file.prm'
diet.file   <- 'your_DietCheck.txt'
food.web(diet.file, grp.file)

```

### Growth of primary producers and limiting factors
```R
nc.initial  <- 'your_initial_conditions.nc'
nc.current  <- 'your_current_output.nc'
grp.csv     <- 'your_groups_definition_file.csv'
prm.file    <- 'your_prm_file.prm'
growth.pp(nc.initial, grp.csv, prm.file, nc.current)
```

### Analysis of recruitment and primary production
```R
nc.initial  <- 'your_initial_conditions.nc'
nc.current  <- 'your_current_output.nc'
yoy.file    <- 'your_yoy_file.txt'
grp.csv     <- 'your_groups_definition_file.csv'
prm.file    <- 'your_prm_file.prm'
recruitment.cal(nc.initial, nc.current, yoy.file, grp.file, prm.file)
```

### Harvest outputs and model skill assessment
```R

catch.nc    <- 'your_output_CATCH.nc'
ext.catch   <- 'external_catch_time_serie.csv'
cum.depths  <- c(0, 20, 50, 150, 250, 400, 650, 1000, 4300) ## This should be the cummulative depth of your model
fsh.csv     <- 'your_fisheries_definition_file.csv'
bgm.file    <- 'your_spatial_configuration_file.bgm'
grp.csv     <- 'your_groups_definition_file.csv'
catch(grp.csv, fsh.csv, catch.nc, ext.catch)
```
## Authors

* **Javier Porobic**

## License

This project is licensed under [GPL3](https://www.gnu.org/licenses/gpl-3.0.en.html)
