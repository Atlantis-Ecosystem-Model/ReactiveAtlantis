# ReactiveAtlantis
*ReactiveAtlantis}* is a R-package builded using the Shiny
package as the main platform for the reactive programming approach.

*ReactiveAtlantis* has several tools that were created to help in the tuning,
parameterization and analysis of the processes and parameters most often modified
during the calibration of Atlantis (e.g. growth rate, predation, recruitment,
Audzijonyte *et al.* 2017. Among the processes performed by this
package are:
\begin{itemize}
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
# install library
install.packages('devtools')   ## you need to do this step just once
# running
library("devtools")
install_github('jporobicg/ReactiveAtlantis','jporobicg', force=TRUE)
```

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Javier Porobic**

## License

This project is licensed under [GPL3](https://www.gnu.org/licenses/gpl-3.0.en.html)
