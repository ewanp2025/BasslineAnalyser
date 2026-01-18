WIP

# Bassline Analyser
A Qt-based application for analysing low frequency audio and detecting / transcribing bass notes to enter into a 16 step sequencer.

## Features
* Visualises waveform data.
* Performs FFT analysis to detect pitch (C2, F1, etc.).
* Not working with all .wav formats yet.

## Future Features
* Low pass filter.

## Third-Party Libraries
This project uses the following open-source libraries:

### QCustomPlot
* **Description:** Used for plotting the audio waveform and grid.
* **Author:** Emanuel Eichhammer
* **License:** GPLv3
* **Source:** [https://www.qcustomplot.com/](https://www.qcustomplot.com/)

### KissFFT
* **Description:** A Fast Fourier Transform library used for frequency analysis.
* **Author:** Mark Borgerding
* **License:** BSD-3-Clause
* **Source:** [https://github.com/mborgerding/kissfft](https://github.com/mborgerding/kissfft)
