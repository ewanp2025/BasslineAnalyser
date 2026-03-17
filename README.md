WIP

# Bassline & Chord Analyser
A Qt-based application for analysing low frequency audio and detecting / transcribing bass notes to enter into a 16 step sequencer.

## Features
* Visualises waveform data.
* Performs FFT analysis to detect pitch (C2, F1, etc.).
* Low pass and notch filter.
* Export as LMMS pattern.

## Current state
* Improving on bass note detection, needs further work on chord detection.
.
## Third-Party Libraries
This project uses the following open-source libraries:

### QCustomPlot
* **Description:** Used for plotting the audio waveform and grid.
* **Author:** Emanuel Eichhammer
* **License:** GPLv3
* **Source:** [https://www.qcustomplot.com/](https://www.qcustomplot.com/)

### KissFFT
* **Description:** A lightweight, single-file C library used for decoding and loading audio files (WAV, MP3, FLAC) into raw PCM float data for DSP analysis.
* **Author:** David Reid (mackron)
* **License:** Public Domain (Unlicense) / MIT-0
* **Source:** [https://github.com/mborgerding/kissfft](https://github.com/mborgerding/kissfft)](https://github.com/mackron/miniaudio)
