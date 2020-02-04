# Dual EEG

This repository contains code written for dual EEG project supported by SNF Sinergia grant CRSII5_170873.


It aims at data collection, processing and analysis for simultaneous scalp (high density) / intracranial (depth) EEG, combined with single pulse electrical stimulations (direct electrical stimulation).


Programming language is mainly Matlab, but there is also some use of Python, bash, Windows batch, ... etc. Most of this code focuses on the interoperability between [Cartool](https://sites.google.com/site/cartoolcommunity) and Matlab, as well as preprocessing of simultaneous EEG data.


Once the data is acquired, the pipeline follows these steps:

* Split .mff files using NetStation if necessary
* Open .mff and .TRC in Cartool and export them to .sef
* If sampling frequencies in both modalities do not resolve to a common divisor (e.g. 2048 & 1000 Hz), use Matlab script to resample (interpolate) one modality, e.g. icEEG
* Use Matlab script to realign markers
* Process CT & MRI data using FreeSurfer, iELVis, and Connectome Mapper 3
* Generate bipolar montages (.mtg) and electrode setup files (.els) for icEEG based on .grid files
* Mark events (e.g. IED) in both modalities (icEEG & hdEEG)
* Use Matlab scripts to preprocess the simultaneous data (this code)


NB: *Most of* general usage functions filenames follow underscore and lowercase, whereas scripts calling the latter functions *generally* follow Pascal case.


This is experimental work in progress: use at your own risks.


Some code in this repository comes from other (open) sources (mentioned in function help section). Please check licensing for relevant code and cite it as required.


[About me and the project](https://sinergiaconsortium.bitbucket.io/resources/html/people/renaudmarquis.html)
