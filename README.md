E2S : Electron to Sample
this project is meant to bridge the gap between AP codes (elegant, AT, ...) used to describe a LATTICE and the
radiation seen in a BEAMLINE

In this Python version
- we use an external call to ELEGANT to determine the linear optics of a machine
- then we isolate the Twiss parameter ata sepecific ID position
- use them to initialize the SRW code to generate wavefronts up to the sample point

Notes: need to define otpical beamlines, start with I13d

