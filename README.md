# oldage-paper-code-released

Code for [Feehan, Dennis M. "Separating the Signal From the Noise: Evidence for Deceleration in Old-Age Death Rates." Demography 55.6 (2018): 2025-2044.](https://link.springer.com/article/10.1007/s13524-018-0728-x); earlier [arxiv versions](https://arxiv.org/abs/1707.09433).

This code makes use of the [mortfit R package](https://github.com/dfeehan/mortfit).

Notes
-----

* you will have to configure the code to run on your machine; in particular, please be sure to set the various directory-related variables correctly. These are at the top of the file
* in particular, you will have to register for the [Kannisto-Thatcher database](https://www.demogr.mpg.de/databases/ktdb/) and download the data from their website; please set the variable `kt.dir` the directory that has all of the raw datafiles
