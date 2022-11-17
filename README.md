# hqr
This repository works on transcribing hqr.f and hqr2.f from eispack into C to convert a Hessenberg matrix into the quasi-Schur form using orthogonal similarity transformation.

# NOTE
* hqr2.f and hqr.f sometimes produce different eigenvalues. An example of this is when the seed is 28 and n is 100. Maybe worth investigating?

# TODO
* Create a port of hqr2.f (Seems to be just a few more lines compared to hqr.f)
	* Implement cdiv.f
    * Implement new blocks inside the check for flag
* Remove goto statements (and preserve readibility)
* Create a single file hqr.c (keep the subroutines structure too)
* Run timing comparisons against the Fortran77 subroutine
* Clean up existing Code
* Make sure comments on all of my stuff is clear

# License
This repository is to be released under GPLv3 (or later). More information can be found in the LICENSE file and at https://www.gnu.org/licenses/gpl.txt
