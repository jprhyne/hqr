# hqr
This repository works on transcribing hqr.f and hqr2.f from eispack into C to convert a Hessenberg matrix into the quasi-Schur form using orthogonal similarity transformation.

# TODO
* Clean up existing Code
* Make sure comments on all of my stuff is clear
* Make more robust testing cases. (Maybe automate? Definitely compare vector norms.)
* Create a single file that contains all the current iteration
* Run timing comparisons against the Fortran77 subroutine
* Remove goto statements (and preserve readibility)
* Create a port of hqr2.f (Seems to be just a few more lines compared to hqr.f
	* Maybe have a flag that will run this extra functionality?

# License
This repository is to be released under GPLv3 (or later). More information can be found in the LICENSE file and at https://www.gnu.org/licenses/gpl.txt
