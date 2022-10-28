# hqr
This repository works on transcribing hqr.f and hqr2.f from eispack into C to convert a Hessenberg matrix into the quasi-Schur form using orthogonal similarity transformation.

# TODO
* Remove goto statements (and preserve readibility)
* Create a single file hqr.c (keep the subroutines structure too)
* Create a port of hqr2.f (Seems to be just a few more lines compared to hqr.f)
	* Maybe have a flag that will run this extra functionality?
* Run timing comparisons against the Fortran77 subroutine
* Clean up existing Code
* Make sure comments on all of my stuff is clear

# Analysis of testing results 
After running the testing cases for n=1,...,500 there are meaningful differences between the real parts of the eigenvalues.
This is likely due to an ordering issue, maybe look into this if worthwhile? The differences in the imaginary part are always
less than 10^{-16}

# License
This repository is to be released under GPLv3 (or later). More information can be found in the LICENSE file and at https://www.gnu.org/licenses/gpl.txt
