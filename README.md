# hqr
This repository works on transcribing hqr.f and hqr2.f from eispack into C to convert a Hessenberg matrix into the quasi-Schur form using orthogonal similarity transformation.

# TODO
* Create a port of hqr2.f (Seems to be just a few more lines compared to hqr.f)
	* Maybe have a flag that will run this extra functionality?
	* Underestimated the time it would take to get this done. Will be done next week
* Remove goto statements (and preserve readibility)
* Create a single file hqr.c (keep the subroutines structure too)
* Run timing comparisons against the Fortran77 subroutine
* Clean up existing Code
* Make sure comments on all of my stuff is clear

# Analysis of testing results 
After running the testing cases for n=1,...,500 there are meaningful differences between the real parts of the eigenvalues.
This is likely due to an ordering issue, maybe look into this if worthwhile? The differences in the imaginary part are always
less than 10^{-16}

Fixed issue that was causing errors in the k=5 case. This was an edge case of dsign not being properly handled. Still looking into other issues
For the larger norm error it takes place at the computation of the second eigenvalues. Differs because on hqr.f, t has a value while in 
the c version it is zero. Issue seems to be resolved. Will run much larger test cases tomorrow, but up to n=135 there are no errors. Forgot
to add `t+=x` to my code inside the form exceptional shift.

# License
This repository is to be released under GPLv3 (or later). More information can be found in the LICENSE file and at https://www.gnu.org/licenses/gpl.txt
