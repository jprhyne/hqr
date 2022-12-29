# hqr
This repository works on transcribing hqr.f and hqr2.f from eispack into C to convert a Hessenberg matrix into the quasi-Schur form using orthogonal similarity transformation.

# Compilation Instructions
There have been issues in compilation with Linux and Mac devices, the current workaround is
* ```cp make.inc.example make.inc```
* run ```make```
* edit ```make.inc``` from ```FC``` to ```CC``` if you have compilation issues

# TODO
* Create a test suite (python or matlab probably)
    * Several matrices with known Quasi-Schur forms
    * 
* Create a port of hqr2.f (Seems to be just a few more lines compared to hqr.f)
    * Implement cdiv.f
    * Implement new blocks inside the check for flag
* Remove goto statements (and preserve readibility)
* Create a single file hqr.c (keep the subroutines structure too)
* Clean up existing Code
* Make sure comments on all of my stuff is clear

# License
This repository is to be released under GPLv3 (or later). More information can be found in the LICENSE file and at https://www.gnu.org/licenses/gpl.txt
