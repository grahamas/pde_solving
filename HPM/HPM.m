(* Wolfram Language Package *)

(* Created by the Wolfram Workbench Sep 19, 2018 *)

BeginPackage["HPM`"]
(* Exported symbols added here with SymbolName::usage *) 

SolutionHomotopy::usage = 
	"SolutionHomotopy[L,A,v0] gives a homotopy H[v,q] of v from the 
	initial guess, v0, at q=0, to the solution of A, at q=1. 
	L is a well-chosen linear operator.
	A is the nonlinear differential operator we wish to solve.
	v0 is the initial guess."

HPMTerms::usage =
	"HPMTerms[L, A, v0, n] gives the first n terms of the Taylor
	expansion of the solution of A, using the homotopy perturbation method
	with linear operator L and initial guess v0."

Begin["`Private`"]
(* Implementation of the package *)

SolutionHomotopy[L_, A_, v0_][v_, q_] := (1-q) L[v - v0] + q A[v]

(* TODO: Not currently maintaining old derivatives, which would be more efficient *)
HPMFromHomotopy[H_, solver_, v_, n_] :=
	solver[v, HPMFromHomotopy[H, solver, v, n-1]]

HPMTerms[Lsym_, Asym_, v0_, n_] :=
	HPMFromHomotopy[SolutionHomotopy[Lsym,Asym,v0], LinearSolver[Lsym, NVars[Asym]], DeformingFunction[NVars[Asym]], n]



End[]

EndPackage[]

