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

HPMTerms[Lsym_, Asym_, guessfn_, n_] :=
	ReleaseHold@ReplaceAll[Hold[SolveDeformation[ 
		LinearSolver[Lsym, p], 
		SolutionHomotopy[Lsym,Asym,v0],
		v, 
		v0,
		p,
		n]],{v0 -> guessfn @@ SpatioTemporalParams[NSpaceDims[Asym]], 
			v -> DeformingFunction[NSpaceDims[Asym], p]}]

Operator[Burgers1] = (D[#,t] + # D[#,x1] + D[#,{x1,2}])&
NSpaceDims[Burgers1] = 1
InitialGuess[Burgers1][x_,t_] = (\[Alpha] + \[Beta] + (\[Beta] - \[Alpha]) e^\[Gamma])/(
 1 + e^\[Gamma]) /. {\[Gamma] -> \[Alpha]/\[Epsilon] (x - \[Lambda])}
RelevantDerivatives[Burgers1, pastsolns_List] := Catenate[RelevantDerivatives /@ pastsolns]
RelevantDerivatives[Burgers1, pastsoln_] := {pastsoln, D[pastsoln, x], D[pastsoln, {x,2}]}

Deformation[H_,v_,n_] := Derivative[0,n][H][v,0]
DerivP[v_, i_,nspacedims_] := ((Derivative @@ {Sequence @@ ConstantArray[0, nspacedims+1], i}) @ v) @@ {SpatioTemporalParams[nspacedims],0} 

Operator[TimeDerivative] = (Derivative[#,t])&
LinearSolver[TimeDerivative, Asym_, v_, H_][pastsolns_] :=
	OneVariableSolve[UsePastSolns[Asym, Deformation[H,v,n], pastsolns], DerivP[v,n,nspacedims], t]

Begin["`Private`"]
(* Implementation of the package *)

UsePastSolns[Asym_, deformation_, pastsolns_] := deformation /. RelevantDerivatives[Asym, pastsolns]

OneVariableSolve[eqn_, target_, var_] := MySimplify[Integrate[Solve[eqn, D[target,var]], var]]

MySimplify[expr_] := FullSimplify[expr /. {Log[e] -> 1}]

DeformingFunction[n_, p_] := Phi @@ DeformingFunctionParams[n, p]

SpatioTemporalParams[1] := {x, t}
SpatioTemporalParams[n_] := {Sequence[Array[x,n]], t}
DeformingFunctionParams[n_, p_] := {Sequence @@ SpatioTemporalParams[n], p}

SolutionHomotopy[Lsym_, Asym_, v0_][v_, q_] := (1-q) Operator[Lsym][v - v0] + q Operator[Asym][v]

(* TODO: Untested for multi-variable *)
SolveDeformation[solver_, H_, v_, v0_List, p_, 0] :=
	Sequence @@ {v, MapThread[((# /. p -> 0) -> #2)&, {v, v0}], H}
	
(* TODO: Untested for one-variable *)
SolveDeformation[solver_, H_, v_, v0_, p_, 0] :=
	Sequence @@ {v, (v /. p -> 0) -> v0, H}

(* TODO: Not currently maintaining old derivatives, which would be more efficient *)
SolveDeformation[solver_, H_, v_, v0_, p_, n_] :=
	solver[SolveDeformation[solver, H, v, v0, p, n-1]]

End[]

EndPackage[]

