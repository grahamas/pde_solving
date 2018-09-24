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
		LinearSolver[Lsym, Asym, v, v0], 
		Asym,
		v, 
		v0,
		n]],{v0 -> guessfn[SpatioTemporalParams[Asym]], 
			v -> Phi}]

(* Burgers1-specific definitions *)
Operator[Burgers1] = (D[#,t] + # D[#,x] + D[#,{x,2}])&
NSpace[Burgers1] = 1
InitialGuess[Burgers1][x_,t_] = (\[Alpha] + \[Beta] + (\[Beta] - \[Alpha]) e^\[Gamma])/(
 1 + e^\[Gamma]) /. {\[Gamma] -> \[Alpha]/\[Epsilon] (x - \[Lambda])}
RelevantDerivatives[Burgers1, pastsolns_List] := Catenate[(RelevantDerivatives[Burgers1,#])& /@ pastsolns]
RelevantDerivatives[Burgers1, pastsoln_] := {pastsoln, D[pastsoln, x], D[pastsoln, {x,2}]}

(* TimeDerivative linear operator definitions *) 
Operator[TimeDerivative] = (D[#,t])&  (* WARNING: Is this right??? Might need to be Derivative and specify by position *)
LinearSolver[TimeDerivative, Asym_, v_, v0_][pastsolns_, n_] :=
	OneVariableSolveZeroBoundary[
		UsePastSolns[
			Asym, 
			Deformation[SolutionHomotopy[TimeDerivative,Asym,v0],v,n], 
			pastsolns], 
		DerivP[Asym,v,n], t][[1]]

(* Implementation of the package *)

Begin["`Private`"]
(* General definitions *)
UsePastSolns[Asym_, deformation_, pastsolns_] := deformation /. RelevantDerivatives[Asym, pastsolns]
Deformation[H_,v_,n_] := Derivative[0,n][H][v,0]
DerivP[Asym_, v_, i_] := ((Derivative @@ {Sequence @@ ConstantArray[0, NSpace[Asym]+1], i}) @ v)[P0Params[NSpace[Asym]]]
OneVariableSolveZeroBoundary[eqn_, target_, var_] := 
	MySimplify[Solve[Integrate[eqn, var] == 0, target]]
MySimplify[expr_] := FullSimplify[expr /. {Log[e] -> 1}]

SpatioTemporalParams[1] := Sequence[x, t]
SpatioTemporalParams[n_Integer] := Sequence[Sequence @@ Array[x,n], t]
SpatioTemporalParams[sym_Symbol] := SpatioTemporalParams[NSpace[sym]]
DeformationParams[norsym_, p_] := Sequence @@ {SpatioTemporalParams[norsym], p}
P0Params[norsym_] := Sequence @@ {SpatioTemporalParams[norsym], 0}


SolutionHomotopy[Lsym_, Asym_, v0_][v_, q_] := (1-q) Operator[Lsym][v[DeformationParams[Asym,q]] - v0] + q Operator[Asym][v[DeformationParams[Asym,q]]]

(* TODO: Untested for multi-variable *)
SolveDeformation[solver_, Asym_, v_, v0_List, 0] :=
	MapThread[(#[P0Params[NSpace[Asym]]] -> #2)&, {v, v0}]
	
(* TODO: Untested for one-variable *)
SolveDeformation[solver_, Asym_, v_, v0_, 0] :=
	 {v[P0Params[NSpace[Asym]]] -> v0}

(* TODO: Not currently maintaining old derivatives, which would be more efficient *)
SolveDeformation[solver_, Asym_, v_, v0_, n_] := 
	ReleaseHold@ReplaceAll[Hold[Catenate[{pastsoln, solver[pastsoln, n]}]], 
		pastsoln -> SolveDeformation[solver, Asym, v, v0, n-1]]

End[]

EndPackage[]

