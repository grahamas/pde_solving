(* Wolfram Language Package *)

(* Created by the Wolfram Workbench Sep 19, 2018 *)

BeginPackage["HPM`"]
(* Exported symbols added here with SymbolName::usage *) 

HPMTerms[Lsym_, Asym_, n_] :=
	ReleaseHold@ReplaceAll[Hold[SolveDeformation[ 
		LinearSolver[Lsym, Asym, v, v0], 
		Asym,
		v, 
		v0,
		n]], {v -> DeformingFunction[Asym], v0 -> InitialGuess[Asym]}]

(* Burgers1-specific definitions *)
Operator[Burgers1,v_,args__] := D[v[args],t] + v[args] D[v[args],x] - \[Epsilon] D[v[args],{x,2}]
NSpace[Burgers1] = 1
NSystem[Burgers1] = 1
InitialGuess[Burgers1] = ((\[Alpha] + \[Beta] + (\[Beta] - \[Alpha]) e^\[Gamma])/(
 1 + e^\[Gamma]) /. {\[Gamma] -> \[Alpha]/\[Epsilon] (# - \[Lambda])})&[x,t]
RelevantDerivatives[Burgers1, pastsolns_List] := Catenate[(RelevantDerivatives[Burgers1,#])& /@ pastsolns]
RelevantDerivatives[Burgers1, pastsoln_] := {pastsoln, D[pastsoln, x], D[pastsoln, {x,2}]}

(* Burgers2-specific definitions *)
Operator[Burgers2,{v1_, v2_},args__] := {
	D[v1[args],t] - D[v1[args], {x,2}] - 2 v1[args] D[v1[args], x] + D[v1[args] v2[args], x],
	D[v2[args],t] - D[v2[args], {x,2}] - 2 v2[args] D[v2[args], x] + D[v2[args] v1[args], x]}
NSpace[Burgers2] = 1
NSystem[Burgers2] = 2
InitialGuess[Burgers2] = ({Sin[#], Sin[#]})&[x,t]
RelevantDerivatives[Burgers2, pastsoln_] := RelevantDerivatives[Burgers1, pastsoln]

(* TimeDerivative linear operator definitions *) 
Operator[TimeDerivative,Asym_Symbol,args__] = Operator[TimeDerivative,InitialGuess[Asym],NSystem[Asym],args]
Operator[TimeDerivative,v0_,1,v_,args__] := D[v[args] - v0,t]  (* WARNING: Is this right??? Might need to be Derivative and specify by position *)
Operator[TimeDerivative,v0_,n_,vl_List,args__] := Function[v, D[v[args] - v0, t]] /@ vl
LinearSolver[TimeDerivative, Asym_, v_, v0_][pastsolns_, n_] :=
	OneVariableSolveZeroBoundary[
		UsePastSolns[
			Asym, 
			Deformation[SolutionHomotopy[TimeDerivative,Asym],v,n], 
			pastsolns], 
		DerivP[Asym,v,n], t][[1]]

(* General definitions *)
UsePastSolns[Asym_, deformation_, pastsolns_] := deformation /. RelevantDerivatives[Asym, pastsolns]
Deformation[H_,v_,n_] := Derivative[0,n][H][v,0]
DerivP[Asym_, v_, i_] := ((Derivative @@ {Sequence @@ ConstantArray[0, NSpace[Asym]+1], i}) @ v)[P0Params[NSpace[Asym]]]
OneVariableSolveZeroBoundary[eqn_, target_, var_] := 
	MySimplify[Solve[Integrate[eqn, var] == 0, target]]
MySimplify[expr_] := FullSimplify[expr /. {Log[e] -> 1}]

DeformingFunction[Asym_Symbol] := DeformingFunction[NSystem[Asym]]
DeformingFunction[1] = Phi 
DeformingFunction[n_Integer] := Array[Phi,n]

SpatioTemporalParams[1] := Sequence[x, t]
SpatioTemporalParams[n_Integer] := Sequence[Sequence @@ Array[x,n], t]
SpatioTemporalParams[sym_Symbol] := SpatioTemporalParams[NSpace[sym]]
DeformationParams[norsym_, p_] := Sequence @@ {SpatioTemporalParams[norsym], p}
P0Params[norsym_] := Sequence @@ {SpatioTemporalParams[norsym], 0}

SolutionHomotopy[Lsym_, Asym_][v_, q_] := (1-q) Operator[Lsym,Asym,v,DeformationParams[Asym,q]] + q Operator[Asym,v,DeformationParams[Asym,q]] 

(*TODO: Untested for multi-variable *)
SolveDeformation[solver_, Asym_, v_, v0_List, 0] :=
	MapThread[(#[P0Params[NSpace[Asym]]] -> #2)&, {v, v0}]
	
(* TODO: Untested for one-variable *)
SolveDeformation[solver_, Asym_, v_, v0_, 0] :=
	 {v[P0Params[NSpace[Asym]]] -> v0}

(* TODO: Not currently maintaining old derivatives, which would be more efficient *)
SolveDeformation[solver_, Asym_, v_, v0_, n_] := 
	ReleaseHold@ReplaceAll[Hold[Catenate[{pastsoln, solver[pastsoln, n]}]], 
		pastsoln -> SolveDeformation[solver, Asym, v, v0, n-1]]
		
Begin["`Private`"]
(* No need for private functions when there's only one module... *)
End[]

EndPackage[]

