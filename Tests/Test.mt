(* Wolfram Language Test file *)

NSpace[TestSym1] = 1
NSpace[TestSym3] = 3

(* Test parameter functions. *)
Test[
	{SpatioTemporalParams[1]}
	,
	{x,t}
	,
	TestID->"Test-20180924-L6P4A3"
]
Test[
	{SpatioTemporalParams[TestSym1]}
	,
	{x,t}
	,
	TestID->"Test-20180924-L6P4A4"
]
Test[
	{SpatioTemporalParams[3]}
	,
	{x[1], x[2], x[3], t}
	,
	TestID->"Test-20180924-C5V5C2"
]
Test[
	{SpatioTemporalParams[TestSym3]}
	,
	{x[1], x[2], x[3], t}
	,
	TestID->"Test-20180924-C5V5C3"
]
Test[
	{DeformationParams[TestSym1,p]}
	,
	{x,t,p}
	,
	TestID->"Test-20180924-O5S9Y7"
]
Test[
	{DeformationParams[TestSym3,p]}
	,
	{x[1],x[2],x[3],t,p}
	,
	TestID->"Test-20180924-O5S9Y8"
]
Test[
	{P0Params[TestSym1]}
	,
	{x,t,0}
	,
	TestID->"Test-20180924-S8E5V0"
]
Test[
	{P0Params[TestSym3]}
	,
	{x[1],x[2],x[3],t,0}
	,
	TestID->"Test-20180924-H7Z6D7"
]

(* Test RelevantDerivatives for Burgers examples *)
Test[
	RelevantDerivatives[Burgers1, {Phi[x,t,0] -> Sin[x]}]
	,
	{Phi[x,t,0] -> Sin[x], D[Phi[x,t,0],x] -> Cos[x], D[Phi[x,t,0],{x,2}] -> -Sin[x]}
	,
	TestID->"Test-20180924-Q5D9X6"
]
Test[
	RelevantDerivatives[Burgers2, {Phi[x,t,0] -> {Sin[x], Sin[x]}}]
	,
	{Phi[x,t,0] -> {Sin[x], Sin[x]}, D[Phi[x,t,0],x] -> {Cos[x], Cos[x]}, D[Phi[x,t,0],{x,2}] -> {-Sin[x],-Sin[x]}}
	,
	TestID->"Test-20180924-S0R4A0"
]

(* Test SolutionHomotopy for Burgers examples *)
Test[
	SolutionHomotopy[TimeDerivative, Burgers1][Phi, p]
	,
	(1 - p) 
\!\(\*SuperscriptBox[\(Phi\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, p] + p (
\!\(\*SuperscriptBox[\(Phi\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, p] + Phi[x, t, p] 
\!\(\*SuperscriptBox[\(Phi\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, p] - \[Epsilon] 
\!\(\*SuperscriptBox[\(Phi\), 
TagBox[
RowBox[{"(", 
RowBox[{"2", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, p])
	,
	TestID->"Test-20180924-Y0U8R5"
]

Test[
	SolutionHomotopy[TimeDerivative, Burgers2][Phi, p]
	,
	result
	,
	TestID->"Test-20180924-M8L8B1"
]


Test[
	HPMTerms[TimeDerivative, Burgers2, 4]
	,
	{Phi[1][x, t, 0] -> Sin[x], Phi[2][x, t, 0] -> Sin[x], 
\!\(\*SuperscriptBox[\(Phi[1]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, 0] -> -t Sin[x], 
\!\(\*SuperscriptBox[\(Phi[2]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, 0] -> -t Sin[x], 
\!\(\*SuperscriptBox[\(Phi[1]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, 0] -> t^2 Sin[x], 
\!\(\*SuperscriptBox[\(Phi[2]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "2"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, 0] -> t^2 Sin[x], 
\!\(\*SuperscriptBox[\(Phi[1]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "3"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, 0] -> -t^3 Sin[x], 
\!\(\*SuperscriptBox[\(Phi[2]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "3"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, 0] -> -t^3 Sin[x], 
\!\(\*SuperscriptBox[\(Phi[1]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "4"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, 0] -> t^4 Sin[x], 
\!\(\*SuperscriptBox[\(Phi[2]\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "4"}], ")"}],
Derivative],
MultilineFunction->None]\)[x, t, 0] -> t^4 Sin[x]}
	,
	TestID->"Test-20181001-G2M8Y7"
]
	