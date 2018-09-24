(* Wolfram Language Test file *)

NSpace[TestSym1] = 1
NSpace[TestSym3] = 3

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