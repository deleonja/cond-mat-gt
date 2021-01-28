(* ::Package:: *)

(* ::Section:: *)
(*Component erasing maps package*)


(* ::Input::Initialization:: *)
(*Author: Jos\[EAcute] Alfredo de Le\[OAcute]n*)
(*Date: August 07, 2020*)
BeginPackage["quantumJA`"]

Reshuffle::usage=
"Reshuffle[\[ScriptCapitalE]_List] reshuffles the superoperator of a qubit quantum channel \[ScriptCapitalE]."
Pauli::usage=
"Pauli[Indices_List] gives the tensor product of Pauli Matrices with indices in Indices_List."
PauliToComp::usage=
"PauliToComp[qubits] calculates the change of basis matrix for a qubit-system map from computational to tensor product of Pauli matrices basis."
PCE::usage=
"PCE[diagElements, qubitsNum] calculates the matrix representation of a map in the tensor product of Pauli matrices given the
diagonal elements of the matrix in computational basis."
CubePositions::usage=
"CubePositions[diagElPos] gives the positions of the painted cubes given the positions of the 1's in the diagonal of the quantum channel."
Cube3q::usage=
"Cube3q[indices] graphs the 3-qubit board given the indices of the painted cubes."
CPtest::usage=
"CPtest[points] returns True if CP, False if not."
PtestM::usage=
"Ptest[A] evaluates the positive-semidefiniteness of A using the principal minors criterion."
Ptest::usage=
"Ptest[A] evaluates the positive-semidefiniteness of A with its eigenvalues."
Dirac::usage=
"Dirac[vector] returns vector in Dirac notation in computational basis."
TwoQBoard::usage=
"TwoQBoard[diagonalPCE] returns a two qubits board."
HSInnerP::usage=
"HSInnerP[A,B] returns the Hilbert-Schmidt inner product between matrices A and B."
KetBra::usage=
"KetBra[A_] returns ket-bra notation of a matrix in computational basis."

Begin["`Private`"]
Reshuffle[\[ScriptCapitalE]_List]:=Module[{dimSubMatrix,subMatrices},
dimSubMatrix=Sqrt[Length[\[ScriptCapitalE]]];
subMatrices=ArrayReshape[#,{dimSubMatrix,dimSubMatrix}]&/@\[ScriptCapitalE];
ArrayFlatten[Partition[subMatrices,dimSubMatrix]]
]

Pauli[0]=Pauli[{0}]={{1,0},{0,1}};
Pauli[1]=Pauli[{1}]={{0,1},{1,0}};
Pauli[2]=Pauli[{2}]={{0,-I},{I,0}};
Pauli[3]=Pauli[{3}]={{1,0},{0,-1}};
Pauli[Indices_List]:=KroneckerProduct@@(Pauli/@Indices);

PauliToComp[qubits_Integer]:=Module[{indices},
indices=Tuples[Range[0,3],qubits];
Transpose[Map[Flatten[Pauli[#]]&,indices]]
]

PCE[diagElements_, qubitsNum_]:=
PauliToComp[qubitsNum].
DiagonalMatrix[diagElements].
Inverse[PauliToComp[qubitsNum]]

CubePositions[diagElPos_]:=
Position[ArrayReshape[SparseArray[diagElPos->ConstantArray[1,Length[diagElPos]],{64}]//Normal,{4,4,4}],1]-1

Cube3q[indices_]:=Graphics3D[{If[Count[#,0]==3,{Black,Cube[#]},
If[Count[#,0]==2,{RGBColor["#CC0000"],Cube[#]},
If[Count[#,0]==1,{RGBColor["#004C99"],Cube[#]},
If[Count[#,0]==0,{RGBColor["#99FF33"],Cube[#]}]]]]&/@indices,
{Thickness[0.012],Line[{{{-0.5,-0.5,-0.5},{-0.5,-0.5,3.5}},{{-0.5,-0.5,-0.5},{-0.5,3.5,-0.5}},{{-0.5,-0.5,-0.5},{3.5,-0.5,-0.5}},
{{3.5,-0.5,-0.5},{3.5,-0.5,3.5}},
{{-0.5,-0.5,3.5},{3.5,-0.5,3.5}},
{{-0.5,3.5,-0.5},{3.5,3.5,-0.5}},
{{3.5,3.5,-0.5},{3.5,3.5,3.5}},
{{3.5,3.5,3.5},{-0.5,3.5,3.5}},
{{-0.5,3.5,3.5},{-0.5,3.5,-0.5}},
{{-0.5,3.5,3.5},{-0.5,-0.5,3.5}},
{{3.5,3.5,3.5},{3.5,-0.5,3.5}},
{{3.5,3.5,-0.5},{3.5,-0.5,-0.5}}}]}},
Axes->False,AxesLabel->{"x","y","z"},LabelStyle->Directive[Bold,Medium,Black],PlotRange->{{-0.5,3.5},{-0.5,3.5},{-0.5,3.5}},AxesOrigin->{0.5,0.5,0.5},AxesStyle->Thickness[0.005],ImageSize->Medium,ImagePadding->45]

CPtest[points_]:=If[(PCE[SparseArray[points+1->ConstantArray[1,{points//Length}],{4,4,4}]//Normal//Flatten,3]//Reshuffle//Eigenvalues//Min)>=0,True,False]

PtestM[A_]:=AllTrue[(Diagonal[Map[Reverse,Minors[A,#],{0,1}]]&/@Range[Length[A]]),#>=0&,2]

Ptest[A_]:=(A//Eigenvalues//Min)>=0

Dirac[vector_List]:=(vector[[#]]Ket[IntegerString[(#-1),2,Log[2,Length[vector]]]])&/@Delete[Range[Length[vector]],Position[vector,0]]//Total

TwoQBoard[diagonalPCE_List]:=ArrayPlot[ArrayReshape[diagonalPCE,{4,4}]]

HSInnerP[A_List,B_List]:=Tr[ConjugateTranspose[A].B]

KetBra[A_]:=Flatten[A].(Flatten[Table[Ket[IntegerString[#-1,2,Log[2,Length[A]]]],Length[A]]&/@Range[Length[A]]][[#]].Flatten[Transpose[Table[Bra[IntegerString[#-1,2,Log[2,Length[A]]]],Length[A]]&/@Range[Length[A]]]][[#]]&/@Range[Length[A]^2])
End[];
EndPackage[]



