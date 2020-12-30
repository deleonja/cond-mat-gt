(* ::Package:: *)

(* ::Section:: *)
(*Component erasing maps package*)


(* ::Input::Initialization:: *)
(*Author: Jos\[EAcute] Alfredo de Le\[OAcute]n*)
(*Date: August 07, 2020*)
BeginPackage["quantumJA`"]

Reshuffle::usage=
"Reshuffle[SqMatrix] reshuffles the matrix SqMatrix."
Pauli::usage=
"Pauli[Indices_List] gives the tensor product of Pauli Matrices with indices in Indices_List."
ChangeOfBasisMatrix::usage=
"ChangeOfBasisMatrix[qbitsNum] calculates the change of basis matrix for a qubit-system map from computational to tensor product of Pauli matrices basis."
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

Begin["`Private`"]
Reshuffle[SqMatrix_]:=ArrayFlatten[ArrayFlatten/@Partition[Partition[ArrayReshape[#,{Sqrt[Dimensions[SqMatrix][[1]]],Sqrt[Dimensions[SqMatrix][[1]]]}]&/@SqMatrix,Sqrt[Dimensions[SqMatrix][[1]]]],Sqrt[Dimensions[SqMatrix][[1]]]],1];

Pauli[0]=Pauli[{0}]={{1,0},{0,1}};
Pauli[1]=Pauli[{1}]={{0,1},{1,0}};
Pauli[2]=Pauli[{2}]={{0,-I},{I,0}};
Pauli[3]=Pauli[{3}]={{1,0},{0,-1}};
Pauli[Indices_List]:=KroneckerProduct@@(Pauli/@Indices);

ChangeOfBasisMatrix[qbitsNum_]:=
Flatten/@(Pauli[#]&/@Tuples[Range[0,3],qbitsNum])//Transpose

PCE[diagElements_, qubitsNum_]:=
ChangeOfBasisMatrix[qubitsNum].
DiagonalMatrix[diagElements].
Inverse[ChangeOfBasisMatrix[qubitsNum]]

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

Dirac[vector_List]:=(vector[[#]]Ket[IntegerString[(#-1),2,Sqrt[Length[vector]]]])&/@Delete[Range[Length[vector]],Position[vector,0]]//Total

TwoQBoard[diagonalPCE_List]:=ArrayPlot[ArrayReshape[diagonalPCE,{4,4}]]
End[];
EndPackage[]









