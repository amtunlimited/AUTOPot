(*
Potlib API
written by Aaron Tagliaboschi and Dr. Jeremy Maddox
To be used with automatically generated MathLink executable given with this package
*)

BeginPackage["potlib`", {"Units`","PhysicalConstants`"}];

(*Electron mass, used throughout*)
me=Convert[ElectronMass,AtomicMassUnit][[1]];

(*
RegFunc:
	Inputs: 
		potlib: Function from Potlib MathLink executable, 
		VBig: a suffeciently lagre number (for finding the energy when atoms are not connected)
	Output:
		a function that takes R space coordintes and outputs the potencial energy of the system at the point 
*)
RegFunc[potFunc_,Vbig_]:=Module[{outFunc,anuzero},
	anuzero=potFunc[
		{0.,0.,0.,0.,0.,Vbig,0.,0.,2*Vbig},
		0.,{},{1,1,0,0,0},
		Prepend[ConstantArray[0,35],1],1][[1]];
	outFunc=potFunc[
		{0.,0.,0.,0.,0.,#1,0.,0.,#1+#2},
		anuzero,{},{1,1,0,0,0},
		Prepend[ConstantArray[0,35],1],1][[ 1]] &;
	Return[outFunc]
];

(*
dFunc:
	Inputs: 
		potlib: Function from Potlib MathLink executable, 
		VBig: a suffeciently lagre number (for finding the energy of not connected atoms)
	Output:
		a function that takes R space coordintes and outputs the derivitave of the potencial energy of the system at the point
			as derivite of the X component of each atom
*)
dFunc[potFunc_,Vbig_]:=Module[{outFunc,anuzero},
	anuzero=potFunc[
		{0.,0.,0.,0.,0.,Vbig,0.,0.,2*Vbig},
		0.,{},{1,1,0,0,0},
		Prepend[ConstantArray[0,35],1],1][[1]];
	outFunc=potFunc[
		{0.,0.,0.,0.,0.,#1,0.,0.,#1+#2},
		anuzero,{},{1,1,0,0,0},
		Prepend[ConstantArray[0,35],1],1][[ {4,7,10}]] &;
	Return[outFunc]
];

(*
MWFunc:
	Input:
		Atoms in the system,
		R space potencial function,
		A suffeciently large number (for finding energy at not connected atoms)
	Output:
		a function that takes R space coordinates and outputs mass weighted Q space coordinates
		(R to Q space transformation)
*)
MWFunc[atoms_,RegFunc_,Rbig_]:=Module[
	{m,a,b,Beta,SigmaM,outFunc},
	
	m=Table[ElementData[atom,"AtomicMass"]/me,{atom,atoms}];
	SigmaM=Sum[i,{i,m}];
	a=(m[[ 1]] *(m[[ 2]] +m[[ 3]] )/SigmaM)^(1/2);
	b=(m[[ 3]] *(m[[ 2]] +m[[ 1]] )/SigmaM)^(1/2);
	Beta=ArcCos[Sqrt[m[[ 1]] *m[[ 3]] /((m[[ 2]] +m[[ 3]] )(m[[ 1]] +m[[ 2]] ))]];
	outFunc=Compile[{{q1,_Real},{q2,_Real}},{(1/a)*q1-(1/a)*Cot[Beta]*q2,(1/(b Sin[Beta]))*q2}];
	Return[outFunc]
];

(*
MWAngle: Inputs atoms in system and outputs angle of Q space transformation
*)
MWAngle[atoms_]:=Module[{m,SigmaM},
	m=Table[ElementData[atom,"AtomicMass"]/me,{atom,atoms}];
	SigmaM=Sum[i,{i,m}];
	Return[ArcCos[Sqrt[m[[ 1]] *m[[ 3]] /((m[[ 2]] +m[[ 3]] )(m[[ 1]] +m[[ 2]] ))]]];
];

(*
AtomMass: Takes atoms of a system and outputs a list of the sum of mass and a list of each each mass (in electron mass)
It is suggested to use this function to precomupte these variables for other functions
*)
AtomMass[atoms_]:=Module[{m,SigmaM},
	m=Table[ElementData[a,"AtomicMass"]/me,{a,atoms}];
	SigmaM=Sum[i,{i,m}];
	Return[{m,SigmaM}]
];

(*
dX2dR: takes the mass (with the same pattern as AtomMass output) 
	and X compontents of the derivitaves for each atom and converts to R space deriviates
*)
dX2dR[{m_,SigmaM_},{x1_,x2_,x3_}]:={(m[[ 1]] )/SigmaM(x1+x2+x3)-x1,x3-(m[[ 3]] )/SigmaM(x1+x2+x3)};

(*
dR2dQ: takes the mass (with the same pattern as AtomMass output) and R space deriviates and converts to Q space
*)
dR2dQ[{m_,SigmaM_},{r1_,r2_}]:=Module[
	{a=((m[[ 1]] *(m[[ 2]] +m[[ 3]] ))/SigmaM)^(1/2),
	b=(((m[[ 1]] +m[[ 2]] )*m[[ 3]] )/SigmaM)^(1/2),
	Beta=ArcCos[Sqrt[(m[[ 1]] *m[[ 3]] )/((m[[ 2]] +m[[ 3]] )(m[[ 1]] +m[[ 2]] ))]]},
	Return[{r1/a,-((b r1 Cot[Beta]-a r2 Csc[Beta])/(a b))}]
];



EndPackage[];
