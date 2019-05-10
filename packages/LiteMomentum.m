(* ::Package:: *)

BeginPackage[ "LiteMomentum`"]


mm::usage = "\
mm[p] represents a momentum p.
mm[p][mu] is an indexed momentum p, with free index mu."
mp::usage = "mp[x,y] is the momentum product between x and y.  Both x and y should be linear combinations of momenta mm[p]."
mp2::usage = "Momentum square: same as mp[#,#]."

mmp::usage = "mmp[a,b] = mp[mm[a],mm[b]]"
mmp2::usage = "mmp2[a] = mp2[mm[a]]"
ExpandMP::usage = "ExpandMP[expr] expands the mp products in expr."

MetricD::usage = "Number of dimensions for metric tensor gD"
MetricEps::usage = "The dimensional regulator.  MetricD -> 4 - 2 MetricEps"
gD::usage = "D-dimensional metric tensor (the dimension is MetricD)"
g4::usage = "4-dimensional metric tensor (the dimension is 4)"
gEps::usage = "Extra-dimensional metric tensor (its dimension is -2 MetricEps)"
mpEps::usage = "Represents (-2 epsilon)-dimensional momentum product."
mp2Eps::usage = "Same as mpEps[#,#]."
mmpEps::usage = "mmpEps[a,b] = mpEps[mm[a],mm[b]]."
mmp2Eps::usage = "mmp2Eps[a] = mp2Eps[mm[a]]."
LContract::usage = "\
LContract[expr] contracts the indexed expression expr.
LContract[expr,mpD,mp4] contracts the indexed expression expr, using the function mpD as D-dimensional scalar product, and mp4 as four-dimensional scalar product."
LEpsContract::usage = "\
CLEpsContract[expr] contracts the extra-dimensions in the indexed expression expr.
CLEpsContract[expr,mpEps] contracts the extra-dimensions in the indexed expression expr, using the function mpEps as the extra-dimensional scalar product."
MomDerivative::usage = "MomDerivative[expr,q,mu] returns the derivative of expr with respect to mm[q][mu]."
IndexedMomExpr::usage = "IndexedMomExpr[epxr,mu] transforms a linear combination of momenta in an indexed expression with free index mu."


Begin[ "`Private`"]


Clear[mm];


(* scalar product *)
SetAttributes[mp,Orderless];
mp[p1_Plus, p2_] := mp[#, p2] & /@ p1;
mp[p2_, coeff__ mm[p1__]] := coeff mp[mm[p1], p2];
mp[0,p_]:=0;
mp2[p_] := mp[p, p];

ExpandMP[expr_]:= expr/. {mp[a_,b_]:>mp[Expand[a,mm],Expand[b,mm]]};


(* shorthands *)
mmp[a_,b_]:=mp[mm[a],mm[b]];
mmp2[a_]:=mp[mm[a],mm[a]];


(* metric tensors in D and 4 dimensions *)
SetAttributes[gD, Orderless];
SetAttributes[g4, Orderless];


SetAttributes[gEps, Orderless];


LContract[expr_,mpD_,mp4_]:= (Expand[expr]/. {gD[mm[k1_],mm[k2_]]:>mp[mm[k1],mm[k2]],
											  gD[mm[k_],mm[k_]]:>mp[mm[k],mm[k]],
											  g4[mm[k1_],mm[k2_]]:>mp[mm[k1],mm[k2]],
	                                           g4[mm[k_],mm[k_]]:>mp[mm[k],mm[k]]}) //. {
	gD[mm[k_],mu_]:>mm[k][mu],
	gD[mu_,mu_] :> MetricD,
	g4[mm[k_],mu_]:>mm[k][mu],
	g4[mu_,mu_] :> 4,
	gD[mu_,nu_]^2 :> MetricD,
    gD[mu_,nu_]gD[nu_,sigma_] :> gD[mu,sigma],
    g4[mu_,nu_]g4[nu_,sigma_] :> g4[mu,sigma],
    gD[mu_,nu_]g4[nu_,sigma_] :> g4[mu,sigma],
	g4[mu_,nu_]^2 :> 4,
	mm[k_][mu_]^2 :> mp[mm[k],mm[k]],
	mm[k1_][mu_] mm[k2_][mu_] :> mp[mm[k1],mm[k2]],
	gD[mu_,nu_] mm[k1_][mu_] mm[k2_][nu_] :> mpD[mm[k1],mm[k2]],
	g4[mu_,nu_] mm[k1_][mu_] mm[k2_][nu_] :> mp4[mm[k1],mm[k2]]}
LContract[expr_]:=LContract[expr,mp,mp] //. {gD[mu_,nu_] mm[k1_][mu_] :> mm[k1][nu],
											 g4[mu_,nu_] mm[k1_][mu_] :> mm[k1][nu]}


(* (-2 epsilon)-dim. scalar product *)
SetAttributes[mpEps,Orderless];
mpEps[p1_Plus, p2_] := mpEps[#, p2] & /@ p1;
mpEps[p2_, coeff__ mm[p1__]] := coeff mpEps[mm[p1], p2];
mpEps[0,p_]:=0;
mp2Eps[p_] := mpEps[p, p];


(* shorthands *)
mmpEps[a_,b_]:=mpEps[mm[a],mm[b]];
mmp2Eps[a_]:=mpEps[mm[a],mm[a]];


LEpsContract[expr_,mpEps_:mpEps]:= Expand[expr]/.{gEps[mm[a_],mm[b_]]:>mpEps[mm[a],mm[b]]} //. {
	gEps[mu_,nu_]^2 :> -2 MetricEps,
	gD[mu_,nu_]g4[nu_,rho_]:> g4[mu,rho],
	gD[mu_,nu_]gEps[nu_,rho_]:> gEps[mu,rho],
	g4[mu_,nu_]gEps[nu_,rho_]:> 0,
	gEps[mu_,mu_] :> -2 MetricEps,
	gEps[mu_,nu_] mm[k1_][mu_] mm[k2_][nu_] :> mpEps[mm[k1],mm[k2]]}


(* Derivatives w.r.t. mmenta *)
MomDerivative[expr_,mm[q_][mu_],dim_,gmetric_,finalcontraction_]:=
 ((D[expr /. {mm[q][nu_] :> IndexedMomentum[mm[q], nu]},mm[q]]/.
	{Derivative[1, 0][mp][a_, b_] :> b[mu],
	 Derivative[0, 1][mp][a_, b_] :> a[mu],
     Derivative[1, 0][IndexedMomentum][mm[q], mu] :> dim,  
     Derivative[1, 0][IndexedMomentum][mm[q], nu_] :> gmetric[mu,nu]}) /.
    {IndexedMomentum[mm[q],nu_] :> mm[q][nu]} ) // finalcontraction
MomDerivative[expr_,mm[q_][mu_]]:=MomDerivative[expr,mm[q][mu],MetricD,gD,LContract]


IndexedMomExpr[expr_Plus,mu_]:=IndexedMomExpr[#,mu]&/@expr;
IndexedMomExpr[coeff__ mm[a__],mu_]:=coeff IndexedMomExpr[mm[a],mu];
IndexedMomExpr[mm[a__],mu_]:=mm[a][mu];


mm/:MakeBoxes[mm[a_][mu_],StandardForm]:=SuperscriptBox[MakeBoxes[a,StandardForm],MakeBoxes[mu,StandardForm]];
mm/:MakeBoxes[mm[a_Integer],StandardForm]:=OverscriptBox[MakeBoxes[a,StandardForm],"\[RightVector]"];
mm/:MakeBoxes[mm[a_],StandardForm]:=RowBox[{MakeBoxes[a,StandardForm]}];
mp/:MakeBoxes[mp[mm[a_],mm[b_]],StandardForm]:=RowBox[{"(",MakeBoxes[a,StandardForm],"\[CenterDot]",
															MakeBoxes[b,StandardForm],")"}];
gD/:MakeBoxes[gD[mu_,nu_],StandardForm]:=SuperscriptBox[SubscriptBox["g","D"],RowBox[{MakeBoxes[mu,StandardForm],MakeBoxes[nu,StandardForm]}]];
gEps/:MakeBoxes[gEps[mu_,nu_],StandardForm]:=SuperscriptBox[SubscriptBox["g","-2\[Epsilon]"],RowBox[{MakeBoxes[mu,StandardForm],MakeBoxes[nu,StandardForm]}]];
g4/:MakeBoxes[g4[mu_,nu_],StandardForm]:=SuperscriptBox[SubscriptBox["g","4"],RowBox[{MakeBoxes[mu,StandardForm],MakeBoxes[nu,StandardForm]}]];


End[] (* `Private` *)


EndPackage[]
