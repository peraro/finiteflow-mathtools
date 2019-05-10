(* ::Package:: *)

(* ::Text:: *)
(*We build a set of integrable symbols using the 2 - loop 5 - point planar alphabet introduced in arXiv : 1712.09610, as well as the first - entry condition, up to weight four.*)
(*  *)
(*We do this with the help of the Symbols` package.  Here we use a parametrization of the letters using Mandelstam variables.  This parametrization contains a square root which we do not rationalize (although we could).*)
(*   *)
(*The example can be easily adapted to the non - planar alphabets or to higher weights. *)


<<Symbols`


(* the planar alphabet *)
AP5={W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,W14,W15,W16,W17,W18,W19,W20,W26,W27,W28,W29,W30,W31};


(* The planar alphabet converted to Mandelstam invarinats.
   It conatains a square root which we do not rationalize. *)
sijvars = {s12,s23,s34,s45,s15};
w2tosij = {W1->s12,W2->s23,W3->s34,W4->s45,W5->s15,W6->s34+s45,W7->s15+s45,W8->s12+s15,W9->s12+s23,W10->s23+s34,W11->s12-s45,W12->-s15+s23,W13->-s12+s34,W14->-s23+s45,W15->s15-s34,W16->s12+s23-s45,W17->-s15+s23+s34,W18->-s12+s34+s45,W19->s15-s23+s45,W20->s12+s15-s34,W26->(-s12 s15+s12 s23-s23 s34-s15 s45+s34 s45-\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2))/(-s12 s15+s12 s23-s23 s34-s15 s45+s34 s45+\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2)),W27->(-s12 s15-s12 s23+s23 s34+s15 s45-s34 s45-\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2))/(-s12 s15-s12 s23+s23 s34+s15 s45-s34 s45+\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2)),W28->(s12 s15-s12 s23-s23 s34-s15 s45+s34 s45-\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2))/(s12 s15-s12 s23-s23 s34-s15 s45+s34 s45+\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2)),W29->(-s12 s15+s12 s23-s23 s34+s15 s45-s34 s45-\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2))/(-s12 s15+s12 s23-s23 s34+s15 s45-s34 s45+\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2)),W30->(s12 s15-s12 s23+s23 s34-s15 s45-s34 s45-\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2))/(s12 s15-s12 s23+s23 s34-s15 s45-s34 s45+\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2)),W31->\[Sqrt](s12^2 s15^2-2 s12^2 s15 s23+s12^2 s23^2+2 s12 s15 s23 s34-2 s12 s23^2 s34+s23^2 s34^2-2 s12 s15^2 s45+2 s12 s15 s23 s45+2 s12 s15 s34 s45+2 s12 s23 s34 s45+2 s15 s23 s34 s45-2 s23 s34^2 s45+s15^2 s45^2-2 s15 s34 s45^2+s34^2 s45^2)};


(* letters can be even or odd *)
AP5even = {W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,W14,W15,W16,W17,W18,W19,W20,W31};
AP5odd = {W26,W27,W28,W29,W30};


(* These functions decide if a symbol is even or odd respectively *)
selecteven = TrueQ[#==((#/.SS_Sym:>(odd^Mod[Length[Select[SS,MemberQ[AP5odd,#]&]],2]) SS)/.odd->0)]&;
selectodd = TrueQ[#==((#/.SS_Sym:>(odd^Mod[Length[Select[SS,MemberQ[AP5odd,#]&]]+1,2]) SS)/.odd->0)]&;


(* First, we compute the crossed derivatives WDlogC of letters *)
dw = WLogDerivs[AP5,sijvars,w2tosij];


(* Then we find all the linear relations between them.  This only needs to be
   done once per alphabet, independently of the weight of the symbols. *)
dwls = WDlogCLinearRelations[AP5,dw,sijvars];


(* In order to compute integrable symbols, only the computed linear relations
 * in dwls are relevant.  The explicit expressions of the letters and the crossed
 * derivates are no longer needed. *)
Clear[w2tosij,dw];


(* At weight one, all letters are integrable,
   but we only select those satisfying the first-entry condition *)
w1even = (Sym/@AP5)[[Range[1,5]]];
w1odd = {};


(* Check against the expected numbers *)
Length[w1even]==5 && Length[w1odd]==0


(* Even and odd symbols at weight 2.  We do not apply the second entry condition
   but in principle we could by modifying the "Select" option. *)
w2even=IntegrableSymbolsWeight[2,AP5,w1even,sijvars,dwls,"Select"->selecteven];
w2odd=IntegrableSymbolsWeight[2,AP5,w1even,sijvars,dwls,"Select"->selectodd];
w2all = Join[w2even,w2odd];


(* Check against the expected numbers *)
Length[w2even]==25 && Length[w2odd]==0


(* Even and odd symbols at weight 3. *)
w3even=IntegrableSymbolsWeight[3,AP5,w2all,sijvars,dwls,"Select"->selecteven];
w3odd=IntegrableSymbolsWeight[3,AP5,w2all,sijvars,dwls,"Select"->selectodd];
w3all = Join[w3even,w3odd];


(* Check against the expected numbers *)
Length[w3even]==125 && Length[w3odd]==1


(* Even and odd symbols at weight 4. *)
w4even=IntegrableSymbolsWeight[4,AP5,w3all,sijvars,dwls,"Select"->selecteven];
w4odd=IntegrableSymbolsWeight[4,AP5,w3all,sijvars,dwls,"Select"->selectodd];
w4all = Join[w4even,w4odd];


(* Check against the expected numbers *)
Length[w4even]==635 && Length[w4odd]==16
