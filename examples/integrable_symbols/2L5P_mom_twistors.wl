(* ::Package:: *)

(* ::Text:: *)
(*We build a set of integrable symbols using the 2-loop 5-point planar alphabet introduced in arXiv:1712.09610, as well as the first-entry condition, up to weight four.*)
(**)
(*We do this with the help of the Symbols` package.  Here we use a rational parametrization of the letters using momentum twistor variables.  For an analogous example where letters have square roots see symbols_mandelstam.wl*)
(* *)
(*The example can be easily adapted to the non-planar alphabets or to higher weights.*)


<<Symbols`


(* the planar alphabet *)
AP5={W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,W14,W15,W16,W17,W18,W19,W20,W26,W27,W28,W29,W30,W31};


(* the planar alphabet converted to mom. twistors *)
w2tomomtw = {W1->ex[1],W2->ex[1] ex[4],W3->(1/ex[2])(-ex[1] ex[2] ex[3]+ex[1] ex[4]+ex[1] ex[3] ex[4]+ex[1] ex[2] ex[3] ex[5]),W4->ex[1] ex[5],W5->ex[1] ex[2] ex[3]-ex[1] ex[3] ex[4]+ex[1] ex[3] ex[5],W6->(1/ex[2])(-ex[1] ex[2] ex[3]+ex[1] ex[4]+ex[1] ex[3] ex[4]+ex[1] ex[2] ex[5]+ex[1] ex[2] ex[3] ex[5]),W7->ex[1] ex[2] ex[3]-ex[1] ex[3] ex[4]+ex[1] ex[5]+ex[1] ex[3] ex[5],W8->ex[1]+ex[1] ex[2] ex[3]-ex[1] ex[3] ex[4]+ex[1] ex[3] ex[5],W9->ex[1]+ex[1] ex[4],W10->(1/ex[2])(-ex[1] ex[2] ex[3]+ex[1] ex[4]+ex[1] ex[2] ex[4]+ex[1] ex[3] ex[4]+ex[1] ex[2] ex[3] ex[5]),W11->ex[1]-ex[1] ex[5],W12->-ex[1] ex[2] ex[3]+ex[1] ex[4]+ex[1] ex[3] ex[4]-ex[1] ex[3] ex[5],W13->(1/ex[2])(-ex[1] ex[2]-ex[1] ex[2] ex[3]+ex[1] ex[4]+ex[1] ex[3] ex[4]+ex[1] ex[2] ex[3] ex[5]),W14->-ex[1] ex[4]+ex[1] ex[5],W15->(1/ex[2])(ex[1] ex[2] ex[3]+ex[1] ex[2]^2 ex[3]-ex[1] ex[4]-ex[1] ex[3] ex[4]-ex[1] ex[2] ex[3] ex[4]),W16->ex[1]+ex[1] ex[4]-ex[1] ex[5],W17->(1/ex[2])(-ex[1] ex[2] ex[3]-ex[1] ex[2]^2 ex[3]+ex[1] ex[4]+ex[1] ex[2] ex[4]+ex[1] ex[3] ex[4]+ex[1] ex[2] ex[3] ex[4]),W18->(1/ex[2])(-ex[1] ex[2]-ex[1] ex[2] ex[3]+ex[1] ex[4]+ex[1] ex[3] ex[4]+ex[1] ex[2] ex[5]+ex[1] ex[2] ex[3] ex[5]),W19->ex[1] ex[2] ex[3]-ex[1] ex[4]-ex[1] ex[3] ex[4]+ex[1] ex[5]+ex[1] ex[3] ex[5],W20->(1/ex[2])(ex[1] ex[2]+ex[1] ex[2] ex[3]+ex[1] ex[2]^2 ex[3]-ex[1] ex[4]-ex[1] ex[3] ex[4]-ex[1] ex[2] ex[3] ex[4]),W26->(ex[2]^2 ex[3]-ex[2] ex[4]-2 ex[2] ex[3] ex[4]+ex[4]^2+ex[3] ex[4]^2+ex[2] ex[3] ex[5]-ex[4] ex[5]-ex[3] ex[4] ex[5])/(ex[2] (1+ex[2]) ex[3] ex[5]),W27->(ex[2] (-ex[2] ex[3]+ex[3] ex[4]+ex[2] ex[3] ex[5]))/((1+ex[3]) ex[4] (-ex[2]+ex[4]-ex[5])),W28->(-ex[2] ex[3] ex[4]+ex[4]^2+ex[3] ex[4]^2-ex[4] ex[5]-ex[3] ex[4] ex[5])/(ex[2] (-ex[2] ex[3]+ex[4]+ex[3] ex[4]+ex[2] ex[3] ex[5])),W29->(ex[2]^2 ex[3]-ex[2] ex[4]-2 ex[2] ex[3] ex[4]+ex[4]^2+ex[3] ex[4]^2-ex[2]^2 ex[3] ex[5]+ex[2] ex[3] ex[4] ex[5])/((1+ex[3]+ex[2] ex[3]) ex[4] ex[5]),W30->(ex[2] (ex[3] ex[5]+ex[3] ex[4] ex[5]-ex[3] ex[5]^2))/((-ex[2]+ex[4]-ex[5]) (-ex[2] ex[3]+ex[4]+ex[3] ex[4]+ex[2] ex[3] ex[5])),W31->-(1/ex[2])ex[1]^2 (-ex[2]^2 ex[3]+ex[2] ex[4]+2 ex[2] ex[3] ex[4]-ex[4]^2-ex[3] ex[4]^2+ex[2]^2 ex[3] ex[5]+ex[4] ex[5]+ex[3] ex[4] ex[5])};


(* letters can be even or odd *)
AP5even = {W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,W14,W15,W16,W17,W18,W19,W20,W31};
AP5odd = {W26,W27,W28,W29,W30};


(* These functions decide if a symbol is even or odd respectively *)
selecteven = TrueQ[#==((#/.SS_Sym:>(odd^Mod[Length[Select[SS,MemberQ[AP5odd,#]&]],2]) SS)/.odd->0)]&;
selectodd = TrueQ[#==((#/.SS_Sym:>(odd^Mod[Length[Select[SS,MemberQ[AP5odd,#]&]]+1,2]) SS)/.odd->0)]&;


(* First, we compute the crossed derivatives WDlogC of letters *)
dw = WLogDerivs[AP5,ex/@Range[5],w2tomomtw];


(* Then we find all the linear relations between them.  This only needs to be
   done once per alphabet, independently of the weight of the symbols. *)
dwls = WDlogCLinearRelations[AP5,dw,ex/@Range[5]];


(* In order to compute integrable symbols, only the computed linear relations
 * in dwls are relevant.  The explicit expressions of the letters and the crossed
 * derivates are no longer needed. *)
Clear[w2tomomtw,dw];


(* At weight one, all letters are integrable,
   but we only select those satisfying the first-entry condition *)
w1even = (Sym/@AP5)[[Range[1,5]]];
w1odd = {};


(* Check against the expected numbers *)
Length[w1even]==5 && Length[w1odd]==0


(* Even and odd symbols at weight 2.  We do not apply the second entry condition
   but in principle we could by modifying the "Select" option. *)
w2even=IntegrableSymbolsWeight[2,AP5,w1even,ex/@Range[5],dwls,"Select"->selecteven];
w2odd=IntegrableSymbolsWeight[2,AP5,w1even,ex/@Range[5],dwls,"Select"->selectodd];
w2all = Join[w2even,w2odd];


(* Check against the expected numbers *)
Length[w2even]==25 && Length[w2odd]==0


(* Even and odd symbols at weight 3. *)
w3even=IntegrableSymbolsWeight[3,AP5,w2all,ex/@Range[5],dwls,"Select"->selecteven];
w3odd=IntegrableSymbolsWeight[3,AP5,w2all,ex/@Range[5],dwls,"Select"->selectodd];
w3all = Join[w3even,w3odd];


(* Check against the expected numbers *)
Length[w3even]==125 && Length[w3odd]==1


(* Even and odd symbols at weight 4. *)
w4even=IntegrableSymbolsWeight[4,AP5,w3all,ex/@Range[5],dwls,"Select"->selecteven];
w4odd=IntegrableSymbolsWeight[4,AP5,w3all,ex/@Range[5],dwls,"Select"->selectodd];
w4all = Join[w4even,w4odd];


(* Check against the expected numbers *)
Length[w4even]==635 && Length[w4odd]==16
