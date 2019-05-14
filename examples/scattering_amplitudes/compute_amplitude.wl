(* ::Package:: *)

(* ::Text:: *)
(*In this file, we compute 4-photon helicity amplitudes.  The input is an integrand in terms of loop-momenta, external momenta and polarization vectors. The calculation is performed using integrand reduction, transverse integration, and IBP identities. We directly reconstruct the coefficients of the polylogarithmic functions appearing in the result.*)
(**)
(*Although this is just a one-loop amplitude, the same setup can be used for multi-loop amplitudes as well.*)
(**)
(*A much simpler (albeit less complete) version of this calculation can be found in the file compute_amplitude_simpler.wl.  We recommend looking at that file instead if you are not interested in integrand reduction or transverse integration methods for computing amplitudes.*)


<<FiniteFlow`


<<LiteMomentum`


<<"topology.wl"


(* ::Text:: *)
(*We select the helicity configuration of the external photons. This may be changed to any different helicity configuration, simply by editing this line.*)


helicity = {"-","+","-","+"};


(* ::Text:: *)
(*The numerator is stored in an external file.*)


numerator = Import["numerator.m"] /. mf -> mf2^(1/2);


(* ::Text:: *)
(*We express momenta, polarization vectors, and spinors in terms of a momentum twistor parametrization using s and t.  The following functions explicitly define their components.  It should be noted that, by doing this, we are neglecting an helicity-dependent overall prefactor, which is however easy to recover (see e.g. arXiv:1605.02172 or arXiv:1608.01902).*)


(* spinor components, {lambda, lambdatilde}
   using momentum twistor parametrization *)
spinor[p1]={{1,0},{1,1}};
spinor[p2]={{0,1},{0,s}};
spinor[p3]={{1/s,1},{t, -s}};
spinor[p4]={{1/s + 1/t,1},{-t,0}};


(* spinor products <a,b> and [a,b] from mom. twistors *)
mtspa[a_,b_]:= spinor[a][[1,1]] spinor[b][[1,2]]-spinor[a][[1,2]] spinor[b][[1,1]];
mtspb[a_,b_]:= spinor[a][[2,2]] spinor[b][[2,1]]-spinor[a][[2,1]] spinor[b][[2,2]];


(* momenta in light cone components from spinors *)
mtlcmom[a_,b_]:={spinor[a][[1,1]]spinor[b][[2,1]],
                 spinor[a][[1,2]]spinor[b][[2,2]],
                 spinor[a][[1,2]]spinor[b][[2,1]],
                 spinor[a][[1,1]]spinor[b][[2,2]]};
mtlcmom[a_]:=mtlcmom[a,a];


(* polarization vectors, without the 1/Sqrt[2] prefactor *)
(* a 1/4 should be included in the overall missing prefactor *)
mtpolvec["+"][p_,ref_]:=mtlcmom[ref,p]/mtspa[ref,p];
mtpolvec["-"][p_,ref_]:=mtlcmom[p,ref]/mtspb[p,ref];


(* external momenta *)
Table[lcmom[pp] = mtlcmom[pp],{pp,{p1,p2,p3,p4}}]//Factor


(* polarization vectors *)
lcmom[eps1] = mtpolvec[helicity[[1]]][p1,p2]//Factor
lcmom[eps2] = mtpolvec[helicity[[2]]][p2,p1]//Factor
lcmom[eps3] = mtpolvec[helicity[[3]]][p3,p1]//Factor
lcmom[eps4] = mtpolvec[helicity[[4]]][p4,p1]//Factor


(* scalar products in light cone components *)
mtlcmp[a_,b_]:=1/2 (a[[1]] b[[2]] + a[[2]] b[[1]] - a[[3]] b[[4]] - a[[4]] b[[3]]);


(* ::Text:: *)
(*We simplify our implementation by defining a global basis of momenta, {p1, p2, v3, v4} and expanding the 4-dim. part of the loop momentum k as*)
(**)
(*  k = y1 p1 + y2 p2 + y3 v3 + y4 v4*)
(**)
(*where v3,v4 are massless and orthogonal to p1,p2.  The extra-dimensional part of the integrand is encoded by the variable mu2, that is (minus) the -2*eps-dimensional part of the scalar product k.k*)
(**)
(*Here we explicitly write the integrand in terms of y1,...,y4 and mu2.  Alternatively, we could also leave it as it is and pass rules for converting scalar products involving k and polarization vectors into these variables to the FFAlgLinearFit procedure used below.*)
(**)
(*We also define vort, which is orthogonal to p1, p2, p3, p4*)


v3 = mtlcmom[p1,p2] mtspa[p2,p3] / mtspa[p1,p3];
v4 = mtlcmom[p2,p1] mtspa[p1,p3] / mtspa[p2,p3];
vort = 1/s mtlcmp[v3,lcmom[p3]] v4 - 1/s mtlcmp[v4,lcmom[p3]] v3 // Factor;


lcmom[k] = y1 lcmom[p1] + y2 lcmom[p2] + y3 v3 + y4 v4 // Together


(* k.vort *)
kvortexpr = Collect[mtlcmp[lcmom[k],vort],y3|y4,Together]


invariants = {t,mf2}


lcmom[vperp]=vort


(* ::Text:: *)
(*Here we use the easiest integrand basis to build (which is also the most IBP friendly).  This is often a good idea, especially at higher loops.  The basis uses, as irreducible scalar products, monomials in the loop denominators of the master topology, and k.vort.  The scalar product k.vort will then be integrated out, and everything else will be reduced via IBPs.*)


(* ::Text:: *)
(*For each of the three permutations of the external legs, we define a complete integrand reduction of the corresponding diagram.*)


FFNewGraph[gired];
FFGraphInputVars[gired,in,invariants];
FFAlgRatNumEval[gired,one,{1}];

Do[
  Module[{applypermutation,loopmoms,num,loopdenoms,allmps,
          cutprops,uncutprops,lincuteqs,cutsol,uncutpropsyms,kvortcut,
          supercuts,subtractions,nonzerosubterms,
          integrand,isps,maxrank,isppowers,coeffs,zerocoeffs,learn},
    
    (* build list of substitution rules for the permutation *)
    applypermutation = Inner[Rule,
                             Join[{p1,p2,p3,p4},{eps1,eps2,eps3,eps4}],
                             Join[{p1,p2,p3,p4}[[permutation]],{eps1,eps2,eps3,eps4}[[permutation]]],
                             List];
  
    (* loop momenta *)
    loopmoms = {
      mm[l1] -> mm[k],
      mm[l2] -> mm[k]-mm[p1],
      mm[l3] -> mm[k]-mm[p1]-mm[p2],
      mm[l4] -> mm[k]+mm[p4]
    } /. applypermutation;
    
    (* numerator *)
    num = numerator /. applypermutation;
    
    (* loop denominators *)
    loopdenoms = (mmp2[#]-mf2)&/@{l1,l2,l3,l4};
    
    (* compute all scalar products... *)
    allmps = Inner[Rule,#,#/.loopmoms/.Replacements/.mmp[k,k]->mmp[k,k]-mu2/.mp[mm[a_],mm[b_]]:>mtlcmp[lcmom[a],lcmom[b]],List]&[Union[Cases[{num,loopdenoms},_mp,Infinity]]]//Factor;
    (* ... and substitute them in the numerator and denominators *)
    num = Together[num /. allmps];
    loopdenoms = Collect[loopdenoms /. allmps,y1|y2|y3|y4|mu2,Together];
    
    (* list cuts which give non-vanishing contributions *)
    cuts = If[TrueQ[mf2==0],
      {{1,2,3,4},{1,2,3},{1,2,4},{1,3,4},{2,3,4},{1,3},{2,4}},
      {{1,2,3,4},{1,2,3},{1,2,4},{1,3,4},{2,3,4},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4},{1},{2},{3},{4}}
     ];
     
    Do[
    
      (* find the cut solution *)
      cutprops = loopdenoms[[thiscut]];
      lincuteqs = Collect[(#-cutprops[[1]])&/@cutprops[[2;;]],y1|y2|y3|y4|mu2,Together];
      cutsol = Together[FFDenseSolve[(#==0)&/@lincuteqs,{y1,y2,y3,y4}]];
      cutsol = Join[cutsol,Together[Solve[Collect[cutprops[[1]]/.cutsol,mu2,Together]==0,mu2][[1]]]];
      
      (* uncut propagators and ortogonal scalar products on the cut solutions *)
      uncutpropsyms = {D1,D2,D3,D4}[[Complement[Range[4],thiscut]]];
      uncutprops = Together[loopdenoms[[Complement[Range[4],thiscut]]]/.cutsol];
      kvortcut = Together[kvortexpr /. cutsol];
      
      (* subtractions *)
      supercuts = Select[cuts,(Length[#]>Length[thiscut]&&SubsetQ[#,thiscut])&];
      subtractions = Join@@Table[-(Times@@({D1,D2,D3,D4}[[thiscut]]))/(Times@@({D1,D2,D3,D4}[[subcut]]))delta[permutation,subcut] ,{subcut,supercuts}];

      (* integrand and delta *)
      integrand = num/(Times@@uncutpropsyms);
      isps = Join[uncutpropsyms,{kvort}];
      maxrank = Length[thiscut];
      isppowers = (Join@@(LIBPPowers[#,Length[isps]]&/@Range[maxrank,0,-1]));
      delta[permutation,thiscut] = (Times@@(isps^#))&/@isppowers;
      coeffs = cc[permutation,thiscut,#]&/@isppowers;
      
      (* define node in the graph *)
      FFAlgLinearFit[gired,cut[permutation,thiscut],
                        Join[{in,one},cut[permutation,#]&/@supercuts],
                        invariants,
                        delta[permutation,thiscut],
                        Join[{integrand},subtractions]/.s->1/.(#->0&/@({D1,D2,D3,D4}[[thiscut]])),
                        Complement[{y1,y2,y3,y4,mu2},First/@cutsol],
                        coeffs,
                        "Substitutions"->(Join[Inner[Rule,uncutpropsyms,uncutprops,List],{kvort->kvortcut},cutsol]/.s->1)];
      FFGraphOutput[gired,cut[permutation,thiscut]];
      learn=FFDenseSolverLearn[gired,coeffs];
      If[!TrueQ[learn[[0]]==List && (("IndepVars"/.learn)=={})],Abort[]];  (* <-- checking system is determined *)
      depcoeffs[permutation,thiscut] = "DepVars" /. learn;
      delta[permutation,thiscut] = (Times@@(isps^#[[-1]]))&/@depcoeffs[permutation,thiscut]; (* <-- update delta removing non-zero coeffs *)

    ,{thiscut,cuts}];
  ];
,{permutation, permutations}];


allcuts = Join@@Table[cut[permutation,thiscut],{permutation, permutations},{thiscut,cuts}];


alldepcoeffs  = Table[depcoeffs@@thiscut,{thiscut,allcuts}];


(* map permutations of external legs to the right master topology *)
permutationtopo = Table[permutations[[ii]]->{box1,box2,box3}[[ii]],{ii,Length[permutations]}]


(* ::Text:: *)
(*Replace denominators and k.vort by standard Feynman integrals, of the form j[family,a1,a2,a3,a4,a5], where the last entry, to be removed later via transverse integration, represents the (inverse of the) k.vort scalar product.*)


permutationsubst = (#[[1]]-> {(D1 -> j[#[[2]],-1,0,0,0,0]),
                              (D2 -> j[#[[2]],0,-1,0,0,0]),
                              (D3 -> j[#[[2]],0,0,-1,0,0]),
                              (D4 -> j[#[[2]],0,0,0,-1,0]),
                              kvort->j[#[[2]],0,0,0,0,-1]})&/@permutationtopo;


(* ::Text:: *)
(*The complete integrand basis:*)


fullintegrandbasis = Join@@Table[(delta@@thiscut)/(Times@@({D1,D2,D3,D4}[[thiscut[[2]]]]))/.(thiscut[[1]]/.permutationsubst),{thiscut,allcuts}];


(* ::Text:: *)
(*Join all the coefficients in the output of the graph gired*)


FFAlgTake[gired,takeall,allcuts,alldepcoeffs->(Join@@alldepcoeffs)];
FFGraphOutput[gired,takeall];


(* ::Text:: *)
(*The graph now looks awfully complicated, but it can still be evaluated very efficiently*)


FFGraphDraw[gired]


(* ::Text:: *)
(*In principle, we can easily reconstruct the full integrand by un-commenting the two lines below.  While in this example that would be very fast and easy, in general this is not a good idea, since it is only an intermediate result, which in more complex processes is very complicated and often not worth reconstructing.*)


(*rec = FFReconstructFunction[gired,invariants];;
fullintegrand = rec.fullintegrandbasis;*)


(* ::Text:: *)
(*Now, before IBP reduction, we want to get rid of scalar products k.vort^a, which are now represented by the 5th power a=a5 of the integrals j[family,a1,a2,a3,a4,a5].  This is easily done, at any loop order, using transverse integration techniques. After this, we will only have integrals of the form j[family,a1,a2,a3,a4], which we will reduced via IBPs.*)


(* ::Text:: *)
(*We compute, for each family, lamda2, defined as*)
(**)
(*     k^2 = kparallel^2 - lambda2*)
(**)
(* where kparallel is the projection of the loop momentum k in the parallel space spanned by the external momenta.*)
(**)
(*Scalar products of the form k.vort^a can be integrated out in terms of lambda2.*)


kparallel = y1 mm[p1] + y2 mm[p2] + y3 mm[p3];
lambda2expr = mp2[kparallel]-mmp2[k] /. Replacements;

Do[Module[{eqs,ysol},
   (* write equations for y1, y2, y3 *)
   eqs = Collect[mp[kparallel,mm[#]]==mmp[k,#]&/@{p1,p2,p3} /. Replacements /. LIBPSpsToJ[fam],_j|y1|y2|y3,Together];
   ysol = FFDenseSolve[eqs,Join[{y1,y2,y3},Union[Cases[eqs,_j,Infinity]]]];
   lamda2[fam] = Collect[lambda2expr /. LIBPSpsToJ[fam] /. ysol,_j,Together];
 ]
,{fam,{box1,box2,box3}}]


(* ::Text:: *)
(*These are the generic results of transverse integration, already encoded in LiteIBP.  These are independent of the topology.*)


dperp = 4 - 2 eps - 3; (* <-- dimension of the orthogonal space *)
gentransverseintegration = {
  kvort->0,
  kvort^2->(LIBPgtInv[dperp,2].LIBPgt[k,k] /. mp[k,k]->-lamda2).LIBPgt[v,v] /. mp[v,v]->vort2,
  kvort^3->0,
  kvort^4->(LIBPgtInv[dperp,4].LIBPgt[k,k,k,k] /. mp[k,k]->-lamda2).LIBPgt[v,v,v,v] /. mp[v,v]->vort2//Together
} /. vort2 -> mtlcmp[vort,vort]//Together


(* ::Text:: *)
(*Now we apply the integration to the individual topologies.*)


Do[
  transverseintegration[fam] = Collect[Together[gentransverseintegration /. lamda2 -> lamda2[fam]],_j,Together];
,{fam,{box1,box2,box3}}]


(* ::Text:: *)
(*This is the integrand basis after transverse integration*)


transverseintegratedbasis = LIBPEliminateZeroSectors[Collect[fullintegrandbasis/.(j[fam_,a__,a5_]:>(j[fam,a]kvort^-a5/.transverseintegration[fam])),_j,Together]];


(* ::Text:: *)
(*The integrand basis after transverse integration depends on epsilon, therefore we need to make a new graph for the full reduction.*)
(**)
(* Because eventually we want to expand in epsilon, epsilon needs to be the first variable.*)


FFNewGraph[gredfull];
FFGraphInputVars[gredfull,in,Join[{eps},invariants]]


(* ::Text:: *)
(*We load the IBP system of equations.  More comments of IBPs can be found in the example "examples/differential_equations/differential_equations.wl".*)


allintegrals = LIBPGetAllInts[FileNames["ibps/ints_*.mx"]];
neededintegrals = LIBPIntegralsIn[transverseintegratedbasis];

(* create system.json, which collects info on the whole IBP system *)
LIBPWriteSystemJSON[
  FileNames["ibps/sids_*.json"],
  allintegrals,neededintegrals,
  {eps,t,mf2},
  "FileName"->"system.json"];
                         
FFAlgJSONSparseSolver[gredfull,ibps,{in},"system.json"];
FFSolverOnlyHomogeneous[gredfull,ibps];
FFGraphOutput[gredfull,ibps];


(* ::Text:: *)
(*Learning phase for IBP system.  This yields the list of master integrals:*)


ibplearn = FFSparseSolverLearn[gredfull,allintegrals];
nonmis = "DepVars"  /. ibplearn;
mis = "IndepVars" /. ibplearn;
Print["Master integrals = ",mis];


(* ::Text:: *)
(*We  now run the mark-and-sweep algorithm (highly recommended).*)


FFSparseSolverMarkAndSweepEqs[gredfull,ibps]
FFSparseSolverDeleteUnneededEqs[gredfull,ibps];


(* ::Text:: *)
(*Full list of relevant IBP integrals, sorted with the masters on the right.*)


ibpintegrals = Join[nonmis,mis];


(* ::Text:: *)
(*We want to implement the IBP reduction as a matrix multiplication.  Because in this example the master integrals also appear in the unreduced amplitude, we must complete the reduction returned by the IBP node with the trivial reduction of the masters to themselves, the latter defined by an Identity matrix.*)


FFAlgRatNumEval[gredfull,misred,Join@@IdentityMatrix[Length[mis]]];
FFAlgChain[gredfull,ibpsfull,{ibps,misred}]


(* ::Text:: *)
(*The graph so far*)


FFGraphDraw[gredfull]


(* ::Text:: *)
(*Some parts of the graph do not depend on epsilon, hence we define  their inputs.*)


FFAlgSlice[gredfull,invonly,{in},2]


(* ::Text:: *)
(*An equivalent alternative to the former is*)


(*FFAlgTake[gredfull,invonly,{in},{{eps,t,mf2}}\[Rule]{t,mf2}]*)


(* ::Text:: *)
(*We take the coefficients of the integrand reduction as a subgraph.  Since it is independent of epsilon, we make it a Memoized Subgraph, so that it is not re-evaluated if only epsilon has changed since the previous evaluation.   This is especially relevant during the Laurent expansion in epsilon.*)


FFAlgMemoizedSubgraph[gredfull,ired,{invonly},gired]


(* ::Text:: *)
(*We want to multiply the output coefficients of the integrand reduction by the matrix of the transverse integration rules, defined from transverseintegratedbasis.  For a complex process the conversion matrix may be huge and very sparse, so we use a sparse matrix multiplication.  In order define the sparse matrix with the orthogonalintegration, we use*)
(**)
(*  FFSparseRowRules[epxression,ibpintegrals]*)
(**)
(*which gives us, for each row, the list of non-zero columns and their analytic coefficients.*)


integralscoeffs = Table[
     FFSparseRowRules[transverseintegratedbasis[[ii]],ibpintegrals],
   {ii,Length[transverseintegratedbasis]}
];


(* ::Text:: *)
(*The analytic coefficients are joined together and evaluated in a node.*)


FFAlgRatFunEval[gredfull,icoeff,{in},Join[{eps},invariants],(Join@@Map[#[[2]]&/@#&,integralscoeffs])/.s->1]


(* ::Text:: *)
(*We now define the sparse matrix multiplication:*)


transverseintcolumns = (First/@#)&/@integralscoeffs;
FFAlgSparseMatMul[gredfull,unreduced,{ired,icoeff},
                     1,FFNParsOut[gredfull,ired],Length[ibpintegrals],
                     {Range[FFNParsOut[gredfull,ired]]},transverseintcolumns]


(* ::Text:: *)
(*The following reconstructs the unreduced amplitude but, again, it's an intermediate result we are not interested in.*)


(*FFGraphOutput[gredfull,unreduced];
unreducedamp = FFReconstructFunction[gredfull,Join[{eps},invariants]].ibpintegrals;*)


(* ::Text:: *)
(*One more matrix multiplication applies the IBP solutions to obtain the reduced result*)


FFAlgMatMul[gredfull,reduced,{unreduced,ibpsfull},1,Length[ibpintegrals],Length[mis]]


(* ::Text:: *)
(*The graph gredfull now gives the amplitude as a linear combination of master integrals.*)


FFGraphDraw[gredfull]


(* ::Text:: *)
(*The following reconstructs the amplitude as a linear combination of master integrals.  This is however not our goal here, since we want to expand the amplitude into polylogarithms.*)


(*FFGraphOutput[gredfull,reduced];
reducedamp = FFReconstructFunction[gredfull,Join[{eps},invariants]].mis;*)


(* ::Text:: *)
(*Now we substitute analytic expressions for the master integrals, which are valid up to O (eps^0)*)
(**)
(*Analytic expressions for the master integrals, taken from https://arxiv.org/pdf/0804.0749.pdf*)
(*In particular:*)
(* - boxfun[X,Y,Z] is H[X,Y]+H[X,Z], with H defined in Eq. (C.28)*)
(* - trianglefun is Log[(1-rho)/(1+rho)], with rho = Sqrt[1-4 mf2/s]*)
(*  - bubblefun is f defined in Eq. (C.15-17)*)


(* the usual prefactor *)
rgamma[eps_] := Gamma[1+eps]Gamma[1-eps]^2/Gamma[1 - 2 eps];


boxint[s_,t_,mf2_]:=-1/(s t) ( boxfun[(s+t) mf2/(s t), mf2/s, mf2/t] );
trint[s_,mf2_]:=-1/(2 s) trianglefun[mf2/s];
bubble[s_,mf2_]:=Normal[Series[1/rgamma[eps](mf2^(- eps) Gamma[1+eps]/eps+2 + bubblefun[mf2/s]),{eps,0,0}]];
tadpole[mf2_]:=Normal[Series[-1/rgamma[eps] Gamma[1-(4-2 eps)/2]/Gamma[1] 1/(mf2)^(1-(4-2 eps)/2),{eps,0,0}]];


box1mis = {
 j[box1,1,1,1,1] -> boxint[s,t,mf2],
 j[box1,0,1,1,1] -> trint[t,mf2],
 j[box1,1,0,1,1] -> trint[s,mf2],
 j[box1,0,1,0,1] -> bubble[t,mf2],
 j[box1,1,0,1,0] -> bubble[s,mf2],
 j[box1,0,0,0,1] -> tadpole[mf2]
};


box1mis


(* All analytic master integrals, including permutations *)
analyticmis = Join[
     box1mis,
     box1mis /. box1->box2 /. InvariantsPermutations[[2]],
     box1mis /. box1->box3 /. InvariantsPermutations[[3]]];
analyticmis = mis /. analyticmis;


(* get list of special functions appearing in the master integrals *)
funcspattern = _Log|_boxfun|_trianglefun|_bubblefun;
allfuncs = Union[Cases[analyticmis,funcspattern,Infinity]]


(* ::Text:: *)
(*We create a list of rules for expressing master integrals in terms of special functions*)


misrules  = Together[(Times@@(allfuncs^#[[1]]))->#[[2]]&/@CoefficientRules[#,allfuncs]&/@analyticmis];


(* ::Text:: *)
(*and obtain a function basis from the rules *)


functionbasis = Union@@(First/@#&/@misrules)


(* ::Text:: *)
(*We build a matrix converting the master integrals into the function basis *)


mis2basis = (functionbasis/.Join[#,(#->0)&/@Complement[functionbasis,First/@#]])&/@misrules


(* ::Text:: *)
(*put it in a graph*)


FFAlgRatFunEval[gredfull,mis2fs,{in},Join[{eps},invariants],Join@@mis2basis /. s -> 1]


(* ::Text:: *)
(*and multiply the reduced amplitude with it*)


FFAlgMatMul[gredfull,amp2fs,{reduced,mis2fs},1,Length[mis],Length[functionbasis]]


FFGraphOutput[gredfull,amp2fs]


(* ::Text:: *)
(*The final graph, which computes the coefficients of the polylog functions.  Now we want to expand these coefficients in epsilon*)


FFGraphDraw[gredfull]


(* ::Text:: *)
(*This graph defines the Laurent expansion*)


FFNewGraph[gexpansion,in,invariants];
FFAlgLaurent[gexpansion,laurent,{in},gredfull,0]


(* ::Text:: *)
(*The Laurent expansion has a learning step*)


FFGraphOutput[gexpansion,laurent];
explearn = FFLaurentLearn[gexpansion];


(* ::Text:: *)
(*We finally reconstruct analytically the coefficients of the expansion*)


rec = FFReconstructFunction[gexpansion,invariants];


(* ::Text:: *)
(*From this we build the final amplitude, up to the finite part.  Because the tree-level is vanishing, the result is finite.*)


amplitude=(Normal/@FFLaurentSol[Factor[rec],eps,explearn]).functionbasis


(* ::Text:: *)
(*In a more complex example, before reconstructing the final expression, we recommend finding linear relations between the output elements of the final graph.  This allows to reconstruct only a subset of linearly independent coefficients of the Laurent expansion of the output, which can simplify both the reconstruction and the result itself.  An explicit example of this can be found in the "examples/linear_relations" directory of this repository.*)
