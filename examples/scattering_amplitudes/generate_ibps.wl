(* ::Package:: *)

<<LiteIBP`


<<"topology.wl"


(* GetSeeds generates the seed integrals for each sector, and for each type of identity *)
rank=4;
maxdots=0;
GetSeeds["IBP"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots},0,rank];
GetSeeds["LI"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots},0,rank];
GetSeeds["SR"][sector_]:=LIBPGenerateSeeds[sector,{0,0},0,rank];
GetSeeds["Map"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots+1},0,rank];
GetSeeds["ExtMap"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots+1},0,rank];


(* You can decide here how many parallel kernels you want to use to generate IBPs,
   but it does not work anymore in Mathematica 11.3 and later versions.  *)
If[$VersionNumber<11.3,
  nkernels = 4;
  LaunchKernels[4];
  ParallelEvaluate[<<"init.m"];
  ParallelNeeds["LiteIBP`"];
  DistributeDefinitions[GetSeeds];
]


(* for each family generate identities and put them in the "ibps/" directory *)
LIBPFastGenIds[#,GetSeeds,"Directory"->"ibps"]&/@{box1,box2,box3};


(* serialize in JSON format *)
LIBPSerializeFastIds[
  FileNames["ibps/ids_box*.mx"],
  s->1,
  {eps,t,mf2}]



