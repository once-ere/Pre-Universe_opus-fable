#!/usr/bin/env wolframscript
(* Evaluate one Wolfram source file under message capture and gate on it.

   Usage:  wolframscript -file tests/verify_wolfram_source.wls wolfram/<file>.wls

   The three sources in wolfram/ each build `report = TestReport[tests]` and
   print their own summaries. Those summaries say nothing about Wolfram
   MESSAGES: a source could pass every assertion while a Simplify::time or a
   Part::partw scrolled past. This wrapper loads the file with $MessageList
   captured, prints a uniform four-line verdict, and exits nonzero on any
   message, any failed assertion, or a missing report. *)

Begin["WolframSourceVerifier`"];

args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1,
  Print["usage: wolframscript -file tests/verify_wolfram_source.wls <source.wls>"];
  Exit[2];
];

repositoryRoot = ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]];
sourcePath = If[
  FileExistsQ[First[args]],
  ExpandFileName[First[args]],
  FileNameJoin[{repositoryRoot, First[args]}]
];
label = FileNameTake[sourcePath];

End[];

If[!FileExistsQ[WolframSourceVerifier`sourcePath],
  Print["MISSING FILE: ", WolframSourceVerifier`sourcePath];
  Exit[2];
];

{WolframSourceVerifier`seconds, WolframSourceVerifier`messages} = Block[
  {$MessageList = {}},
  {
    First[AbsoluteTiming[Get[WolframSourceVerifier`sourcePath]]],
    DeleteDuplicates[ToString[#, InputForm] & /@ $MessageList]
  }
];

WolframSourceVerifier`failed = Quiet[
  Check[Global`report["TestsFailedCount"], $Failed]];
WolframSourceVerifier`succeeded = Quiet[
  Check[Global`report["TestsSucceededCount"], $Failed]];

Print[];
Print["Source verifier file     : ", WolframSourceVerifier`label];
Print["Source verifier seconds  : ", N[WolframSourceVerifier`seconds]];
Print["Source verifier messages : ", Length[WolframSourceVerifier`messages]];
Print["Source verifier succeeded: ", WolframSourceVerifier`succeeded];
Print["Source verifier failed   : ", WolframSourceVerifier`failed];
If[WolframSourceVerifier`messages =!= {},
  Print[Column[WolframSourceVerifier`messages]]];

If[
  WolframSourceVerifier`messages =!= {} ||
  WolframSourceVerifier`failed === $Failed ||
  WolframSourceVerifier`failed =!= 0,
  Print["FAILED: ", WolframSourceVerifier`label];
  Exit[1]
];
Print["SUCCESS: ", WolframSourceVerifier`label,
  " passed every assertion with no messages."];
Exit[0];
