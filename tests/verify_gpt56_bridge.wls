#!/usr/bin/env wolframscript

Begin["GPT56Verifier`"];
repositoryRoot = ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]];
sourcePath = FileNameJoin[{repositoryRoot, "wolfram", "gpt56_bridge.wls"}];
End[];

If[!FileExistsQ[GPT56Verifier`sourcePath],
  Print["MISSING FILE: ", GPT56Verifier`sourcePath];
  Exit[2];
];

{GPT56Verifier`seconds, GPT56Verifier`messages} = Block[
  {$MessageList = {}},
  {
    First[AbsoluteTiming[Get[GPT56Verifier`sourcePath]]],
    DeleteDuplicates[ToString[#, InputForm] & /@ $MessageList]
  }
];

GPT56Verifier`failed = Quiet[Check[Global`report["TestsFailedCount"], $Failed]];
GPT56Verifier`succeeded = Quiet[Check[Global`report["TestsSucceededCount"], $Failed]];
Print["Focused verifier seconds: ", N[GPT56Verifier`seconds]];
Print["Focused verifier messages: ", Length[GPT56Verifier`messages]];
Print["Focused verifier succeeded: ", GPT56Verifier`succeeded];
Print["Focused verifier failed: ", GPT56Verifier`failed];
If[GPT56Verifier`messages =!= {}, Print[Column[GPT56Verifier`messages]]];

If[
  GPT56Verifier`messages =!= {} ||
  GPT56Verifier`failed === $Failed ||
  GPT56Verifier`failed =!= 0,
  Exit[1]
];
Exit[0];
