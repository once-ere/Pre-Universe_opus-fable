#!/usr/bin/env wolframscript
(* Execute notebooks/gpt-5.6_bridge.nb IN PLACE and store the outputs.

   scripts/run_gpt56_notebook.wls evaluates the notebook's Input cells but
   discards the results, so the file on disk stayed unevaluated: 5 Input cells
   and 0 Output cells. This script evaluates each Input cell, captures both the
   printed text and the returned value, and writes the notebook back with an
   Output cell after every Input cell.

   Run:  wolframscript -file scripts/execute_gpt56_notebook.wls 2>&1 | tee logs/execute_gpt56_notebook.log
*)

Begin["GPT56Executor`"];

repositoryRoot = ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]];
notebookPath = FileNameJoin[{repositoryRoot, "notebooks", "gpt-5.6_bridge.nb"}];

asOne[held_HoldComplete] := Replace[
  held,
  HoldComplete[expressions___] :> HoldComplete[CompoundExpression[expressions]]
];

(* Evaluate one cell body, returning {result, printedText, seconds, messages}. *)
evaluateSource[source_String] := Module[
  {stream, file, result, printed, timing, messages},
  stream = OpenWrite[];
  file = First[stream];
  {timing, messages} = Block[
    {$MessageList = {}},
    {
      First[AbsoluteTiming[
        Block[{$Output = {stream}},
          result = ReleaseHold[
            asOne[ToExpression[source, InputForm, HoldComplete]]
          ]
        ]
      ]],
      DeleteDuplicates[ToString[#, InputForm] & /@ $MessageList]
    }
  ];
  Close[stream];
  printed = If[FileExistsQ[file], Import[file, "Text"], ""];
  If[FileExistsQ[file], DeleteFile[file]];
  {result, printed, timing, messages}
];

End[];

If[!FileExistsQ[GPT56Executor`notebookPath],
  Print["MISSING FILE: ", GPT56Executor`notebookPath];
  Print["Run: python3 scripts/build_gpt56_notebook.py"];
  Exit[2];
];

GPT56Executor`notebook = Get[GPT56Executor`notebookPath];
GPT56Executor`cells = First[GPT56Executor`notebook];
GPT56Executor`options = Rest[List @@ GPT56Executor`notebook];

GPT56Executor`newCells = {};
GPT56Executor`inputCount = 0;
GPT56Executor`outputCount = 0;
GPT56Executor`messageCount = 0;
GPT56Executor`totalSeconds = 0;

(* The iterator must NOT live in the Global context: cell one evaluates
   ClearAll["Global`*"], which would try to clear a Block-scoped Global`i and
   raise ClearAll::clloc. *)
Do[
  Module[{cell = GPT56Executor`cells[[GPT56Executor`index]], source, evaluated,
          result, printed, seconds, messages},
    AppendTo[GPT56Executor`newCells, cell];
    If[MatchQ[cell, Cell[BoxData[_String], "Input", ___]],
      GPT56Executor`inputCount++;
      source = First[First[cell]];
      evaluated = GPT56Executor`evaluateSource[source];
      {result, printed, seconds, messages} = evaluated;
      GPT56Executor`totalSeconds += seconds;
      GPT56Executor`messageCount += Length[messages];
      If[StringLength[StringTrim[printed]] > 0,
        AppendTo[GPT56Executor`newCells,
          Cell[StringTrim[printed], "Print", CellLabel -> "During evaluation of In[" <> ToString[GPT56Executor`inputCount] <> "]:="]];
        GPT56Executor`outputCount++;
      ];
      If[result =!= Null,
        AppendTo[GPT56Executor`newCells,
          Cell[BoxData[ToString[Short[result, 8], InputForm]], "Output",
            CellLabel -> "Out[" <> ToString[GPT56Executor`inputCount] <> "]="]];
        GPT56Executor`outputCount++;
      ];
      Print["CELL ", GPT56Executor`inputCount,
        " seconds=", N[seconds],
        " printed=", StringLength[StringTrim[printed]],
        " messages=", Length[messages]];
    ]
  ],
  {GPT56Executor`index, Length[GPT56Executor`cells]}
];

GPT56Executor`header = StringRiffle[
  {
    "(* Content-type: application/vnd.wolfram.mathematica *)",
    "",
    "(*** Wolfram Notebook File ***)",
    "(* http://www.wolfram.com/nb *)",
    "",
    "(* Built by scripts/build_gpt56_notebook.py, then evaluated in place by",
    "   scripts/execute_gpt56_notebook.wls. Output cells below are captured",
    "   from a real evaluation; they are not hand-written. *)",
    ""
  },
  "\n"
];

GPT56Executor`stream = OpenWrite[GPT56Executor`notebookPath];
WriteString[GPT56Executor`stream, GPT56Executor`header];
WriteString[GPT56Executor`stream, "\nNotebook[{\n"];
WriteString[
  GPT56Executor`stream,
  StringRiffle[ToString[#, InputForm] & /@ GPT56Executor`newCells, ",\n"]
];
WriteString[GPT56Executor`stream, "},\n"];
WriteString[
  GPT56Executor`stream,
  StringRiffle[ToString[#, InputForm] & /@ GPT56Executor`options, ",\n"]
];
WriteString[GPT56Executor`stream, "\n]\n"];
Close[GPT56Executor`stream];

GPT56Executor`failed = Quiet[
  Check[Global`report["TestsFailedCount"], $Failed]
];
GPT56Executor`succeeded = Quiet[
  Check[Global`report["TestsSucceededCount"], $Failed]
];

Print[];
Print["Notebook path            : notebooks/gpt-5.6_bridge.nb"];
Print["Input cells evaluated    : ", GPT56Executor`inputCount];
Print["Output cells stored      : ", GPT56Executor`outputCount];
Print["Total evaluation seconds : ", GPT56Executor`totalSeconds];
Print["Cells with messages      : ", GPT56Executor`messageCount];
Print["Notebook tests succeeded : ", GPT56Executor`succeeded];
Print["Notebook tests failed    : ", GPT56Executor`failed];

If[
  GPT56Executor`messageCount =!= 0 ||
  GPT56Executor`outputCount === 0 ||
  GPT56Executor`failed === $Failed ||
  GPT56Executor`failed =!= 0,
  Print["FAILED: the notebook did not evaluate cleanly."];
  Exit[1]
];
Print["SUCCESS: the notebook is evaluated and its outputs are stored."];
Exit[0];
