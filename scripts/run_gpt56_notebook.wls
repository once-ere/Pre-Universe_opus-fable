#!/usr/bin/env wolframscript

Begin["GPT56NotebookRunner`"];
repositoryRoot = ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]];
notebookPath = FileNameJoin[{repositoryRoot, "notebooks", "gpt-5.6_bridge.nb"}];
asOne[held_HoldComplete] := Replace[
  held,
  HoldComplete[expressions___] :> HoldComplete[CompoundExpression[expressions]]
];
End[];

If[!FileExistsQ[GPT56NotebookRunner`notebookPath],
  Print["MISSING FILE: ", GPT56NotebookRunner`notebookPath];
  Print["Run: python3 scripts/build_gpt56_notebook.py"];
  Exit[2];
];

GPT56NotebookRunner`notebook = Import[
  GPT56NotebookRunner`notebookPath, "Notebook"
];
GPT56NotebookRunner`inputs = Cases[
  GPT56NotebookRunner`notebook,
  Cell[BoxData[source_String], "Input", ___] :> source,
  Infinity
];
Print["Evaluating ", Length[GPT56NotebookRunner`inputs],
  " Input cells from notebooks/gpt-5.6_bridge.nb"];

GPT56NotebookRunner`messageRecords = {};
GPT56NotebookRunner`timings = {};
Do[
  {GPT56NotebookRunner`timedResult, GPT56NotebookRunner`cellMessages} = Block[
    {$MessageList = {}},
    {
      AbsoluteTiming[
        ReleaseHold[
          GPT56NotebookRunner`asOne[
            ToExpression[
              GPT56NotebookRunner`inputs[[GPT56NotebookRunner`cell]],
              InputForm,
              HoldComplete
            ]
          ]
        ]
      ],
      DeleteDuplicates[ToString[#, InputForm] & /@ $MessageList]
    }
  ];
  AppendTo[
    GPT56NotebookRunner`timings,
    GPT56NotebookRunner`timedResult[[1]]
  ];
  If[GPT56NotebookRunner`cellMessages =!= {},
    AppendTo[
      GPT56NotebookRunner`messageRecords,
      {GPT56NotebookRunner`cell, GPT56NotebookRunner`cellMessages}
    ]
  ];
  Print[
    "CELL ", GPT56NotebookRunner`cell,
    " seconds=", N[GPT56NotebookRunner`timedResult[[1]]],
    " messages=", Length[GPT56NotebookRunner`cellMessages]
  ],
  {GPT56NotebookRunner`cell, Length[GPT56NotebookRunner`inputs]}
];

GPT56NotebookRunner`failed = Quiet[
  Check[Global`report["TestsFailedCount"], $Failed]
];
GPT56NotebookRunner`succeeded = Quiet[
  Check[Global`report["TestsSucceededCount"], $Failed]
];
Print["Notebook input cells: ", Length[GPT56NotebookRunner`inputs]];
Print["Notebook total seconds: ", Total[GPT56NotebookRunner`timings]];
Print["Notebook cells with messages: ",
  Length[GPT56NotebookRunner`messageRecords]];
Print["Notebook tests succeeded: ", GPT56NotebookRunner`succeeded];
Print["Notebook tests failed: ", GPT56NotebookRunner`failed];
If[GPT56NotebookRunner`messageRecords =!= {},
  Print[Column[GPT56NotebookRunner`messageRecords]]
];

If[
  GPT56NotebookRunner`messageRecords =!= {} ||
  GPT56NotebookRunner`failed === $Failed ||
  GPT56NotebookRunner`failed =!= 0,
  Exit[1]
];
Exit[0];
