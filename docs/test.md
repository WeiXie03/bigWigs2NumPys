### Specifying coordinates
-> Repeat `sequential with missings` toy data, but just following subset of intervals (0-based, half-open):

>[!note]
>'`-`' here denotes a missing/unmappable value in the data

chr | start | end | expected result
--- | --- | --- | ---
1 | 0 | 4 | 0 - 0 1
1 | 5 | 7 | 3 -
1 | 9 | 10 | 1
2 | 1 | 2 | 1
2 | 3 | 4 | 3
3 | 3 | 6 | 0 1 2

_Jan 18_
_Why `bbGetOverlappingEntries()` returning no overlapping entries when test program in libBigWig test directory prints correct ones __on same bigBed file?!___