# P*のMPS親監査：全八 rational-jet 証明書と全四十IDの残候補

2026-11-08、PRIMARY E_multiple_hyperplane。
producerのONE symbolic finite batch全八例 PASS、whole-source proper-incidence論理の独立査読待ち。
新component数は零。全parentを排除したとは主張しない。

## 1. 算術の正確な範囲

job `jobs/unimodular_pole_eight_mps_fixed_quadratic_jets.json` はfresh gate PASS。
60秒／1GiBの有限契約下、SymPyのexact有理関数体上で元source全parameterを保持した。
実測0.684708283秒、最大RSS50,336KiB、prescribed八例全てがcertifiedとなった。
新Poisson source点、周囲Jacobi operator/M1、rank/minor、normal deformation jet又はprofileは無い。
全true sourceとgeometric coordinate零点が事前の有限リストに固定され、別点を選び直していない。

各例について全六affine括弧の0/1jet恒等零、二次leading Q、cubic alpha、PfQ、factorization、
全log-residue係数恒等式、固定rational residue比、integer-exponent F/G、
全六 `alpha wedge (G*dF-F*dG)` 係数零と非定数性を保存した。
executed source全文とhashはreportに埋込み、producer/helper import無し。

## 2. 全parameterの実first integral

v0,…,v3はsource元順からmarked coordinateを除いた四変数。
下表のF/Gは全source parameterをsymbolicのまま保持する。

| source | marked order二zero | (v0,v1,v2,v3) | F/G | degree d |
|---|---|---|---|---:|
| T24 | x2 | (x0,x1,x3,x4) | v3³/(2p0v0v1v2+6t0v2v3²+3t2v0²v3) | 3 |
| T26 | x0 | (x1,x2,x3,x4) | v1v3/(p0v0v2+4t2v3²) | 2 |
| T27 | x0 | 同上 | v1v3³/(p0v0v2+7t2v3²)² | 4 |
| T31 | x0 | 同上 | v1²/(p0v0v2+2t2v1v3) | 2 |
| T32 | x0 | 同上 | v1/v0 | 1 |
| T33 | x3 | (x0,x1,x2,x4) | v3³/(p0v0v1v2+3t0v2v3²-t1v1³) | 3 |
| T34 | x0 | (x1,x2,x3,x4) | v1²(4p0v2v3+t2v1²)/v0⁴ | 4 |
| T35 | x0 | 同上 | v1⁴(3t2v0²-4p0v2v3)/v0⁶ | 6 |

各F/Gは同じdegreeの独立polynomial二つから成り、元source open上で非定数。
constant coefficient parameterを数字へ特殊化した結果ではない。

## 3. 全component閉包への最弱矢印

純method `UNIMODULAR_POLE_MPS_BOUNDED_RATIONAL_JET_METHOD.md` により、
固定degree dのGr(2,Sym^d V*)とalphaのwedge-incidenceはclosedでprojectionがproper。
従って「あるorder二zero qでそのquadratic-jet foliationにdegree dのrational first integralがある」
という条件はwhole projective tensor space上のclosed条件になる。
sourceのgeneric rational F/Gが特殊parameterでdependent又は共通因子付きになっても、
Grassmannian上のproperな二次元spanのlimitを用いるので全零／critical補集合を失わない。
alpha零もclosed coneへ含める。

P*はorder二zeroが一つだけで、その非零alphaは非定数rational first integralを持たない。
よって表の八つの**whole closed MPS componentは全てP*を含まない**。
この結論はlocal dimension又は一つのgeneric numerical pointからの推測ではない。

## 4. 全四十MPS IDの監査合成

現在のpure/whole source証明の役割を分ける。

- T01,T39,T40：projective dimension25で、P*の同一点全周囲tangent上界24に反する。
- T02,T03,T04,T06,T07,T08,T09,T13,T15,T16,T17,T20,T21,T22,T23,T25,T28,T29,T30：
  full H-Casimir sourceがC+Q、marked rank[dQ|L]≤4のproper閉条件で排除。
- T04,T05,T06,T07,T08,T09,T10,T11,T12,T13,T14,T15：
  coordinate order二zeroのtoric四linear-log jet。triple-pole最高quadratic exact numeratorのrank≤1条件に反して排除。
  上行との重複を新成分／別排除数に二重計数しない。
- T24,T26,T27,T31,T32,T33,T34,T35：本全parameter bounded-RFI certificateで排除。
- T36：既whole order二zero/quadric-pencil条件とP*の無rational-first-integralで排除。
- T38：既whole rankzero-surface条件とP*のzero support一次元で排除。

残るnamed MPS親は厳密に **T18、T19、T37**。
T18はH-Casimir helperに真のodd terms tD,uBを持ち、pure-evenと取り違えない。
T19は別作者がactual曲線を研究中、secant-conic jetをtangentと仮定しない。
T37はgenericにorder二zeroが無いので、本point-marked jet条件を勝手に転用しない。

他のnonMPS成分の除外はrootの別prior auditに依存する。
P*の既知親、成分性、smoothness又は全h5 E完了はまだ結論しない。
今の局所D101/M47 gap二に新点／higher jet／minorを足さない。
