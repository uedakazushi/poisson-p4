# L-frame structured transcript / cubic11 generic rank10 の独立監査

2026-11-04。対象は `S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63`。
入力は先に独立受理した L-frame quad10×26（317 core terms）と genuine16×31
（50085 terms）。original source、q3 行、列順、C row scaling はその REVIEW で全係数を照合済み。
今回は新しい消去法を実行せず、保存された全 polynomial row identity を別実装で再生する。

## literal transcript の照合

行列の元行番号は、pair の二行、他の八 quad 行、genuine G0,…,G15 の順に固定する。
従って G6 は全26行の row16。保存済み十 quad pivot と五 heavy pivot は既知のものを読み、
新 pivot、minor、core 点、標数、GB、gcd を探さない。
レビューの polynomial context は author の DEGREVLEX と異なる LEX を使う。

各 `ELIMINATE t p c` の保存 u,v を読み、単に指定された column を零に置くのでなく、
最初に u*A[t,c]−v*A[p,c]=0 を検証し、**全31列**を u*A[t]−v*A[p] として計算する。
各 `NORMALIZE` の保存 polynomial g と有理数 c について、g≠0、c≠0 を確認し、
全31係数の除算が exact であることと、商*(g*c)=除算前の係数を確認する。
u≠0 と pivot coefficient≠0 も毎回確認する。これらは QQ(a,e,d) 上の可逆 row operation である。

十 quad 処理後の16×21、五 heavy 処理後の **全11×16**、十 quad echelon の全10×31、
最終 checkpoint の全26×31を、それぞれ保存済み係数と比較する。
source-orderや一部 entry の一致のみから complete-R11 の正しさを推論しない。
全項目が PASS であることは同名 report/native.log に記録する。

## 零行と generic rank10

11×16 の relative row4 は G6 であり、十 quad と五 heavy を処理し終えた時点で、
残りの16係数だけでなく元の31係数すべてが既に零である。
その後の最初の cubic pivot に到達しただけで job が零を検出して終了したため、
歴史的な名前 `six-remainders.txt` は実際には **15個**の零係数を保存している。
完成した十 cubic pivot や六 determinant 展開を実行したという主張ではない。

元 G6 の係数は initial normalization と target-row operations における非零の係数の積/商。
G6=row16 は全 transcript で pivot として使われず、他行の元 G6 係数は零のままなので、
この係数は QQ(a,e,d) 上で非零である。従って零行は自明な source relation ではない。
実 canonical11 map はこの generic chart 上で rank≤10。
先に独立受理した実 saved-QQ minor から generic rank≥10 なので、
**generic canonical11 rank は正確に10**である。

既存 lower10 は元の保存有理点の nonzero34minor に基づき、L-frame change と九 monomial
unit quotient を経た同一の intrinsic map に適用される。新しい point evaluation はない。

## chart と未証明範囲

入力 L-frame の row normalization だけは既知 aC unit による。
今回の transcript の非定数 pivot、normalization factor、elimination multiplier は一般に
追加 chart factor である。非零 polynomial と検証しただけで FULLU の全域 unit にはしない。
従って全点 lower10、全域 source splitting、uniform kernel38 は主張しない。

この REVIEW は別途抽出中の compact relation の小ささや、新方向が既知37から独立であることを
証明したものではない。colon による lift と既知 Bianchi submodule の比較には別の監査が必要。
全点 cubic57、全 quartic239、C122 の全域全射も未証明のまま保持する。

再実行（空 output directory）：

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_ell_Lframe_borders_review_20261104_v63.py --output work/marked-square-v63/Lframe-borders-review-REPLAY
```
