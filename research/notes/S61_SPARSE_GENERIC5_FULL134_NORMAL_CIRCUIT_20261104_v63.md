# 五変数源上の全134式のfactored normal42回路

元のq45に関する全134式を保持し、受理済みのアフィン変換q=center+E*z+std_nonpivot*n42を有限多項式回路として代入する。元の係数は五変数QQ多項式のscalar atom、Eは既存の未展開行列式・adjugate回路から読み込む。Gammaや全parameter多項式を展開しない。

直接代入した134出力、q45のforward回路、z3/n42のinverse回路を保存する。E_pivot=I3を使い、形式的なE係数上で45×45線形変換の両方向の積がIになることを4050成分で確認する。さらにq次数2以下の1081単項式をすべて扱い、係数block組立ての普遍的恒等式を確認する。

係数queryはn²、z*n、linear-n、plane上のz次数2以下だけに限る。n²は元の非pivot二次係数をそのまま使い、厳密なzero/nonzeroが分かる。linear-nとplane定数は小さいcenter多項式だけで正確に求め、元への逆center移動で全係数を確認する。Eを含む係数は未展開scalar回路として保持し、既知の零・既知の非零・field加算の相殺が未判定の部分を区別する。構文上の非零を数学的な非零と呼ばない。

すべての逆元はparameterだけの式で、未知q,z,nを反転しない。元134式は非零で、アフィン変換が可逆なので、代入後も全134式が同じアフィン問題を表す。center自体がPoisson点だとは仮定しない。

## 既知五点との適合

元の平面の三方向V=[N*T0,N*T1,Q52-center]とW=P Vを用いると、sigma=P*row2(adj B_W)/GammaW*zとなる。普遍的な逆行列恒等式からtoric三点ではsigma=0、Q51では-1/2、Q52では1である。既知五点quotient上の三冪等元もsigmaの回路として保存する。ただし、未知の平面外解を含む全ファイバー上の冪等元とは主張しない。

以前の単一β標本のquintic係数はコピーしない。既知3+1+1の積代数は座標変換で保存されるが、新しいprimitive separatorやquinticの明示化には別の有限計算が必要である。

## 次のcost

保存blockからplane制限134×10の実際の非零5minorを選んで引けば129残差が得られ、z*nの129×126係数が構成できる。cross126rankの生成的な事実と、以前の固定126行minorがこの源で非零であることは分ける。今回、新しい126rank、Euler/J0行列、Macaulay閉包、Groebner基底、一般残余の空性は計算しない。

再現：

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_sparse_generic5_full134_normal_circuit_20261104_v63.py --output work/NEW-output
```

一回60秒・2GB・32MiBで、直接回路を先に保存してから有界な係数queryを行う。

## 今回の確定範囲と数値

2026-11-04 の一回の実行は **5.8409秒で完了**した。直接回路27163 nodeと追加係数回路201487 node、計228650 node、4296 scalar atomを保存した。成果物の圧縮合計は約1.5MiB、展開したJSONの合計は7246948 byteであり、32MiBを下回る。45未知数は `z0,z1,z2,n0,...,n41` で、元の全134非零式を保持する。

| 係数block | 総成分 | 厳密な零 | 厳密な非零 | field和の相殺が未判定 |
| --- | ---: | ---: | ---: | ---: |
| n² (i≤j) | 121002 | 108480 | 12522 | 0 |
| linear-n | 5628 | 1065 | 4563 | 0 |
| z*n | 16884 | 5571 | 145 | 11168 |
| plane (1,z,z²) | 1340 | 0 | 134 | 1206 |

最後の134個の非零は各式のplane定数項であり、plane上の他の係数がすべて非零という意味ではない。非零判定は非零QQ多項式atomと整域の積だけを使う。未判定のfield和は回路のまま保存する。

元の式を `F=c+Σ l_i q_i+Σ_{i≤j}Q_ij q_i q_j`、`q_i=C_i+y_i` と書くと、保存する平行移動係数は

```
c'=c+Σ l_i C_i+Σ Q_ij C_i C_j,
l'_i=l_i+Σ_{j≠i}Q_{min(i,j),max(i,j)} C_j+2Q_ii C_i.
```

この小さいscalar計算だけをQQ[H,L,s0,q0,t0]で行い、逆平行移動による全6164恒等式を検証した。続いて `y_i=Σ_a E_ia z_a+n_i` (pivot行ではn_i=0) を使う。全1081基底単項式に対する係数組立ての形式的恒等式と、両方向の45×45逆変換4050成分が成立する。この二段階により、保存した係数blockは直接代入回路と同一の134多項式を表す。E自体の大きいparameter分数の展開は行わない。

係数体は、前段で証明済みのrank3 chart上の `K=QQ(H,L,s0,q0,t0)` である。Γ、P、前段pivotの分母等はparameterだけであり、Kの非零元として反転する。この表現はK上の全アフィン解を保存する。chartの分母が零になる特殊parameter点まで同じ表示が有効とは主張しない。

## 保存形式

`full134-direct-circuit.json.gz` は全134出力、forward45とinverse3+42、E45×3、center、sigmaを含む。scalar atomは `scalar-atoms-direct.json.gz` のprimitive整数多項式とQQ contentで表す。外部葉 `original_W_poly(r,c)` は前段 `plane-W-original.json.gz` のW45×3、`rank_pivot_rf(i)` と `selected_row_factor_poly(i)` は `rank-certificate.json.gz` の同名配列を参照する。葉の入力ファイルと全入力SHA256はreportに固定してある。

`coefficient-block-circuit.json.gz` は直接回路をprefixとして参照する追加nodeと、全134行のplane10、linear-n42、z*n126、非零n²の疎配列を含む。atomの追加後の全表は `scalar-atoms-full.json.gz` に保存する。演算はQQ定数、parameter多項式atom、外部parameter葉、和、積、parameter逆元、未知数葉だけであり、有限個の厳密な体演算として再生できる。

`structural-proof.json`、`next-cost.json`、`report.json` は、それぞれ普遍的逆変換検証、次段の規模・未確認事項、全入力と証拠ファイルのhashを記録する。元134式のcanonical疎表に重複係数がなく、q次数が2以下で全134行が非零であることも保存入力の直接確認で確かめた。

## 未完了

新しいplane制限のrank5と実際の5minor、129×126 cross行列のrankおよび非零126minorは未計算である。上表の疎性だけからrankを推定しない。以前の一般的なcross-rank証明の座標不変性を利用する場合も、制限した源で使う式空間とknown-plane quotientを対応させる必要がある。既知五点のcubic×linear×linearの積構造は保存されるが、primitive quinticの新係数や全未知ファイバーの空性はここでは証明していない。
