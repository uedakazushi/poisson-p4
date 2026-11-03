# Cubic 9行判定と quartic 319行縮約の独立監査

2026-11-03。著者 packet `S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63` と `S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63` を、既存の三次 row-basis 独立監査を依存として検証する。著者 producer は実行しない。

## Cubic: 追加 chart 上の厳密な判定

圧縮された元の57×91行列を streaming で再読し、全5187成分、164333係数項、最大 core 次数16と非圧縮 hash を確認する。5190という旧 prose の誤記は独立した COUNT_CLARIFICATION note/report が訂正しており、header は行列成分ではない。

七 heavy 列と84 pure列、最初の41 product 行の全 heavy 成分が恒等的に零であること、残る16 genuine 行との分割を元 label・全係数から確認する。固定 product41列と genuine7行の選択により48角 pivot D が block 三角形となり、`det D=DeltaP*DeltaH` である。残る9行・43 pure列の Schur 行列を T とすると、各49角 bordered minor は `Sij=det(D)*Tij`。従って **追加 chart `D(DeltaP*DeltaH)` 上**で `rank A=48+rank S` であり、full57 と rank S=9 は同値である。

QQ 点(1,2,3)で元全行列を評価し、独自 Fraction Gauss 消去により D の逆行列・両 block の determinant・全387 Schur 成分・9小行列式・border scale を確認する。一つの49角 bordered determinant も直接計算して符号を照合する。387個の parametric determinant 自体は展開しない。この点は chart が非空である証拠であり、全点 rank57、defect 有限性を示さない。

全 product/heavy minor の有限 chart 被覆という構成は正しい。ただし各 chart の DeltaP/DeltaH は追加の単元であり、全 aFBC-open の既知単元ではない。product rank41 と heavy rank7 と抽象商 free16 だけから rank57 は従わない。著者の module countermodel は、固定 pure41、heavy7、別 pure8 と最後の行 `x*e49+y*e50` により全体 rank56への曲線上の低下を実現する。この論理的反例を actual Jacobi 解と解釈しない。

## Quartic: 全 aFBC-open で row span が同一

元の full weight4 表示は6二次式×29重み二単項式、21 genuine cubic×7 light、33 genuine quartic の354×239である。元 chart の equation weights と selected34 から33本を独立列挙し、全て保持する。

既存の独立 row-basis 監査は63行の既知単元正規化と全6三次恒等式を全91列で確認している。その最初5本を light7通りで掛け、35本の有理係数恒等式を独立構成する。product `Ri=qi*vj/ui` は quadratic×weight2 の行へ写し、genuine row はその light 倍の行へ写す。保存された全728非零係数を numerator/denominator の多項式交差積で比較する。全63 normalizer と各保存分母が a,F,B,C と非零QQ定数のみの積であることを確認する。

削除する cubic 元行は48,55,59,57,61で、各 light に対する pivot は1,1,1,a,−a。後の恒等式に既削除行が現れないことを全35行で確認し、35角三角 determinant は `−a^14` になる。従って35行は残る319行で解け、逆包含は部分集合であることから自明である。全239係数の恒等性は、監査済みの全91列三次恒等式の七つの単項式写像（各91→239が単射）によって得られる。元の重い四次展開を新しく実行する必要はない。

この row span 等号は **全 aFBC-open** で有効であり、d、J、product pointwise lower41、heavy lower7、full cubic lower57 を仮定しない。319×239 の surjectivity、長65、quartic pivot の全域可逆性、ambient removal は未証明である。

再現（新しい出力先を指定）:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_structural_reductions_review_20261103_v63.py \
  --output work/marked-square-v63/structural-reductions-independent-new
```

公開 report の `all_passed` が実行結果である。GB、新しい体、巨大な parametric 小行列式、四次全射の探索は行わない。
