# 新成分のgeneric二法線quadratic obstructionは正確にrank1

2026年11月3日。係数体は正しいsheared dense normalization n1=1のQQ(n0,n2)である。
点(1,1)だけの二次rank1をgenericへ延長するのではなく、全関数体の上界を証明する。

前稿で構成した二つの普遍degree1 polynomial tangent Tを使う。
全350×126 Jacobi行列JとJ*T=0を再構成し、原物理二次Jacobiから
三つのmonomial e1²,e1e2,e2²に対応する全350×3列Q(T)を直接計算する。
既存の26個の実際の群・scalar・族方向とTは、generic kernel28を張る。

## 実族を法として二次obstructionを計算してよい理由

任意の実族方向P'_sに沿って普遍恒等式J*T=0を微分すると

    B(P'_s,T)+J*T'_s=0

となるので、mixed Hessian項はimJに入る。
実族自身の二次項も、全次数Poisson恒等式を二度微分すればimJに入る。
群・scalar方向についても同じことが有限作用の同変性から従う。
従ってQ(T)をimJで割ったrankが、26 smooth族方向を除いたintrinsic二法線obstructionのrankである。
一つの点で任意に選んだkernel vectorだけを使う議論ではない。

## full augmented matrixの全関数体rank

350×129行列 `[J | Q(T)]` を作る。
QQ[n0,n2]で79個の非零有理定数pivotを検査・適用すると271×50になる。
全entryが恒等的に零である98行と8列だけを除き、173×42の残差を得る。
全41327項を保存したnative inputをFLINT有理関数体rank engineへ渡す。
実際の完了runは77.9066秒、正確な残差rank20で、全最終零行を検査した。
よって全generic augmented rankは79+20=99である。

前稿の普遍28 kernelと98 QQ-minorによりrankJ=98なので、

    rank(coker quadratic map)=99−98=1

が全係数体で厳密に成り立つ。
別に(1,1)での99×99非零minorと、rank1 cokerのbinary quadratic formの
非零discriminantも保存する。後者はgeneric二次形式が非退化となる開条件の下界証拠である。

degree1/2 polynomial syzygy ansatzで二つ目のHessian relationが見つからなかった事実は
generic rankの上界にも下界にも使わない。
また六つの有理標本が全てrank99だったこともgeneric上界には使わない。
上界は全二変数関数体の行消去による。

## Artin長へのscope

本稿はgeneric接方向二つとquadratic関係一つを証明するので、
generic normal ringのdegree≤2のHilbert値は少なくとも(1,2,2)、長さは少なくとも5である。
別packetの点(1,1)でのexact Artin長5と合わせてgeneric長5を結論するには、
その点のformal sourceが実際の26変数branch baseに対してfiniteな5生成moduleであることが必要である。
その有限性・free性・全完成環の同定は別のflatness packetで証明する。
本稿単体ではgeneric長や全local ringを分類したとは記載しない。

* `code/q_new_family_generic_quadratic_rank_20261103_v63.py`
* `reports/Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63/{report.json,certificate.json.gz,native-input.txt,native-log.txt}`

保存時は既に完了した実際の77.9066秒のnative runを再利用し、
再構成した全native inputとのbyte一致を検査する。
`--native` を指定すれば別の新規native replayも実行できる。
全source・証明書・実行入力・完了logのSHA256を保存する。
