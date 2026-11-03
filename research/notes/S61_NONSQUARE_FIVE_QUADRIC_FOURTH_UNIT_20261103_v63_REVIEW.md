# 補助五二次式の第四単元pivot：独立監査

2026-11-03。`S61_NONSQUARE_FIVE_QUADRIC_FOURTH_UNIT_20261103_v63`を独立に確認した。
対象は補助5×21係数行列だけであり、主三次/product行列のconstructive pivot34へ加算しない。
source依存を含めて、全a*P*J≠0で補助rank≥4を受理する。

## 最初の三pivotまでのsource依存

V49 structural、third-pivot frontier、fourth-unitの三archiveと全input SHAを照合した。
凍結inputを変更せず、copyだけを新しいdirectoryに展開した。
804個の実四quadric係数から

    Q2prime=(Q2−eQ1)/(a²J),
    D=(40*x0*Q1−A*Q2prime)/W

を再計算した。DのW除算はliteral polynomial divisionであり、Wをinvertしていない。
全5×21係数が元のnative sourceへ一致する。

第一pivotは元行列のrow0,col2のa/4。
既知a,P,J因子とQQ定数だけで行をnormalizeし、残るQ2prime行について独立な小さいQQ係数方程式を解いた。
残存20列のlocal index(2,9,16,18)のcofactor

    (7/60,−1/10,7/360,−1/60)

がP*Jを与える。
local col2（元col3）の係数7/60は非零QQ定数なので、列置換は可逆である。
この第二pivot後の全3×19行は、保存V49出力と非零QQ定数倍で一致した。
最初の二pivotを未監査のlogだけから数えていない。

第三pivotは元3×19のrow1に対する

    (55000/3)*M[1,5]−(11000/3)*M[1,7]=a³*P*J²

をliteral QQ積で確認した。
col5の置換係数は非零QQ定数で、残存二行の全Schur係数は保存2×18に対し共通known unit
−1800*a²*P*J倍で一致する。
従ってsourceから2×18まで全て既知単元だけの安全な操作である。

## 第四の175項恒等式

2×18元行列の第0行とsource-row0ファイルが全係数で一致する。
18cofactorは全て非零、合計175項、最大次数4である。
行とcofactorのQQ分母をそれぞれ一度だけclearし、Pythonの独立な整数polynomial積で

    sum A[0,j]*c[j]=a³*P³*J²

を正確に確認した。
authorのFLINT row-dot/helperは呼び出していない。
modp候補選択の結果も証明の前提としていない。

この組合せを新列としてappendしてもimageは変わらない。
新列の第0成分g=a³P³J²は既知単元なので、残存行は

    R_j=g*A[1,j]−(sum A[1,k]*c[k])*A[0,j]

というfraction-free Schur行で表せる。
全18多項式を独立に計算し、a,P,Jだけの共通因子とQQ定数を除いた後、
保存された全8475項・次数19の1×18出力と一致した。
その正確な共通ratioは、非零QQ定数*a³*P*J³である。
未知polynomial factorを剥がしたりinvertしたりしていない。
元出力のSHAもarchive/reportと一致する。

## 全18個のprimitive column normalization

保存された各columnについて

    last_j = QQ_content_j * a^i_j * P^k_j * J^l_j * primitive_j

の全係数恒等式を、再度整数分母clearした積で検査した。
全powerは非負、QQ_contentは非零定数であり、許可された単元以外を除かない。
primitive出力は全7865項・次数18、全係数が整数であることも確認した。
これは列の単元倍なのでlocalized row ideal/imageを保つ。

## 受理する範囲

補助行列の最初の四pivotをsourceから追跡でき、全a*P*J≠0で補助rank≥4となる。
元主三次行列の34pivotという記録は変更しない。
補助rank5、pair gate、全QQ raw24/product41、quartic、ambient removalは未完了のままである。
残る行の有限次数探索が失敗したことから、ideal非単元性や全次数不可能性を推論しない。
新しいGB、CRT、cofactor探索は行っていない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/s61_nonsquare_five_quadric_fourth_unit_review_20261103_v63.py --output <new-output-directory>

全known-unit ratioと独立証明書:
`reports/S61_NONSQUARE_FIVE_QUADRIC_FOURTH_UNIT_20261103_v63_REVIEW/`。
実行約3秒、制限120秒・2.2GB。copy snapshotは検査後に削除し、frozen sourceは不変。
