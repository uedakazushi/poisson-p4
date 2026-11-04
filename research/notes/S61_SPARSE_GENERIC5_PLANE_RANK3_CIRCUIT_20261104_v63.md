# 平面rank3の証明と未展開normal42アフィン変換

保存された共同分母Pの平面をV=[N*T0,N*T1,Q52-center]とする。多項式行列W=[n0,n1,P*(Q52-center)]はW=P Vを満たす。Pは非零なパラメータ多項式であり、列全体への共通scalarなのでWとVは一般係数体で同じ平面を張る。

一回のジョブは最初に、各行で三列共通の最大lex単項式を選び、その係数45×3行列の階数を調べる。これは係数を平坦化した行列ではない。rank3なら選択行の行列式の最高単項式係数が非零となり、点を使わず多項式minorの非零が従う。

結論が出なければ各行をexactな多項式gcdでprimitive化して、三回の一般係数体pivotを行う。因子を除いた各成分は再乗算し、すべての非零有理pivotの分子・分母、row/column swapと消去係数を保存する。この行scalar変更はrank判定だけに使う。選択三行の元Wのminorは、primitive行列の非零pivot積に除いた三つの行因子と列置換符号を掛けたものである。

normal42のアダプタは**元のW**から構成する。primitive行列を45個のqの平面と取り違えない。選択三行blockをBとすると

```
Gamma=det(B), E=W*adj(B)/Gamma,
q=center+E*z+std_nonpivot*n42,
z=q_pivot-center_pivot
```

である。det、adj、E、全45出力は有限SLPで保存し、大きな三角式やadjugateは展開しない。普遍的な可換環の恒等式B*adj(B)=adj(B)*B=det(B)IとGamma非零からE_pivot=I3、従って45未知数の可逆なアフィン変換が従う。

単元は一般係数体QQ(H,L,s0,q0,t0)のパラメータ式だけで、未知qを反転しない。全特殊パラメータ上の正則性、全域有限性、一般ファイバーの網羅は主張しない。134式やJ0への代入、Groebner基底、点や素数の試行も行わない。

再現：

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_sparse_generic5_plane_rank3_circuit_20261104_v63.py --output work/NEW-output --work work/NEW-work
```

総60秒・2GB・32MiB、中間150k項を保ち、先頭係数証明書・primitive行・field操作の完成分を順次保存する。

## 完了結果と正確な係数chart

一回のジョブは7.032秒で完了した。先頭単項式の係数行列はrank1で、これだけでは結論できなかった。規定のfallbackを実行し、QQ多項式gcdによる行primitive化と三つの非零一般係数体pivotでrank3を証明した。選択した元のq行は3,6,13（0-based）、pivot列順は2,1,0で符号は-1である。

nativeは1.740秒、258回のfield操作で、最大分子12,611項・最大分母343項だった。三つの非零pivotはそれぞれ分子/分母支持290/1、343/1、5772/343である。選択した三つの行因子はすべて単項式であり、元のWの3minorはそれらの積とpivot積に符号-1を掛けたものとなる。全除算と再乗算、row/column swap、正規化行、消去係数を保存した。

normal42は元のWから作った約1600nodeの有限SLPである。Gammaの六つの行列式項も保存するが、その逆数は上の非零行因子とpivotの逆数の積として保持し、大きな多項式を掛け合わせない。E_selected=I3は普遍的adjugate恒等式と正確なGamma因子恒等式から従う。残る42座標はq行3,6,13以外を順に使う。

選択された歴史S61座標はu6,u9,u19（1-based）である。全45個のq-to-u対応とSLP leafの係数表現も保存し、源の五変数入力から既存の61変数モデルへ戻せるようにした。最終SLPは1608nodeである。

この結果により一般係数体で45未知数の可逆なアフィンnormal42座標が得られた。新たな単元条件は保存された非零parameter式と未展開Gammaだけで、未知qに依存しない。元134式のnormal42への代入・残余の空性はまだ計算していない。大きいliteral証明書は、nativeが未圧縮32MiB制限を守ったことを確認した後にgzip化して保存する。

centerは既知五点を含むアフィン平面の中心であり、center自体がPoisson解だとは仮定しない。この変換は全アフィン45変数の変換であって、中心点における局所完備化の指定でもない。
