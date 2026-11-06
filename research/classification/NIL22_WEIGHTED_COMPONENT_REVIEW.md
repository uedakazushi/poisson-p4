# nil22 A=−1：22次元全成分の独立検証

2026-11-07。主課題 E_multiple_hyperplane。
`MODULAR_NIL22_TOP_SUPPORT.md` §5–7 の主張を受理する。
全七変数源の rank四部分は一つの22次元射影軌道であり、その閉包は全 Poisson 空間の
既約成分となる。従前70成分とは相異なる。全分類の完成・全成分数の上界は主張しない。

## 1. 今回独立に再構成した数値証拠

`classification_nil22_weighted_slice_review.py` は標準 Python の Fraction と AST のみを使う。
producer、producer の Jacobi engine、SymPy、NumPy、旧数学 helper は import しない。
公開正規形の係数から多項式積・微分・Schouten 括弧・Euler 外積を組み立て、
全350本の四次係数と全150方向を再構成した。
元の全90係数から七変数源への reduction は、この独立数値チェックの対象外である。

| 検証 | 独立 QQ 結果 |
|---|---:|
| 全 Jacobi 微分 D の形 | 350×150 |
| D の階数 | 102 |
| 実 PGL・鉛直・scale の接空間 G | 47 |
| G と公開法線を合わせた階数 | 48 |
| 真の安定線形切片 S | 103 |
| S 内の陰関数消去座標・方程式 | 102 |
| 法線の相対重み | −1 |

D の階数は、102次の非零 QQ 行列式と、48個の独立 kernel 方向の全積零による
上下の比較で確定した。保存された IFT minor の正確な行列式とも一致した。
G の独立な47 pivot 行と補う103座標は全150座標を分割し、
G の pivot 行への射影は階数47。従って保存された S は真の線形補空間である。
法線は S に入り、102個の陰関数座標以外の自由座標の係数は1である。
全群方向の重み安定性、全 ambient 相対重みの範囲 [−8,6] も確認した。

公開された有限曲線を別算術で展開すると、全350式の τ次数0,…,5の係数は零となる。
次数6の全残差は行228,318の二箇所だけで、いずれも 455/497664。
全102選択方程式の次数6係数も零である。
公開された四項の functional

\[
 \ell=\tfrac{321}{80}e_{15}^*+\tfrac{351}{40}e_{39}^*
       +\tfrac{39}{40}e_{86}^*+e_{228}^*
\]

について、**全150列**に対する `ℓD=0` と

\[
 \ell(F_6)=455/497664\ne0
\]

を有理数で直接確認した。全保存補正が真の S の陰関数座標にだけ存在し、
次数kの係数が相対重み−kを持つことも確認した。
従って任意の部分行列や法線だけの方程式の検査ではない。

実測2.452秒。出力は
`reports/CLASSIFICATION_NIL22_WEIGHTED_SLICE_REVIEW_20261107.json`。
producer・正規形源・証明書の hash、独立 checker の hash、正確な minor・functional・残差を含む。
有限契約は `jobs/nil22_weighted_slice_review.json`。

## 2. 真の形式支持から成分性への論理

G は実際の GL5、Euler 鉛直加算、scale の軌道の接空間である。
47個の独立な作用パラメータを局所的に選び、S と作用写像を組み合わせると、
原点で微分同型の150次元形式座標系を得る。
Poisson 条件はこの作用を保つので、S 上の支持次元は全周囲の支持次元を支配する。

102次の可逆 minor により、S の102座標は一つの法線パラメータ τ の一意な形式級数となる。
相対トーラス作用の同変性から各級数は一つの単項式または零で、次数は8以下。
この有限化は、任意に計算次数を増やして真の支持を推定する操作とは異なる。
独立検証した次数0,…,5の消去は一意な IFT 写像のその次数の係数に一致する。
次数6で現れた非零余核 pairing は、次数6の未知補正や以後の補正では消せない。
従って一変数の真の形式切片の ideal は非零で、被約支持は原点だけである。
この被約支持結論には、正確な scheme の厚さの決定は不要である。

全周囲の局所次元上界47と実軌道の下界47が一致する。
鉛直24次元と射影 scale 一次元を除けば局所次元22である。
22次元既約軌道の閉包は全空間の成分となる。
さらに連結群は各既約成分を保つので、代表を通る任意の成分はこの軌道を含み、
次元が22に限られるため同じ軌道閉包に等しい。
全七変数源の rank四部分がこの一軌道であることから、全源が被覆される。

## 3. 全源と相違の適用範囲

`MODULAR_NIL22_TOP_SUPPORT.md` §5 の `w−u,z−uy` と中心化群 shear は全パラメータで定義される。
その後の `δ=−3c/4` は Pxy の `z−yw` を消し、
`e0=−3(b−δf)/5,c0=−3(d−δe)/2` は残る w 項と定数項を消す。
`T=z+3fw−3e/2` は Pyz の affine 項を消す。
これらの変換は modular 場を保ち、Pfaffian は tw、rank四は t≠0 と同値である。
二つの Jordan block の相対 scale で t は1へ移る。
この手検算と、既存の全係数 producer 恒等式を区別して入力として記録する。
今回新たに全七変数の polynomial identity を別実装で再実行したわけではない。

今回の軌道の実因子は `4H+L`、H≠L。
旧70の同次元22成分の一般因子は `H+2Q,3H+Q,2H+C3,2H+L+Q,3H+2L` であり、
どれも `4H+L` ではない。異なる次元の成分、一般 rank二成分とも一致しない。
したがって既知の大きな graph に対する非所属検査は今回の相違証明に必要ない。
四重 H の canonical cofactor が Jordan 型(3,2)、h∈im(N*)²という補足も正しいが、
現行70成分との相違には因子型だけで足りる。

よって新たな個別成分を一つ追加する数学的条件を満たす。
台帳・主原稿・下界の更新は統合担当者が行う。

## 4. 共通runner更新後の歴史証明書の再現

独立検証後、共通 producer は複数法線と CLI に対応するため更新された。
上記 producer 証明書と旧 job は不変であり、その `producer_bytes` が検証当時の正確な入口である。
現在の共通 runner の hash を旧 job の hash に読み替えない。
fresh clone から旧 bytes と最小入力だけを一時ディレクトリへ置けば、現行ファイルを変更せず独立検証を再現できる。

```python
from pathlib import Path
import hashlib, json, shutil, subprocess, tempfile

root = Path.cwd()  # repository root
scratch = Path(tempfile.mkdtemp(prefix="poisson-nil22-review-"))
(scratch / "code").mkdir()
(scratch / "reports").mkdir()
packet_name = "CLASSIFICATION_MODULAR_WEIGHTED_SLICE_20261107.json"
packet = json.loads((root / "reports" / packet_name).read_text())
old_bytes = packet["producer_bytes"].encode()
assert hashlib.sha256(old_bytes).hexdigest() == packet["producer_sha256"]
(scratch / "code/classification_modular_weighted_slice.py").write_bytes(old_bytes)
for name in (packet_name, "CLASSIFICATION_NIL22_EQUIVARIANT_20261107.json"):
    shutil.copy2(root / "reports" / name, scratch / "reports" / name)
checker = "classification_nil22_weighted_slice_review.py"
shutil.copy2(root / "code" / checker, scratch / "code" / checker)
subprocess.run(["python3", str(scratch / "code" / checker)], check=True)
print(scratch / "reports/CLASSIFICATION_NIL22_WEIGHTED_SLICE_REVIEW_20261107.json")
```

この独立 checker は producer を import しない。旧 producer bytes は依存 hash の固定のために配置し、
数値は公開正規形と曲線の文字式から独立に再構成する。
