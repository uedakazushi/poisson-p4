# Whole QQ pair gate の独立監査（2026-11-03）

対象は `S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63`。実際の四二次式から得た pair の rank は、全 `a*P*J!=0` 上で2である。先行 finite theorem の bytes は変更せず、今回の空性を別の監査として保存する。

独立 verifier は元の804係数、`A=x0=W=0` かつ **B 自由** の三変数制限、全12係数、既知 unit のみの行正規化を再構成する既存の独立 finite review を別の出力先で再生する。特殊値 `v=-3,-4/3` の全小行列式と QQ Bezout 恒等式も再確認する。著者の source helper や今回の verifier は実行しない。

今回追加された全恒等式は `Fraction` 辞書の係数和・積と一変数の長除法で検査する。resultant、因数分解、GB、有限体計算は生成しない。

自然座標 `t=a/P, v=3e+20t, H=P+8t` は既知 unit の開集合で可逆である。`f=(v+3)(3v+4)` が非零の残りでは、元の小行列式から得た三つの線形 t 方程式と三つの eliminant が必要条件になる。保存された射影多項式について全係数で

`R(v)=k*f(v)^2*C(v)^2*Q23(v)`, `k∈QQ*`, `C=8v³+30v²+71v+58`

を検査した。従って残る幾何学的点は C=0 または Q23=0 に属する。C の重複は点の被覆を変えない。

C=0 では `beta3=-f*C` と t の可逆性から alpha3=0。必要条件 N01=0 に対する保存 membership `U*alpha3+V*N01=Rc(v)` と `Rc*inverse=1 mod C` を独立に検査し、矛盾を得る。

Q23=0 では、保存された `UA*A+UB*B=L(v)H+K(v)` と L の明示 inverse により H=Hbar が強制される。さらに alpha3(Hbar) の明示 inverse により t=tbar が強制される。ここでは QQ[v]/Q23 を有限代数としてのみ扱い、既約性、squarefree 性、体であることを仮定しない。全三本の線形 t 方程式と全三つの eliminant がこの評価で零になることも検査した。

しかし元の全12係数から **第0列と u²列の actual 2-minor** を再構成し `(tbar,Hbar,v)` で評価すると、その剰余 m に対して保存 inverse が `m*inverse=1 mod Q23` を満たす。rank failure なら minor=0 なので矛盾する。この inverse は追加の未知因子を開条件にしていない。

これで QQbar 上の rank-failure 点は全て排除された。係数環は QQ 上有限型の `aPJ` 局所化であり、maximal-minor ideal が全ての閉点で unit となる。Nullstellensatz と局所 Nakayama により pair の cokernel は零である。幾何学的点での判定から非被約 base の module も扱うこの段階と、Q23 の有限代数内での明示 inverse の段階を区別する。

独立監査済みの binary-D 縮約により、補助五二次式は全 aPJ-open 上で独立となる。さらに同じ独立 binary review の線形 syzygy 論証、および既に監査済みの高 rank proportionality 排除と有限 Σ 全域の QQ raw module surjectivity を合わせると、**FULLU=`a*P*F*B*C*J!=0` 上の raw24×56 module は全域で surjective** である。有限 Σ の証明に限って用いた proper finite Z_(32003) / Nakayama と、今回の純 QQ な全域 pair 排除を混同しない。

元の product41 への17-kernel bridge は別の独立 receipt が担当する。本監査はその適用、full57、quartic gate、ambient removal、全 Poisson 成分の尽くしを昇格しない。

再現:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_pair_empty_QQ_review_20261103_v63.py \
  --output work/marked-square-v63/pair-empty-independent-new
```

新しい出力先が必要である。保存済み author/review 入力は一切変更しない。公開 report の `all_passed` が実行結果である。
