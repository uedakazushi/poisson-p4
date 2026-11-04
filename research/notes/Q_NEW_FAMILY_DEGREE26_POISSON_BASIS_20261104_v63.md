# Poisson括弧による全Inv₂₆の45本のQQ基底

2026-11-04。**普遍的なSp₄不変式空間Inv₂₆の45本のQQ基底を構成した。** 全16方向のambient gradientからPoisson括弧を作り、同じ固定形式ambient二平面上の係数minorとcharacter上界で証明した。NEW25のsourceへの26次制限は計算しておらず、定義次数は26以上のままである。

## 入力と微分の規約

4×4 C,Dの全136 quadratic monomial係数、元のΩ、旧probeの最初の三つのambient β vectorを保持した。平面は前の24次検査と同じβ₀+sβ₁+tβ₂、素数は32003。actual NEW25の点は使わない。

**C,Dを周囲の16個のβ座標で微分してから平面へ代入した。** βᵢβⱼをβₖで微分する係数tensorを全k=0,…,15について構成し、対角βₖ²の係数2も保持した。s,tだけの微分で代用しない。全32 quadratic成分と全8 trace invariantについて、Σβₖ∂ₖf=deg(f)fの全係数恒等式も確認した。

Poisson括弧の符号は元の規約

\[
 \{f,g\}=\nabla f^{\mathsf T}\Omega^{-1}\nabla g
\]

で、Ω⁻¹を転置・符号反転せず入力した。traceのgradientはcyclic product ruleで構成した。

## 選択した14次と18次の基底

次の5本をInv₈の基底とする。

\[
 (T_0,T_1,T_2,T_3,T_4)
 =(\operatorname{tr}C^4,\operatorname{tr}C^3D,
 \operatorname{tr}C^2D^2,\operatorname{tr}CDCD,
 \operatorname{tr}CD^3).
\]

全10個の括弧を係数比較し、3本

\[
 Q_0=\{T_0,T_1\},\quad Q_1=\{T_0,T_2\},\quad
 Q_2=\{T_1,T_2\}
\]

を選んだ。14次係数120行の選択3×3 minorは **7584 mod32003**。

12次の三つのC/D wordを

\[
 (U_0,U_1,U_2)=(\operatorname{tr}C^3DCD,
 \operatorname{tr}C^3D^3,\operatorname{tr}CD^2CD^2)
\]

とする。N₁₂との差はF₄Inv₈に入り、F₄は不変式のPoisson代数の中心に属する。従ってこれらを使った候補は、N₁₂を使った候補とF₄Inv₁₄を法として同じspanである。

18次では3本F₄Qᵢを先に置き、15個の{Tᵢ,Uⱼ}を追加した。全rankは9で、最初の3本に加え

\[
 (P_0,\ldots,P_5)=
 (\{T_0,U_0\},\{T_0,U_1\},\{T_0,U_2\},
 \{T_1,U_0\},\{T_1,U_1\},\{T_1,U_2\})
\]

を選んだ。190係数行の選択9×9 minorは **26750 mod32003**。従ってP₀,…,P₅はF₄Inv₁₄に対する6次元の補空間である。

H₁₆のgradientは一切使っていない。旧θ⁵H₁₆は後述のInv₁₂基底に、微分しない因子としてのみ残る。

## 57個の26次候補と45次元基底

受理済み全Inv₁₂の9本をB₀,…,B₈とする。順序は`Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS`で固定したF₄Inv₈の5本、N₁₂、C/Dの二word、θ⁵H₁₆である。次の57列を作った。

\[
 Q_iB_j\quad(i=0,1,2;\ j=0,\ldots,8)\quad(27列),
\]
\[
 P_iT_j\quad(i=0,\ldots,5;\ j=0,\ldots,4)\quad(30列).
\]

26次のsⁱtʲ係数は378行。選択45行は0,…,44、選択45列は0-basedで

```text
0..35, 37,40,42,43,44,45,47,52,55
```

であり、その45×45 minorは **7998 mod32003** で非零。全57列・semantic product定義・全8×16 gradient配列を保存した。Python側の全係数消去でも選択行・列・minorが一致した。

同じ16 weightsのSym冪とWeyl8項を26次までexact integer arithmeticで計算し、14,18,26次の次元3,9,45を得た。全Weyl項はcertificateに同梱している。候補の全分母は32003で単元なので、非零有限体係数minorはQQ上の独立性を証明する。この下界がcharacter上界と一致し、**選択45本は全Inv₂₆のQQ基底**となる。同様にQᵢ、F₄QᵢとPᵢは全Inv₁₄、全Inv₁₈の基底である。

これは指定次数の全空間の基底構成であり、不変式環全体の生成・表示やNEW25の方程式の構成ではない。57列には普遍関係が12本ある。将来のsource関係検査では、非零57-vectorだけを周囲で非零な方程式と扱わず、固定45本を使う必要がある。

## 実行・保存・残るsource問題

native一回約0.0213秒、wrapper約2.89秒。polynomial products8210回、最大support378。60秒/2GB/certificate32MiBの上限内で完了した。保存全体約138KB。

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_degree26_poisson_basis_20261104_v63.py \
  --output work/degree26-poisson-replay-report --work work/degree26-poisson-replay-native
```

output/workは新しいdirectoryを指定する。一次保存は
`reports/Q_NEW_FAMILY_DEGREE26_POISSON_BASIS_20261104_v63/`。
全semantic定義、全係数、選択minor、全character記録、入力・証拠hashを保持する。producerは同名Python/C++。C++は凍結済みnmod helperをincludeし、旧mainは実行しない。26次の積に限り二変数378項のguardを使う。

sourceには新たにQ₀,Q₁,Q₂、P₀,…,P₅の9 mapsが必要だが、未展開である。選択結果によりsource上で必要な8次gradientはT₀,T₁,T₂の三本だけになった。全16方向の微分は依然必要。B¹³H²⁶、c¹³、(u,X)奇parityとu=aX後の共通X除去というpreflightの分母・奇偶性は維持する。

18次mapの理論support上界124730は出力100k項capを超え、全9 mapsのraw保存上界も32MiB内とは証明できていない。別のbounded budgetと逐次保存・停止条件を設定するまで、source26展開や45列のrelation検査を自動開始しない。NEW25の定義次数は依然26以上、方程式・正確な次数は未解決である。
