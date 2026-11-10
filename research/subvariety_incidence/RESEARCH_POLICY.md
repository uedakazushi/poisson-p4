# Poisson不変部分多様体による147成分の幾何学的理解

2026-11-10 / v7。再開順序：ルート RESEARCH_POLICY.md → REPRODUCE_GITHUB.md → 本書 → NOTES.md → papers/poisson_subvariety_incidence_ja.tex。

## 上位原則：研究価値の最大化

数学的重要性、新規性、理解の深化、他の成分への波及、正確性、実現可能性を総合して現在の目標を柔軟に見直す。重要な方向転換は元の目標の達成点と残存点を短く整理して現行方針を書き直す。作業履歴を蓄積しない。

## 現行の数学的成果

- **既存の一般Q3分類からの延長**：一般Rat(1,2)の全延長は被約一点（rank2）、一般Log(1,1,1)は被約A1（全rank2）、一般Affは被約二点（rank2とrank4）。特別なintrinsic Poisson元の全分類は未完成。
- **modular因子重複度定理**：strict q-Casimir lift Pのmodular場vがQ上のPoisson葉に接するとき、かつそのときに限って、一般rank4のPfaffian五次DにQが二重因子となる。通常の一般証明は現行主稿。
- **共鳴Logの非退化ブロック**：q=qU+w0²+w1²でqUが非退化なら、rank H=2/1の全源を直交性で判定。追加rank4枝の全JacobiイデアルとPfaffianを論文共有稿に保持。
- **今回の新しい境界：smooth Qでもq|Wが退化する**。非退化qUを直交分解したまま退化させればQも特異になるため、研究対象をq|W（Wは対数三線形形式の共通核）の退化に変更。三つの非退化qによる正確な全延長例を得た。rank(q|W)=1でA²と交叉二直線の両方が現れ、rank(q|W)=0では非被約Spec C[a,b,c]/(ac,c²)となり、被約台のA²には一般rank4延長が存在する。通常の制限空間・発散と全線形stabilizer計算を使う証明は papers/quadric_log_isotropic_boundary.tex（主稿にinput）。10変数の正確な係数行列・Jacobi ideal・Pfaffianはcode/subvariety_incidence/check_isotropic_log_boundaries.pyと同名JSONで再現。これらは三明示例の全延長を証明するもので、全特異制限層の完全分類ではない。
- 価値のある既存成果（cubic一意性、一般CI22 del Pezzo剛性、rational型次元36/28、単純H対称性、既存147成分定理）は現行TeXと正確な証明資産に保持。新規性の一次文献による全面確認は未完了。

## 現在の優先順位

1. **延長スキームのscheme構造とambient147成分の接続**。新しい被約A²/交叉二直線/非被約厚みのrank4開部分が、既存のどの成分に属するかを実Poisson変形・退化因子の分割と比較し、非被約のembedded線がどのような一次・二次障害を表すか解明する。
2. smooth Qでrank(q|W)=1,0の全特殊Log源を、重みの二次形式とJacobiの小さな可解Lie代数により層別する。rankだけではファイバーは決まらないので、rankのみの表は作らない。
3. special Aff、4次del Pezzoの真の例外延長、P3横断Poisson対称性へ波及する一般定理を優先する。複雑な有限場合分けが理解を深めなければ、目標を変更する。
4. 最終的に147成分の幾何的特徴づけへ統合。新規性・先行研究の精査は並行する。

次の着手はSpec C[a,b,c]/(ac,c²)のembedded一次方向の幾何と、被約なrank4開部分の既存成分への帰属を調べる。Poissonスキーム全体の非被約性や新しい全空間成分を直ちに主張しない。

## 再現・保存

GitHub mainが唯一の正本。現行日本語LuaLaTeX（主稿・共有TeX）、方針・ノート、8つのSymPy producerと8つの正確なJSON証明書をmainに保持。PDFは追跡せずlatexmkで生成。GitHub Actionsは使わず、Libraryは実push失敗時の緊急避難だけに使用。数学の通常証明、今回の有限検算、外部査読、新規性を区別する。
