# Poisson不変部分多様体と147成分分類の幾何

2026-11-10 / v8。研究再開：ルート RESEARCH_POLICY.md → REPRODUCE_GITHUB.md → 本方針 → NOTES.md → papers/poisson_subvariety_incidence_ja.tex。

## 上位原則：研究価値の最大化

数学的重要性・新規性・概念的説明力・他成分への波及・正確性・実現可能性を総合して研究目標を柔軟に再評価する。重要な変更は現在有効な方針として書き直し、元の目標の達成点・残る核心のみ保持する。研究の作業履歴は蓄積しない。

## 現在の数学的成果と状態

- **一般の滑らかなQ3のintrinsic三型**：Rat(1,2)は被約一点（rank2）、Log(1,1,1)は被約A¹（全rank2）、Affは被約二点（rank2とrank4）。一般元の定理であり全特殊点の分類ではない。
- **modular因子重複度**：strict q-Casimir lift Pのdivergence vがQ上の葉方向へ接するとき、かつそのときに限り、一般rank4の実PfaffianにQが二重因子となる。Q+C3対2Q+Hの幾何学的判定。
- **共鳴Logの非退化ブロックと直交性wall**：rank1/2の全qU非退化配置を現行共有TeXに証明。rank4枝の正確なJacobi idealとPfaffianを保持。
- **滑らかなQ上の等方的境界**：制限q|Wが退化した三つの明示例の全延長はA²、二直線、非被約Spec C[a,b,c]/(ac,c²)。これは特定三例の計算であり、全退化配置の分類ではない。
- **今回の新定理：ambient Schouten二次障害**。上の完全等方例のa=b=c=0においてc方向は一次Poisson変形を与えるが、固定Qやσの制限を全部外しても二次Poisson変形へ持ち上がらない。全150係数の射影Jacobi線形化はrank95、二次障害でaugmented rank96。五係数汎関数が全ambient線形化を消し、二次障害に-9を与える。対応するambient Poissonスキームの点は特異（接空間dim30）。**ambient非被約性は未証明**。共有 `papers/quadric_ambient_quadratic_obstruction.tex` に通常の証明と5×6行列を収録し、`code/subvariety_incidence/check_ambient_quadratic_obstruction.py` と証明書で全係数を正確検証。
- 以前のsmooth cubic一意延長、一般CI(2,2) del Pezzo剛性、次元36/28、P3の単純H横断対称性、147成分の全分類と有限証明書は現行数学資産として保持。既存研究と本研究の新規性を区別し、外部先行研究の網羅調査と外部査読は未完了。

## 次の優先順位と具体的着手

## 二次接錐の現行結果と残る問題

既存の五係数消滅汎関数はambient Poissonスキームの30次元Zariski接空間上に非零二次式を与え、完全等方Log点の局所次元を29以下に抑える。固定quadric延長の三方向(a,b,c)への制限は正確に -9c²。36次元のR2_RAT113と39次元のR2_LPBはこの点を含まない。証明はpapers/quadric_ambient_quadratic_obstruction.texの系cor:ambient-isotropic-tangent-cone-bound。三方向の直接展開は今回独立に検算した。全接錐の成分分解、局所次元の確定、固定ファイバーのacのambient障害、ambient非被約性は未解決。次はカタログの次元29以下の族と退化因子Q+2H+H'を照合し、帰属を特定する。


1. **ambient二次障害と既存147成分の局所幾何を比較**：P0を含む既約成分を現行カタログの実源・退化因子型から特定し、全ambient tangent coneのquadratic idealと既知成分の局所切片を比較する。非被約性の証明には局所環のnilradicalを独立に示す必要がある。単に二次障害が非零・接空間が過大という理由で非被約と主張しない。
2. Q3の特殊Log/Aff延長のより普遍的な幾何層別：重み空間とqの直交・等方配置、strict stabilizer、Jacobi二次障害をまとめる。無意味な有限場合分けが増えるだけなら優先度を下げる。
3. 4次del Pezzoの例外的rank4延長、P3の六成分の横断対称性、147成分の概念的特徴づけへ波及する定理を探す。研究価値が高ければ方向変更する。

次の一手はP0における周囲局所二次接錐の幾何学的分解と、既知147成分の一般元への退化族を比較すること。元の目標（特殊Q3延長の全分類）は依然未完成である。

## 正本・再現

GitHub mainが唯一の正本。現行LuaLaTeX・共有TeX、研究方針・ノート、9本の正確なSymPy producerと9件のJSON証明書を追跡し、`bash scripts/verify_github_source.sh incidence`で再生成する。PDFは追跡せず、latexmkで組版可能。GitHub Actions禁止、ChatGPT Libraryは実際のpush失敗時だけ緊急退避。数学的に証明したこと、今回の厳密計算、未実行の全分類再検証、外部査読を分ける。
