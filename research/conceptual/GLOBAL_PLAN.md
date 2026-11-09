# 最終成果までの連続実行計画

2026-11-09 / v72-normalization-residue-integration

研究方針を先に読む。最終定理と終了条件は `FINAL_TARGET.md`、現行の証拠・残課題は
`research/classification/CURRENT.md` に置く。本書の各段階は中間工程であり、ユーザーへの承認待ちや終了点ではない。

## 1. 変更しない主軸

全点の認識 → 標識付き復元 → 実際のPoisson変形 → 極大な相異なる閉像。
147は現行稿の結論で、計画の前提ではない。固定因子の復元と動く因子の持上げを分け、
型の認識と全空間の成分性を分ける。一般論だけ、一覧表だけ、検算だけでは完成としない。

Aは超平面因子あり（5Hを含む）、Bは超平面なしの既約二次＋三次、
Cは被約既約五次（normal/nonnormal）、R2は一般階数二である。
この区別は成分閉包を互いに素に分割するものではない。

## 2. 再利用する現在有効な共通証明

| 現行原稿 | 統合されている内容 | なお区別する依存 |
|---|---|---|
| marked_hyperplane_reconstruction.tex | 一意な標識affine模型、modular原始関数、単純H/重複Hの共通入口、共鳴 | 全源の有限認識と全域の帰属 |
| poisson_lifting_examples.tex | TRIPLE→Ga9の次数補正、可換場による超平面分裂 | 任意方向の持上げや全境界の分類ではない |
| quadric_cubic_reconstruction.tex | rank3/4/5 cofactor、平面五次の全源、因子留数による磁気商 | 31保存源の完全性入力、全有限配置・成分性 |
| conductor_conormal_reconstruction.tex | rank0/2の可換余法線、planeの全weight源・隠れた三係数、leading symbol輸送 | 非plane認識、normalization/降下、全higherの解と閉像 |

| normalization_residue_reconstruction.tex | 被約fixed-D線形開、全正規化留数とB/A降下、scroll交点kernel、部分導手・一意相対lift、rank2 pencil | 因子baseの全rank strata、非被約A、全極大閉像 |

パスは全て `papers/` 内。主稿と構造稿はこの共有本文を参照する。
既存の正しい成果を再発見するために戻らず、実際の仮定・証明を確認して再利用する。
新規性未確認と、証明済み／今回確認した範囲の違いを維持する。

## 3. 実行する工程（全てを一つの課題として続ける）

### 工程A：非plane導手の認識と復元を閉じる

`research/classification/IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md` と
`NONPLANE_CONDUCTOR_POINTWISE_KNOWN_COVERAGE.md` から原証明を辿る。
固定した正規化・導手square、Hamiltonian導分の降下、逆形式の留数を定式化し、
Σδ_i deg S_i=6−gの数値制約、restriction rank0/2、各低次数曲面の全源をつなぐ。
元の逆写像と全zero complementを保持し、少なくとも重複する認識・降下の議論を共通命題へ統合する。
planeの完成源を非planeへそのまま延長せず、必要な追加条件を証明する。
導手の被約性、family正規化のflatness、正規化後留数の単極性を仮定に足さない。
normalの認識は `papers/normal_actual_divisor_recognition_ja.tex` の適用範囲を確認して接続する。
終了時にCの全点が有限閉像へ至る証明を本文に置き、**工程Bへ続行する**。
被約fixed-D fiberは全compatible線形空間からquadricを除いた開集合として復元した。
正規化上の導手pole bound・全B/A降下、scroll全留数kernel、有限MCM部分導手の相対liftを共有TeXへ統合した。
次数五・六と二つの実T40曲線の重複を置換し、平面型の全10係数とrank2 determinant境界を明示した。
一般因子が被約な相異なるambient成分の閉因子像は相異なるが、rank jumpでの像包含は別に比較する。
詳細はpapers/normalization_residue_reconstruction.tex、検算37項目はreports/normalization_residue_checks.json。
これは新しい中間証明であり、最終条件を満たすまでR2・極大閉像・全依存・最終保存へ続行する。


### 工程B：rank 2支持を接続する

`RANK_TWO_SUPPORT.md`、`KNOWN_COMPONENT_AUDIT.md`、独自stable_rank2の元証明を読む。
飽和foliationと零因子（d+e=2）から各全源を復元し、安定分岐、非安定分岐、因子付き分岐を尽くす。
六つの全空間成分と、toric外積・Aff32等のrank4閉包に入る族を区別する。
Corrêa–Munizの代数的可積分性だけを固定Jacobian族への帰属と読まない。
rank4閉包項を全rank4被覆の後に代入する依存順を保持し、**工程Cへ続行する**。

### 工程C：全点被覆から真の全成分へ

`GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md` と各全支持・成分性証明を、
共通復元と実変形の枠組みに接続する。標識重複、特殊底、Galois結合を戻し、
有限閉被覆の各項の極大性・相違、不要な候補の吸収を証明する。
全147×147の交叉表ではなく、極大閉像の決定に必要な有限比較を選ぶ。
主定理の全IDについて構成・次元・階数・幾何・証明参照を一対一に確認する。
共通判定で処理しきれない比較は、有効な有限補題として完全な証明・証明書を残す。
個数を導き、誤りが出たら依存先を訂正してこの工程を閉じ、**工程Dへ続行する**。

### 工程D：公表原稿を完成させる

主定理と幾何的総覧、全入口、共通復元・持上げ、A/B/Cの認識、R2、極大閉像と一般元表の順に統合する。
共通定理を挿入しただけの古い構成を残さず、置換可能な重複は参照にする。
有効な全結果・十分な証明は現行原稿または同時配布論文に保持する。
先行研究の一次本文・仮定・適用範囲を照合し、Introductionで貢献と非主張を具体的に説明する。
H+2Qの0/1/6の局所理論、全軌道、全scheme厚みを完了の必須条件へ格上げしない。
最短・最少を目指す無限改稿に入らず、FINAL_TARGETのC1/C2を満たしたら**工程Eへ続行する**。

### 工程E：必要な証明確認・再現・保存を完了する

新しい主証明の依存を読み、主定理に必要な有限計算を正確な入力とともに確認する。
旧全実験を無差別に反復する必要はないが、旧PASSを見ただけで未確認の推論を受理しない。
変更原稿を組版し表示・参照・引用を確認する。現在の依存閉包からbundleを構成し、
必須証明書の未同梱、旧ZIPへの参照、旧ノートだけにある証明、切れた参照を解消する。
`make_structural_bundle.py` の現行の「証明書依存閉包は不完全」という制約は、実際に解消するまで残す。
manifestの文言だけを完了に変えない。非force main更新を確認し、Library保存を別に確認する。

最後にFINAL_TARGETの全条件を証拠付きで見直す。未達があれば関連工程へ戻る。
**全条件達成時だけ通常の最終報告を出す。**

## 4. 方法変更・中間保存・実行制約

工程順は依存に応じ変更してよい。行き詰まりは別法・反例・より弱い十分命題・独立な必須課題への切替えで処理する。
二回続けて同じ方法が進まないことは方法変更の目安であり、停止条件ではない。
中間のcommit/pushとLibrary保存は行い、その直後に次の未達課題へ進む。
未完の数学的義務を単に範囲外へ移して終了しない。
実際の実行上限・権限等による不可避な中断はFINAL_TARGET §5に従い、未完を明記して再開可能な状態を保存する。
