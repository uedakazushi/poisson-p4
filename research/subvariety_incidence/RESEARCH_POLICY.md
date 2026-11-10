# Poisson部分多様体と射影四空間の成分交叉：研究方針

2026-11-10 / v11。再開順：ルートRESEARCH_POLICY.md → REPRODUCE_GITHUB.md → 本方針 → NOTES.md → papers/poisson_subvariety_incidence_ja.tex。

## 恒常的原則と保存許可

数学的重要性・新規性・理解の深化・正確性・実現可能性・他成分への応用価値を総合し、研究目標を柔軟に調整する。主要結果の仮定、証明、計算検証、未完成箇所、先行研究との差を分ける。公表原稿は日本語LuaLaTeX、ノートは非重複の研究知識のみとし、価値のある内容・再取得困難な入力・証明書は軽量化で失わない。履歴や旧版を蓄積せず現在の一式だけから再開できるようにする。ユーザーの継続保存許可の現行運用は**GitHub mainを唯一の正本**とし、非force push後にremote先端を確認する。ChatGPT Libraryへの保存は実際のpush失敗時の緊急退避に限る。GitHub Actionsは使わず、生成PDFは原則Git追跡しない。Library上の既存ファイルは明示的指示なしには削除しない。

## 現在の中心成果とその確実性

- **局所次元28（今回の独立証明）**：完全等方Log点p=[P0]、q=u0²+2u1w0+2zw1について全ambient Poissonスキームの接空間は30次元、局所Krull次元は正確に28。二つの互いに素な二次障害式9*t21*t29と-4*t25*t27が接錐を拘束し、Log(1,1,1,2)源の28次minor（-57177414000000）で下界を得る。全147成分の完全性はこの次元決定に用いない。
- **25次元rank4既約成分の明示帰属（今回の独立同定＋既存成分性）**：P_a=P0+aK1 (a≠0) の六つのアフィンPoisson括弧は、実際のPGL5座標変更によって既存rank一長鎖・非定数g源に完全に一致する。g=6Z、U=60Z/a² ∂W、K=-ZW/a²+50Z²/a⁴。従ってpはカタログのR4_MODULAR_RANK1_LONG_NONCONSTANT_25に属する。modular冪零Jordan型はa≠0で(3,1,1)、a=0で(2,2,1)。
- **非等次元・非Cohen--Macaulay（既存成分性に依存する系）**：pはrank2 Log1112（28次元）と上記rank4（25次元）の相異なる既約成分の交点。局所環には異なる次元の極小成分が存在するので非等次元、非CM、非Gorenstein、非局所完全交叉。ambient局所環の非被約性は証明していない。
- **保持する以前の有効な数学**：固定smooth Qとσの非被約延長ファイバーSpec C[a,b,c]/(ac,c²)、全ambientでも消えないc方向の二次障害、一般Q3 Rat/Log/Affの延長（一点／A¹／二点）、共鳴Log直交性wallと等方的境界、smooth cubic一意延長、一般CI(2,2) del Pezzo剛性、P3横断対称性、全147既約成分分類の現行証明と計算入力。新規性・優先権の全面確認と外部査読は未完了。

## 次の核心と具体的着手

第一に**pを通るrank4成分が長鎖25次元成分だけか**を判定し、二つの既知成分の接錐・交叉スキームを同定する。接空間上の二つの単項式二次障害を出発点として追加二次関係の線形代数と、カタログの閉包源による実退化の両側から検査する。第二に局所環の深さ、embeddedな部分、固定Qファイバーのac方向がambientで真の障害かを独立に検証する。ambient非被約性は二次障害・非CM性だけから推論しない。研究価値が局所場合分けによって低下すれば、del Pezzo例外・P3対称性やより普遍的なPoissonスキームの非CM機構へ切り替える。

## 論文・再現

main原稿 papers/poisson_subvariety_incidence_ja.tex に二つの主証明 papers/quadric_log_component_intersection.tex、papers/rank_four_long_component_membership.tex をinclude。従来の障害は papers/quadric_ambient_quadratic_obstruction.tex、固定ファイバーは papers/quadric_log_isotropic_boundary.tex。今回の正確なproducer code/subvariety_incidence/check_rank4_long_nonconstant_membership.py、証明書 research/subvariety_incidence/certificates/rank4_long_nonconstant_membership.json。手動全再現：bash scripts/verify_github_source.sh incidence（11 producer/11 JSON、LuaLaTeX）。GitHub正本の確認と計算・組版の実行状況は分けて報告する。
