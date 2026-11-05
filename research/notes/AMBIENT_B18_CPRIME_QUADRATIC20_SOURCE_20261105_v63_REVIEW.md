# Cprime の完全二次 initial 写像：独立 SOURCE 監査

2026-11-05。本稿は SOURCE gate の独立監査であり、canonical gauge・新 source 値・二次 rank の計算を含まない。Primary `5c3bf716…`、選択入力 `bb05969e…`、未実行 worker V2 `8c0c791a…` と保存証明書 `5b2086b7…` を固定し、下記の必須 runtime 検査を伴う一回の値 job に対する SOURCE gate を受理する。

受理済み full-nine 行列 V は原点で零であり、一次係数行列 `V1:QQ^44 -> QQ^27` は rank24 を持つ。保存済み24次 minor inverse を用い、非 pivot 二十列を補正して全定数 kernel の基底 K20 を作る。この二十列は全二十七行を消し、非 pivot 座標への射影が I20 となる。実際の非零係数は四十九個だけである。

任意の formal source 係数を `a0+a1+...` と書く。出力の一次項が零となる条件は `V1*a0=0`、二次項は `V2*a0+V1*a1` である。逆にこの形の二次項はすべて実現する。従って既知 `S1*im(V1)` を割った後、全二次 initial 関係は厳密に `pi4*V2*K20` で尽くされる。これは補助的な conductor 八十四列だけの上界とは異なり、実際の HF2 を決める完全写像である。ただしその値は本 SOURCE 監査では計算しない。

固定 pi4 は target0/2 の `rA²,rA*ell,ell²` だけに支えを持つ。係数函数を `B27=QQ[z0,...,z6,rA,ell]/(m³,z_i*z_j)` に取れば、定数一個・全九一次・十七二次の二十七基底となる。`B26=B27/(ell²)` には既存 actual Artin33 から algebra homomorphism があり、全1089積を独立に確認する。pure-plane 六基底は `1,rA,ell,rA²,rA*ell,ell²` である。既存33の五係数に欠けるのは ell² のみだが、full-nine CB adapter の二次曲率と全一次係数は捨てない。

共通105は両族各57 source の合併で、共有九原 physical 座標 ID は 2..10 である。各族57式の最後九本は F9 target 自体なので、全114検査ではそれぞれ full YA/YB を引く。CB zero7 の二次曲率は adapter 全九行から保持する。source DAG 閉包は20383 node、単元逆演算の原点値はすべて既存証明書に固定され、last-use peak は13689 QQ である。

未実行 V1 `1c3d6584…` は旧164 root との検査で CB f63=28539 だけを保持せず、163 本を通る構成だった。この狭い custody 不一致を値計算前に検出し、V1 と snapshot を保存した。V2 は28539の保持と `overlap==164` だけを追加する。既に閉包内にある node なので20383/13689は変わらず、全164の literal overlap を要求できる。凍結済み V1 を書き換えない。

FULL9 全 JSON の310866個の数値列を保持しない。独立 raw-decode は一つの選択 block 29310 numeric leaves を読み、source2 tensor6270個とラベル6733個へ直ちに縮める。抽出 phase と後段 arithmetic phase は別の生存期間を持つ。同じ J30 の両側 inverse 全1800成分、det=4/9、group15・fixed fiber10・free fiber9・proper kappa9 frame を固定する。旧 B26 gauge を ell²=0 の初期 seed にするだけであり、未知の実 ell² 値を零とはしない。三十個を補正後に全30制約と全 B26 proper overlap を確認する。

full proper polynomial/rational formula を monic trace algebra で評価し、nu の全二次 carry を保持してから成分を読む。base 最大 ideal の三乗が零で、nuA,nuB がその ideal にあるため trace 最大 ideal の七乗も零である。D=1-d²、P=5+d² は既存同 frame で residue1,5の単元であり、未知因子の逆数は使わない。

whole 生存上界118472 QQ は、600 physical root・nu、共通105、canonical30、43200の有限 GL scratch、保存値と metadata を含む。2048 bit・60秒・2GB・32MiB guard と phase 保存を必須とする。worker は入力部分だけを選択し、旧 source2 大配列、旧 source33/adapter と一時値を使用後に解放する。新しい physical source の implicit solve、source3、旧 r/g 除算を含まない。

この受理は値計算の SOURCE gate に限る。実際の4×20 rank、HF2、長い ell-torsion 指数、nilcube は未計算である。optional cubic はこの worker に含めない。後続で実 h2<=3 が判明した場合にも、その時点の新 SOURCE gate を別に確認する。
