# K22の20段partial split：独立固定trace監査

2026-11-04。対象は frozen `S61_NONSQUARE_KERNEL22_KNOWN_UNIT_SPLIT_20261104_v63`。
受理済みの実K174×22とsource304を入力とする。追加pivot・minor・Bézoutの探索、
元四次行列・GB・新点の計算は行わない。

元Kの全係数とnative inputを独立に比較し、native stateのA20,L20,U20と全174×22
Schur係数を読む。別のreview用native codeで **保存20 pivotだけ**を指定順に
再演算する。演算utilityは既受理のFLINT QQ(a,e,d) headerを明示依存とし、
候補探索loop・pivot選択heuristicは実行しない。各pivotはlogの非零QQ定数と
a,d,F,B,C,Jの整数冪を直接乗除した値と等しいことを確認する。
全3828 Schur entry、全A/L/U、およびA=L Uを厳密に比較する。

Lはunit lower triangular、Uはupper triangularで対角は既知単元。
保存された2044-node/73-atom回路について、Lのforward solveとUのback solve
が全20列の指定nodeを与えることを検査する。従って

    P20=U^(-1)L^(-1)π_sigma,  P20 K[:,tau]=I20.

これは有限exact回路の同定であり、P20の全core係数を新たに展開したものではない。
特にP20が全22列のleft inverseであるとはしない。

選択20 source座標のA20は単元minorなので、対応20関係を引いてそのsource座標を
残りへ解ける。以前の15削除と合わせ、指定284個の元source行だけを残す
**固定M284が元と同じimage/cokernel**を持つ。全indexを照合する。

残る二列はKoszul(3,4)のcolumn12と新22本目のcolumn21。
`Krest−Kselected P20 Krest` を残る154座標へ写したものが保存Schur residualであり、
全154×2 entryをfrontier recordと比較する。元Kのpointwise rank22とこの単元
部分splitから、residualはFULLU上でpointwise rank2となる。
これは新しいrank実験ではない。二行minorの非零開が覆うが、固定minorや
最後のP2/P22を計算済みとはしない。全195非零entryと統計をread-onlyで確認し、
提案されたcandidate pairのdeterminantは一つも評価しない。

記録の「個別known-unit entryが無くなった」はauthorの停止理由である。
本レビューは追加候補を探索せず、entry combinationやunit minorの不存在を
結論しない。現在の282 source bundleは
`coker(R2→R154)⊕R130` とも表示できる。双対D45は既受理codim2定理と同じ
pd≤1/Ext¹(D45,R)=actual Cという論理を持つが、これでCのcurveや点は消えない。

fullP22、単一global282 submatrix、full239、C122 rank64、残る43方向、
affine exhaustionはいずれも未取得のまま。作者のimmutable packetを変更しない。

再現には既存 `/workspace/cloud-setup/activate.sh` をsourceする。最初の未activation
入口はpkg-configがflint.pcを見つけられず、compile/native算術の前で停止した。
既存activationを使った正しい入口で一回のnative固定traceを実行し、環境の再構築はしない。
