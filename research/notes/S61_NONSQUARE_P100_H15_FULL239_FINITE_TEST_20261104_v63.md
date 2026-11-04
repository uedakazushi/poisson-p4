# P100 の H15 全有限候補における実四次入力の一回限定検査

本 packet は凍結 preflight `958a02292e9a88ea4151eb3c077f4198b8e5aaf84dcc041fdcabb1cfea691109` に従う一回限定の実計算である。結果は `report.json` を参照し、途中停止・欠損を標数零の欠損と解釈しない。

入力 Hbar の全七因子（次数4,1,1,1,2,3,3）を順に扱う。新しい点や素数は使わない。各体は既存の形式的有限代数 F_32003[z]/f である。全319源は元の174二次積、112実三次積、33原quarticを保持し、受理済みの同一 cokernel 操作だけで pure98 の Gamma に移す。

曲線は a=-10(z+20)/(9z²), e=(z+100)/(9z), d=-3z(z+20)/100。全 literal QQ 係数の分母が32003単元であることを入力時に検査する。各因子で五つの曲線極を検査し、以後の全ての逆元はその体で積を照合する。逆元の値と逆元を binary に逐次記録する。非単元が一つでも現れれば停止し、その因子を除去しない。

Oracle の物理変数へ diag(Tden) を入れ、原 quartic を D_H² で割った Qhat を用いる。その後に物理 light から joint light へ移し、NF117 を適用する。この順序を逆にしない。H2 源 chi=uQhat77+vQhat131 を一つずつ計算し、実係数 U2=(50aBC)² との一致を確認する。

三次は元の正規化57×91を一係数ずつ読み、全16×91の体係数を計算する。112積の NF122 を使って

    Section=Lactual^t Ddelta7 PiOrig

を形成し、全23×23の H1 section 恒等式を各体で確認する。delta7 は行6を -64000 C²*(行6-行1) に変える。全 pure tail を保持する。112積を Gamma へ補正して書き出した後、H1 lift と112行キャッシュを解放する。

174二次積を処理し、続いて原33 quarticを一根ずつ評価する。全33根を同時に保持する方式は禁止する。各行は

    Gamma=f0-(f2/U2)chi0-(f1-(f2/U2)chi1)Section

であり、元の319源番号を保って一時ファイルに書き出す。全係数キャッシュが破棄され live field object 数が零になってから階数処理を開始する。

係数段階は107743 prime-field slots、階数段階は81712 slotsを上限とし、さらに実行全体60秒・2GB・成果物32MiBを守る。多項式は疎な係数配列で扱い、異なる物理/joint monomial frame の整数索引は通常次数4以下の共通 metadata のみを使う。これは係数を495個ずつ稠密に保持する指示ではない。

正の証明書では選択した98源行、元の98×98体行列、固定pivot列・pivot値・元列順での行列式を保存する。これは独立検査が探索なしで固定pivotを追える形式である。全因子で階数98なら CRT と有限自由 Z_(32003)[z]/H モデルへの Nakayama により、H15 上の標数零の全有限候補で full239 が従う。任意の欠損や上限停止では、その結論を出さない。H15 外の全曲線 full239 や一般アフィン尽くしも結論しない。

実行:

    source /workspace/cloud-setup/activate.sh
    .venv/bin/python code/s61_nonsquare_P100_H15_full239_finite_test_20261104_v63.py \
      --output reports/S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63 \
      --work work/P100-H15-full239-once-20261104

型付き入力・実係数処理・固定pivot行列の source binding は同名のコードと報告書に固定する。自動再試行・別の素数・別の点・候補の拡大は行わない。
