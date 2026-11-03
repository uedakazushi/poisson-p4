# Q_b actual Dの残るcubic socle：全150 QQ(b)恒等零

2026-11-03。point b=1での零と、今回のgeneric QQ(b)恒等零を区別する。
full-D applicationは別のproof packetへ保存する。

## 実sourceとimplicit corrections

同じ実29 target関数をbaseにする三つのactual smooth source（CB1、CB2、CA）を使用する。
全150物理source Jacobianを確認し、CB sourceは常にcanonical(rawCB(f,mu*f,H))でmu*fを保持する。
point helperと新しいQQ(b) helperは、有限nilpotent exp(sum flag matrices)と実逆座標による係数演算を行い、別のLie法とも全150行で一致する。

58次blockは二CB sourceの同一base29、20base cuts、8個のlinear conductor座標、rhoを同時に解く。
87次blockは三sourceの同一base58式、20base cuts、7個のactual CB1−CA regular式、actual N equation f、rhoを同時に解く。
両square inverseをQQ(b)上で厳密に確認し、point matrixへのspecializationも一致する。
20baseをlinear tangentだけで置き換えたり、7個のnonlinear regular constraintを落としたりしない。
sourceのorder-n offsetはoriginal source Jacobianに厳密にlinearに入るため、最終係数の補正はその全150×29 Jacobianで更新する。

normal Taylor次数≤3の全係数は四つのrational slopeによる厳密なpolarizationで抽出する。
これは既知の四monomialの係数演算であり、**bは常にsymbolic**である。
有限b標本による恒等零推論は一切行わない。

## 全ambient odd vectorのregseq division

実base29を除く全121 physical座標を含む全150行を保持する。
CB二sheetの差を2で割ったb-vector、CAのodd sectionとの差から得るc-vectorを使用し、actual f-row18によってc0=1とする。
constant row normalization後のlinear conductor generatorはrows

    18,19,28,29,46,60,67,68

で、quadratic conductor generatorはrow63のlinear部分を除いたものである。
8linearと1quadraticはactual double-S conductorのregular sequenceである。
各ambient b_jをこの9generatorでtelescoping形式に分解する。

8linear generatorのcofactorをdegree2まで求めるにはb_jをdegree3まで保持する。
quadratic generatorのcofactorはdegree1までで十分だが、その理由を省略しない：
対応するactual multiplier c_quadはconstant termが**厳密に零**であり、publisherも明示的に確認する。
したがって未計算のdegree2 cofactorはc_quadと掛けてもdegree3以上で、目的のdegree2 defectへ寄与しない。
残り7linear generatorもc_i(0)=0なので、そのcofactorは目的の評価上でdegree1までで足りる。
この次数boundはambient redundant generatorを無視するものではない。

全150 b-syzygyを、actual CA sectionへ代入したdegree3までの全恒等式として確認する。
全150のc-defectはconstant/linearが零であり、全quadraticがactual Z4の二つのquadrics（rows18,63）のidealへ入る。
従って追加し得る唯一のsocle class f²rhoの係数は**QQ(b)上で全150行について恒等的に零**である。

## 計算の範囲と保存

point actual source計算は4.27秒、generic setupは約7秒、完成generic cubic workerは22.07秒で完了した。
最初のgeneric adapterにあったFracElementの0^0 exceptionは、zeroth powerを1と明記して修正した。
数学的な方程式やsource modelを変更せず、完成したrecovery runは保守的90秒cap内で完了した。
全source・inverse・cofactor・defectを公開gzipへ保存し、実行済みworkerをbyte同一でcodeへ保存する。
通常のpublisherは完成済み計算を再利用して全150恒等式を再監査し、22秒workerを再実行したと主張しない。

全GL×b像は実smooth20 stratumをdominantに覆う。
quadzero後の残差moduleはそのbase上のrank1 free moduleのsubmoduleであり、generic QQ(b)上の全ambient class零はGL-equivarianceによりgeneric stratumで零となる。
domain上のfree moduleへ入るsubmoduleにgeneric torsionを残す余地はない。
これで別packetのfinite-colon/full-D lemmaへの最後のgateが満たされた。
高次normal冪は既知core idealへ既に入るため、blind fourth-jet runを要しない。

再現（保存済みcompleted evidenceの再監査）:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CB_actual_D_cubic_syzygy_20261103_v63.py --output <copied-evidence-directory>

全物理workerの明示的再実行は`--recompute`を付ける。
各workerは2.2GB/120秒以内の有限計算で、GBは使用しない。
